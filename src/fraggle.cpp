#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <cmath>
#include <deque>
#include <docopt/docopt.h>
#include <iostream>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <set>

#include "scindo/bam.hpp"
#include "scindo/fastq.hpp"
#include "scindo/files.hpp"
#include "scindo/kmers.hpp"

using namespace scindo;
using namespace boost::math;

using counts = std::vector<size_t>;
using position_counts = std::vector<counts>;
using distr = std::vector<double>;
using position_distr = std::vector<distr>;

namespace // anonymous
{
const char usage[] =
    R"(fraggle - detecting k-mer bias in reads.

    Usage:
      fraggle count [options] (<fastq1> <fastq2>)...
      fraggle score [options] <counts>...

    Options:
      -h, --help                    Show this help message
      -k SIZE                       k-mer size [default: 6]
      --sample ID                   Label for the sample.
      --raw-counts FILE             Output the raw k-mer/position/strand count matrix.
      --raw-distributions FILE      Output the raw k-mer/position/strand fraction matrix.
      --save-counts FILE            Save the count data as a JSON object.
      --save-distributions FILE     Save the distribution data as a JSON object.
      -o FILE, --output-file FILE   File to write the final output to [default: -]
)";

std::string longest_common_prefix(const std::string &a, const std::string &b) {
  const size_t n = std::min(a.size(), b.size());
  std::string res;
  res.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (a[i] == b[i]) {
      res.push_back(a[i]);
    } else {
      break;
    }
  }
  return res;
}

std::string longest_common_prefix(const std::vector<std::string> &p_strings) {
  if (p_strings.size() == 1) {
    return p_strings.front();
  }
  if (p_strings.size() == 2) {
    return longest_common_prefix(p_strings[0], p_strings[1]);
  }
  std::vector<std::string> tmp;
  const size_t n = p_strings.size();
  size_t i = 0;
  while (i + 1 < n) {
    tmp.push_back(longest_common_prefix(p_strings[i], p_strings[i + 1]));
    i += 2;
  }
  if (i < n) {
    tmp.push_back(p_strings[i]);
  }
  return longest_common_prefix(tmp);
}

std::string chop_left(const std::string &p_orig, char p_ch) {
  size_t lhs = 0;
  size_t pos = 0;
  while (lhs < p_orig.size() &&
         (pos = p_orig.find(p_ch, lhs)) != std::string::npos) {
    lhs = pos + 1;
  }
  return std::string(p_orig.begin() + lhs, p_orig.end());
}

using gamma_distribution_ptr = std::shared_ptr<gamma_distribution<>>;

struct gamma_estimator {
  size_t N;
  double sx;
  double slx;
  double sxlx;

  gamma_estimator() {
    N = 0;
    sx = 0;
    slx = 0;
    sxlx = 0;
  }

  void add(double x) {
    if (x <= 0) {
      return;
    }
    double lx = std::log(x);
    N += 1;
    sx += x;
    slx += lx;
    sxlx += x * lx;
  }

  double kHat() const {
    double kHat0 = (N * sx) / (N * sxlx - slx * sx);
    return kHat0 - 1.0 / N *
                       (3 * kHat0 - (2.0 / 3.0) * (kHat0 / (1 + kHat0)) -
                        (4.0 / 5) * (kHat0 / ((1 + kHat0) * (1 + kHat0))));
  }

  double thetaHat() const {
    double thetaHat0 = 1.0 / (N * N) * (N * sxlx - slx * sx);
    return double(N) / double(N - 1) * thetaHat0;
  }

  gamma_distribution_ptr distribution() const {
    if (N <= 1) {
      return gamma_distribution_ptr();
    }
    return gamma_distribution_ptr(new gamma_distribution<>(kHat(), thetaHat()));
  }
};

double klDivergence(const std::vector<double> &p_P,
                    const std::vector<double> &p_Q) {
  double d = 0;
  for (size_t i = 0; i < p_P.size(); ++i) {
    const double p = p_P[i];
    const double q = p_Q[i];
    if (p == 0) {
      continue;
    }
    d += p * std::log(p / q);
  }
  return d;
}

void add_counts(counts &lhs, const counts &rhs) {
  for (size_t i = 0; i < rhs.size(); ++i) {
    lhs[i] += rhs[i];
  }
}

counts counts_add(const counts &lhs, const counts &rhs) {
  counts res(lhs.begin(), lhs.end());
  add_counts(res, rhs);
  return res;
}

void add_position_counts(position_counts &cts, const position_counts &oth) {
  for (size_t i = 0; i < oth.size(); ++i) {
    if (i < cts.size()) {
      add_counts(cts[i], oth[i]);
    } else {
      cts.push_back(oth[i]);
    }
  }
}

counts position_counts_sum(const position_counts &cts) {
  const size_t J = cts[0].size();
  counts res(J, 0);
  for (size_t i = 0; i < cts.size(); ++i) {
    add_counts(res, cts[i]);
  }
  return res;
}

distr compute_distr(const counts &cts) {
  double t = 0;
  const size_t J = cts.size();
  for (size_t j = 0; j < J; ++j) {
    t += cts[j];
  }
  distr res(J, 0);
  for (size_t j = 0; j < J; ++j) {
    res[j] = cts[j] / t;
  }
  return res;
}

position_distr compute_position_distr(const position_counts &cts) {
  position_distr res;
  for (size_t i = 0; i < cts.size(); ++i) {
    res.push_back(compute_distr(cts[i]));
  }
  return res;
}

struct fwd_and_rev_counts {
  position_counts fwd;
  position_counts rev;

  fwd_and_rev_counts &operator+=(const fwd_and_rev_counts &other) {
    add_position_counts(fwd, other.fwd);
    add_position_counts(rev, other.rev);
    return *this;
  }

  counts sum() const {
    return counts_add(position_counts_sum(fwd), position_counts_sum(rev));
  }

  nlohmann::json to_json() const {
    nlohmann::json res = nlohmann::json::object();
    res["fwd"] = fwd;
    res["rev"] = rev;
    return res;
  }

  static fwd_and_rev_counts from_json(const nlohmann::json &p_json) {
    fwd_and_rev_counts res;
    res.fwd = p_json["fwd"];
    res.rev = p_json["rev"];
    return res;
  }
};

struct counts_state {
  size_t K;
  std::string sample;
  fwd_and_rev_counts read1;
  fwd_and_rev_counts read2;

  counts_state &operator+=(const counts_state &other) {
    if (other.K != K) {
      throw std::runtime_error("mismatched values of K");
    }
    read1 += other.read1;
    read2 += other.read2;
    return *this;
  }

  nlohmann::json to_json() const {
    nlohmann::json res = nlohmann::json::object();
    res["K"] = K;
    res["sample"] = sample;
    res["read1"] = read1.to_json();
    res["read2"] = read2.to_json();
    return res;
  }

  static counts_state from_json(const nlohmann::json &p_json) {
    counts_state res;
    res.K = p_json["K"];
    res.sample = p_json["sample"];
    res.read1 = fwd_and_rev_counts::from_json(p_json["read1"]);
    res.read2 = fwd_and_rev_counts::from_json(p_json["read2"]);
    return res;
  }

  static counts_state load(const std::string &p_filename) {
    input_file_holder_ptr inp = files::in(p_filename);
    std::istream &in = **inp;
    nlohmann::json ctsJson = nlohmann::json::parse(in);
    return counts_state::from_json(ctsJson);
  }
};

struct fwd_and_rev_distr {
  position_distr fwd;
  position_distr rev;

  static fwd_and_rev_distr from_counts(const fwd_and_rev_counts &cts) {
    fwd_and_rev_distr res;
    res.fwd = compute_position_distr(cts.fwd);
    res.rev = compute_position_distr(cts.rev);
    return res;
  }

  nlohmann::json to_json() const {
    nlohmann::json res = nlohmann::json::object();
    res["fwd"] = fwd;
    res["rev"] = rev;
    return res;
  }

  static fwd_and_rev_distr from_json(const nlohmann::json &p_json) {
    fwd_and_rev_distr res;
    res.fwd = p_json["fwd"];
    res.rev = p_json["rev"];
    return res;
  }
};

struct distr_state {
  size_t K;
  std::string sample;
  distr global;
  fwd_and_rev_distr read1;
  fwd_and_rev_distr read2;

  static distr_state from_counts(const counts_state &cts) {
    distr_state res;
    res.K = cts.K;
    res.sample = cts.sample;
    res.global = compute_distr(counts_add(cts.read1.sum(), cts.read2.sum()));
    res.read1 = fwd_and_rev_distr::from_counts(cts.read1);
    res.read2 = fwd_and_rev_distr::from_counts(cts.read2);
    return res;
  }

  nlohmann::json to_json() const {
    nlohmann::json res = nlohmann::json::object();
    res["K"] = K;
    res["sample"] = sample;
    res["global"] = global;
    res["read1"] = read1.to_json();
    res["read2"] = read2.to_json();
    return res;
  }

  static distr_state from_json(const nlohmann::json &p_json) {
    distr_state res;
    res.K = p_json["K"];
    res.sample = p_json["sample"];
    res.global = p_json["global"].get<std::vector<double>>();
    res.read1 = fwd_and_rev_distr::from_json(p_json["read1"]);
    res.read2 = fwd_and_rev_distr::from_json(p_json["read2"]);
    return res;
  }
};


struct gamma_estimator_state {
  gamma_estimator global;
  std::vector<gamma_estimator> read1;
  std::vector<gamma_estimator> read2;

  void add_sample(const distr_state& p_global, const distr_state& p_sample) {
    for (size_t i = 0; i < p_sample.read1.fwd.size(); ++i) {
      double kldFwd = klDivergence(p_sample.read1.fwd[i], p_global.global);
      global.add(kldFwd);
      while (i >= read1.size()) {
        read1.push_back(gamma_estimator());
      }
      double kldFwdI = klDivergence(p_sample.read1.fwd[i], p_global.read1.fwd[i]);
      read1[i].add(kldFwdI);
    }
    for (size_t i = 0; i < p_sample.read2.fwd.size(); ++i) {
      double kldFwd = klDivergence(p_sample.read2.fwd[i], p_global.global);
      global.add(kldFwd);
      while (i >= read2.size()) {
        read2.push_back(gamma_estimator());
      }
      double kldFwdI = klDivergence(p_sample.read2.fwd[i], p_global.read2.fwd[i]);
      read2[i].add(kldFwdI);
    }
  }
};

struct gamma_state {
  gamma_distribution_ptr global;
  std::vector<gamma_distribution_ptr> read1;
  std::vector<gamma_distribution_ptr> read2;

  static gamma_state from_estimator(const gamma_estimator_state& p_estimator) {
    gamma_state res;
    res.global = p_estimator.global.distribution();
    for (size_t i = 0; i < p_estimator.read1.size(); ++i) {
      res.read1.push_back(p_estimator.read1[i].distribution());
    }
    for (size_t i = 0; i < p_estimator.read2.size(); ++i) {
      res.read2.push_back(p_estimator.read2[i].distribution());
    }
    return res;
  }

  template <typename X>
  typename std::enable_if<std::is_convertible<X,std::function<void(bool,size_t,double,double,double,double)>>::value,void>::type
  score_sample(const distr_state& p_global, const distr_state& p_sample,
                    X p_acceptor) const {
    for (size_t i = 0; i < p_sample.read1.fwd.size(); ++i) {
      double kldFwd = klDivergence(p_sample.read1.fwd[i], p_global.global);
      double pvalFwd = (kldFwd > 0 ? cdf(complement(*global, kldFwd)) : 1);
      double kldFwdI = klDivergence(p_sample.read1.fwd[i], p_global.read1.fwd[i]);
      double pvalFwdI = (kldFwdI > 0 ? cdf(complement(*read1[i], kldFwdI)) : 1);
      p_acceptor(true, i, kldFwd, pvalFwd, kldFwdI, pvalFwdI);
    }
    for (size_t i = 0; i < p_sample.read2.fwd.size(); ++i) {
      double kldFwd = klDivergence(p_sample.read2.fwd[i], p_global.global);
      double pvalFwd = (kldFwd > 0 ? cdf(complement(*global, kldFwd)) : 1);
      double kldFwdI = klDivergence(p_sample.read2.fwd[i], p_global.read2.fwd[i]);
      double pvalFwdI = (kldFwdI > 0 ? cdf(complement(*read2[i], kldFwdI)): 1);
      p_acceptor(false, i, kldFwd, pvalFwd, kldFwdI, pvalFwdI);
    }
  }
};

void score_position_distr(gamma_estimator &p_gam, const distr &p_global,
                          const distr_state &p_sample) {
  for (size_t i = 0; i < p_sample.read1.fwd.size(); ++i) {
    double kldFwd = klDivergence(p_sample.read1.fwd[i], p_global);
    p_gam.add(kldFwd);
  }
  for (size_t i = 0; i < p_sample.read2.fwd.size(); ++i) {
    double kldFwd = klDivergence(p_sample.read2.fwd[i], p_global);
    p_gam.add(kldFwd);
  }
}

void print_position_distr(const gamma_state &p_gamma,
                          const distr_state &p_global, const distr_state &p_sample,
                          std::ostream &p_out, bool p_printHeader) {

  if (p_printHeader) {
    p_out << "sample\tread\tpos\tkldGlobal\tpvalGlobal\tkldPosition\tpvalPosition" << std::endl;
  }

  p_gamma.score_sample(p_global, p_sample, [&](bool isFirst, size_t pos, double kldGlobal, double pvalGlobal, double kldPosition, double pvalPosition) {
    p_out << p_sample.sample << '\t' << (isFirst ? "R1" : "R2") << '\t' << pos << '\t' << kldGlobal << '\t'
    << pvalGlobal << '\t' << kldPosition << '\t' << pvalPosition << std::endl;
  });
}

void score_and_print(const distr_state &p_global,
                     const std::vector<std::string> &p_names,
                     const std::string &p_filename) {

  if (p_filename == "/dev/null") {
    return;
  }

  gamma_estimator_state estimator;

  for (size_t n = 0; n < p_names.size(); ++n) {
    BOOST_LOG_TRIVIAL(info) << "reloading " << p_names[n];
    counts_state cts = counts_state::load(p_names[n]);
    distr_state dst = distr_state::from_counts(cts);
    estimator.add_sample(p_global, dst);
  }

  gamma_state gamma = gamma_state::from_estimator(estimator);

  output_file_holder_ptr outp = files::out(p_filename);
  std::ostream &out = **outp;

  for (size_t n = 0; n < p_names.size(); ++n) {
    BOOST_LOG_TRIVIAL(info) << "scoring " << p_names[n];
    counts_state cts = counts_state::load(p_names[n]);
    distr_state dst = distr_state::from_counts(cts);
    print_position_distr(gamma, p_global, dst, out, n==0);
  }
}

void print_position_distr(const gamma_distribution<> &p_gamma,
                          const distr &p_global, const distr_state &p_sample,
                          std::ostream &p_out) {

  p_out << "sample\tread\tpos\tkld\tpval" << std::endl;

  for (size_t i = 0; i < p_sample.read1.fwd.size(); ++i) {
    double kldFwd = klDivergence(p_sample.read1.fwd[i], p_global);
    double pvalFwd = cdf(complement(p_gamma, kldFwd));
    p_out << p_sample.sample << '\t' << "R1" << '\t' << i << '\t' << kldFwd << '\t'
        << pvalFwd << std::endl;
  }
  for (size_t i = 0; i < p_sample.read2.fwd.size(); ++i) {
    double kldFwd = klDivergence(p_sample.read2.fwd[i], p_global);
    double pvalFwd = cdf(complement(p_gamma, kldFwd));
    p_out << p_sample.sample << '\t' << "R2" << '\t' << i << '\t' << kldFwd << '\t'
        << pvalFwd << std::endl;
  }
}

void score_and_print(const distr_state& p_sample, const std::string& p_filename) {
  if (p_filename == "/dev/null") {
    return;
  }

  BOOST_LOG_TRIVIAL(info) << "estimating gamma.";

  gamma_estimator gam;

  score_position_distr(gam, p_sample.global, p_sample);

  const double kHat = gam.kHat();
  const double thetaHat = gam.thetaHat();
  BOOST_LOG_TRIVIAL(info) << "gamma parameter k = " << kHat;
  BOOST_LOG_TRIVIAL(info) << "gamma parameter theta = " << thetaHat;

  gamma_distribution<> GammaDist(kHat, thetaHat);

  output_file_holder_ptr outp = files::out(p_filename);
  std::ostream &out = **outp;

  print_position_distr(GammaDist, p_sample.global, p_sample, out);
}

int main_merge(std::map<std::string, docopt::value> &opts) {
  const size_t K = opts["-k"].asLong();
  const size_t J = 1ULL << (2 * K);

  const std::vector<std::string> names = opts["<counts>"].asStringList();

  counts_state agg;
  agg.K = K;

  for (size_t i = 0; i < names.size(); ++i) {
    BOOST_LOG_TRIVIAL(info) << "loading " << names[i];
    counts_state cts = counts_state::load(names[i]);
    agg += cts;
  }
  distr_state global = distr_state::from_counts(agg);

  score_and_print(global, names, opts["--output-file"].asString());

  return 0;
}

int main_count(std::map<std::string, docopt::value> &opts) {
  const size_t K = opts["-k"].asLong();
  const size_t J = 1ULL << (2 * K);

  const std::vector<std::string> fns1 = opts["<fastq1>"].asStringList();
  const std::vector<std::string> fns2 = opts["<fastq2>"].asStringList();

  counts_state cts;
  cts.K = K;
  if (opts["--sample"]) {
    cts.sample = opts["--sample"].asString();
  } else {
    std::vector<std::string> labs(fns1.begin(), fns1.end());
    labs.insert(labs.end(), fns2.begin(), fns2.end());
    std::string sample = longest_common_prefix(labs);
    const std::set<char> trim{'_', '.', '-', '/'};
    while (sample.size() > 0 && trim.contains(sample.back())) {
      sample.pop_back();
    }
    sample = chop_left(sample, '/');
    cts.sample = sample;
  }
  BOOST_LOG_TRIVIAL(info) << "sample = " << cts.sample;

  size_t rn = 0;
  size_t r1Len = 0;
  size_t r2Len = 0;
  for (size_t i = 0; i < fns1.size(); ++i) {
    BOOST_LOG_TRIVIAL(info) << "scanning " << fns1[i] << " & " << fns2[i];
    fastq_reader::with(
        fns1[i], fns2[i],
        [&](const fastq_tuple &lhsRead, const fastq_tuple &rhsRead,
            bool &stop) {
          rn += 1;
          if ((rn & ((1ULL << 20) - 1)) == 0) {
            BOOST_LOG_TRIVIAL(info) << "reads processed: " << rn;
          }
          const fastq_text &lhsSeq = std::get<1>(lhsRead);
          size_t curR1Len = (lhsSeq.second - lhsSeq.first) - K + 1;
          if (curR1Len > r1Len) {
            while (cts.read1.fwd.size() < curR1Len) {
              cts.read1.fwd.push_back(counts(J, 0));
              cts.read1.rev.push_back(counts(J, 0));
            }
            r1Len = curR1Len;
          }
          kmers::make(lhsSeq, K, [&](size_t pos, kmer x, kmer xb) {
            cts.read1.fwd[pos][x] += 1;
            cts.read1.rev[pos][xb] += 1;
          });

          const fastq_text &rhsSeq = std::get<1>(rhsRead);
          size_t curR2Len = (rhsSeq.second - rhsSeq.first) - K + 1;
          if (curR2Len > r2Len) {
            while (cts.read2.fwd.size() < curR2Len) {
              cts.read2.fwd.push_back(counts(J, 0));
              cts.read2.rev.push_back(counts(J, 0));
            }
            r2Len = curR2Len;
          }
          kmers::make(rhsSeq, K, [&](size_t pos, kmer x, kmer xb) {
            cts.read2.fwd[pos][x] += 1;
            cts.read2.rev[pos][xb] += 1;
          });
        });
  }

  if (opts["--raw-counts"]) {
    output_file_holder_ptr outp = files::out(opts["--raw-counts"].asString());
    std::ostream &out = **outp;
    out << "kmer";
    for (size_t i = 0; i < cts.read1.fwd.size(); ++i) {
      out << "\tL" << i << 'f' << "\tL" << i << 'r';
    }
    for (size_t i = 0; i < cts.read2.fwd.size(); ++i) {
      out << "\tR" << i << 'f' << "\tR" << i << 'r';
    }
    out << std::endl;

    for (size_t j = 0; j < J; ++j) {
      out << kmers::render(K, j);
      for (size_t i = 0; i < cts.read1.fwd.size(); ++i) {
        out << '\t' << cts.read1.fwd[i][j];
        out << '\t' << cts.read1.rev[i][j];
      }
      for (size_t i = 0; i < cts.read2.fwd.size(); ++i) {
        out << '\t' << cts.read2.fwd[i][j];
        out << '\t' << cts.read2.rev[i][j];
      }
      out << std::endl;
    }
  }

  if (opts["--save-counts"]) {
    nlohmann::json res = cts.to_json();

    output_file_holder_ptr outp = files::out(opts["--save-counts"].asString());
    std::ostream &out = **outp;
    out << res << std::endl;
  }

  BOOST_LOG_TRIVIAL(info) << "computing distributions.";

  distr_state dist = distr_state::from_counts(cts);

  if (opts["--raw-distributions"]) {
    output_file_holder_ptr outp =
        files::out(opts["--raw-distributions"].asString());
    std::ostream &out = **outp;

    out << "kmer\tglobal";
    for (size_t i = 0; i < dist.read1.fwd.size(); ++i) {
      out << "\tL" << i << 'f' << "\tL" << i << 'r';
    }
    for (size_t i = 0; i < dist.read2.fwd.size(); ++i) {
      out << "\tR" << i << 'f' << "\tR" << i << 'r';
    }
    out << std::endl;

    for (size_t j = 0; j < J; ++j) {
      out << kmers::render(K, j) << '\t' << dist.global[j];
      for (size_t i = 0; i < dist.read1.fwd.size(); ++i) {
        out << '\t' << dist.read1.fwd[i][j];
        out << '\t' << dist.read1.rev[i][j];
      }
      for (size_t i = 0; i < dist.read2.fwd.size(); ++i) {
        out << '\t' << dist.read2.fwd[i][j];
        out << '\t' << dist.read2.rev[i][j];
      }
      out << std::endl;
    }
  }

  if (opts["--save-distributions"]) {
    nlohmann::json res = dist.to_json();

    output_file_holder_ptr outp =
        files::out(opts["--save-distributions"].asString());
    std::ostream &out = **outp;
    out << res << std::endl;
  }

  score_and_print(dist, opts["--output-file"].asString());

  return 0;
}

int main0(int argc, const char *argv[]) {
  boost::log::add_console_log(
      std::cerr, boost::log::keywords::format =
                     "[%TimeStamp%] [%ThreadID%] [%Severity%] %Message%");
  boost::log::add_common_attributes();

  std::map<std::string, docopt::value> opts =
      docopt::docopt(usage, {argv + 1, argv + argc}, true, "fraggle 0.1");

  for (auto itr = opts.begin(); itr != opts.end(); ++itr) {
    BOOST_LOG_TRIVIAL(info) << itr->first << '\t' << itr->second;
  }

  if (opts["score"].asBool()) {
    return main_merge(opts);
  }
  if (opts["count"].asBool()) {
    return main_count(opts);
  }

  return 1;
}

} // namespace
// namespace anonymous

int main(int argc, const char *argv[]) {
  try {
    return main0(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    ;
    return -1;
  }
}
