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
      fraggle [options] (<fastq1> <fastq2>)...
      fraggle --combine [options] <counts>...

    Options:
      -h, --help                Show this help message
      -k SIZE                   k-mer size [default: 6]
      --raw-counts FILE         Output the raw k-mer/position/strand count matrix.
      --raw-distributions FILE  Output the raw k-mer/position/strand fraction matrix.
      --save-counts FILE        Save the count data as a JSON object.
)";

struct gamma_estimator_state {
  size_t N;
  double sx;
  double slx;
  double sxlx;
};

struct gamma_estimator : gamma_estimator_state {

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
};

double klDivergence(const std::vector<double> &p_Q,
                    const std::vector<double> &p_P) {
  double d = 0;
  for (size_t i = 0; i < p_P.size(); ++i)
  {
    const double p = p_P[i];
    const double q = p_Q[i];
    if (p == 0) {
      continue;
    }
    d += p * std::log(p / q);
  }
  return d;
}

void add_counts(counts& lhs, const counts& rhs) {
  for (size_t i = 0; i < rhs.size(); ++i) {
    lhs[i] += rhs[i];
  }
}

counts counts_add(const counts& lhs, const counts& rhs) {
  counts res(lhs.begin(), lhs.end());
  add_counts(res, rhs);
  return res;
}

void add_position_counts(position_counts& cts, const position_counts& oth) {
  for (size_t i = 0; i < oth.size(); ++i) {
    if (i < cts.size()) {
      add_counts(cts[i], oth[i]);
    } else {
      cts.push_back(oth[i]);
    }
  }
}

counts position_counts_sum(const position_counts& cts) {
  const size_t J = cts[0].size();
  counts res(J, 0);
  for (size_t i = 0; i < cts.size(); ++i) {
    add_counts(res, cts[i]);
  }
  return res;
}

distr compute_distr(const counts& cts) {
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

position_distr compute_position_distr(const position_counts& cts) {
  position_distr res;
  for (size_t i = 0; i < cts.size(); ++i) {
    res.push_back(compute_distr(cts[i]));
  }
  return res;
}

struct fwd_and_rev_counts {
  position_counts fwd;
  position_counts rev;


  fwd_and_rev_counts& operator+=(const fwd_and_rev_counts& other) {
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

  static fwd_and_rev_counts from_json(const nlohmann::json& p_json) {
    fwd_and_rev_counts res;
    res.fwd = p_json["fwd"];
    res.rev = p_json["rev"];
    return res;
  }
};

struct counts_state {
  size_t K;
  fwd_and_rev_counts read1;
  fwd_and_rev_counts read2;

  counts_state& operator+=(const counts_state& other) {
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
    res["read1"] = read1.to_json();
    res["read2"] = read2.to_json();
    return res;
  }

  static counts_state from_json(const nlohmann::json& p_json) {
    counts_state res;
    res.K = p_json["K"];
    res.read1 = fwd_and_rev_counts::from_json(p_json["read1"]);
    res.read2 = fwd_and_rev_counts::from_json(p_json["read2"]);
    return res;
  }

  static counts_state load(const std::string& p_filename) {
    input_file_holder_ptr inp = files::in(p_filename);
    std::istream& in = **inp;
    nlohmann::json ctsJson = nlohmann::json::parse(in);
    return counts_state::from_json(ctsJson);
  }
};

struct fwd_and_rev_distr {
  position_distr fwd;
  position_distr rev;

  static fwd_and_rev_distr from_counts(const fwd_and_rev_counts& cts) {
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

  static fwd_and_rev_counts from_json(const nlohmann::json& p_json) {
    fwd_and_rev_counts res;
    res.fwd = p_json["fwd"];
    res.rev = p_json["rev"];
    return res;
  }
};

struct distr_state {
  distr global;
  fwd_and_rev_distr read1;
  fwd_and_rev_distr read2;

  static distr_state from_counts(const counts_state& cts) {
    distr_state res;
    res.global = compute_distr(counts_add(cts.read1.sum(), cts.read2.sum()));
    res.read1 = fwd_and_rev_distr::from_counts(cts.read1);
    res.read2 = fwd_and_rev_distr::from_counts(cts.read2);
    return res;
  }
};


int main_merge(std::map<std::string, docopt::value>& opts)
{
  const size_t K = opts["-k"].asLong();
  const size_t J = 1ULL << (2*K);

  const std::vector<std::string> names = opts["<counts>"].asStringList();

  counts_state agg;
  agg.K = K;

  for (size_t i = 0; i < names.size(); ++i) {
    BOOST_LOG_TRIVIAL(info) << "loading " << names[i];
    counts_state cts = counts_state::load(names[i]);
    agg += cts;
  }
  distr_state global = distr_state::from_counts(agg);

  gamma_estimator gam;

  for (size_t n = 0; n < names.size(); ++n) {
    BOOST_LOG_TRIVIAL(info) << "reloading " << names[n];
    counts_state cts = counts_state::load(names[n]);
    distr_state dst = distr_state::from_counts(cts);
    for (size_t i = 0; i < dst.read1.fwd.size(); ++i) {
      double kldFwd = klDivergence(dst.read1.fwd[i], global.global);
      gam.add(kldFwd);
      double kldRev = klDivergence(dst.read1.rev[i], global.global);
      gam.add(kldRev);
    }
    for (size_t i = 0; i < dst.read2.fwd.size(); ++i) {
      double kldFwd = klDivergence(dst.read2.fwd[i], global.global);
      gam.add(kldFwd);
      double kldRev = klDivergence(dst.read2.rev[i], global.global);
      gam.add(kldRev);
    }
  }

  const double kHat = gam.kHat();
  const double thetaHat = gam.thetaHat();
  gamma_distribution<> GammaDist(kHat, thetaHat);

  for (size_t n = 0; n < names.size(); ++n) {
    BOOST_LOG_TRIVIAL(info) << "scoring " << names[n];
    counts_state cts = counts_state::load(names[n]);
    distr_state dst = distr_state::from_counts(cts);
    for (size_t i = 0; i < dst.read1.fwd.size(); ++i) {
      double kldFwd = klDivergence(dst.read1.fwd[i], global.global);
      double pvalFwd = cdf(complement(GammaDist, kldFwd));
      double kldRev = klDivergence(dst.read1.rev[i], global.global);
      double pvalRev = cdf(complement(GammaDist, kldRev));
      std::cout << names[n]
        << '\t' << "R1"
        << '\t' << i
        << '\t' << kldFwd
        << '\t' << pvalFwd
        << '\t' << kldRev
        << '\t' << pvalRev
        << std::endl;
    }
    for (size_t i = 0; i < dst.read2.fwd.size(); ++i) {
      double kldFwd = klDivergence(dst.read2.fwd[i], global.global);
      double pvalFwd = cdf(complement(GammaDist, kldFwd));
      double kldRev = klDivergence(dst.read2.rev[i], global.global);
      double pvalRev = cdf(complement(GammaDist, kldRev));
      std::cout << names[n]
        << '\t' << "R2"
        << '\t' << i
        << '\t' << kldFwd
        << '\t' << pvalFwd
        << '\t' << kldRev
        << '\t' << pvalRev
        << std::endl;
    }
  }

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

  if (opts["--combine"].asBool()) {
    return main_merge(opts);
  }

  const size_t K = opts["-k"].asLong();
  const size_t J = 1ULL << (2*K);

  const std::vector<std::string> fns1 = opts["<fastq1>"].asStringList();
  const std::vector<std::string> fns2 = opts["<fastq2>"].asStringList();

  counts_state cts;
  cts.K = K;

  size_t rn = 0;
  for (size_t i = 0; i < fns1.size(); ++i) {
    BOOST_LOG_TRIVIAL(info) << "scanning " << fns1[i] << " & " << fns2[i];
    fastq_reader::with(fns1[i], fns2[i], [&](const fastq_tuple& lhsRead, const fastq_tuple& rhsRead, bool& stop) {
      rn += 1;
      if ((rn & ((1ULL << 20) - 1)) == 0) {
        BOOST_LOG_TRIVIAL(info) << "reads processed: " << rn;
      }
      const fastq_text& lhsSeq = std::get<1>(lhsRead);
      kmers::make(lhsSeq, K, [&](size_t pos, kmer x, kmer xb) {
        while (cts.read1.fwd.size() <= pos) {
          cts.read1.fwd.push_back(counts(J, 0));
        }
        cts.read1.fwd[pos][x] += 1;
        while (cts.read1.rev.size() <= pos) {
          cts.read1.rev.push_back(counts(J, 0));
        }
        cts.read1.rev[pos][xb] += 1;
      });
      const fastq_text& rhsSeq = std::get<1>(rhsRead);
      kmers::make(rhsSeq, K, [&](size_t pos, kmer x, kmer xb) {
        while (cts.read2.fwd.size() <= pos) {
          cts.read2.fwd.push_back(counts(J, 0));
        }
        cts.read2.fwd[pos][x] += 1;
        while (cts.read2.rev.size() <= pos) {
          cts.read2.rev.push_back(counts(J, 0));
        }
        cts.read2.rev[pos][xb] += 1;
      });
    });
  }

  if (opts["--raw-counts"]) {
    output_file_holder_ptr outp = files::out(opts["--raw-counts"].asString());
    std::ostream& out = **outp;
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
    std::ostream& out = **outp;
    out << res << std::endl;
  }

  
  BOOST_LOG_TRIVIAL(info) << "computing distributions.";

  distr_state dist = distr_state::from_counts(cts);
  
  if (opts["--raw-distributions"]) {
    output_file_holder_ptr outp = files::out(opts["--raw-distributions"].asString());
    std::ostream& out = **outp;

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

  BOOST_LOG_TRIVIAL(info) << "estimating gamma.";

  gamma_estimator gam;
  for (size_t i = 0; i < dist.read1.fwd.size(); ++i) {
    double kldFwd = klDivergence(dist.read1.fwd[i], dist.global);
    gam.add(kldFwd);
    double kldRev = klDivergence(dist.read1.rev[i], dist.global);
    gam.add(kldRev);
  }
  for (size_t i = 0; i < dist.read2.fwd.size(); ++i) {
    double kldFwd = klDivergence(dist.read2.fwd[i], dist.global);
    gam.add(kldFwd);
    double kldRev = klDivergence(dist.read2.rev[i], dist.global);
    gam.add(kldRev);
  }

  const double kHat = gam.kHat();
  const double thetaHat = gam.thetaHat();

  gamma_distribution<> GammaDist(kHat, thetaHat);

  for (size_t i = 0; i < dist.read1.fwd.size(); ++i) {
    double kldFwd = klDivergence(dist.read1.fwd[i], dist.global);
    double pvalFwd = cdf(complement(GammaDist, kldFwd));
    double kldRev = klDivergence(dist.read1.rev[i], dist.global);
    double pvalRev = cdf(complement(GammaDist, kldRev));
    std::cout
      << '\t' << "R1"
      << '\t' << i
      << '\t' << kldFwd
      << '\t' << pvalFwd
      << '\t' << kldRev
      << '\t' << pvalRev
      << std::endl;
  }
  for (size_t i = 0; i < dist.read2.fwd.size(); ++i) {
    double kldFwd = klDivergence(dist.read2.fwd[i], dist.global);
    double pvalFwd = cdf(complement(GammaDist, kldFwd));
    double kldRev = klDivergence(dist.read2.rev[i], dist.global);
    double pvalRev = cdf(complement(GammaDist, kldRev));
    std::cout
      << '\t' << "R2"
      << '\t' << i
      << '\t' << kldFwd
      << '\t' << pvalFwd
      << '\t' << kldRev
      << '\t' << pvalRev
      << std::endl;
  }

  return 0;
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
