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
#include "scindo/fasta.hpp"
#include "scindo/files.hpp"
#include "scindo/kmers.hpp"
#include "scindo/okeefe.hpp"
#include "scindo/tsv.hpp"
#include "scindo/vcf.hpp"

using namespace scindo;

namespace // anonymous
{
const char usage[] =
    R"(crispin - detecting differences between BAMs.

    Usage:
      crispin [options] <reference> <BAM>...

    Options:
      -h, --help                        Show this help message
      --index-only                      Exit after indexing step.
      -S INDEX, --save-index INDEX      Save the index.
      -X INDEX, --index INDEX           Load a previously computed index.
      -c COV, --min-coverage COV        Minimum coverage to consider [default: 10]
      -R BED, --regions BED             Report against only the regions in the given BED file.
)";

struct gamma_estimator {
  size_t N;
  double sx;
  double slx;
  double sxlx;

  gamma_estimator() : N(0), sx(0), slx(0), sxlx(0) {}

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

struct summarizer {
  size_t N;
  double sx;
  double sx2;

  summarizer() : N(0), sx(0), sx2(0) {}

  void add(double x) {
    N += 1;
    sx += x;
    sx2 += x * x;
  }

  double mean() const { return sx / N; }

  double var() const {
    double m = mean();
    return sx2 / N - m * m;
  }

  double sd() const { return std::sqrt(var()); }
};

double sqr(double x) { return x * x; }

double fit_homozygous(const std::vector<double> &p_probs) {
  if (p_probs.size() < 2) {
    return 0;
  }
  size_t n = p_probs.size();
  const double eps = 0.05;
  const double p0 = 1 - eps;
  const double p1 = eps / n;
  double chi2 = sqr(p_probs[0] - p0) / p0;
  for (size_t i = 1; i < n; ++i) {
    chi2 += sqr(p_probs[i] - p1) / p1;
  }
  boost::math::chi_squared D(n - 1);
  return boost::math::cdf(D, chi2);
}

double fit_heterozygous(const std::vector<double> &p_probs) {
  if (p_probs.size() < 2) {
    return 1;
  }
  size_t n = p_probs.size();
  const double eps = 0.05;
  const double p0 = 0.5 - eps;
  const double p1 = eps / n;
  double chi2 = 0;
  chi2 += sqr(p_probs[0] - p0) / p0;
  chi2 += sqr(p_probs[1] - p0) / p0;
  for (size_t i = 2; i < n; ++i) {
    chi2 += sqr(p_probs[i] - p1) / p1;
  }
  boost::math::chi_squared D(n - 1);
  return boost::math::cdf(D, chi2);
}

int chromOrd(const std::string &p_name) {
  static std::map<std::string, int> chroms = {
      {"chr1", 1},   {"chr2", 2},   {"chr3", 3},   {"chr4", 4},   {"chr5", 5},
      {"chr6", 6},   {"chr7", 7},   {"chr8", 8},   {"chr9", 9},   {"chr10", 10},
      {"chr11", 11}, {"chr12", 12}, {"chr13", 13}, {"chr14", 14}, {"chr15", 15},
      {"chr16", 16}, {"chr17", 17}, {"chr18", 18}, {"chr19", 19}, {"chr20", 20},
      {"chr21", 21}, {"chr22", 22}, {"chrX", 23},  {"chrY", 24},  {"chrM", 25}};
  auto itr = chroms.find(p_name);
  if (itr != chroms.end()) {
    return itr->second;
  }
  return -1;
}

bool chromLess(const std::string &p_lhs, const std::string &p_rhs) {
  int lhsOrd = chromOrd(p_lhs);
  int rhsOrd = chromOrd(p_rhs);
  if (lhsOrd > 0 && rhsOrd > 0) {
    return lhsOrd < rhsOrd;
  }
  if (lhsOrd > 0) {
    return true;
  }
  if (rhsOrd > 0) {
    return false;
  }
  return p_lhs < p_rhs;
}

using vcf_ptr = std::shared_ptr<vcf::vcf_file_reader>;
using vcf_ptrs = std::vector<vcf_ptr>;
using vcf_itr = vcf_ptrs::const_iterator;

bool cmpVcf(const vcf_ptr &a, const vcf_ptr &b) {

  if (a->chromId() != b->chromId()) {
    return a->chromId() < b->chromId();
  }
  return a->pos0() < b->pos0();
}

void mergeVcfs(const vcf_ptrs &p_vcfs,
               std::function<void(vcf_itr, vcf_itr)> p_consumer) {

  vcf_ptrs vcfs;
  for (auto itr = p_vcfs.begin(); itr != p_vcfs.end(); ++itr) {
    if ((*itr)->next()) {
      vcfs.push_back(*itr);
    }
  }

  while (vcfs.size() > 0) {
    std::sort(vcfs.begin(), vcfs.end(), cmpVcf);

    for (auto itr = vcfs.begin(); itr != vcfs.end(); ++itr) {
      BOOST_LOG_TRIVIAL(info)
          << "pre  " << nlohmann::json((*itr)->samples()) << " : "
          << (*itr)->chrom() << "\t" << (*itr)->pos1();
    }

    size_t i = 1;
    while (i < vcfs.size()) {
      if (vcfs[i]->chromId() != vcfs[0]->chromId() ||
          vcfs[i]->pos0() != vcfs[0]->pos0()) {
        break;
      }
      i += 1;
    }

    p_consumer(vcfs.begin(), vcfs.begin() + i);

    vcf_ptrs tmp;
    for (size_t j = 0; j < i; ++j) {
      if (vcfs[j]->next()) {
        tmp.push_back(vcfs[j]);
      }
    }
    tmp.insert(tmp.end(), vcfs.begin() + i, vcfs.end());
    vcfs.swap(tmp);

    for (auto itr = vcfs.begin(); itr != vcfs.end(); ++itr) {
      BOOST_LOG_TRIVIAL(info)
          << "post " << nlohmann::json((*itr)->samples()) << " : "
          << (*itr)->chrom() << "\t" << (*itr)->pos1();
    }
  }
}

using range_vector = std::vector<std::pair<uint32_t, uint32_t>>;

void index_bed_file(const std::string &p_filename,
                    std::map<std::string, range_vector> &p_regions) {

  // Read it all in.
  //
  scindo::input_file_holder_ptr inp = scindo::files::in(p_filename);
  scindo::tsv::with(**inp, [&](const std::vector<std::string> &p_row) {
    const std::string &chrom = p_row[0];
    const size_t start = boost::lexical_cast<uint64_t>(p_row[1]);
    const size_t stop = boost::lexical_cast<uint64_t>(p_row[2]);
    p_regions[chrom].push_back(std::make_pair(start, stop));
  });

  // Sort, and merge overlapping regions.
  //
  for (auto itr = p_regions.begin(); itr != p_regions.end(); ++itr) {
    range_vector &rs = itr->second;
    std::sort(rs.begin(), rs.end());
    rs.erase(std::unique(rs.begin(), rs.end()), rs.end());
    range_vector tmp;
    for (auto jtr = rs.begin(); jtr != rs.end(); ++jtr) {
      if (tmp.size() == 0 || jtr->first > tmp.back().second) {
        tmp.push_back(*jtr);
      } else {
        tmp.back().second = std::max(tmp.back().second, jtr->second);
      }
    }
    rs.swap(tmp);
  }
}

struct ranges_cursor {
  const std::map<std::string, range_vector> &m_ranges;
  std::string m_chrom;
  range_vector::const_iterator m_cur;
  range_vector::const_iterator m_end;

  ranges_cursor(const std::map<std::string, range_vector> &p_ranges)
      : m_ranges(p_ranges) {
    m_cur = m_end;
  }

  bool operator()(const std::string &p_chrom, int64_t p_pos) {
    if (p_chrom != m_chrom) {
      m_chrom = p_chrom;
      auto itr = m_ranges.find(m_chrom);
      if (itr == m_ranges.end()) {
        m_cur = m_end;
        return false;
      }
      m_cur = itr->second.begin();
      m_end = itr->second.end();
    }
    while (m_cur != m_end && p_pos >= m_cur->second) {
      ++m_cur;
    }
    if (m_cur == m_end || p_pos < m_cur->first) {
      return false;
    }
    return true;
  }
};

using bam_ptr = std::shared_ptr<bam::bam_file_reader>;
using bam_ptrs = std::vector<bam_ptr>;
using bam_itr = bam_ptrs::const_iterator;

struct coverage {
  size_t bases[5];
  std::map<std::string, size_t> ins;
  std::map<size_t, size_t> del;

  coverage() {
    bases[0] = 0;
    bases[1] = 0;
    bases[2] = 0;
    bases[3] = 0;
    bases[4] = 0;
  }

  size_t total() const {
    size_t t = 0;
    t += bases[0];
    t += bases[1];
    t += bases[2];
    t += bases[3];
    t += bases[4];
    for (auto itr = ins.begin(); itr != ins.end(); ++itr) {
      t += itr->second;
    }
    for (auto itr = del.begin(); itr != del.end(); ++itr) {
      t += itr->second;
    }
    return t;
  }

  coverage &operator+=(const coverage &p_rhs) {
    bases[0] += p_rhs.bases[0];
    bases[1] += p_rhs.bases[1];
    bases[2] += p_rhs.bases[2];
    bases[3] += p_rhs.bases[3];
    bases[4] += p_rhs.bases[4];
    for (auto itr = p_rhs.ins.begin(); itr != p_rhs.ins.end(); ++itr) {
      ins[itr->first] += itr->second;
    }
    for (auto itr = p_rhs.del.begin(); itr != p_rhs.del.end(); ++itr) {
      del[itr->first] += itr->second;
    }
    return *this;
  }

  void add(const char &p_ch) {
    switch (p_ch) {
    case 'a':
    case 'A': {
      bases[0] += 1;
      return;
    }
    case 'c':
    case 'C': {
      bases[1] += 1;
      return;
    }
    case 'g':
    case 'G': {
      bases[2] += 1;
      return;
    }
    case 't':
    case 'T':
    case 'u':
    case 'U': {
      bases[3] += 1;
      return;
    }
    default: {
      bases[4] += 1;
    }
    }
  }

  void addIns(const std::string &p_seq) { ins[p_seq] += 1; }

  void addDel(const size_t &p_len) { del[p_len] += 1; }

  std::map<std::string, double> distribution() {
    double t = total();
    std::map<std::string, double> res;
    if (bases[0] > 0)
      res["A"] = double(bases[0]) / t;
    if (bases[1] > 0)
      res["C"] = double(bases[1]) / t;
    if (bases[2] > 0)
      res["G"] = double(bases[2]) / t;
    if (bases[3] > 0)
      res["T"] = double(bases[3]) / t;
    if (bases[4] > 0)
      res["N"] = double(bases[4]) / t;
    for (auto itr = ins.begin(); itr != ins.end(); ++itr) {
      std::string lab = std::string("+") + itr->first;
      res[lab] = double(itr->second) / t;
    }
    for (auto itr = del.begin(); itr != del.end(); ++itr) {
      std::string lab = std::string("-") + std::to_string(itr->first);
      res[lab] = double(itr->second) / t;
    }
    return res;
  }

  std::string indels() const {
    if (ins.size() == 0 && del.size() == 0) {
      return ".";
    }
    std::ostringstream out;
    auto itr = ins.begin();
    while (itr != ins.end()) {
      out << "+" << itr->first << ":" << itr->second;
      ++itr;
      if (itr != ins.end() || del.size() > 0) {
        out << ",";
      }
    }
    auto jtr = del.begin();
    while (jtr != del.end()) {
      out << "-" << jtr->first << ":" << jtr->second;
      ++jtr;
      if (jtr != del.end()) {
        out << ",";
      }
    }
    return out.str();
  }

  size_t pseudo_allele_count() const {
    size_t n = 0;
    for (size_t i = 0; i < 5; ++i) {
      if (bases[i] > 0) {
        n += 1;
      }
    }
    n += ins.size();
    n += del.size();
    return n;
  }
};

struct pileup {
  // using queue_type = std::deque<coverage>;
  using queue_type = okeefe<coverage>;
  bam::bam_file_reader &bam;
  bool hasNext;
  std::string chrom;
  int64_t pos;
  queue_type cov;

  pileup(bam::bam_file_reader &p_bam) : bam(p_bam), pos(-1) {
    hasNext = bam.next();
    while (hasNext && skipThisRead()) {
      hasNext = bam.next();
    }
  }

  bool next(std::string &p_chrom, int64_t &p_pos, coverage &p_cov) {
    if (cov.size() == 0) {
      if (hasNext) {
        addRead();
        hasNext = bam.next();
        while (hasNext && skipThisRead()) {
          hasNext = bam.next();
        }
      } else {
        return false;
      }
    }
    while (hasNext && bam.chrom() == chrom && bam.pos0() == pos) {
      addRead();
      hasNext = bam.next();
      while (hasNext && skipThisRead()) {
        hasNext = bam.next();
      }
    }
    if (chrom != p_chrom) {
      p_chrom = chrom;
    }
    p_pos = pos;
    p_cov = cov.front();
    pos += 1;
    cov.pop_front();
    return true;
  }

  bool skipThisRead() const {
    return hasNext && bam::flag(bam.flag()).is<bam::flag::unmapped>() &&
           bam.mapq() < 55;
  }

  bool addRead() {
    if (bam.chrom() != chrom) {
      if (cov.size() > 0) {
        return false;
      } else {
        chrom = bam.chrom();
        pos = bam.pos0();
      }
    }

    if (cov.size() > 0 && pos + cov.size() < bam.pos0()) {
      return false;
    }

    if (cov.size() == 0) {
      pos = bam.pos0();
    }

    const auto cig = bam.cigar();
    const auto seq = bam.seq();
    int64_t p = bam.pos0() - pos; // ref relative to pos.
    int64_t q = 0;                // read

    size_t z = cov.size();
    for (size_t i = 0; i < cig.size(); ++i) {
      const size_t n = cig.len(i);
      switch (cig.op(i)) {
      case 'M': {
        if (n > 200) {
          BOOST_LOG_TRIVIAL(info) << "big match " << (p + n);
        }
        cov.reserve(p + n);
        while (cov.size() <= p + n) {
          cov.push_back(coverage());
        }
        for (size_t j = 0; j < n; ++j) {
          cov[p + j].add(seq[q + j]);
        }
        p += n;
        q += n;
        break;
      }
      case 'I': {
        if (n > 200) {
          BOOST_LOG_TRIVIAL(info) << "big ins " << n;
        }
        cov.reserve(p);
        while (cov.size() <= p) {
          cov.push_back(coverage());
        }
        const auto ins = seq.range(q, q + n).asString();
        cov[p].addIns(ins);
        q += n;
        break;
      }
      case 'D': {
        if (n > 200) {
          BOOST_LOG_TRIVIAL(info) << "big del " << n;
        }
        cov.reserve(p);
        while (cov.size() <= p) {
          cov.push_back(coverage());
        }
        cov[p].addDel(n);
        p += n;
        break;
      }
      case 'N': {
        p += n;
        break;
      }
      case 'S': {
        q += n;
        break;
      }
      case 'H': {
        break;
      }
      case 'P': {
        break;
      }
      }
    }
    if (z < 1000 && cov.size() > 1000) {
      std::ostringstream msg;
      for (size_t i = 0; i < cig.size(); ++i) {
        msg << cig.op(i);
        msg << cig.len(i);
      }
      BOOST_LOG_TRIVIAL(info)
          << "z = " << z << "\tand cov.size() = " << cov.size() << "\t" << pos
          << "\t" << bam.pos0() << "\t" << p << "\t" << msg.str();
    }
    return true;
  }
};
using pileup_ptr = std::shared_ptr<pileup>;
using pileup_ptrs = std::vector<pileup_ptr>;

struct pileup_cursor {
  pileup &pile;
  std::string name;
  bool valid;
  std::string chrom;
  int64_t pos;
  coverage cov;

  pileup_cursor(pileup &p_pile) : pile(p_pile), pos(0) {
    name = pile.bam.find_header("RG", "SM");
    if (name.size() == 0) {
      name = pile.bam.filename;
    }
    valid = pile.next(chrom, pos, cov);
  }

  void next() { valid = pile.next(chrom, pos, cov); }
};
using pileup_cursor_ptr = std::shared_ptr<pileup_cursor>;
using pileup_cursor_ptrs = std::vector<pileup_cursor_ptr>;

bool cursorLess(const pileup_cursor_ptr &p_lhs,
                const pileup_cursor_ptr &p_rhs) {
  if (!p_lhs->valid || !p_rhs->valid) {
    return p_lhs->valid;
  }
  if (p_lhs->chrom != p_rhs->chrom) {
    return chromLess(p_lhs->chrom, p_rhs->chrom);
  }
  return p_lhs->pos < p_rhs->pos;
}

double klDivergence(const std::map<std::string, double> &p_Q,
                    const std::map<std::string, double> &p_P) {
  double d = 0;
  for (auto itr = p_P.begin(); itr != p_P.end(); ++itr) {
    if (itr->second == 0) {
      continue;
    }
    d += itr->second * std::log(itr->second / p_Q.find(itr->first)->second);
  }
  if (d < 0) {
    BOOST_LOG_TRIVIAL(info)
        << d << '\t' << nlohmann::json(p_Q) << '\t' << nlohmann::json(p_P);
  }
  return d;
}

struct fasta_cursor {
  const std::string filename;
  scindo::input_file_holder_ptr inp;
  scindo::seq::fasta_reader src;
  std::string name;

  fasta_cursor(const std::string &p_filename)
      : filename(p_filename), inp(scindo::files::in(filename)), src(**inp) {}

  bool seek(const std::string &p_name) {
    if (p_name == name) {
      return true;
    }

    std::vector<std::string> parts;
    while (true) {
      name_parts((*src).first, parts);
      for (auto itr = parts.begin(); itr != parts.end(); ++itr) {
        if (*itr == p_name) {
          name = *itr;
          return true;
        }
      }
      if (!src.more()) {
        break;
      }
      ++src;
    }
    return false;
  }

  const std::string &seq() const { return (*src).second; }

  scindo::kmer context(size_t p_pos) const {
    static constexpr size_t K = 3;
    scindo::kmer res = 0;
    if (p_pos < K) {
      const auto itr = seq().begin();
      scindo::kmers::make(std::make_pair(itr, itr + 2 * K), 2 * K,
                          [&](scindo::kmer x) { res = x; });
      return res;
    } else {
      const auto itr = seq().begin() + p_pos - K;
      scindo::kmers::make(std::make_pair(itr, itr + 2 * K), 2 * K,
                          [&](scindo::kmer x) { res = x; });
      return res;
    }
  }

  static void name_parts(const std::string &p_str,
                         std::vector<std::string> &p_parts) {
    p_parts.clear();
    size_t n = 0;
    while (n < p_str.size()) {
      size_t i = p_str.find(' ', n);
      if (i == std::string::npos) {
        p_parts.push_back(p_str.substr(n));
        if (p_parts.back().size() == 0) {
          p_parts.pop_back();
        }
        return;
      }
      p_parts.push_back(p_str.substr(n, i - n));
      if (p_parts.back().size() == 0) {
        p_parts.pop_back();
      }
      n = i + 1;
    }
  }
};

struct gamma {
  size_t N;
  double k;
  double theta;

  gamma() : N(0), k(0), theta(0) {}

  gamma(size_t p_N, double p_k, double p_theta)
      : N(p_N), k(p_k), theta(p_theta) {}

  gamma(const nlohmann::json &p_json)
      : N(p_json["N"]), k(p_json["k"]), theta(p_json["theta"]) {}

  nlohmann::json json() const {
    nlohmann::json res;
    res["N"] = N;
    res["k"] = k;
    res["theta"] = theta;
    return res;
  }
};

void build_estimate(const std::string &p_reference,
                    const std::vector<std::string> &p_filenames,
                    const size_t &p_C,
                    const std::map<std::string, range_vector> &p_ranges,
                    std::map<scindo::kmer, gamma> &p_estimators) {
  bam_ptrs bams;
  for (auto bamName : p_filenames) {
    bams.push_back(bam_ptr(new bam::bam_file_reader(bamName)));
  }
  pileup_ptrs piles;
  for (auto itr = bams.begin(); itr != bams.end(); ++itr) {
    piles.push_back(pileup_ptr(new pileup(**itr)));
  }
  pileup_cursor_ptrs cursors;
  for (auto itr = piles.begin(); itr != piles.end(); ++itr) {
    cursors.push_back(pileup_cursor_ptr(new pileup_cursor(**itr)));
    if (!cursors.back()->valid) {
      cursors.pop_back();
    }
  }
  std::sort(cursors.begin(), cursors.end(), cursorLess);

  fasta_cursor ref(p_reference);

  std::map<scindo::kmer, gamma_estimator> estimators;

  size_t N = 0;
  std::string chrom;
  std::vector<coverage> covs;

  const bool useRanges = (p_ranges.size() > 0);
  ranges_cursor inRange(p_ranges);

  while (cursors.size()) {
    chrom = cursors.front()->chrom;
    int64_t pos = cursors.front()->pos;

    bool seekRes = ref.seek(chrom);

    scindo::kmer ctxt = ref.context(pos);

    coverage cov;
    covs.clear();
    for (size_t i = 0; i < cursors.size() && cursors[i]->chrom == chrom &&
                       cursors[i]->pos == pos;
         ++i) {
      cov += cursors[i]->cov;
      covs.push_back(cursors[i]->cov);
      cursors[i]->next();
    }
    std::sort(cursors.begin(), cursors.end(), cursorLess);
    while (cursors.size() > 0 && !cursors.back()->valid) {
      cursors.pop_back();
    }

    if (useRanges && !inRange(chrom, pos)) {
      continue;
    }

    const size_t t = cov.total();
    if (t <= p_C) {
      continue;
    }

    if (cov.pseudo_allele_count() <= 1) {
      continue;
    }

    auto dist = cov.distribution();

    gamma_estimator &E = estimators[ctxt];
    for (auto itr = covs.begin(); itr != covs.end(); ++itr) {
      if (itr->total() < p_C) {
        continue;
      }
      auto dist1 = itr->distribution();
      auto kld = klDivergence(dist, dist1);

      E.add(kld);
    }

    N += 1;
    if ((N & 0xfffff) == 0) {
      BOOST_LOG_TRIVIAL(info)
          << "updated " << N << " sites (chrom = " << chrom << ").";
    }
  }
  for (auto itr = estimators.begin(); itr != estimators.end(); ++itr) {
    const gamma_estimator &e = itr->second;
    p_estimators[itr->first] = gamma(e.N, e.kHat(), e.thetaHat());
  }
}

void compute_significance(const std::string &p_reference,
                          const std::vector<std::string> &p_filenames,
                          const size_t &p_C,
                          const std::map<std::string, range_vector> &p_ranges,
                          const std::map<scindo::kmer, gamma> &p_estimators,
                          std::ostream &p_out) {
  bam_ptrs bams;
  for (auto bamName : p_filenames) {
    bams.push_back(bam_ptr(new bam::bam_file_reader(bamName)));
  }
  pileup_ptrs piles;
  for (auto itr = bams.begin(); itr != bams.end(); ++itr) {
    piles.push_back(pileup_ptr(new pileup(**itr)));
  }
  pileup_cursor_ptrs cursors;
  for (auto itr = piles.begin(); itr != piles.end(); ++itr) {
    cursors.push_back(pileup_cursor_ptr(new pileup_cursor(**itr)));
    if (!cursors.back()->valid) {
      cursors.pop_back();
    }
  }
  std::sort(cursors.begin(), cursors.end(), cursorLess);

  using namespace boost::math;

  fasta_cursor ref(p_reference);

  using item = std::pair<std::string, coverage>;

  std::string chrom;
  std::vector<item> covs;
  std::vector<double> ps;

  const bool useRanges = (p_ranges.size() > 0);
  ranges_cursor inRange(p_ranges);

  p_out << "locus" << '\t' << "num.samples" << '\t' << "context" << '\t'
        << "sample" << '\t' << "divergence" << '\t' << "pvalue" << '\t'
        << "pseudo.allele.count" << '\t' << "gt" << '\t' << "fit" << '\t' << "A" << '\t' << "C"
        << '\t' << "G" << '\t' << "T" << '\t' << "N" << '\t' << "indels"
        << std::endl;

  while (cursors.size()) {
    chrom = cursors.front()->chrom;
    int64_t pos = cursors.front()->pos;

    bool seekRes = ref.seek(chrom);

    scindo::kmer ctxt = ref.context(pos);

    coverage cov;
    covs.clear();
    for (size_t i = 0; i < cursors.size() && cursors[i]->chrom == chrom &&
                       cursors[i]->pos == pos;
         ++i) {
      cov += cursors[i]->cov;
      if (cursors[i]->cov.total() >= p_C) {
        covs.push_back({cursors[i]->name, cursors[i]->cov});
      }
      cursors[i]->next();
    }
    std::sort(cursors.begin(), cursors.end(), cursorLess);
    while (cursors.size() > 0 && !cursors.back()->valid) {
      cursors.pop_back();
    }

    if (useRanges && !inRange(chrom, pos)) {
      continue;
    }

    if (covs.size() <= 1) {
      continue;
    }

    const size_t t = cov.total();
    if (t < p_C) {
      continue;
    }

    if (cov.pseudo_allele_count() <= 1) {
      continue;
    }

    auto dist = cov.distribution();

    double k = 0;
    double theta = 0;
    if (true) {
      auto itr = p_estimators.find(ctxt);
      if (itr == p_estimators.end()) {
        BOOST_LOG_TRIVIAL(fatal)
            << "unable to find context '" << scindo::kmers::render(6, ctxt)
            << "'. Index is inconsistent with data.";
        return;
      }
      k = itr->second.k;
      theta = itr->second.theta;
    }
    gamma_distribution<> G(k, theta);

    std::string ctxtStr = scindo::kmers::render(6, ctxt);

    for (auto itr = covs.begin(); itr != covs.end(); ++itr) {
      if (itr->second.total() < p_C) {
        continue;
      }
      auto dist1 = itr->second.distribution();
      auto kld = klDivergence(dist, dist1);
      double pv = cdf(complement(G, kld));
      if (pv > 1e-5) {
        continue;
      }

      summarizer S;
      ps.clear();
      for (auto jtr = dist1.begin(); jtr != dist1.end(); ++jtr) {
        S.add(jtr->second);
        ps.push_back(jtr->second);
      }
      std::sort(ps.rbegin(), ps.rend());
      double v = S.var();
      double h1 = fit_homozygous(ps);
      double h2 = fit_heterozygous(ps);
      double h = h1;
      bool isHet = false;
      if (h2 < h1) {
        h = h2;
        isHet = true;
      }
      if (h > 0.1) {
        continue;
      }
      const auto &c = itr->second;
      p_out << chrom << ':' << (pos + 1) << '\t' << covs.size() << '\t'
            << ctxtStr << '\t' << itr->first << '\t' << kld << '\t' << pv
            << '\t' << c.pseudo_allele_count() << '\t'
            << (isHet ? "het" : "hom") << '\t' << h << '\t' << c.bases[0]
            << '\t' << c.bases[1] << '\t' << c.bases[2] << '\t' << c.bases[3]
            << '\t' << c.bases[4] << '\t' << c.indels() << std::endl;
    }
  }
}

int main0(int argc, const char *argv[]) {
  boost::log::add_console_log(
      std::cerr, boost::log::keywords::format =
                     "[%TimeStamp%] [%ThreadID%] [%Severity%] %Message%");
  boost::log::add_common_attributes();

  std::map<std::string, docopt::value> opts =
      docopt::docopt(usage, {argv + 1, argv + argc}, true, "perya 0.1");

  for (auto itr = opts.begin(); itr != opts.end(); ++itr) {
    BOOST_LOG_TRIVIAL(info) << itr->first << '\t' << itr->second;
  }

  const size_t C = std::stoi(opts["--min-coverage"].asString());

  std::map<std::string, range_vector> R;
  if (opts["--regions"]) {
    BOOST_LOG_TRIVIAL(info)
        << "Loading regions from: " << opts["--regions"].asString();
    index_bed_file(opts["--regions"].asString(), R);
    BOOST_LOG_TRIVIAL(info) << "done.";
  }

  std::map<scindo::kmer, gamma> estimators;

  if (opts["--index"]) {
    scindo::input_file_holder_ptr inp =
        scindo::files::in(opts["--index"].asString());
    nlohmann::json x = nlohmann::json::parse(**inp);
    for (nlohmann::json::iterator itr = x.begin(); itr != x.end(); ++itr) {
      estimators[scindo::kmers::make(itr.key())] = gamma(itr.value());
    }
    BOOST_LOG_TRIVIAL(info) << "loaded " << estimators.size() << " contexts";
  } else {
    build_estimate(opts["<reference>"].asString(), opts["<BAM>"].asStringList(),
                   C, R, estimators);
  }

  if (opts["--save-index"]) {
    nlohmann::json res;
    for (auto itr = estimators.begin(); itr != estimators.end(); ++itr) {
      res[scindo::kmers::render(6, itr->first)] = itr->second.json();
    }
    scindo::output_file_holder_ptr outp =
        scindo::files::out(opts["--save-index"].asString());
    (**outp) << res << std::endl;
  }

  if (opts["--index-only"].asBool()) {
    BOOST_LOG_TRIVIAL(info) << "index information processed. Exiting.";
    return 0;
  }

  compute_significance(opts["<reference>"].asString(),
                       opts["<BAM>"].asStringList(), C, R, estimators,
                       std::cout);

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
