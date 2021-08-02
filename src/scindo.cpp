#include <iostream>
#include <sstream>
#include <unordered_set>
#include <boost/flyweight.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include "scindo/bam.hpp"
#include "scindo/benjamini_hochberg.hpp"
#include "scindo/gtf.hpp"
#include "scindo/stabby.hpp"
#include "scindo/tsv.hpp"
#include "scindo/vcf.hpp"
#include "scindo/summerizer.hpp"
#include "scindo/table.hpp"
#include "scindo/profile.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(scindo - find allele specific expression

    Usage:
      scindo [options] <annotation-gtf> <vcf-file> <bam-file>

    Options:
      -h --help                         Show this screen
      -e FRAC,      --faction FRAC      Minimum fraction for calling allele specific expression [default: 0.6]
      -f FEATS,     --features FEATS    Comma separated list GTF features to capture [default: gene]
      -c COV,       --min-coverage NUM  Minumum coverage for sites [default: 10]
      -m NUM,       --min-hets NUM      Minimum number of heterozygous sites for a gene [default: 5]
      -o FILE,      --output-file FILE  File to write the output to [default: -]
      -q FDR,       --false-discovery-rate FDR False discovery rate to use [default: 0.05]
      -s SAMPLE,    --sample SAMPLE     Select the given sample to detect heterozygous sites.
      -g FILE,      --gene-list FILE    A TSV containing gene names to report against.
      -G FILE,      --exclusion FILE    A TSV containing gene names to exclude.
)";

    using allele_seq = boost::flyweight<std::string>;

    allele_seq allele_seq_1(char p_ch)
    {
        std::string s;
        s.push_back(p_ch);
        return allele_seq(s);
    }

    std::string locus(const std::string& p_chrom, const scindo::stabby::interval& p_ivl)
    {
        std::ostringstream out;
        out << p_chrom << ':' << p_ivl.first << '-' << p_ivl.second;
        return out.str();
    }

    template <typename Sep, typename Itr>
    std::string join(const Sep& p_sep, Itr p_begin, Itr p_end)
    {
        std::ostringstream out;
        for (auto itr = p_begin; itr != p_end; ++itr)
        {
            if (itr != p_begin)
            {
                out << p_sep;
            }
            out << *itr;
        }
        return out.str();
    }

    template <typename T>
    std::string join(const std::vector<T>& p_items, char p_sep = ',')
    {
        return join(p_sep, p_items.begin(), p_items.end());
    }

    std::string joinHap(const std::vector<int>& p_gts)
    {
        std::ostringstream out;
        for (auto itr = p_gts.begin(); itr != p_gts.end(); ++itr)
        {
            if (itr != p_gts.begin())
            {
                out << '/';
            }
            if (*itr >= 0)
            {
                out << *itr;
            }
            else
            {
                out << '.';
            }
        }
        return out.str();
    }

    std::unordered_set<std::string> gene_list(const std::string& p_filename)
    {
        std::unordered_set<std::string> res;
        {
            std::ifstream in(p_filename);
            scindo::tsv T(in);
            size_t cn = 0;
            if (T.idx.contains("Gene Symbol"))
            {
                cn = T.idx.at("Gene Symbol");
            }
            else if (T.idx.contains("gene"))
            {
                cn = T.idx.at("gene");
            }
            for (auto itr = T.begin(); itr != T.end(); ++itr)
            {
                res.insert((*itr)[cn]);
            }
        }
        return res;
    }

    static constexpr bool enabled = false;

    using vcf::allele;
    using gt_vector = vcf::vcf_file_reader::format_data<int32_t>;
    using diploid_vector = std::vector<std::pair<allele,allele>>;

    void repack_diploid(const gt_vector& p_raw_gts, diploid_vector& p_res, int p_idx)
    {
        p_res.clear();
        for (size_t i = 0; i < p_raw_gts.size(); i += 2)
        {
            if (p_idx < 0 || i == 2*p_idx)
            {
                p_res.push_back(std::make_pair(allele(p_raw_gts[i]), allele(p_raw_gts[i+1])));
            }
        }
    }

    struct stuff
    {
        std::string gene_id;
        std::string gene_name;
    };

    using match = std::function<void(uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases)>;
    using ins = std::function<void(uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases)>;
    using del = std::function<void(uint32_t p_ref, uint32_t p_seq, uint32_t p_len)>;
    using skip = std::function<void(uint32_t p_ref, uint32_t p_seq, uint32_t p_len)>;
    using clip = std::function<void(uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases)>;

    template <typename XM, typename XI, typename XD, typename XS, typename XC>
    void with(const std::string& p_chrom, uint32_t p_pos, const scindo::bam::seq& p_seq, const scindo::bam::cigar& p_cig,
              XM p_match, XI p_ins, XD p_del, XS p_skip, XC p_clip)
    {
        static_assert(std::is_convertible<XM,match>::value);
        static_assert(std::is_convertible<XI,ins>::value);
        static_assert(std::is_convertible<XD,del>::value);
        static_assert(std::is_convertible<XS,skip>::value);
        static_assert(std::is_convertible<XC,clip>::value);

        uint32_t r = p_pos;
        uint32_t q = 0;
        for (size_t i = 0; i < p_cig.size(); ++i)
        {
            char op = p_cig.op(i);
            uint32_t len = p_cig.len(i);

            switch (op)
            {
                case 'M':
                case '=':
                case 'X':
                {
                    p_match(r, q, p_seq.range(q, q+len));
                    r += len;
                    q += len;
                    break;
                }
                case 'I':
                {
                    p_ins(r, q, p_seq.range(q, q+len));
                    q += len;
                    break;
                }
                case 'D':
                {
                    p_del(r, q, len);
                    r += len;
                    break;
                }
                case 'N':
                {
                    p_skip(r, q, len);
                    r += len;
                    break;
                }
                case 'S':
                {
                    p_clip(r, q, p_seq.range(q, q+len));
                    q += len;
                    break;
                }
            }
        }
    }

    using ref_and_alts = std::vector<allele_seq>;

    int find_allele(const ref_and_alts& p_alleles, const allele_seq& p_seq,
                    const std::string& p_chrom, const uint32_t& p_pos, const uint32_t p_count)
    {
        for (int i = 0; i < p_alleles.size(); ++i)
        {
            if (p_alleles[i] == p_seq)
            {
                return i;
            }
        }
        //std::cerr << "find_allele: " << p_seq << " in " << nlohmann::json(p_alleles)
        //    << '\t' << p_chrom
        //    << '\t' << p_pos
        //    << '\t' << p_count
        //    << std::endl;
        return -1;
    }

    int main0(int argc, const char* argv[])
    {
        boost::log::add_console_log(std::cerr, boost::log::keywords::format = "[%TimeStamp%] [%ThreadID%] [%Severity%] %Message%");
        boost::log::add_common_attributes();

        std::map<std::string, docopt::value>
            opts = docopt::docopt(usage, { argv + 1, argv + argc }, true, "scindo 0.1");

        for (auto itr = opts.begin(); itr != opts.end(); ++itr)
        {
            BOOST_LOG_TRIVIAL(info) << itr->first << '\t' << itr->second;
        }

        std::unordered_set<std::string> features;
        {
            std::string s = "gene";
            if (opts.at("-f"))
            {
                s = opts.at("-f").asString();
            }
            size_t p = 0;
            size_t i = s.find(',');
            while (i != std::string::npos)
            {
                features.insert(std::string(s.begin() + p, s.begin() + i));
                p = i + 1;
                i = s.find(',', p);
            }
            features.insert(std::string(s.begin() + p, s.end()));
        }

        std::unordered_set<std::string> wanted_genes;
        if (opts.at("-g"))
        {
            wanted_genes = gene_list(opts.at("-g").asString());
        }

        std::unordered_set<std::string> unwanted_genes;
        if (opts.at("-G"))
        {
            unwanted_genes = gene_list(opts.at("-G").asString());
        }

        std::unordered_map<std::string, std::vector<stabby::interval>> ivls;
        std::unordered_map<std::string, std::unordered_map<stabby::interval,stuff>> idx;
        {
            BOOST_LOG_TRIVIAL(info) << "scanning annotation: " << opts.at("<annotation-gtf>").asString();
            profile<enabled> P("scan GTF");

            gtf::gtf_file G(opts.at("<annotation-gtf>").asString());
            G.parse([&](const std::string& p_seqname, const std::string& p_source, const std::string& p_feature,
                        const int64_t& p_start, const int64_t& p_end,
                        const std::optional<double>& p_score, const char& p_strand,
                        const std::optional<int>& p_frame, const std::vector<gtf::attribute>& p_attrs) {
                if (!features.contains(p_feature))
                {
                    return;
                }
                stabby::interval ivl = {p_start, p_end};

                stuff s;
                for (auto itr = p_attrs.begin(); itr != p_attrs.end(); ++itr)
                {
                    if (itr->first == "gene_id")
                    {
                        s.gene_id = std::get<std::string>(itr->second);
                    }
                    if (itr->first == "gene_name")
                    {
                        s.gene_name = std::get<std::string>(itr->second);
                    }
                    if (itr->first == "gene_type")
                    {
                        if (std::get<std::string>(itr->second) != "protein_coding")
                        {
                            return;
                        }
                    }
                }

                if (wanted_genes.size() > 0 && !wanted_genes.contains(s.gene_name))
                {
                    return;
                }
                else if (unwanted_genes.contains(s.gene_name))
                {
                    return;
                }

                idx[p_seqname][ivl] = s;
                ivls[p_seqname].push_back(ivl);
            });
        }

        std::unordered_map<std::string,std::shared_ptr<stabby>> annot;
        {
            profile<false> P("build annot");
            for (auto itr = ivls.begin(); itr != ivls.end(); ++itr)
            {
                const auto& chrom = itr->first;
                std::vector<stabby::interval>& raw = itr->second;
                std::sort(raw.begin(), raw.end());
                raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
                annot[chrom] = std::shared_ptr<stabby>(new stabby(raw));
            }
            //ivls.clear();
        }

        std::optional<std::string> sample_id;
        if (opts.at("-s"))
        {
            sample_id = opts.at("-s").asString();
        }

        std::unordered_map<std::string,std::unordered_map<stabby::interval,std::vector<uint32_t>>> vcf_positions;
        std::unordered_map<std::string,std::unordered_map<uint32_t,ref_and_alts>> position_seqs;
        {
            BOOST_LOG_TRIVIAL(info) << "scanning vcf: " << opts.at("<vcf-file>").asString();
            profile<enabled> P("scan VCF");

            vcf::vcf_file_reader V(opts.at("<vcf-file>").asString());
            std::vector<std::string> samples = V.samples();
            int sample_idx = -1;
            if (sample_id)
            {
                const std::string& wanted = sample_id.value();
                for (int i = 0; i < samples.size(); ++i)
                {
                    if (samples[i] == wanted)
                    {
                        sample_idx = i;
                        break;
                    }
                }
                if (sample_idx < 0)
                {
                        BOOST_LOG_TRIVIAL(warning) << "requested sample '" << wanted << "' not found. Using all samples.";
                }
            }
            gt_vector gt;
            diploid_vector D;
            std::vector<stabby::interval> hits;
            while (V.next())
            {
                if(!V.format_values("GT", gt))
                {
                    std::cerr << "error getting GT" << std::endl;
                    return -1;
                }
                bool hets_found = false;
                {
                    profile<enabled> P1("massaging");
                    repack_diploid(gt, D, sample_idx);
                    for (auto itr = D.begin(); itr != D.end(); ++itr)
                    {
                        if (*(itr->first) != *(itr->second))
                        {
                            hets_found = true;
                            break;
                        }
                    }
                }
                if (hets_found)
                {
                    const std::string& chrom = V.chrom();
                    if (!annot.contains(chrom))
                    {
                        continue;
                    }
                    uint32_t pos = V.pos1();
                    {
                        profile<enabled> P1("stabbing");
                        hits.clear();
                        annot.at(chrom)->find(pos, hits);
                    }
                    if (hits.size() > 0)
                    {
                        ref_and_alts& aa = position_seqs[chrom][pos];
                        const size_t n = V.num_alts();
                        for (size_t i = 0; i <= n; ++i)
                        {
                            std::string x = V.alt(i);
                            bool found = false;
                            for (size_t j = 0; j < aa.size(); ++j)
                            {
                                if (aa[j] == x)
                                {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found)
                            {
                                aa.push_back(allele_seq(x));
                            }
                        }
                        //std::cerr << chrom << '\t' << pos << '\t' << nlohmann::json(aa) << std::endl;
                    }
                    for (auto itr = hits.begin(); itr != hits.end(); ++itr)
                    {
                        if (!(itr->first <= pos && pos <= itr->second))
                        {
                            std::cerr << "failed: " << itr-> first << " <= " << pos << " <= " << itr->second << std::endl;
                            nlohmann::json blob;
                            blob["query"] = pos;
                            if (ivls.contains(chrom))
                            {
                                blob["intervals"] = ivls.at(chrom);
                            }
                            std::cerr << blob << std::endl;
                            throw std::runtime_error("bad stab!");
                        }
                        vcf_positions[chrom][*itr].push_back(pos);
                    }
                }
            }
        }
        size_t min_hets = 5;
        if (opts.at("-m"))
        {
            min_hets = opts.at("-m").asLong();
        }
        size_t min_cov = 10;
        if (opts.at("-c"))
        {
            min_cov = opts.at("-c").asLong();
        }
        double Frac = 0.6;
        if (opts.at("-e"))
        {
            Frac = std::stod(opts.at("-e").asString());
        }
        double Q = 0.05;
        if (opts.at("-q"))
        {
            Q = std::stod(opts.at("-q").asString());
        }
        size_t num_low_hets = 0;
        std::unordered_map<std::string,std::vector<uint32_t>> het_positions;
        {
            profile<enabled> P("filter genes");
            for (auto itr = vcf_positions.begin(); itr != vcf_positions.end(); ++itr)
            {
                const auto& chrom = itr->first;
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    if (jtr->second.size() < min_hets)
                    {
                        num_low_hets += 1;
                        continue;
                    }
                    std::vector<uint32_t>& pos_vec = het_positions[chrom];
                    pos_vec.insert(pos_vec.end(), jtr->second.begin(), jtr->second.end());
                }
                if (het_positions.contains(chrom))
                {
                    std::vector<uint32_t>& tmp = het_positions[chrom];
                    std::sort(tmp.begin(), tmp.end());
                    tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
                }
            }
        }
        BOOST_LOG_TRIVIAL(info) << "number of genes with too few hets: " << num_low_hets;
        std::unordered_map<std::string,std::unordered_map<uint32_t,std::unordered_map<allele_seq,uint32_t>>> counts;
        {
            BOOST_LOG_TRIVIAL(info) << "scanning bam: " << opts.at("<bam-file>").asString();

            using scindo::bam::flag;
            profile<enabled> P("scan BAM");

            bam::bam_file_reader V(opts.at("<bam-file>").asString());

            const std::string* pchrom = NULL;
            uint32_t ppos = 0;
            size_t chrom_hit_count = 0;
            std::vector<uint32_t>::const_iterator cursor;
            std::vector<uint32_t>::const_iterator end;

            while (V.next())
            {
                flag flg(V.flag());
                if (flg.is<flag::unmapped>() || flg.is<flag::duplicate>() || !flg.is<flag::proper_pair>())
                {
                    continue;
                }
                const std::string& chrom = V.chrom();

                if (&chrom != pchrom)
                {
                    if (!het_positions.contains(chrom))
                    {
                        continue;
                    }
                    if (pchrom != NULL)
                    {
                        BOOST_LOG_TRIVIAL(info) << "number of hits: " << chrom_hit_count;
                    }
                    BOOST_LOG_TRIVIAL(info) << "chromosome " << chrom;
                    pchrom = &chrom;
                    ppos = 0;
                    chrom_hit_count = 0;
                    cursor = het_positions[chrom].begin();
                    end = het_positions[chrom].end();
                }

                if (cursor == end)
                {
                    // no more positions of interest on this chromosome
                    continue;
                }

                uint32_t pos = V.pos1();

                if (pos < ppos)
                {
                    std::cerr << "BAM file not sorted" << std:: endl;
                    return -1;
                }

                // Scan forward to the right part of the position vector.
                // Yes, this is linear search, but it's monotonic, so it's
                // once through the list for the whole bam/chromosome.
                //
                while (cursor != end && *cursor < pos)
                {
                    ++cursor;
                }
                if (cursor == end)
                {
                    continue;
                }
                auto tmp = cursor;
                with(chrom, pos, V.seq(), V.cigar(),
                    [&](uint32_t p_ref, uint32_t p_seq, const bam::seq& p_bases) {
                        // Match
                        for (uint32_t i = 0; i < p_bases.size(); ++i)
                        {
                            uint32_t p = p_ref + i;
                            while (tmp != end && *tmp < p)
                            {
                                ++tmp;
                            }
                            if (tmp == end)
                            {
                                break;
                            }
                            if (*tmp == p)
                            {
                                allele_seq b = allele_seq_1(p_bases[i]);
                                counts[chrom][p][b] += 1;
                                chrom_hit_count += 1;
                            }
                        }
                    },
                    [&](uint32_t p_ref, uint32_t p_seq, const bam::seq& p_bases) {
                        // Insertion
                        while (tmp != end && *tmp < p_ref)
                        {
                            ++tmp;
                        }
                        if (*tmp == p_ref)
                        {
                            allele_seq a(p_bases.asString());
                            counts[chrom][p_ref][a] += 1;
                            chrom_hit_count += 1;
                        }
                    },
                    [&](uint32_t p_ref, uint32_t p_seq, uint32_t p_len) {
                        // Deletion
                    },
                    [&](uint32_t p_ref, uint32_t p_seq, uint32_t p_len) {
                        // Split read
                    },
                    [&](uint32_t p_ref, uint32_t p_seq, const bam::seq& p_bases) {
                        // Soft clip
                    });
            }
            if (pchrom != NULL)
            {
                BOOST_LOG_TRIVIAL(info) << "number of hits: " << chrom_hit_count;
            }

            scindo::table<std::string, std::string, std::string, uint32_t,
                          uint32_t, uint32_t, double, double, double, double, double, std::string, std::string, std::string>
                tbl({"geneName", "geneId", "locus", "numHets",
                     "highCount", "lowCount", "avgCov", "covSd", "highFrac", "pValue", "qValue",
                     "highHaplo", "lowHaplo", "positions"});

            std::vector<double> chi2s;
            summerizer majors;
            summerizer minors;
            summerizer totals;
            std::vector<int> maj_alleles;
            std::vector<int> mor_alleles;
            std::vector<std::pair<uint32_t,allele_seq>> qq;
            std::vector<uint32_t> used_positions;
            allele_seq no_such_allele = allele_seq(std::string("N"));
            for (auto itr = vcf_positions.begin(); itr != vcf_positions.end(); ++itr)
            {
                const auto& chrom = itr->first;
                if (!counts.contains(chrom))
                {
                    continue;
                }
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    const auto& ivl = jtr->first;
                    const auto& positions = jtr->second;
                    chi2s.clear();
                    majors.clear();
                    minors.clear();
                    totals.clear();
                    maj_alleles.clear();
                    mor_alleles.clear();
                    used_positions.clear();
                    for (auto ktr = positions.begin(); ktr != positions.end(); ++ktr)
                    {
                        const auto& pos = *ktr;
                        if (!counts.contains(chrom) || !counts.at(chrom).contains(pos))
                        {
                            continue;
                        }
                        const auto& alleleCounts = counts.at(chrom).at(pos);
                        uint32_t tot = 0;
                        const auto& allele_seqs = position_seqs.at(chrom).at(pos);

                        qq.clear();
                        qq.push_back(std::make_pair(0, no_such_allele));
                        qq.push_back(std::make_pair(0, no_such_allele));
                        for (auto ltr = alleleCounts.begin(); ltr != alleleCounts.end(); ++ltr)
                        {
                            qq.push_back(std::make_pair(ltr->second, ltr->first));
                            tot += ltr->second;
                        }

                        std::sort(qq.rbegin(), qq.rend());
                        uint32_t maj = qq[0].first;
                        const auto& majAllele = qq[0].second;
                        uint32_t mor = qq[1].first;
                        const auto& morAllele = qq[1].second;
                        double m = (maj + mor)/2.0;
                        double chi2 = 0;
                        if (maj < min_cov)
                        {
                            continue;
                        }
                        if (maj > 0)
                        {
                            double p = maj/m;
                            chi2 += maj * std::log(p);
                        }
                        if (mor > 0)
                        {
                            double q = mor/m;
                            chi2 += mor * std::log(q);
                        }
                        chi2 *= 2;
                        chi2s.push_back(chi2);

                        majors.push_back(maj);
                        minors.push_back(mor);
                        totals.push_back(maj + mor);

                        int majAlleleNum = find_allele(allele_seqs, majAllele, chrom, pos, maj);
                        int morAlleleNum = find_allele(allele_seqs, morAllele, chrom, pos, mor);

                        maj_alleles.push_back(majAlleleNum);
                        mor_alleles.push_back(morAlleleNum);
                        used_positions.push_back(pos);
                    }
                    uint32_t n = chi2s.size();
                    if (n == 0 || majors.sum() + minors.sum() == 0)
                    {
                        continue;
                    }
                    if (used_positions.size() < min_hets)
                    {
                        continue;
                    }

                    // Compute the p-value as aScan
                    //
                    double chi2sum = 0;
                    for (auto i = 0; i < n; ++i)
                    {
                        chi2sum += chi2s[i];
                    }
                    using namespace boost::math;
                    chi_squared_distribution<long double> dist(n);
                    double pv = cdf(complement(dist, chi2sum));
                    double majorFrac = majors.sum() / (majors.sum() + minors.sum());

                    tbl.push_back({
                            idx.at(chrom).at(ivl).gene_name, idx.at(chrom).at(ivl).gene_id, locus(chrom, ivl), used_positions.size(),
                            majors.sum(), minors.sum(), totals.mean(), totals.sd(), majorFrac, pv, 0,
                            joinHap(maj_alleles), joinHap(mor_alleles), join(used_positions)});
                }
            }
            {
                std::vector<double> pvals;
                pvals.reserve(tbl.size());
                for (size_t i = 0; i < tbl.size(); ++i)
                {
                    pvals.push_back(std::get<9>(tbl[i]));
                }
                std::vector<double> qvals;
                benjamini_hochberg(static_cast<const double&>(Q))(pvals, qvals);
                for (size_t i = 0; i < tbl.size(); ++i)
                {
                    std::get<10>(tbl[i]) = qvals[i];
                }

            }
            tbl.filter([&](const auto& p_row) {
                double frac = std::get<8>(p_row);
                if (frac < Frac)
                {
                    return false;
                }
                double pval = std::get<9>(p_row);
                double qval = std::get<10>(p_row);
                return pval < qval;
            });
            tbl.sort(std::vector<std::string>({"pValue", "-highFrac"}));

            output_file_holder_ptr outp = files::out(opts.at("-o").asString());
            tbl.write(**outp);
        }
        profile<true>::report();

        return 0;
    }

}
// namespace anonymous

int main(int argc, const char* argv[])
{
    try
    {
        return main0(argc, argv);
    }
    catch (std::exception e)
    {
        std::cerr << e.what();
        return -1;
    }
}
