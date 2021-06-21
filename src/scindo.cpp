#include <iostream>
#include <sstream>
#include <unordered_set>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include "scindo/bam.hpp"
#include "scindo/gtf.hpp"
#include "scindo/stabby.hpp"
#include "scindo/vcf.hpp"
#include "scindo/summerizer.hpp"
#include "scindo/table.hpp"
#include "scindo/profile.hpp"

using namespace scindo;

namespace // anonymous
{
    std::string locus(const std::string& p_chrom, const scindo::stabby::interval& p_ivl)
    {
        std::ostringstream out;
        out << p_chrom << ':' << p_ivl.first << '-' << p_ivl.second;
        return out.str();
    }

    template <typename Itr>
    std::string join(const std::string& p_sep, Itr p_begin, Itr p_end)
    {
        std::string res;
        for (auto itr = p_begin; itr != p_end; ++itr)
        {
            if (itr != p_begin)
            {
                res += p_sep;
            }
            res += *itr;
        }
        return res;
    }

    const char usage[] =
R"(scindo - find allele specific expression

    Usage:
      scindo [options] <annotation-gtf> <vcf-file> <bam-file>

    Options:
      -h --help                     Show this screen
      -f FEATS  --features FEATS    Comma separated list GTF features to capture [default: gene]
      -c COV    --min-coverage NUM  Minumum coverage for sites [default: 10]
      -m NUM    --min-hets NUM      Minimum number of heterozygous sites for a gene [default: 5]
      -s SAMPLE --sample SAMPLE     Select the given sample to detect heterozygous sites.
)";

    static constexpr bool enabled = false;

    using vcf::allele;
    using gt_vector = vcf::vcf_file_reader::format_data<int32_t>;
    using diploid_vector = std::vector<std::pair<allele,allele>>;

    void repack_diploid(const gt_vector& p_raw_gts, diploid_vector& p_res, int p_idx)
    {
        p_res.clear();
        for (size_t i = 0; i < p_raw_gts.size(); i += 2)
        {
            if (p_idx < 0 || i == p_idx)
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
    void with(uint32_t p_pos, const scindo::bam::seq& p_seq, const scindo::bam::cigar& p_cig,
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

    int main0(int argc, const char* argv[])
    {
        boost::log::add_console_log(std::cerr);

        std::map<std::string, docopt::value>
            opts = docopt::docopt(usage, { argv + 1, argv + argc }, true, "scindo 0.1");

        for (auto itr = opts.begin(); itr != opts.end(); ++itr)
        {
            std::cerr << itr->first << '\t' << itr->second << std::endl;
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

        BOOST_LOG_TRIVIAL(info) << "features: " << join(",", features.begin(), features.end());

        std::unordered_map<std::string, std::vector<stabby::interval>> ivls;
        std::unordered_map<std::string, std::unordered_map<stabby::interval,stuff>> idx;
        {
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
                idx[p_seqname][ivl] = s;
                ivls[p_seqname].push_back(ivl);
            });
        }
        BOOST_LOG_TRIVIAL(info) << "features: " << join(",", features.begin(), features.end());
        std::unordered_map<std::string,std::shared_ptr<stabby>> annot;
        {
            profile<enabled> P("build annot");
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

        std::unordered_map<std::string,std::unordered_map<stabby::interval,std::vector<uint32_t>>> feature_positions;
        std::unordered_map<std::string,std::unordered_map<uint32_t,char>> feature_reference;
        {
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
                        std::string r = V.ref();
                        if (r.size() == 1)
                        {
                            feature_reference[chrom][pos] = r[0];
                        }
                    }
                    for (auto itr = hits.begin(); itr != hits.end(); ++itr)
                    {
                        if (!(itr->first <= pos && pos <= itr->second))
                        {
                            std::cerr << "failed: " << itr-> first << " <= " << pos << " <= " << itr->second << std::endl;
                            nlohmann::json blob;
                            blob["query"] = pos;
                            blob["intervals"] = ivls.at(chrom);
                            std::cerr << blob << std::endl;
                            throw std::runtime_error("bad stab!");
                        }
                        feature_positions[chrom][*itr].push_back(pos);
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
        std::unordered_map<std::string,std::vector<uint32_t>> positions;
        {
            profile<enabled> P("filter genes");
            for (auto itr = feature_positions.begin(); itr != feature_positions.end(); ++itr)
            {
                const auto& chrom = itr->first;
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    if (jtr->second.size() < min_hets)
                    {
                        continue;
                    }
                    std::vector<uint32_t>& pos_vec = positions[chrom];
                    pos_vec.insert(pos_vec.end(), jtr->second.begin(), jtr->second.end());
                }
                if (positions.contains(chrom))
                {
                    std::vector<uint32_t>& tmp = positions[chrom];
                    std::sort(tmp.begin(), tmp.end());
                    tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
                }
            }
        }
        std::unordered_map<std::string,std::unordered_map<uint32_t,std::unordered_map<char,uint32_t>>> counts;
        {
            profile<enabled> P("scan BAM");

            bam::bam_file_reader V(opts.at("<bam-file>").asString());

            const std::string* pchrom = NULL;
            uint32_t ppos = 0;
            size_t chrom_hit_count = 0;
            std::vector<uint32_t>::const_iterator cursor;
            std::vector<uint32_t>::const_iterator end;

            while (V.next())
            {
                const std::string& chrom = V.chrom();
                if (&chrom != pchrom)
                {
                    if (!positions.contains(chrom))
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
                    cursor = positions[chrom].begin();
                    end = positions[chrom].end();
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
                with(pos, V.seq(), V.cigar(),
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
                                auto b = p_bases[i];
                                counts[chrom][p][b] += 1;
                                chrom_hit_count += 1;
                            }
                        }
                    },
                    [&](uint32_t p_ref, uint32_t p_seq, const bam::seq& p_bases) {
                        // Insertion
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

            scindo::table<std::string, std::string, std::string, double, uint32_t, uint32_t, uint32_t, uint32_t, double>
                tbl({"geneName", "geneId", "locus", "majFrac", "coverage", "major", "minor", "numHets", "pValue"});

            std::vector<double> chi2s;
            summerizer cov_summary;
            summerizer frac_summary;
            std::vector<uint32_t> majors;
            std::vector<uint32_t> minors;
            std::vector<std::pair<uint32_t,char>> qq;
            for (auto itr = feature_positions.begin(); itr != feature_positions.end(); ++itr)
            {
                const auto& chrom = itr->first;
                if (!counts.contains(chrom))
                {
                    continue;
                }
                for (auto jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr)
                {
                    const auto& ivl = jtr->first;
                    chi2s.clear();
                    cov_summary.clear();
                    frac_summary.clear();
                    majors.clear();
                    minors.clear();
                    uint64_t maj_sum = 0;
                    uint64_t mor_sum = 0;
                    for (auto ktr = jtr->second.begin(); ktr != jtr->second.end(); ++ktr)
                    {
                        const auto& pos = *ktr;
                        if (!counts.at(chrom).contains(pos))
                        {
                            continue;
                        }
                        const auto& baseCounts = counts.at(chrom).at(pos);
                        uint32_t tot = 0;
                        uint32_t n = 0;
                        uint32_t ref = 0;
                        char ref_base = 'N';
                        if (feature_reference.at(chrom).contains(pos))
                        {
                            ref_base = feature_reference.at(chrom).at(pos);
                        }
                        qq.clear();
                        for (auto ltr = baseCounts.begin(); ltr != baseCounts.end(); ++ltr, ++n)
                        {
                            if (ltr->first == ref_base)
                            {
                                ref = ltr->second;
                            }
                            else
                            {
                                qq.push_back(std::make_pair(ltr->second, ltr->first));
                            }
                            tot += ltr->second;
                        }
                        if (qq.size() < 1 || tot < min_cov)
                        {
                            continue;
                        }
                        std::sort(qq.rbegin(), qq.rend());
                        uint32_t alt = qq[0].first;
                        double m = (ref + alt)/2.0;
                        double chi2 = 0;
                        if (ref > 0)
                        {
                            double p = ref/m;
                            chi2 += ref * std::log(p);
                        }
                        if (alt > 0)
                        {
                            double q = alt/m;
                            chi2 += alt * std::log(q);
                        }
                        chi2 *= 2;
                        chi2s.push_back(chi2);
                        cov_summary.push_back(ref + alt);
                        majors.push_back(std::max(ref, alt));
                        minors.push_back(std::min(ref, alt));
                        frac_summary.push_back(double(majors.back())/double(ref + alt));
                        maj_sum += majors.back();
                        mor_sum += minors.back();
                    }
                    if (chi2s.size() < min_hets)
                    {
                        continue;
                    }
                    uint32_t n = chi2s.size();

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

                    tbl.push_back({
                            idx.at(chrom).at(ivl).gene_name, idx.at(chrom).at(ivl).gene_id, locus(chrom, ivl),
                            double(maj_sum)/double(maj_sum + mor_sum), maj_sum + mor_sum, maj_sum, mor_sum, n, pv});
                }
            }
            tbl.sort({9, -4, 1});
            tbl.write(std::cout);
        }
        profile<enabled>::report();

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
