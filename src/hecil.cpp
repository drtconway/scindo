#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <boost/flyweight.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include "scindo/gtf.hpp"
#include "scindo/stabby.hpp"
#include "scindo/files.hpp"
#include "scindo/fasta.hpp"
#include "scindo/kmers.hpp"
#include "scindo/annotated_kmer_set.hpp"
#include "scindo/profile.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(hecil - in silico depletion for RNASeq

    Usage:
      hecil index [options] <annotation-file> <reference-fasta>

    Options:
      -h --help                         Show this screen
)";

    struct cursor
    {
        const std::vector<stabby::interval>& ivls;
        const std::vector<uint64_t>& annot;
        size_t curr;
        std::deque<size_t> live;

        cursor(const std::vector<stabby::interval>& p_ivls, const std::vector<uint64_t>& p_annot)
            : ivls(p_ivls), annot(p_annot), curr(0)
        {
        }

        void seek(const uint32_t& p_pos)
        {
            while (curr < ivls.size() && ivls[curr].first <= p_pos)
            {
                live.push_back(curr);
                ++curr;
            }
            std::sort(live.begin(), live.end(), [&](auto i, auto j) {
                if (ivls[i].second != ivls[j].second)
                {
                    return ivls[i].second < ivls[j].second;
                }
                return ivls[i].first < ivls[j].first;

            });
            while (live.size() > 0 && ivls[live.front()].second < p_pos)
            {
                live.pop_front();
            }
        }
    };

    struct stringdex
    {
        std::vector<std::string> things;
        std::unordered_map<std::string,size_t> index;

        const std::string& operator[](size_t p_idx) const
        {
            return things[p_idx];
        }

        size_t operator[](const std::string& p_thing)
        {
            auto itr = index.find(p_thing);
            if (itr != index.end())
            {
                return itr->second;
            }
            size_t n = index.size();
            things.push_back(p_thing);
            index[p_thing] = n;
            return n;
        }

        std::vector<std::string> decode(const uint64_t p_set) const
        {
            std::vector<std::string> res;
            uint64_t x = p_set;
            for (size_t i = 0; i < 64 && x > 0; ++i, x >>= 1)
            {
                if (x & 1)
                {
                    res.push_back(things[i]);
                }
            }
            std::sort(res.begin(), res.end());
            return res;
        }
    };

    using annot_num = size_t;
    using annotation = uint64_t;

    struct chrom_data
    {
        const size_t N;
        std::vector<stabby::interval> ivls;
        std::vector<annotation> annot;

        chrom_data(std::unordered_map<stabby::interval,annotation>& p_annots)
            : N(p_annots.size())
        {
            profile<true> P("constructing chrom_data");

            ivls.reserve(N);
            for (auto itr = p_annots.begin(); itr != p_annots.end(); ++itr)
            {
                ivls.push_back(itr->first);
            }
            std::sort(ivls.begin(), ivls.end());

            annot.reserve(N);
            for (auto itr = ivls.begin(); itr != ivls.end(); ++itr)
            {
                annot.push_back(p_annots.at(*itr));
            }
            p_annots.clear();
        }
    };
    using chrom_data_ptr = std::shared_ptr<chrom_data>;

    std::string first_word(const std::string& p_str)
    {
        auto n = p_str.find(' ');
        if (n == std::string::npos)
        {
            return p_str;
        }
        return p_str.substr(0, n);
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

        const uint64_t K = 25;

        stringdex idx;
        std::unordered_map<std::string,chrom_data_ptr> all_chroms;
        std::string prev_chrom;
        std::unordered_map<stabby::interval,annotation> annots;
        uint64_t gene_types = 0;
        uint64_t feature_types = 0;
        uint64_t transcript_types = 0;
        if (1)
        {
            profile<true> P("scan GTF");

            std::string gtf_name = opts.at("<annotation-file>").asString();
            gtf::gtf_file G(gtf_name);
            G.parse([&](const std::string& p_seqname, const std::string& p_source, const std::string& p_feature,
                        const int64_t& p_start, const int64_t& p_end,
                        const std::optional<double>& p_score, const char& p_strand,
                        const std::optional<int>& p_frame, const std::vector<gtf::attribute>& p_attrs) {

                if (prev_chrom != p_seqname)
                {
                    if (annots.size() > 0)
                    {
                        BOOST_LOG_TRIVIAL(info) << "saving: " << prev_chrom << " (" << annots.size() << ")";
                        all_chroms[prev_chrom] = chrom_data_ptr(new chrom_data(annots));
                    }
                    BOOST_LOG_TRIVIAL(info) << "scanning: " << p_seqname;
                    prev_chrom = p_seqname;
                }
                size_t ft = 1ULL << idx[p_feature];
                size_t gt = 0;
                size_t tt = 0;
                for (auto itr = p_attrs.begin(); itr != p_attrs.end(); ++itr)
                {
                    if (itr->first == "gene_type")
                    {
                        gt |= 1ULL << idx[std::get<std::string>(itr->second)];
                    }
                    if (itr->first == "transcript_type")
                    {
                        tt = 1ULL << idx[std::get<std::string>(itr->second)];
                    }
                }
                uint64_t ann = ft | gt | tt;
                annots[stabby::interval(p_start, p_end)] |= ann;
                feature_types |= ft;
                gene_types |= gt;
                transcript_types |= tt;
            });
            if (annots.size() > 0)
            {
                BOOST_LOG_TRIVIAL(info) << "saving: " << prev_chrom << " (" << annots.size() << ")";
                all_chroms[prev_chrom] = chrom_data_ptr(new chrom_data(annots));
            }
        }
        const uint64_t rRNA_mask = 1ULL << idx["rRNA"];
        if (1)
        {
            std::string ref_name = opts.at("<reference-fasta>").asString();
            input_file_holder_ptr in_ptr = files::in(ref_name);
            size_t xc = 0;
            annotated_kmer_set X;
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                std::string chrom = first_word(r.first);
                if (!all_chroms.contains(chrom))
                {
                    BOOST_LOG_TRIVIAL(warning) << "no annotations for " << chrom;
                    continue;
                }

                cursor cur(all_chroms.at(chrom)->ivls, all_chroms.at(chrom)->annot);

                BOOST_LOG_TRIVIAL(info) << "scanning: " << chrom << " (" << r.second.size() << ")";

                scindo::kmers::make_canonical(r.second, K, [&](const scindo::kmer_and_pos& p_xp) {
                    const auto x = p_xp.first;
                    const auto p = p_xp.second;
                    cur.seek(p);
                    if (cur.live.size() == 0)
                    {
                        return;
                    }
                    xc += 1;
                    uint64_t ann = 0;
                    for (auto itr = cur.live.begin(); itr != cur.live.end(); ++itr)
                    {
                        ann |= cur.annot[*itr];
                    }
                    //std::cout << std::endl;
                    auto fa = ann & feature_types;
                    auto ga = ann & gene_types;
                    X.add(x, fa | ga);
                });
                X.flush();
                BOOST_LOG_TRIVIAL(info) << "found: " << xc << " kmers of interest";
                BOOST_LOG_TRIVIAL(info) << "found: " << X.size() << " distinct kmers of interest";
                if (1)
                {
                    std::cout << ">" << chrom << std::endl;
                    std::unordered_map<uint64_t,size_t> counts;
                    for (size_t i = 0; i < X.seen.size(); ++i)
                    {
                        const auto& anns = X.seen[i]->annot;
                        for (size_t j = 0; j < anns.size(); ++j)
                        {
                            counts[anns[j]] += 1;
                        }
                    }
                    for (auto itr = counts.begin(); itr != counts.end(); ++itr)
                    {
                        std::cout << itr->second << '\t' << nlohmann::json(idx.decode(itr->first)) << std::endl;
                    }
                }
            }
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
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;;
        return -1;
    }
}
