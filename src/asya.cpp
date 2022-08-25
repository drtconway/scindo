#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/flyweight.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include "scindo/fasta.hpp"
#include "scindo/files.hpp"
#include "scindo/gtf.hpp"
#include "scindo/kmers.hpp"
#include "scindo/murmur3.hpp"
#include "scindo/stabby.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(asya - scan genome for bait kmers

    Usage:
      asya [options] <reference-fasta> <bait-fasta> <annotation-gtf>

    Options:
      -h --help                         Show this screen
)";

    using label = boost::flyweight<std::string>;

    std::string first_word(const std::string& p_str)
    {
        auto n = p_str.find(' ');
        if (n == std::string::npos)
        {
            return p_str;
        }
        return p_str.substr(0, n);
    }

    template <size_t I>
    double entropy(size_t K, kmer x)
    {
        constexpr size_t N = 1ULL << (2*I);
        constexpr kmer M = N - 1;
        std::vector<size_t> cs(N, 0);
        for (size_t i = 0; i < K + 1 - I; ++i)
        {
            cs[x&M] += 1;
            x >>= 2;
        }
        double e = 0;
        for (size_t i = 0; i < N; ++i)
        {
            double p = double(cs[i])/K;
            if (p == 0)
            {
                continue;
            }
            e += -p*std::log2(p);
        }
        return e;
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

        constexpr size_t K = 27;

        std::unordered_map<kmer,std::vector<label>> idx;
        std::unordered_map<label,std::string> bait_labels;
        {
            BOOST_LOG_TRIVIAL(info) << "scanning baits: " << opts.at("<bait-fasta>").asString();

            std::string ref_name = opts.at("<bait-fasta>").asString();
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                label lab = label(first_word(r.first));
                bait_labels[lab] = r.first;
                kmers::make(r.second, K, [&](kmer p_x, kmer p_xb) {
                    //if (entropy<3>(K, p_x) < 1.5)
                    //{
                    //    return;
                    //}
                    idx[p_x].push_back(lab);
                    idx[p_xb].push_back(lab);
                });
            }

            for (auto itr = idx.begin(); itr != idx.end(); ++itr)
            {
                auto& vec = itr->second;
                std::sort(vec.begin(), vec.end());
                vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
            }
        }

        std::unordered_map<label,label> gene_names;
        std::unordered_map<label,std::unordered_map<stabby::interval,label>> annot;
        std::unordered_map<label,std::shared_ptr<stabby>> annot_idx;
        if (0)
        {
            BOOST_LOG_TRIVIAL(info) << "scanning annotation: " << opts.at("<annotation-gtf>").asString();
            gtf::gtf_file G(opts.at("<annotation-gtf>").asString());

            G.parse([&](const std::string& p_seqname, const std::string& p_source, const std::string& p_feature,
                        const int64_t& p_start, const int64_t& p_end,
                        const std::optional<double>& p_score, const char& p_strand,
                        const std::optional<int>& p_frame, const std::vector<gtf::attribute>& p_attrs) {
                if (p_feature != "exon")
                {
                    return;
                }

                stabby::interval ivl = {p_start, p_end};

                std::string gene_id;
                std::string gene_name;
                for (auto itr = p_attrs.begin(); itr != p_attrs.end(); ++itr)
                {
                    if (itr->first == "gene_id")
                    {
                        gene_id = std::get<std::string>(itr->second);
                    }
                    if (itr->first == "gene_name")
                    {
                        gene_name = std::get<std::string>(itr->second);
                    }
                }
                if (gene_id.size() == 0)
                {
                    return;
                }

                annot[label(p_seqname)][ivl] = label(gene_id);
                gene_names[label(gene_id)] = label(gene_name);
            });

            for (auto itr = annot.begin(); itr != annot.end(); ++itr)
            {
                const auto& chrom = itr->first;
                const auto& ivl_map = itr->second;

                std::vector<stabby::interval> raw;
                raw.reserve(ivl_map.size());
                for (auto jtr = ivl_map.begin(); jtr != ivl_map.end(); ++jtr)
                {
                    raw.push_back(jtr->first);
                }
                std::sort(raw.begin(), raw.end());
                raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
                annot_idx[chrom] = std::shared_ptr<stabby>(new stabby(raw));
            }
        }

        std::unordered_map<label,std::unordered_map<label,std::vector<std::pair<size_t,size_t>>>> res;
        {
            std::string ref_name = opts.at("<reference-fasta>").asString();
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                label chrom = label(first_word(r.first));
                BOOST_LOG_TRIVIAL(info) << "scanning " << chrom << ".";

                kmers::make(r.second, K, [&](kmer_and_pos p_xp) {
                    kmer x = p_xp.first;
                    size_t p = p_xp.second;

                    auto itr = idx.find(x);
                    if (itr == idx.end())
                    {
                        return;
                    }
                    const auto& labs = itr->second;
                    for (auto jtr = labs.begin(); jtr != labs.end(); ++jtr)
                    {
                        auto& vec = res[*jtr][chrom];
                        if (vec.size() > 0 && vec.back().second + 1 == p)
                        {
                            vec.back().second = p;
                        }
                        else
                        {
                            vec.push_back(std::make_pair(p, p));
                        }
                    }
                });
            }

            for (auto itr = res.begin(); itr != res.end(); ++itr)
            {
                auto& chrom_map = itr->second;
                for (auto jtr = chrom_map.begin(); jtr != chrom_map.end(); ++jtr)
                {
                    auto& vec = jtr->second;
                    for (auto ktr = vec.begin(); ktr != vec.end(); ++ktr)
                    {
                        ktr->second += K - 1;
                    }
                }
            }
        }

        for (auto itr = res.begin(); itr != res.end(); ++itr)
        {
            const auto& bait_id = itr->first;
            const auto& chrom_map = itr->second;

            //std::cout << "#track name=\"" << bait_id << "\""
            //    << " description=\"" << bait_labels[bait_id] << "\""
            //    << std::endl;
            for (auto jtr = chrom_map.begin(); jtr != chrom_map.end(); ++jtr)
            {
                const auto& chrom = jtr->first;
                const auto& vec = jtr->second;
                for (auto ktr = vec.begin(); ktr != vec.end(); ++ktr)
                {
                    const auto& ivl = *ktr;
                    std::cout << chrom
                        << '\t' << ivl.first
                        << '\t' << ivl.second
                        << '\t' << bait_id
                        << std::endl;
                }
            }
        }
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
