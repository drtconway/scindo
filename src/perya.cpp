#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <execution>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include "scindo/fasta.hpp"
#include "scindo/fastq.hpp"
#include "scindo/files.hpp"
#include "scindo/stabby.hpp"
#include "scindo/murmur3.hpp"
#include "scindo/kmers.hpp"
#include "scindo/profile.hpp"
#include "scindo/tsv.hpp"
#include "scindo/summerizer.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(perya - exome panel detection

    Usage:
      perya -X [options] <index-name> <reference-fasta> <bed-file>...
      perya [options] <index-name> <fastq1> <fastq2>

    Options:
      -h, --help                        Show this screen
      -k NUM, --kmer-size NUM           K-mer size to use [default: 25]
      -S FRAC, --sample-kmers FRAC      Sub-sample given fraction of k-mers
)";

    struct succinct_kmer_set
    {
        const sdsl::sd_vector<> X;
        const kmer maxx;

        succinct_kmer_set(const std::vector<kmer>& p_X)
            : X(p_X.begin(), p_X.end()), maxx(p_X.size() > 0 ? p_X.back() : 0)
        {
        }

        bool contains(kmer x) const
        {
            if (x > maxx)
            {
                return false;
            }
            return X[x];
        }

    };

    std::string first_word(const std::string& p_str)
    {
        auto n = p_str.find(' ');
        if (n == std::string::npos)
        {
            return p_str;
        }
        return p_str.substr(0, n);
    }

    void remove_dups(const std::vector<kmer>& p_lhs, std::vector<kmer>& p_rhs)
    {
        auto lhsItr = p_lhs.begin();
        auto rhsItr = p_rhs.begin();
        auto writeItr = p_rhs.begin();
        while (lhsItr != p_lhs.end() && rhsItr != p_rhs.end())
        {
            if (*lhsItr < *rhsItr)
            {
                ++lhsItr;
                continue;
            }
            else if (*rhsItr < *lhsItr)
            {
                if (writeItr != rhsItr)
                {
                    *writeItr = *rhsItr;
                }
                ++writeItr;
                ++rhsItr;
            }
            else // *rhsItr == *lhsItr
            {
                ++rhsItr;
                ++lhsItr;
            }
        }
        while (rhsItr != p_rhs.end())
        {
            if (writeItr != rhsItr)
            {
                *writeItr = *rhsItr;
            }
            ++writeItr;
            ++rhsItr;
        }
        p_rhs.erase(writeItr, p_rhs.end());
    }

    bool sample(const double& prob, const kmer& p_x)
    {
        uint64_t u = murmur3(17).update(p_x)();
        double p = 0.0;
        double q = 0.5;
        while (u > 0)
        {
            if ((u & 1) == 1)
            {
                p += q;
            }
            u >>= 1;
            q /= 2;
        }
        return p < prob;
    }

    bool access_and_rank(const std::vector<kmer>& p_X, const kmer& p_x, size_t& p_rnk)
    {
        auto itr = std::lower_bound(p_X.begin(), p_X.end(), p_x);
        if (itr != p_X.end() && *itr == p_x)
        {
            p_rnk = itr - p_X.begin();
            return true;
        }
        return false;
    }

    void write_index(std::ostream& p_out,
                     const uint64_t& p_K,
                     const std::vector<kmer>& p_kmers,
                     const std::vector<uint16_t>& p_masks,
                     const std::vector<std::string>& p_bed_names)
    {
        p_out.write(reinterpret_cast<const char*>(&p_K), sizeof(p_K));
        const uint64_t N = p_kmers.size();
        p_out.write(reinterpret_cast<const char*>(&N), sizeof(N));
        p_out.write(reinterpret_cast<const char*>(p_kmers.data()), sizeof(kmer)*p_kmers.size());
        p_out.write(reinterpret_cast<const char*>(p_masks.data()), sizeof(uint16_t)*p_masks.size());
        const uint64_t M = p_bed_names.size();
        p_out.write(reinterpret_cast<const char*>(&M), sizeof(M));
        for (auto itr = p_bed_names.begin(); itr != p_bed_names.end(); ++itr)
        {
            const std::string& name = *itr;
            const uint64_t L = name.size();
            p_out.write(reinterpret_cast<const char*>(&L), sizeof(L));
            p_out.write(name.data(), name.size());
        }
    }

    void read_index(std::istream& p_in,
                    uint64_t& p_K,
                    std::vector<kmer>& p_kmers,
                    std::vector<uint16_t>& p_masks,
                    std::vector<std::string>& p_bed_names)
    {
        p_in.read(reinterpret_cast<char*>(&p_K), sizeof(p_K));
        uint64_t N = 0;
        p_in.read(reinterpret_cast<char*>(&N), sizeof(N));
        p_kmers.resize(N);
        p_in.read(reinterpret_cast<char*>(p_kmers.data()), sizeof(kmer)*p_kmers.size());
        p_masks.resize(N);
        p_in.read(reinterpret_cast<char*>(p_masks.data()), sizeof(uint16_t)*p_masks.size());
        uint64_t M = 0;
        p_in.read(reinterpret_cast<char*>(&M), sizeof(M));
        p_bed_names.resize(M);
        for (auto itr = p_bed_names.begin(); itr != p_bed_names.end(); ++itr)
        {
            std::string& name = *itr;
            uint64_t L = 0;
            p_in.read(reinterpret_cast<char*>(&L), sizeof(L));
            name.resize(L);
            p_in.read(name.data(), name.size());
        }
    }

    template <typename X>
    void collect(const std::vector<uint64_t>& xs, X p_acceptor)
    {
        static_assert(std::is_convertible<X, std::function<void(uint64_t,size_t)>>::value);
        uint64_t x = 0;
        uint64_t cnt = 0;
        for(auto itr = xs.begin(); itr != xs.end(); ++itr)
        {
            if (*itr != x)
            {
                if (cnt > 0)
                {
                    p_acceptor(x, cnt);
                }
                x = *itr;
                cnt = 0;
            }
            cnt += 1;
        }
        if (cnt > 0)
        {
            p_acceptor(x, cnt);
        }
    }

    void jaccard(const std::vector<kmer>& p_allKmers, const std::vector<uint16_t>& p_allMasks,
                 const std::vector<kmer>& p_xs, std::vector<std::pair<size_t,size_t>>& p_res)
    {
        size_t i = 0;
        size_t j = 0;

        while (i < p_allKmers.size() && j < p_xs.size())
        {
            kmer y = p_allKmers[i];
            kmer x = p_xs[j];
            uint64_t msk = p_allMasks[i];
            if (y < x)
            {
                for (size_t p = 0; msk > 0 && p < p_res.size(); ++p, msk >>= 1)
                {
                    if (msk & 1)
                    {
                        p_res[p].second += 1;
                    }
                }
                ++i;
                continue;
            }
            if (x < y)
            {
                for (size_t p = 0; p < p_res.size(); ++p)
                {
                    p_res[p].second += 1;
                }
                ++j;
                continue;
            }
            // x == y
            for (size_t p = 0; msk > 0 && p < p_res.size(); ++p, msk >>= 1)
            {
                if (msk & 1)
                {
                    p_res[p].first += 1;
                }
                p_res[p].second += 1;
            }
            ++i;
            ++j;
        }
        while (i < p_allKmers.size())
        {
            uint64_t msk = p_allMasks[i];
            for (size_t p = 0; msk > 0 && p < p_res.size(); ++p, msk >>= 1)
            {
                if (msk & 1)
                {
                    p_res[p].second += 1;
                }
            }
            ++i;
        }
        if (j < p_xs.size())
        {
            size_t n = p_xs.size() - j;
            for (size_t p = 0; p < p_res.size(); ++p)
            {
                p_res[p].second += n;
            }
        }
    }

    int indexMain(const std::map<std::string, docopt::value>& p_opts)
    {
        const size_t K = p_opts.at("--kmer-size").asLong();
        const std::string ref_name = p_opts.at("<reference-fasta>").asString();
        const std::vector<std::string> beds = p_opts.at("<bed-file>").asStringList();
        double frac = 1.0;
        if (p_opts.at("--sample-kmers"))
        {
            frac = boost::lexical_cast<double>(p_opts.at("--sample-kmers").asString());
        }

        using locus = std::pair<uint32_t,uint32_t>;
        using bed_nums = std::vector<size_t>;
        std::unordered_map<std::string, std::map<locus,bed_nums>> loci;

        for (size_t i = 0; i < beds.size(); ++i)
        {
            input_file_holder_ptr inp = files::in(beds[i]);
            tsv::with(**inp, [&](const std::vector<std::string>& p_row) {
                const std::string& chrom = p_row[0];
                uint32_t st = boost::lexical_cast<uint32_t>(p_row[1]);
                uint32_t en = boost::lexical_cast<uint32_t>(p_row[2]);
                loci[chrom][locus(st,en)].push_back(i);
            });
        }

        std::vector<kmer> allKmers;
        {
            std::vector<kmer> buffer;
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                const std::string chrom = first_word(r.first);
                const auto& seq = r.second;
                if (!loci.contains(chrom))
                {
                    continue;
                }
                BOOST_LOG_TRIVIAL(info) << "scanning (1st pass) " << chrom;
                const auto itr = loci.find(chrom);
                const auto& ivls = itr->second;

                buffer.clear();
                for (auto jtr = ivls.begin(); jtr != ivls.end(); ++jtr)
                {
                    const auto& ivl = jtr->first;
                    const auto& bns = jtr->second;
                    std::pair<std::string::const_iterator,std::string::const_iterator>
                        bait{seq.begin()+ivl.first,
                             seq.begin()+ivl.second};
                    //std::cout << std::string(bait.first, bait.second) << std::endl;
                    kmers::make(bait, K, [&](kmer p_x, kmer p_xb) {
                        kmer xc = std::min(p_x, p_xb);
                        //std::cout << kmers::render(K, xc) << '\t' << sample(frac, xc) << std::endl;
                        if (sample(frac, xc))
                        {
                            buffer.push_back(p_x);
                            buffer.push_back(p_xb);
                        }
                    });
                }
                std::sort(buffer.begin(), buffer.end());
                buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());
                //BOOST_LOG_TRIVIAL(info) << "distinct k-mers " << buffer.size();
                remove_dups(allKmers, buffer);
                //BOOST_LOG_TRIVIAL(info) << "new k-mers " << buffer.size();
                allKmers.insert(allKmers.end(), buffer.begin(), buffer.end());
                std::sort(allKmers.begin(), allKmers.end());
            }
        }

        BOOST_LOG_TRIVIAL(info) << "total distinct k-mers " << allKmers.size();
        std::vector<uint16_t> allMasks(allKmers.size());
        {
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                const std::string chrom = first_word(r.first);
                const auto& seq = r.second;
                if (!loci.contains(chrom))
                {
                    continue;
                }
                BOOST_LOG_TRIVIAL(info) << "scanning (2nd pass) " << chrom;
                const auto itr = loci.find(chrom);
                const auto& ivls = itr->second;

                for (auto jtr = ivls.begin(); jtr != ivls.end(); ++jtr)
                {
                    const auto& ivl = jtr->first;
                    const auto& bns = jtr->second;
                    uint64_t msk = 0;
                    for (auto ktr = bns.begin(); ktr != bns.end(); ++ktr)
                    {
                        msk |= (1ULL << (*ktr));
                    }
                    std::pair<std::string::const_iterator,std::string::const_iterator>
                        bait{seq.begin()+ivl.first,
                             seq.begin()+ivl.second};
                    kmers::make(bait, K, [&](kmer p_x, kmer p_xb) {
                        size_t rnk = 0;
                        if (access_and_rank(allKmers, p_x, rnk))
                        {
                            allMasks[rnk] |= msk;
                        }
                        if (access_and_rank(allKmers, p_xb, rnk))
                        {
                            allMasks[rnk] |= msk;
                        }
                    });
                }
            }
            std::map<size_t,size_t> hist;
            for (auto itr = allMasks.begin(); itr != allMasks.end(); ++itr)
            {
                hist[std::popcount(*itr)] += 1;
            }
            for (auto itr = hist.begin(); itr != hist.end(); ++itr)
            {
                BOOST_LOG_TRIVIAL(info) << "mask popcount " << itr->first << "\t" << itr->second;

            }
        }

        output_file_holder_ptr out_ptr = files::out(p_opts.at("<index-name>").asString());
        write_index(**out_ptr, K, allKmers, allMasks, beds);

        return 0;
    }

    int main0(int argc, const char* argv[])
    {
        boost::log::add_console_log(std::cerr, boost::log::keywords::format = "[%TimeStamp%] [%ThreadID%] [%Severity%] %Message%");
        boost::log::add_common_attributes();

        std::map<std::string, docopt::value>
            opts = docopt::docopt(usage, { argv + 1, argv + argc }, true, "perya 0.1");

        for (auto itr = opts.begin(); itr != opts.end(); ++itr)
        {
            BOOST_LOG_TRIVIAL(info) << itr->first << '\t' << itr->second;
        }

        if (opts.at("-X").asBool())
        {
            return indexMain(opts);
        }

        double frac = 1.0;
        if (opts.at("--sample-kmers"))
        {
            frac = boost::lexical_cast<double>(opts.at("--sample-kmers").asString());
        }

        uint64_t K = 0;
        std::vector<kmer> allKmers;
        std::vector<uint16_t> allMasks;
        std::vector<std::string> bedNames;

        {
            input_file_holder_ptr in_ptr = files::in(opts.at("<index-name>").asString());
            read_index(**in_ptr, K, allKmers, allMasks, bedNames);
        }

        std::string fq1_name = opts.at("<fastq1>").asString();
        std::string fq2_name = opts.at("<fastq2>").asString();

        std::vector<kmer> xs;

        size_t rn = 0;
        fastq_reader::with(fq1_name, fq2_name, [&](const fastq_tuple& r1,
                                                  const fastq_tuple& r2, bool& stop) {
            rn += 1;

            kmers::make(std::get<1>(r1), K, [&](kmer p_x, kmer p_xb) {
                kmer xc = std::min(p_x, p_xb);
                //std::cout << kmers::render(K, xc) << '\t' << sample(frac, xc) << std::endl;
                if (sample(frac, xc))
                {
                    xs.push_back(xc);
                }
            });
            kmers::make(std::get<1>(r2), K, [&](kmer p_x, kmer p_xb) {
                kmer xc = std::min(p_x, p_xb);
                //std::cout << kmers::render(K, xc) << '\t' << sample(frac, xc) << std::endl;
                if (sample(frac, xc))
                {
                    xs.push_back(xc);
                }
            });

            constexpr size_t N = (1ULL << 20);
            constexpr size_t M = N - 1;
            if ((rn & M) == 0)
            {
                std::sort(xs.begin(), xs.end());
                xs.erase(std::unique(xs.begin(), xs.end()), xs.end());

                std::vector<std::pair<size_t,size_t>> intersectionAndUnions;
                intersectionAndUnions.resize(bedNames.size());
                jaccard(allKmers, allMasks, xs, intersectionAndUnions);

                double t = 0;
                for (size_t i = 0; i < bedNames.size(); ++i)
                {
                    size_t jacI = intersectionAndUnions[i].first;
                    size_t jacU = intersectionAndUnions[i].second;
                    double jac = double(jacI)/double(jacU);
                    t += jac;
                }

                std::cout << "# " << rn << std::endl;
                for (size_t i = 0; i < bedNames.size(); ++i)
                {
                    size_t jacI = intersectionAndUnions[i].first;
                    size_t jacU = intersectionAndUnions[i].second;
                    double jac = double(jacI)/double(jacU);
                    std::cout << i
                        << '\t' << jacI
                        << '\t' << jacU
                        << '\t' << jac
                        << '\t' << (jac/t)
                        << '\t' << bedNames[i]
                        << std::endl;
                }
                xs.clear();
            }
        });

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
