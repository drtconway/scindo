#include <iostream>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <execution>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include <sdsl/sd_vector.hpp>
#include "scindo/fasta.hpp"
#include "scindo/files.hpp"
#include "scindo/gtf.hpp"
#include "scindo/kmer_set.hpp"
#include "scindo/kmers.hpp"
#include "scindo/profile.hpp"
#include "scindo/stabby.hpp"
#include "scindo/summerizer.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(perya - in silico binning of reads

    Usage:
      perya -X [options] <genome-fasta> <annotation>

    Options:
      -h --help                         Show this screen
)";

    using sparse_domain_ptr = std::shared_ptr<sparse_domain>;

    struct item
    {
        sparse_domain_ptr ptr;
        size_t  idx;
        uint64_t cur;

        item(const sparse_domain_ptr& p_ptr)
            : ptr(p_ptr), idx(0)
        {
            if (idx < ptr->count())
            {
                cur = ptr->select(idx);
            }
        }

        bool more() const
        {
            return idx < ptr->count();
        }

        void next()
        {
            ++idx;
            if (idx < ptr->count())
            {
                cur = ptr->select(idx);
            }
        }

        bool operator<(const item& other) const
        {
            return cur > other.cur;
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

    template <typename Vec, typename Itm>
    size_t rank(const Vec& X, const Itm& x)
    {
        auto itr = std::lower_bound(X.begin(), X.end(), x);
        return itr - X.begin();
    }

    template <typename Vec, typename Itm>
    bool contains(const Vec& X, const Itm& x)
    {
        auto itr = std::lower_bound(X.begin(), X.end(), x);
        return itr != X.end() && *itr == x;
    }

    template <typename Vec, typename Itm>
    bool contains(const Vec& X, const Itm& x, size_t& rnk)
    {
        auto itr = std::lower_bound(X.begin(), X.end(), x);
        rnk = itr - X.begin();
        return itr != X.end() && *itr == x;
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

        const size_t K = 23;
        const size_t S = 64 - 2*K;
        const uint64_t M = (1ULL << S) - 1;

        std::unordered_map<std::string,std::unordered_map<std::string,std::vector<stabby::interval>>> exons;
        std::vector<std::string> toc;
        {
            gtf::gtf_file G(opts.at("<annotation>").asString());
            G.parse([&](const std::string& p_seqname, const std::string& p_source, const std::string& p_feature,
                        const int64_t& p_start, const int64_t& p_end,
                        const std::optional<double>& p_score, const char& p_strand,
                        const std::optional<int>& p_frame, const std::vector<gtf::attribute>& p_attrs) {
                if (p_feature != "exon")
                {
                    return;
                }
                // 0-based, half-open
                stabby::interval ivl = {p_start - 1, p_end};

                std::string gene_id;
                for (auto itr = p_attrs.begin(); itr != p_attrs.end(); ++itr)
                {
                    if (itr->first == "gene_id")
                    {
                        gene_id = std::get<std::string>(itr->second);
                    }
                }
                if (gene_id.size() == 0)
                {
                    return;
                }
                exons[p_seqname][gene_id].push_back(ivl);
            });
            for (auto itr = exons.begin(); itr != exons.end(); ++itr)
            {
                const auto& chrom = itr->first;
                auto& genes = itr->second;
                for (auto jtr = genes.begin(); jtr != genes.end(); ++jtr)
                {
                    const auto& gene_id = jtr->first;
                    auto& ivls = jtr->second;

                    // Sort and remove duplicates
                    //
                    std::sort(ivls.begin(), ivls.end());
                    ivls.erase(std::unique(ivls.begin(), ivls.end()), ivls.end());

                    toc.push_back(gene_id);
                }
            }
            std::sort(toc.begin(), toc.end());
        }

        std::vector<size_t> total(toc.size(), 0);

        std::unordered_map<std::string,std::shared_ptr<sparse_domain>> X;
        {
            using str_range = std::pair<std::string::const_iterator,std::string::const_iterator>;

            std::string ref_name = opts.at("<genome-fasta>").asString();
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                std::string chrom = first_word(r.first);
                if (!exons.contains(chrom))
                {
                    continue;
                }
                BOOST_LOG_TRIVIAL(info) << "scanning chromosome " << chrom;
                const auto& genes = exons.at(chrom);
                std::vector<kmer> xs;
                std::vector<kmer> ys;
                for (auto itr = genes.begin(); itr != genes.end(); ++itr)
                {
                    const auto& gene_id = itr->first;
                    const auto& ivls = itr->second;
                    size_t rnk = 0;
                    if (!contains(toc, gene_id, rnk))
                    {
                        throw std::logic_error("missing gene_id");
                    }

                    ys.clear();
                    for (auto jtr = ivls.begin(); jtr != ivls.end(); ++jtr)
                    {
                        const auto& ivl = *jtr;
                        str_range v{r.second.begin() + ivl.first, r.second.begin() + ivl.second};
                        kmers::make(v, K, [&](kmer x, kmer xb) {
                            ys.push_back((x << S) | rnk);
                            ys.push_back((xb << S) | rnk);
                        });
                    }
                    std::sort(std::execution::unseq, ys.begin(), ys.end());
                    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
                    xs.insert(xs.end(), ys.begin(), ys.end());
                    total[rnk] = ys.size();
                }
                std::sort(std::execution::unseq, xs.begin(), xs.end());
                xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
                X[chrom] = std::shared_ptr<sparse_domain>(new sparse_domain(xs.begin(), xs.end()));
                size_t z = xs.size();
                BOOST_LOG_TRIVIAL(info) << "k-mer count " << z;
                size_t Z = size_in_bytes(X.at(chrom)->array)
                            + size_in_bytes(X.at(chrom)->array_rank)
                            + size_in_bytes(X.at(chrom)->array_select);
                BOOST_LOG_TRIVIAL(info) << "succinct size " << Z;
                BOOST_LOG_TRIVIAL(info) << "bytes/kmer " << (double(Z)/double(z));
            }
        }

        std::vector<size_t> unique(toc.size(), 0);
        size_t unique_count = 0;
        {
            std::priority_queue<item, std::vector<item>> items;
            for (auto itr = X.begin(); itr != X.end(); ++itr)
            {
                items.push(item(itr->second));
            }

            kmer x = 0;
            std::vector<size_t> gs;
            while (items.size())
            {
                item w = items.top();
                items.pop();
                kmer y = w.cur >> S;
                if (y != x)
                {
                    //std::cerr << "merge progress: " << kmers::render(K, x) << '\t' << nlohmann::json(gs) << std::endl;
                    if (gs.size() == 1)
                    {
                        unique[gs.back()] += 1;
                        unique_count += 1;
                    }
                    gs.clear();
                    x = y;
                }
                gs.push_back(w.cur & M);
                w.next();
                if (w.more())
                {
                    items.push(w);
                }
            }
        }

        {
            for (size_t i = 0; i < toc.size(); ++i)
            {
                std::cout << i
                    << '\t' << toc[i]
                    << '\t' << total[i]
                    << '\t' << unique[i]
                    << '\t' << (double(unique[i])/double(total[i]))
                    << std::endl;
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
