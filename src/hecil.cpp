#include <iostream>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include "scindo/files.hpp"
#include "scindo/fasta.hpp"
#include "scindo/fastq.hpp"
#include "scindo/kmers.hpp"
#include "scindo/kmer_set.hpp"
#include "scindo/profile.hpp"
#include "scindo/summerizer.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(hecil - in silico depletion for RNASeq

    Usage:
      hecil [options] <reference-fasta> <fastq1> <fastq2> <keep1> <keep2> <toss1> <toss2>

    Options:
      -h --help                         Show this screen
      --sample-rate PROB                probability to sub-sample reference kmers [default: 1.0].
      --read-buffer SIZE                buffer size for reading reads [default: 1024].
      --write-buffer SIZE               buffer size for writing reads [default: 1024].
)";

    std::string first_word(const std::string& p_str)
    {
        auto n = p_str.find(' ');
        if (n == std::string::npos)
        {
            return p_str;
        }
        return p_str.substr(0, n);
    }

    template <typename X>
    void with(size_t p_readBufferSize, std::istream& p_fq1, std::istream& p_fq2, X p_acceptor)
    {
        fastq_reader fq1(p_fq1, p_readBufferSize);
        fastq_reader fq2(p_fq2, p_readBufferSize);

        bool stop = false;
        while (fq1.more() && fq2.more())
        {
            p_acceptor(*fq1, *fq2, stop);
            if (stop)
            {
                return;
            }
            ++fq1;
            ++fq2;
        }
        if (fq1.more() != fq2.more())
        {
            throw std::runtime_error("fastq files had unequal length.");
        }
    }

    bool contains(const std::vector<kmer>& X, const kmer& x)
    {
        auto itr = std::lower_bound(X.begin(), X.end(), x);
        return itr != X.end() && *itr == x;
    }

    uint64_t roundup(uint64_t v)
    {
        v--;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v |= v >> 32;
        return v + 1;
    }

    template <typename X>
    void with_contains(std::vector<kmer>& Xs, const std::vector<kmer>& xs, X p_acceptor)
    {
        const uint32_t P = roundup(Xs.size()) >> 1;
        Xs.resize(P << 1, uint64_t(-1LL));
        for (size_t i = 0; i < xs.size(); ++i)
        {
            const kmer x = xs[i];
            uint32_t j = 0;
            uint32_t k = P;
            while (k > 0)
            {
                uint32_t r = j | k;
                //std::cerr << "probing " << r << " = " << Xs[r] << " for " << x << std::endl;
                if (x >= Xs[r])
                {
                    j = r;
                }
                k >>= 1;
            }
            if (Xs[j] == x)
            {
                p_acceptor(x);
            }
        }
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

        uint64_t seed = 0;
        double prob = 1.0;
        {
            double d = std::stod(opts.at("--sample-rate").asString());
            seed = uint64_t(d);
            prob = d - seed;
            if (prob == 0.0)
            {
                prob = 1.0;
            }
        }

        std::unordered_map<kmer,std::vector<std::string>> X;
        std::vector<kmer> Y;
        {
            profile<true> P("indexing");
            std::mt19937_64 gen(seed);
            std::uniform_real_distribution<> U;
            std::string ref_name = opts.at("<reference-fasta>").asString();
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                std::string chrom = r.first;
                kmers::make(r.second, K, [&](kmer p_x, kmer p_xb) {
                    if (U(gen) <= prob)
                    {
                        X[p_x].push_back(chrom);
                        X[p_xb].push_back(chrom);
                        Y.push_back(p_x);
                        Y.push_back(p_xb);
                    }
                });
            }
            std::sort(Y.begin(), Y.end());
            Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
        }
        kmer_set Z(Y);

        {
            profile<true> P("processing");

            const size_t readBufferSize = opts.at("--read-buffer").asLong();

            std::string fq1_name = opts.at("<fastq1>").asString();
            input_file_holder_ptr fq1 = files::in(fq1_name);
            std::string fq2_name = opts.at("<fastq2>").asString();
            input_file_holder_ptr fq2 = files::in(fq2_name);

            const size_t writeBufferSize = opts.at("--write-buffer").asLong();

            std::string keep1_name = opts.at("<keep1>").asString();
            output_file_holder_ptr keep1 = files::out(keep1_name);
            fastq_writer keeper1(**keep1, writeBufferSize);
            std::string keep2_name = opts.at("<keep2>").asString();
            output_file_holder_ptr keep2 = files::out(keep2_name);
            fastq_writer keeper2(**keep2, writeBufferSize);

            std::string toss1_name = opts.at("<toss1>").asString();
            output_file_holder_ptr toss1 = files::out(toss1_name);
            fastq_writer tosser1(**toss1, writeBufferSize);
            std::string toss2_name = opts.at("<toss2>").asString();
            output_file_holder_ptr toss2 = files::out(toss2_name);
            fastq_writer tosser2(**toss2, writeBufferSize);

            size_t rn = 0;
            size_t rn_d = 0;
            size_t kept = 0;
            size_t tossed = 0;
            std::vector<kmer> xs;
            const kmer zmax = Y.back();
            std::map<size_t,size_t> hitHist;
            auto start_time = std::chrono::high_resolution_clock::now();
            with(readBufferSize, **fq1, **fq2, [&](const fastq_read& r1, const fastq_read& r2, bool& stop) {
                //profile<true> P("read handling");

                if (rn > 1ULL << 20)
                {
                    stop = true;
                    return;
                }

                rn += 1;
                rn_d += 1;
                if ((rn & ((1ULL << 18) - 1)) == 0)
                {
                    auto end_time = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> delta = (end_time - start_time);
                    double rps = rn_d / delta.count();
                    BOOST_LOG_TRIVIAL(info) << "processed records: " << rn << '\t' << rps << " reads/second";
                    rn_d = 0;
                    start_time = std::chrono::high_resolution_clock::now();
                }
                xs.clear();
                kmers::make(std::get<1>(r1), K, [&](kmer p_x) {
                    xs.push_back(p_x);
                });
                kmers::make(std::get<1>(r2), K, [&](kmer p_x) {
                    xs.push_back(p_x);
                });
                std::sort(xs.begin(), xs.end());

                size_t hits = 0;
                Z.with_contains(xs, [&](const kmer p_x) {
                    hits += 1;
                });

                if (hits > 0)
                {
                    //profile<true> P("writing tosses");
                    tosser1.write(r1);
                    tosser2.write(r2);
                    tossed += 1;
                }
                else
                {
                    //profile<true> P("writing keeps");
                    keeper1.write(r1);
                    keeper2.write(r2);
                    kept += 1;
                }
            });

            if (0)
            {
                BOOST_LOG_TRIVIAL(info) << "kmer hits histogram:";
                for (auto itr = hitHist.begin(); itr != hitHist.end(); ++itr)
                {
                    BOOST_LOG_TRIVIAL(info) << '\t' << itr->first << '\t' << itr->second;
                    
                }
            }
            BOOST_LOG_TRIVIAL(info) << "finished";
            BOOST_LOG_TRIVIAL(info) << "reads kept: " << kept;
            BOOST_LOG_TRIVIAL(info) << "reads tossed: " << tossed;
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
