#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <docopt/docopt.h>
#include <nlohmann/json.hpp>
#include <sdsl/sd_vector.hpp>
#include "scindo/files.hpp"
#include "scindo/fasta.hpp"
#include "scindo/fastq.hpp"
#include "scindo/kmers.hpp"
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
    void with(std::istream& p_fq1, std::istream& p_fq2, X p_acceptor)
    {
        fastq_reader fq1(p_fq1);
        fastq_reader fq2(p_fq2);

        while (fq1.more() && fq2.more())
        {
            p_acceptor(*fq1, *fq2);
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

    template <typename X>
    void with_contains(const std::vector<kmer>& Xs, const std::vector<kmer>& xs, X p_acceptor)
    {
        auto from = Xs.begin();
        for (auto itr = xs.begin(); itr != xs.end(); ++itr)
        {
            auto jtr = std::lower_bound(from, Xs.end(), *itr);
            if (jtr == Xs.end())
            {
                return;
            }
            if (*jtr == *itr)
            {
                p_acceptor(*itr);
            }
            from = jtr;
        }
    }

    template <typename X>
    void with_contains(const sdsl::sd_vector<>& Xs, const kmer p_max, const std::vector<kmer>& xs, X p_acceptor)
    {
        for (auto itr = xs.begin(); itr != xs.end(); ++itr)
        {
            kmer x = *itr;
            if (x > p_max)
            {
                return;
            }

            if (Xs[x])
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

        std::unordered_map<kmer,std::vector<std::string>> X;
        std::vector<kmer> Y;
        if (1)
        {
            std::string ref_name = opts.at("<reference-fasta>").asString();
            input_file_holder_ptr in_ptr = files::in(ref_name);
            for (scindo::seq::fasta_reader R(**in_ptr); R.more(); ++R)
            {
                const auto& r = *R;
                std::string chrom = r.first;
                kmers::make(r.second, K, [&](kmer p_x, kmer p_xb) {
                    X[p_x].push_back(chrom);
                    X[p_xb].push_back(chrom);
                    Y.push_back(p_x);
                    Y.push_back(p_xb);
                });
            }
            std::sort(Y.begin(), Y.end());
            Y.erase(std::unique(Y.begin(), Y.end()), Y.end());
        }
        sdsl::sd_vector<> Z(Y.begin(), Y.end());

        std::string fq1_name = opts.at("<fastq1>").asString();
        input_file_holder_ptr fq1 = files::in(fq1_name);
        std::string fq2_name = opts.at("<fastq2>").asString();
        input_file_holder_ptr fq2 = files::in(fq2_name);

        std::string keep1_name = opts.at("<keep1>").asString();
        output_file_holder_ptr keep1 = files::out(keep1_name);
        fastq_writer keeper1(**keep1);
        std::string keep2_name = opts.at("<keep2>").asString();
        output_file_holder_ptr keep2 = files::out(keep2_name);
        fastq_writer keeper2(**keep2);

        std::string toss1_name = opts.at("<toss1>").asString();
        output_file_holder_ptr toss1 = files::out(toss1_name);
        fastq_writer tosser1(**toss1);
        std::string toss2_name = opts.at("<toss2>").asString();
        output_file_holder_ptr toss2 = files::out(toss2_name);
        fastq_writer tosser2(**toss2);

        size_t rn = 0;
        size_t rn_d = 0;
        std::vector<kmer> xs;
        const kmer zmax = Y.back();
        auto start_time = std::chrono::high_resolution_clock::now();
        with(**fq1, **fq2, [&](const fastq_read& r1, const fastq_read& r2) {
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
            with_contains(Z, zmax, xs, [&](const kmer p_x) {
                hits += 1;
            });

            if (hits > 0)
            {
                tosser1.write(r1);
                tosser2.write(r2);
                //fastq_writer::write(**toss1, r1);
                //fastq_writer::write(**toss2, r2);
            }
            else
            {
                keeper1.write(r1);
                keeper1.write(r2);
                //fastq_writer::write(**keep1, r1);
                //fastq_writer::write(**keep2, r2);
            }
        });

        BOOST_LOG_TRIVIAL(info) << "finished";

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
