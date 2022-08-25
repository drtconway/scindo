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
#include "scindo/nmf.hpp"

using namespace scindo;

namespace // anonymous
{
    const char usage[] =
R"(crispin - detecting differences between GVCFs

    Usage:
      crispin [options] <VCF>...

    Options:
      -h, --help                        Show this help message
)";

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
