#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/gzip.hpp"

namespace // anonymous
{
}
// namespace anonymous

TEST_CASE("test reading blocks", "[gzip]") {
    std::string str;
    scindo::gzip::with_blocks("data/tale.txt.gz", [&](std::string::const_iterator beg, std::string::const_iterator end) {
        str.insert(str.end(), beg, end);
    });
    std::cout << str << std::endl;
    REQUIRE( str.size() == 771503 );
}

TEST_CASE("test reading lines", "[gzip]") {
    std::vector<std::string> lines;
    scindo::gzip::with_lines("data/tale.txt.gz", [&](std::string::const_iterator beg, std::string::const_iterator end) {
        lines.push_back(std::string(beg, end));
    });
    REQUIRE( lines.size() == 15900 );
}

TEST_CASE("reader lines", "[gzip]") {
    std::vector<std::string> lines;
    scindo::gzip::with_lines("data/tale.txt.gz", [&](const scindo::gzip::reader::span& p_span) {
        lines.push_back(std::string(p_span.first, p_span.second));
    });
    REQUIRE( lines.size() == 15900 );
}

TEST_CASE("lines iterator", "[gzip]") {
    std::vector<std::string> ls;
    for (scindo::gzip::lines ltr("data/tale.txt.gz"); ltr.more(); ++ltr)
    {
        const auto& s = *ltr;
        ls.push_back(std::string(s.first, s.second));
    }
    REQUIRE( ls.size() == 15900 );
}

TEST_CASE("lines paired", "[gzip]") {
    size_t n = 0;
    scindo::gzip::with_lines("data/SRR5364317-few_1.fastq.gz", "data/SRR5364317-few_2.fastq.gz", [&](const scindo::gzip::reader::span& p_lhs, const scindo::gzip::reader::span& p_rhs, bool& stop) {
        n += 1;
    });
    REQUIRE( n == 20000 );
}
