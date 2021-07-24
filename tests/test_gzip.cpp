#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/gzip.hpp"

namespace // anonymous
{
}
// namespace anonymous

TEST_CASE("test reading blocks", "[gzip]") {
    std::string str;
    scindo::gzip::with_blocks("data/tale.txt.gz", [&](auto beg, auto end) {
        str.insert(str.end(), beg, end);
    });
    std::cout << str << std::endl;
    REQUIRE( str.size() == 771503 );
}

TEST_CASE("test reading lines", "[gzip]") {
    std::vector<std::string> lines;
    scindo::gzip::with_lines("data/tale.txt.gz", [&](auto beg, auto end) {
        lines.push_back(std::string(beg, end));
    });
    REQUIRE( lines.size() == 15900 );
}
