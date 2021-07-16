#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/annotated_kmer_set.hpp"

namespace // anonymous
{
}
// namespace anonymous

TEST_CASE("test reordering", "[annotated_kmer_set]") {
    std::vector<scindo::kmer> xs({5, 6, 1, 2, 5, 2, 3});
    std::vector<uint64_t>     ys({1, 1, 1, 1, 2, 2, 1});

    scindo::annotated_kmer_set X;
    for (size_t i = 0; i < xs.size(); ++i)
    {
        X.add(xs[i], ys[i]);
    }
    REQUIRE( X.kmer_buf.size() == X.annot_buf.size() );
    REQUIRE( X.kmer_buf.size() == 7 );
    X.prepare();
    REQUIRE( X.kmer_buf.size() == X.annot_buf.size() );
    REQUIRE( X.kmer_buf.size() == 5 );
    REQUIRE( X.annot_buf[0] == 1 );
    REQUIRE( X.annot_buf[1] == 3 );
    REQUIRE( X.annot_buf[2] == 1 );
    REQUIRE( X.annot_buf[3] == 3 );
    REQUIRE( X.annot_buf[4] == 1 );
}
