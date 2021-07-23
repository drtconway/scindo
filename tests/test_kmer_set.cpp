#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/kmer_set.hpp"

#include <random>
#include <set>

namespace // anonymous
{
}
// namespace anonymous

TEST_CASE("test 2", "[kmer_set]") {
    const size_t K = 8;
    const size_t N = 100;

    std::mt19937_64 gen(19);
    std::uniform_int_distribution<uint64_t> U(0, (1ULL << (2*K)) - 1);

    std::vector<scindo::kmer> xs;

    {
        for (size_t i = 0; i < N; ++i)
        {
            xs.push_back(U(gen));
        }
        std::sort(xs.begin(), xs.end());
    }

    std::set<scindo::kmer> Y(xs.begin(), xs.end());

    scindo::kmer_set X(xs);

    for (size_t i = 0; i < N; ++i)
    {
        const scindo::kmer x = xs[i];
        REQUIRE( X.contains(x) == true );
    }

    for (size_t i = 0; i < 100; ++i)
    {
        const scindo::kmer x = U(gen);
        REQUIRE( X.contains(x) == Y.contains(x) );
    }
}
