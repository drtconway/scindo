#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/stabby.hpp"

#include <fstream>
#include <random>
#include <set>
#include <nlohmann/json.hpp>

namespace // anonymous
{
    nlohmann::json data(const std::string& p_filename)
    {
        std::ifstream in(p_filename);
        nlohmann::json dat;
        in >> dat;
        return dat;
    }
}
// namespace anonymous

TEST_CASE("test regression 1", "[stabby]") {
    nlohmann::json blob = data("stabby/case001.json");
    std::vector<std::pair<uint32_t,uint32_t>> ivls = blob["intervals"];
    scindo::stabby X(ivls);
    std::vector<scindo::stabby::interval> ans;
    uint32_t q = blob["query"];
    X.find(q, ans);
    for (auto itr = ans.begin(); itr != ans.end(); ++itr)
    {
        REQUIRE( itr->first <= q );
        REQUIRE( q <= itr->second );
    }
}
