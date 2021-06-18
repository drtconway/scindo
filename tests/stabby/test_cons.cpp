#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/stabby.hpp"

#include <random>
#include <set>

namespace // anonymous
{
    std::vector<scindo::stabby::interval> trappc4_genomic = {
        {119018763,119025454},
        {119018766,119025454},
        {119018766,119024134},
        {119018766,119023672},
        {119018766,119023665},
        {119018766,119023642},
        {119018766,119023421},
        {119018766,119018970},
        {119019143,119019317},
        {119019143,119019210},
        {119019143,119019188},
        {119020150,119020253},
        {119020205,119020253},
        {119021760,119021886},
        {119021874,119021886},
        {119023321,119024134},
        {119023321,119023672},
        {119023321,119023665},
        {119023321,119023642},
        {119023321,119023421},
        {119024857,119025454}
    };

    std::vector<scindo::stabby::interval> trappc4_ranks = {
        {0,20},
        {1,20},
        {1,18},
        {1,17},
        {1,16},
        {1,15},
        {1,14},
        {1,2},
        {3,6},
        {3,5},
        {3,4},
        {7,9},
        {8,9},
        {10,12},
        {11,12},
        {13,18},
        {13,17},
        {13,16},
        {13,15},
        {13,14},
        {19,20}
    };

    std::vector<scindo::stabby::interval> trappc4_aug = {
        {1,41},
        {3,41},
        {3,37},
        {3,35},
        {3,33},
        {3,31},
        {3,29},
        {3,5},
        {7,13},
        {7,11},
        {7,9},
        {15,19},
        {17,19},
        {21,25},
        {23,25},
        {27,37},
        {27,35},
        {27,33},
        {27,31},
        {27,29},
        {39,41}
    };

    void make_intervals(const uint64_t& S, const uint64_t& N, const uint64_t& M, const uint64_t& D,
                        std::vector<scindo::stabby::interval>& p_res)
    {
        std::mt19937_64 rdev(S);
        std::uniform_int_distribution<int64_t> U(0, M);
        std::normal_distribution<> W(0, D);

        std::set<scindo::stabby::interval> res;
        while (res.size() < N)
        {
            uint32_t a = U(rdev);
            int32_t b = W(rdev);
            uint32_t ab = (b < 0 ? -b : b);
            res.insert(std::make_pair(a, a+ab));
        }
        p_res.insert(p_res.end(), res.begin(), res.end());
    }

#if 0
    void verify_schmidt(const std::vector<scindo::stabby::interval>& p_intervals,
                        const scindo::schmidt& X)
    {
        for (auto itr = p_intervals.begin(); itr != p_intervals.end(); ++itr)
        {
            REQUIRE( X.parent.contains(*itr) );
        }

        for (auto itr = X.parent.begin(); itr != X.parent.end(); ++itr)
        {
            const auto& a = itr->first;
            const auto& b = itr->second;

            REQUIRE( a.first <= a.second );
            REQUIRE( b.first <= b.second );

            if (b == scindo::schmidt::root())
            {
                continue;
            }

            REQUIRE( a.first > b.first );
            REQUIRE( a.second <= b.second );
        }

        for (auto itr = X.last.begin(); itr != X.last.end(); ++itr)
        {
            const auto& p = itr->first;
            scindo::schmidt::interval a = itr->second;
            while (X.left.contains(a))
            {
                a = X.left.at(a);
                REQUIRE( X.parent.at(a) == p );
            }
        }
    }

    void check_answers(const scindo::stabby& X, const scindo::schmidt& Y)
    {
        std::unordered_map<scindo::schmidt::position,std::set<scindo::schmidt::interval>> answers;
        for (auto itr = X.dense.begin(); itr != X.dense.end(); ++itr)
        {
            scindo::schmidt::interval a = *itr;
            //std::cerr << '[' << a.first << ',' << a.second << ']' << std::endl;
            for (auto q = a.first; q <= a.second; ++q)
            {
                answers[q].insert(a);
            }
        }

        std::vector<scindo::schmidt::interval> r;
        for (scindo::schmidt::position q = 0; q < Y.Q; ++q)
        {
            r.clear();
            Y.stab(q, r);

            if (0)
            {
                std::cerr << "q = " << q << std::endl;
                std::cerr << "A:";
                for (auto itr = answers[q].begin(); itr != answers[q].end(); ++itr)
                {
                    std::cerr << ' ' << '[' << itr->first << ',' << itr->second << ']';
                }
                std::cerr << std::endl;
                std::cerr << "S:";
                for (auto itr = r.begin(); itr != r.end(); ++itr)
                {
                    std::cerr << ' ' << '[' << itr->first << ',' << itr->second << ']';
                }
                std::cerr << std::endl;
            }

            REQUIRE( r.size() == answers[q].size() );
            
            auto r_itr = r.begin();
            auto s_itr = answers[q].begin();
            while (r_itr != r.end() && s_itr != answers[q].end())
            {
                REQUIRE( *r_itr == *s_itr );
                ++r_itr;
                ++s_itr;
            }

            REQUIRE( r_itr == r.end() );
            REQUIRE( s_itr == answers[q].end() );
        }
    }
#endif
}
// namespace anonymous

TEST_CASE("test construction TRAPPC4", "[stabby]") {
    scindo::stabby X(trappc4_genomic);

    std::set<uint32_t> xs;
    for (size_t i = 0; i < trappc4_genomic.size(); ++i)
    {
        xs.insert(trappc4_genomic[i].first);
        xs.insert(trappc4_genomic[i].second);
    }
    std::vector<uint32_t> ys(xs.begin(), xs.end());
    for (uint32_t i = 0; i < ys.size(); ++i)
    {
        REQUIRE( X.domain->rank(ys[i]) == i );
        REQUIRE( X.domain->select(i) == ys[i] );
    }

    REQUIRE( trappc4_genomic.size() == trappc4_ranks.size() );
    for (size_t i = 0; i < trappc4_genomic.size(); ++i)
    {
        const auto& g = trappc4_genomic[i];
        const auto& r = trappc4_ranks[i];
        std::cerr << i
                  << '\t' << '[' << g.first << ',' << g.second << ']'
                  << '\t' << '[' << r.first << ',' << r.second << ']'
                  << std::endl;
        REQUIRE( X.domain->rank(g.first) == r.first );
        REQUIRE( X.domain->rank(g.second) == r.second );
        REQUIRE( X.domain->select(r.first) == g.first );
        REQUIRE( X.domain->select(r.second) == g.second );
    }

    REQUIRE( trappc4_ranks.size() == trappc4_aug.size() );
    for (size_t i = 0; i < trappc4_ranks.size(); ++i)
    {
        const auto& r = trappc4_ranks[i];
        const auto& a = trappc4_aug[i];
        std::cerr << i
                  << '\t' << '[' << r.first << ',' << r.second << ']'
                  << '\t' << '[' << a.first << ',' << a.second << ']'
                  << std::endl;
        REQUIRE( X.gapped->rank(a.first) == r.first );
        REQUIRE( X.gapped->rank(a.second) == r.second );
        REQUIRE( X.gapped->select(r.first) == a.first );
        REQUIRE( X.gapped->select(r.second) == a.second );
    }

    for (size_t i = 0; i < trappc4_genomic.size(); ++i)
    {
        const auto& g = trappc4_genomic[i];
        const auto& a = trappc4_aug[i];
        std::cerr << i
                  << '\t' << '[' << g.first << ',' << g.second << ']'
                  << '\t' << '[' << a.first << ',' << a.second << ']'
                  << std::endl;
        REQUIRE( X.sparse_to_dense(g.first) == a.first );
        REQUIRE( X.sparse_to_dense(g.second) == a.second );
        REQUIRE( X.dense_to_sparse(a.first) == g.first );
        REQUIRE( X.dense_to_sparse(a.second) == g.second );
    }

    uint32_t min_pos = trappc4_genomic[0].first;
    uint32_t max_pos = trappc4_genomic[0].second;
    for (auto itr = trappc4_genomic.begin(); itr != trappc4_genomic.end(); ++itr)
    {
        min_pos = std::min(min_pos, itr->first);
        max_pos = std::max(max_pos, itr->second);
    }

    std::unordered_map<uint32_t,std::set<scindo::stabby::interval>> answers;
    for (size_t q = 0; q <= X.Q; ++q)
    {
        for (size_t j = 0; j < trappc4_aug.size(); ++j)
        {
            const auto a = trappc4_aug[j];
            if (a.first <= q && q <= a.second)
            {
                answers[q].insert(trappc4_genomic[j]);
            }
        }
    }

    for (uint32_t p = min_pos - 1; p <= max_pos + 1; ++p)
    {
        auto q = X.sparse_to_dense(p);

        std::vector<scindo::stabby::interval> ans;
        X.find(p, ans);
        if (answers.contains(q))
        {
            REQUIRE( ans.size() == answers.at(q).size() );
            for (auto itr = ans.begin(); itr != ans.end(); ++itr)
            {
                REQUIRE( answers.at(q).contains(*itr) );
            }
        }
        else
        {
            REQUIRE( ans.size() == 0 );
        }
    }
}

TEST_CASE("test construction", "[stabby]") {
    const uint64_t N = 30;
    std::vector<scindo::stabby::interval> intervals;

    make_intervals(19, N, 100, 10, intervals);

    scindo::stabby X(intervals);
}
