#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/bam.hpp"

#include <map>
#include <set>
#include <sstream>

namespace // anonymous
{
    using match = std::function<void(uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases)>;
    using ins = std::function<void(uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases)>;
    using del = std::function<void(uint32_t p_ref, uint32_t p_seq, uint32_t p_len)>;
    using skip = std::function<void(uint32_t p_ref, uint32_t p_seq, uint32_t p_len)>;
    using clip = std::function<void(uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases)>;

    std::map<std::string,std::map<uint32_t,uint32_t>> ans;

    void make_ans()
    {
        std::ifstream in("data/tiny-bam.pileup");
        std::string chr;
        uint32_t pos;
        std::string ref;
        uint32_t cov;
        std::string bas;
        std::string qal;

        while (in >> chr >> pos >> ref >> cov >> bas >> qal)
        {
            ans[chr][pos] = cov;
        }
    }

    template <typename XM, typename XI, typename XD, typename XS, typename XC>
    void with(uint32_t p_pos, const scindo::bam::seq& p_seq, const scindo::bam::cigar& p_cig,
              XM p_match, XI p_ins, XD p_del, XS p_skip, XC p_clip)
    {
        static_assert(std::is_convertible<XM,match>::value);
        static_assert(std::is_convertible<XI,ins>::value);
        static_assert(std::is_convertible<XD,del>::value);
        static_assert(std::is_convertible<XS,skip>::value);
        static_assert(std::is_convertible<XC,clip>::value);

        uint32_t r = p_pos;
        uint32_t q = 0;
        for (size_t i = 0; i < p_cig.size(); ++i)
        {
            char op = p_cig.op(i);
            uint32_t len = p_cig.len(i);

            switch (op)
            {
                case 'M':
                case '=':
                case 'X':
                {
                    p_match(r, q, p_seq.range(q, q+len));
                    r += len;
                    q += len;
                    break;
                }
                case 'I':
                {
                    p_ins(r, q, p_seq.range(q, q+len));
                    q += len;
                    break;
                }
                case 'D':
                {
                    p_del(r, q, len);
                    r += len;
                    break;
                }
                case 'N':
                {
                    p_skip(r, q, len);
                    r += len;
                    break;
                }
                case 'S':
                {
                    p_clip(r, q, p_seq.range(q, q+len));
                    q += len;
                    break;
                }
            }
        }
    }
}
// namespace anonymous

TEST_CASE("test bam parsing", "[bam]") {
    using scindo::bam::flag;

    make_ans();

    scindo::bam::bam_file_reader V("data/tiny-bam.bam");

    std::map<std::string,std::map<uint32_t,uint32_t>> cov;
    while (V.next())
    {
        flag flg(V.flag());
        if (flg.is<flag::unmapped>() || flg.is<flag::duplicate>() || !flg.is<flag::proper_pair>())
        {
            continue;
        }

        const std::string& chrom = V.chrom();
        int64_t pos = V.pos1();
        scindo::bam::seq s = V.seq();
        scindo::bam::cigar c = V.cigar();
        with(pos, s, c,
            // match
            //
            [&](uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases) {
                for (uint32_t i = 0; i < p_bases.size(); ++i)
                {
                    cov[chrom][p_ref+i] += 1;
                }
            },
            // ins
            //
            [&](uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases) {
                // nothing
            },
            // del
            //
            [&](uint32_t p_ref, uint32_t p_seq, uint32_t p_len) {
                for (uint32_t i = 0; i < p_len; ++i)
                {
                    cov[chrom][p_ref+i] += 1;
                }
            },
            // skip
            //
            [&](uint32_t p_ref, uint32_t p_seq, uint32_t p_len) {
                // nothing
            },
            // clip
            //
            [&](uint32_t p_ref, uint32_t p_seq, const scindo::bam::seq& p_bases) {
                // nothing
            });
    }

    REQUIRE( cov.size() == ans.size() );
    for (auto itr = cov.begin(); itr != cov.end(); ++itr)
    {
        REQUIRE( ans.contains(itr->first) );
        const auto& xs = ans.at(itr->first);
        const auto& ys = itr->second;
        if (ys.size() != xs.size())
        {
            std::set<uint32_t> xSet;
            for (auto jtr = xs.begin(); jtr != xs.end(); ++jtr)
            {
                xSet.insert(jtr->first);
            }
            std::set<uint32_t> ySet;
            for (auto jtr = ys.begin(); jtr != ys.end(); ++jtr)
            {
                ySet.insert(jtr->first);
            }
            std::cerr << "missing from xs:";
            for (auto jtr = ySet.begin(); jtr != ySet.end(); ++jtr)
            {
                if (!xSet.contains(*jtr))
                {
                    std::cerr << ' ' << *jtr;
                }
            }
            std::cerr << std::endl;
            std::cerr << "missing from ys:";
            for (auto jtr = xSet.begin(); jtr != xSet.end(); ++jtr)
            {
                if (!ySet.contains(*jtr))
                {
                    std::cerr << ' ' << *jtr;
                }
            }
            std::cerr << std::endl;
        }
        REQUIRE( ys.size() == xs.size() );
        for (auto jtr = ys.begin(); jtr != ys.end(); ++jtr)
        {
            REQUIRE( xs.contains(jtr->first) );
            if (jtr->second != xs.at(jtr->first))
            {
                std::cerr << itr->first << '\t' << jtr->first << '\t' << jtr->second << '\t' << xs.at(jtr->first) << std::endl;
            }
            REQUIRE( jtr->second == xs.at(jtr->first) );
        }
    }
}
