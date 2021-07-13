#ifndef SCINDO_STABBY_HPP
#define SCINDO_STABBY_HPP

#include <bit>
#include <memory>
#include "scindo/stabby/schmidt.hpp"
#include "scindo/stabby/dense_domain.hpp"
#include "scindo/stabby/sparse_domain.hpp"

namespace scindo
{
    struct stabby
    {
        using position = uint32_t;
        using interval = std::pair<position,position>;

        stabby(const std::vector<interval>& p_intervals)
        {
            // First, collect all the positions in the intervals:
            //
            std::vector<position> tmp;
            tmp.reserve(2*p_intervals.size());
            for (auto itr = p_intervals.begin(); itr != p_intervals.end(); ++itr)
            {
                tmp.push_back(itr->first);
                tmp.push_back(itr->second);
            }
            std::sort(tmp.begin(), tmp.end());
            tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());

            std::vector<uint32_t> ump;
            uint32_t q = 1;
            for (size_t i = 0; i < tmp.size(); ++i, ++q)
            {
                if (i > 0 && tmp[i - 1] + 1 < tmp[i])
                {
                    ++q;
                }
                ump.push_back(q);
            }
            ++q;
            Q = q;

            domain = std::shared_ptr<sparse_domain>(new sparse_domain(tmp.begin(), tmp.end()));
            
            sdsl::bit_vector vmp(ump.back() + 1);
            for (auto itr = ump.begin(); itr != ump.end(); ++itr)
            {
                vmp[*itr] = 1;
            }
            gapped = std::shared_ptr<dense_domain>(new dense_domain(vmp));

            std::vector<schmidt::interval> dense;
            dense.reserve(p_intervals.size());
            for (auto itr = p_intervals.begin(); itr != p_intervals.end(); ++itr)
            {
                schmidt::interval x = {sparse_to_dense(itr->first), sparse_to_dense(itr->second)};
                dense.push_back(x);
            }

            impl = std::shared_ptr<schmidt>(new schmidt(Q, dense));
        }

        size_t find(const position& p_pos, std::vector<interval>& p_res) const
        {
            uint32_t q = sparse_to_dense(p_pos);

            std::vector<schmidt::interval> res;
            impl->stab(q, res);

            for (auto itr = res.begin(); itr != res.end(); ++itr)
            {
                interval r = {dense_to_sparse(itr->first), dense_to_sparse(itr->second)};
                p_res.push_back(r);
            }
            return res.size();
        }

        uint32_t sparse_to_dense(uint32_t p_pos) const
        {
            uint32_t r0 = domain->rank(p_pos);
            uint32_t r1 = domain->rank(p_pos+1);
            if (r0 == gapped->count())
            {
                return gapped->size();
            }
            uint32_t q = gapped->select(r0);

            return (r1 > r0 ? q : q - 1);
        }

        uint32_t dense_to_sparse(uint32_t p_r) const
        {
            uint32_t r = gapped->rank(p_r);
            return domain->select(r);
        }

        uint32_t Q;
        std::shared_ptr<sparse_domain> domain;
        std::shared_ptr<dense_domain> gapped;
        std::shared_ptr<schmidt> impl;
    };
}
// namespace scindo

#endif // SCINDO_STABBY_HPP
