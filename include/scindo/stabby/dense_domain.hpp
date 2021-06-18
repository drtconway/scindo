#ifndef SCINDO_STABBY_DENSE_DOMAIN_HPP
#define SCINDO_STABBY_DENSE_DOMAIN_HPP

#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_mcl.hpp>

namespace scindo
{
    struct dense_domain
    {
        using sparse_type = uint32_t;
        using dense_type = uint32_t;

        dense_domain(const sdsl::bit_vector& p_array)
            : array(p_array),
              array_rank(&array),
              array_select(&array)
        {
        }

        bool contains(const sparse_type& p_x) const
        {
            uint64_t r0 = array_rank.rank(p_x);
            uint64_t r1 = array_rank.rank(p_x+1);
            return r1 == r0 + 1;
        }

        sparse_type select(const dense_type& p_x) const
        {
            return array_select.select(p_x + 1);
        }

        dense_type rank(const sparse_type& p_x) const
        {
            return array_rank.rank(p_x);
        }

        sdsl::bit_vector array;
        typename sdsl::rank_support_v5<> array_rank;
        typename sdsl::select_support_mcl<> array_select;
    };
}
// namespace scindo

#endif // SCINDO_STABBY_DENSE_DOMAIN_HPP
