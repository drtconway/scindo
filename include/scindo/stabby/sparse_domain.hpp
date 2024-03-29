#ifndef SCINDO_STABBY_SPARSE_DOMAIN_HPP
#define SCINDO_STABBY_SPARSE_DOMAIN_HPP

#include <sdsl/sd_vector.hpp>

namespace scindo
{
    struct sparse_domain
    {
        using sparse_type = uint32_t;
        using dense_type = uint32_t;

        template <typename Itr>
        sparse_domain(Itr p_begin, Itr p_end)
            : array(p_begin, p_end),
              array_rank(&array),
              array_select(&array),
              N(p_end - p_begin),
              Z(p_end != p_begin ? *std::prev(p_end) + 1: 0)
        {
        }

        size_t size() const
        {
            return Z;
        }

        size_t count() const
        {
            return N;
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
            if (p_x >= size())
            {
                return count();
            }
            return array_rank.rank(p_x);
        }

        sdsl::sd_vector<> array;
        typename sdsl::sd_vector<>::rank_1_type array_rank;
        typename sdsl::sd_vector<>::select_1_type array_select;
        const size_t N;
        const uint32_t Z;
    };
}
// namespace scindo

#endif // SCINDO_STABBY_SPARSE_DOMAIN_HPP
