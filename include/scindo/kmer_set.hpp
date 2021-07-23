#ifndef SCINDO_KMER_SET_HPP
#define SCINDO_KMER_SET_HPP

#include <vector>
#include <bit>
#include <sdsl/sd_vector.hpp>

#include "scindo/kmers.hpp"

namespace scindo
{
    namespace detail
    {
        struct heap_kmer_set
        {
            const uint64_t B;
            std::vector<kmer> X;

            heap_kmer_set(const std::vector<kmer>& p_X)
                : B(roundup(p_X.size()+1))
            {
                std::vector<kmer> tmp(p_X);
                tmp.resize(B - 1, uint64_t(-1LL));
                X.resize(B, uint64_t(-1LL));
                eytzinger(tmp, X);
            }

            bool contains(const kmer& p_x) const
            {
                size_t k = 1;
                const size_t N = X.size();
                while (k < N)
                {
                    k = (k << 1) + (X[k] < p_x ? 1 : 0);
                }
                //k >>= __builtin_ffs(~k);
                k >>= std::countr_zero(~k) + 1;
                return X[k] == p_x;
            }

            static uint64_t roundup(uint64_t v)
            {
                v--;
                v |= v >> 1;
                v |= v >> 2;
                v |= v >> 4;
                v |= v >> 8;
                v |= v >> 16;
                v |= v >> 32;
                return v + 1;
            }

            static size_t eytzinger(const std::vector<kmer>& p_orig, std::vector<kmer>& p_permuted, size_t p_i = 0, size_t p_k = 1)
            {
                if (p_k <= p_orig.size())
                {
                    p_i = eytzinger(p_orig, p_permuted, p_i, 2 * p_k);
                    p_permuted[p_k] = p_orig[p_i++];
                    p_i = eytzinger(p_orig, p_permuted, p_i, 2 * p_k + 1);
                }
                return p_i;
            }
        };

        struct succinct_kmer_set
        {
            const sdsl::sd_vector<> X;
            const kmer maxx;

            succinct_kmer_set(const std::vector<kmer>& p_X)
                : X(p_X.begin(), p_X.end()), maxx(p_X.size() > 0 ? p_X.back() : 0)
            {
            }

            bool contains(kmer x) const
            {
                if (x > maxx)
                {
                    return false;
                }
                return X[x];
            }

        };

        struct vector_kmer_set
        {
            const size_t P;
            std::vector<kmer> X;

            vector_kmer_set(const std::vector<kmer>& p_X)
                : P(roundup(p_X.size()) >> 1), X(P << 1, uint64_t(-1LL))
            {
                for (size_t i = 0; i < p_X.size(); ++i)
                {
                    X[i] = p_X[i];
                }
            }

            bool contains(kmer x) const
            {
                uint32_t j = 0;
                uint32_t k = P;
                while (k > 0)
                {
                    uint32_t r = j | k;
                    if (x >= X[r])
                    {
                        j = r;
                    }
                    k >>= 1;
                }
                return X[j] == x;
            }

            static uint64_t roundup(uint64_t v)
            {
                v--;
                v |= v >> 1;
                v |= v >> 2;
                v |= v >> 4;
                v |= v >> 8;
                v |= v >> 16;
                v |= v >> 32;
                return v + 1;
            }


        };
    }
    // namespace detail

    struct kmer_set : detail::vector_kmer_set
    {
        using detail::vector_kmer_set::vector_kmer_set;

        template <typename X>
        void with_contains(const std::vector<kmer>& p_xs, X p_acceptor) const
        {
            for (size_t i = 0; i < p_xs.size(); ++i)
            {
                kmer x = p_xs[i];
                if (contains(x))
                {
                    p_acceptor(x);
                }
            }
        }

    };

}
// namespace scindo

#endif // SCINDO_KMER_SET_HPP
