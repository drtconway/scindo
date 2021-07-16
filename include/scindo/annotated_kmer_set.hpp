#ifndef SCINDO_ANNOTATED_KMER_SET_HPP
#define SCINDO_ANNOTATED_KMER_SET_HPP

#include <deque>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include <sdsl/sd_vector.hpp>
#include "scindo/kmers.hpp"

namespace scindo
{
    namespace detail
    {
        struct seen_annotated_kmer_set
        {
            seen_annotated_kmer_set(const std::vector<kmer>& p_kmers, std::vector<uint64_t>&& p_annot)
                : kmers(p_kmers.begin(), p_kmers.end()),
                  kmers_rank(&kmers),
                  kmers_select(&kmers),
                  annot(std::move(p_annot)),
                  N(p_kmers.size()),
                  xMax(p_kmers.size() > 0 ? p_kmers.back() : 0)
            {
            }

            size_t size() const
            {
                return N;
            }

            bool contains(const kmer& p_x) const
            {
                uint64_t r0 = rank(p_x);
                uint64_t r1 = rank(p_x+1);
                return r1 == r0 + 1;
            }

            bool update(const kmer& p_x, const uint64_t& p_a)
            {
                uint64_t r0 = rank(p_x);
                uint64_t r1 = rank(p_x+1);
                if ( r1 == r0 + 1)
                {
                    annot[r0] |= p_a;
                    return true;
                }
                return false;
            }

            size_t rank(const kmer& p_x) const
            {
                if (p_x > xMax)
                {
                    return size();
                }
                return kmers_rank.rank(p_x);
            }

            kmer select(const size_t& p_i) const
            {
                return kmers_select.select(p_i + 1);
            }

            static std::shared_ptr<seen_annotated_kmer_set> merge(
                    std::shared_ptr<seen_annotated_kmer_set>& p_lhs,
                    std::shared_ptr<seen_annotated_kmer_set>& p_rhs)
            {
                std::vector<kmer> xs;
                std::vector<uint64_t> as;

                {
                    const seen_annotated_kmer_set& lhs = *p_lhs;
                    const seen_annotated_kmer_set& rhs = *p_rhs;

                    const size_t n = lhs.size() + rhs.size();
                    xs.reserve(n);
                    as.reserve(n);

                    const size_t I = lhs.size();
                    size_t i = 0;
                    kmer xi = 0;
                    uint64_t ai = 0;

                    const size_t J = rhs.size();
                    size_t j = 0;
                    kmer xj = 0;
                    uint64_t aj = 0;

                    if (i < I)
                    {
                        xi = lhs.select(i);
                        ai = lhs.annot[i];
                    }
                    if (j < J)
                    {
                        xj = rhs.select(j);
                        aj = rhs.annot[j];
                    }

                    while (i < I && j < J)
                    {
                        if (xi < xj)
                        {
                            xs.push_back(xi);
                            as.push_back(ai);
                            i += 1;
                            if (i < I)
                            {
                                xi = lhs.select(i);
                                ai = lhs.annot[i];
                            }
                            continue;
                        }
                        if (xi > xj)
                        {
                            xs.push_back(xj);
                            as.push_back(aj);
                            j += 1;
                            if (j < J)
                            {
                                xj = rhs.select(j);
                                aj = rhs.annot[j];
                            }
                            continue;
                        }
                        std::cerr << i << '\t' << kmers::render(25, xi)
                            << '\t' << j << '\t' << kmers::render(25, xj)
                            << std::endl;
                        throw std::logic_error("succinct sets should not overlap");
                    }
                    while (i < I)
                    {
                        xs.push_back(xi);
                        as.push_back(ai);
                        i += 1;
                        if (i < I)
                        {
                            xi = lhs.select(i);
                            ai = lhs.annot[i];
                        }
                    }
                    while (j < J)
                    {
                        xs.push_back(xj);
                        as.push_back(aj);
                        j += 1;
                        if (j < J)
                        {
                            xj = rhs.select(j);
                            aj = rhs.annot[j];
                        }
                    }
                }
                p_lhs = std::shared_ptr<seen_annotated_kmer_set>();
                p_rhs = std::shared_ptr<seen_annotated_kmer_set>();

                return std::shared_ptr<seen_annotated_kmer_set>(new seen_annotated_kmer_set(xs, std::move(as)));
            }

            const sdsl::sd_vector<> kmers;
            const typename sdsl::sd_vector<>::rank_1_type kmers_rank;
            const typename sdsl::sd_vector<>::select_1_type kmers_select;
            std::vector<uint64_t> annot;
            const size_t N;
            const kmer xMax;
        };
        using seen_annotated_kmer_set_ptr = std::shared_ptr<seen_annotated_kmer_set>;
    };

    struct annotated_kmer_set
    {
        std::deque<detail::seen_annotated_kmer_set_ptr> seen;
        std::vector<std::pair<kmer,uint64_t>> mini_buf;
        std::vector<kmer> kmer_buf;
        std::vector<uint64_t> annot_buf;

        void add(const kmer& p_x, const uint64_t& p_a)
        {
            mini_buf.push_back(std::make_pair(p_x, p_a));
            if (mini_buf.size() > 1024*1024)
            {
                std::sort(mini_buf.begin(), mini_buf.end());
                for (size_t i = 0; i < mini_buf.size(); ++i)
                {
                    add0(mini_buf[i].first, mini_buf[i].second);
                }
                mini_buf.clear();
            }
        }

        void add0(const kmer& p_x, const uint64_t& p_a)
        {
            for (size_t i = 0; i < seen.size(); ++i)
            {
                if (seen[i]->update(p_x, p_a))
                {
                    return;
                }
            }
            kmer_buf.push_back(p_x);
            annot_buf.push_back(p_a);
        }

        size_t size() const
        {
            size_t z = 0;
            for (size_t i = 0; i < seen.size(); ++i)
            {
                z += seen[i]->size();
            }
            return z + kmer_buf.size();
        }

        void flush()
        {
            if (mini_buf.size() > 0)
            {
                std::sort(mini_buf.begin(), mini_buf.end());
                for (size_t i = 0; i < mini_buf.size(); ++i)
                {
                    add0(mini_buf[i].first, mini_buf[i].second);
                }
                mini_buf.clear();
            }

            if (kmer_buf.size() == 0)
            {
                return;
            }

            prepare();
            seen.push_back(
                detail::seen_annotated_kmer_set_ptr(new detail::seen_annotated_kmer_set(kmer_buf, std::move(annot_buf))));
            kmer_buf.clear();
            std::vector<kmer>().swap(kmer_buf);
            annot_buf.clear();
            std::vector<uint64_t>().swap(annot_buf);

            std::sort(seen.rbegin(), seen.rend(), [](auto l, auto r) { return l->size() < r->size(); });

            while (seen.size() >= 2)
            {
                if (seen.back()->size() > 256*1024*1024)
                {
                    // Don't try and merge really big sets.
                    break;
                }
                auto rhs = seen.back();
                seen.pop_back();
                auto lhs = seen.back();
                seen.pop_back();
                seen.push_back(detail::seen_annotated_kmer_set::merge(lhs, rhs));
            }

        }

        void prepare()
        {
            {
                // Create a permutation vector & sort
                //
                std::vector<uint32_t> perm;
                perm.reserve(kmer_buf.size());
                for (uint32_t i = 0; i < kmer_buf.size(); ++i)
                {
                    perm.push_back(i);
                }
                std::sort(perm.begin(), perm.end(), [&](uint32_t i, uint32_t j) { return kmer_buf[i] < kmer_buf[j]; });

                // Apply the permutation vector
                //
                if (0)
                {
                    std::vector<bool> done(perm.size());
                    for (uint32_t i = 0; i < perm.size(); ++i)
                    {
                        if (done[i])
                        {
                            continue;
                        }
                        done[i] = true;
                        uint32_t k = i;
                        uint32_t j = perm[i];
                        while (i != j)
                        {
                            std::swap(kmer_buf[k], kmer_buf[j]);
                            done[j] = true;
                            k = j;
                            j = perm[j];
                        }
                    }
                }

                for (uint32_t i = 0; i < perm.size(); i++)
                {
                    kmer x = kmer_buf[i];
                    uint64_t a = annot_buf[i];
                    uint32_t j = i;
                    while (i != perm[j])
                    {
                        uint32_t k = perm[j];
                        kmer_buf[j] = kmer_buf[k];
                        annot_buf[j] = annot_buf[k];
                        perm[j] = j;
                        j = k;
                    }
                    kmer_buf[j] = x;
                    annot_buf[j] = a;
                    perm[j] = j;
                }

                // and now we're done with the permutation vector.
            }

            // merge the annotations for duplicate kmers
            //
            uint32_t i = 0;
            uint32_t j = 0;
            while (i < kmer_buf.size())
            {
                uint32_t i0 = i;
                i += 1;
                while (i < kmer_buf.size() && kmer_buf[i] == kmer_buf[i - 1])
                {
                    ++i;
                }

                if (j == i0 && i0 + 1 == i)
                {
                    // the trivial case
                    continue;
                }

                uint64_t a = annot_buf[i0];
                for (uint32_t k = i0 + 1; k < i; ++k)
                {
                    a |= annot_buf[k];
                }
                kmer_buf[j] = kmer_buf[i0];
                annot_buf[j] = a;
                j += 1;
            }
            kmer_buf.resize(j);
            annot_buf.resize(j);
        }
    };
}
// namespace scindo

#endif // scindo_annotated_kmer_set_hpp
