#ifndef SCINDO_MURMER3_HPP
#define SCINDO_MURMER3_HPP

namespace scindo
{
    // Implementation based on Wikipedia page.
    //
    struct murmur3
    {
        static constexpr uint32_t c1 = 0xcc9e2d51U;
        static constexpr uint32_t c2 = 0x1b873593U;
        static constexpr uint32_t r1 = 15;
        static constexpr uint32_t r2 = 13;
        static constexpr uint32_t m = 5;
        static constexpr uint32_t n = 0xe6546b64U;

        murmur3(uint32_t p_seed)
            : h(p_seed), l(0)
        {
        }

        murmur3& update(const std::string& p_str)
        {
            update(p_str.size());
            uint32_t x = 0;
            for (size_t i = 0; i < p_str.size(); ++i)
            {
                x = (x << 8) | uint32_t(p_str[i]);
                if ((i & 3) == 3)
                {
                    update(x);
                    x = 0;
                }
            }
            update(x);
            return *this;
        }

        murmur3& update(uint64_t x)
        {
            uint32_t hi = static_cast<uint32_t>(x >> 32);
            uint32_t lo = static_cast<uint32_t>(x);
            return update(hi).update(lo);
        }

        murmur3& update(uint32_t k)
        {
            k *= c1;
            k = rotl(k, r1);
            k *= c2;

            h ^= k;
            h = rotl(h, r2);
            h = (h * m) + n;

            l += 4;

            return *this;
        }

        uint32_t operator()() const
        {
            uint32_t r = h ^ l;
            r ^= (r >> 16);
            r *= 0x85ebca6bU;
            r ^= (r >> 13);
            r *= 0xc2b2ae35U;
            r ^= (r >> 16);
            return r;
        }

        static inline uint32_t rotl(uint32_t x, uint32_t r)
        {
            return (x << r) ^ (x >> (32 - r));
        }

        static inline uint32_t rotr(uint32_t x, uint32_t r)
        {
            return (x >> r) ^ (x << (32 - r));
        }

        uint32_t h;
        uint64_t l;
    };
}
// namespace scindo

#endif // SCINDO_MURMER3_HPP
