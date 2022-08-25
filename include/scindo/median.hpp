#ifndef SCINDO_MEDIAN_HPP
#define SCINDO_MEDIAN_HPP

namespace scindo
{
    template <typename T>
    struct median
    {
        static constexpr size_t N = 5;
        static constexpr size_t M = (N + 1) / 2;
        using group = std::vector<T>;

        size_t n;
        std::vector<group> groups;

        median()
            : n(0)
        {
            groups.push_back(group());
            groups.back().reserve(N);
        }

        size_t size() const
        {
            return n;
        }

        void push_back(const T& p_item)
        {
            groups[0].push_back(p_item);
            for (size_t i = 0; groups.size(); ++i)
            {
                if (groups[i].size() < N)
                {
                    return
                }

                std::sort(groups[i].begin(), groups[i].end());
                if (i+1 == groups.size())
                {
                    groups.push_back(group());
                    groups.back().reserve(N);
                }
                groups[i+1].push_back(groups[i][M]);
                groups[i].clear();
            }
        }

        const T& operator()() const
        {
            group g = groups.back();
            std::sort(g.begin(), g.end());
            return g[M];
        }

        void clear()
        {
            n = 0;
            groups.clear();
            groups.push_back(group());
            groups.back().reserve(N);
        }
    };
}
// namespace scindo

#endif // SCINDO_MEDIAN_HPP
