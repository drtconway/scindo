#ifndef SCINDO_STABBY_SCHMIDT_HPP
#define SCINDO_STABBY_SCHMIDT_HPP

#include <list>
#include <deque>
#include <unordered_map>
#include <vector>

namespace std
{
    template <typename T1, typename T2>
    struct hash<std::pair<T1,T2>>
    {
        std::size_t operator()(const std::pair<T1,T2>& x) const
        {
            std::size_t h1 = std::hash<T1>()(x.first);
            std::size_t h2 = std::hash<T1>()(x.second);
            return std::rotl(h1, 13) ^ std::rotr(h2, 15);
        }
    };

    template <>
    struct hash<std::pair<uint32_t,uint32_t>>
    {
        std::size_t operator()(const std::pair<uint32_t,uint32_t>& x) const
        {
            uint64_t h = (static_cast<uint64_t>(x.first) << 32) | static_cast<uint64_t>(x.second);
            return std::hash<uint64_t>()(h);
        }
    };

}
// namespace

namespace scindo
{
    struct schmidt
    {
        using position = uint32_t;
        using interval = std::pair<position,position>;

        // This uses the vector of intervals destructively.
        //
        schmidt(const position& p_Q, std::vector<interval>& p_intervals)
            : Q(p_Q)
        {
            std::sort(p_intervals.begin(), p_intervals.end());

            make_smaller(p_intervals);

            std::unordered_map<position,std::vector<interval>> event;

            for (size_t i = 0; i < p_intervals.size(); ++i)
            {
                const interval& a = p_intervals[i];

                event[a.second].push_back(a);
                event[a.first].push_back(a);
            }

            std::list<interval> L;
            std::unordered_map<interval,std::list<interval>::iterator> saved;

            for (size_t q = 0; q < Q; ++q)
            {
                if (0)
                {
                    std::cerr << q;
                    for (auto itr = event[q].begin(); itr != event[q].end(); ++itr)
                    {
                        const interval& a = *itr;
                        std::cerr << ' '
                                  << '[' << a.first << ',' << a.second << ']';
                    }
                    std::cerr << std::endl;
                }

                if (L.size() > 0)
                {
                    const interval& a = L.back();
                    start[q] = a;
                    if (0)
                    {
                        std::cerr << "start(" << q << ") <- "
                                  << '[' << a.first << ',' << a.second << ']'
                                  << std::endl;
                    }
                    
                }

                for (auto itr = event[q].rbegin(); itr != event[q].rend(); ++itr)
                {
                    const interval& a = *itr;

                    if (a.first == q && !saved.contains(a))
                    {
                        start[q] = a;
                        if (0)
                        {
                            std::cerr << "start(" << q << ") <- "
                                      << '[' << a.first << ',' << a.second << ']'
                                      << std::endl;
                        }

                        L.push_back(a);
                        saved[a] = std::prev(L.end());
                    }
                    else // a.second == q
                    {
                        interval p = root();
                        if (saved.at(a) != L.begin())
                        {
                            p = *std::prev(saved[a]);
                        }
                        parent[a] = p;
                        if (last.contains(p))
                        {
                            left[a] = last.at(p);
                        }
                        last[p] = a;
                        if (0)
                        {
                            const auto& b = parent.at(a);
                            std::cerr << "parent("
                                      << '[' << a.first << ',' << a.second << ']'
                                      << ") <- "
                                      << '[' << b.first << ',' << b.second << ']'
                                      << std::endl;
                        }
                        L.erase(saved[a]);
                        saved.erase(a);
                    }
                }
            }
        }

        size_t stab(const position& p_q, std::vector<interval>& p_res) const
        {
            if (!start.contains(p_q))
            {
                return 0;
            }

            p_res.clear();
            std::deque<interval> kew;

            for (auto v = start.at(p_q); v != root(); v = parent.at(v))
            {
                kew.push_front(v);
            }

            while (kew.size() > 0)
            {
                interval a = kew.back();
                kew.pop_back();
                p_res.push_back(a);

                if (smaller.contains(a))
                {
                    const auto& s = smaller.at(a);
                    for (auto itr = s.rbegin(); itr != s.rend() & itr->second >= p_q; ++itr)
                    {
                        p_res.push_back(*itr);

                    }
                }

                if (left.contains(a))
                {
                    interval t = left.at(a);
                    while (true)
                    {
                        if (t.second < p_q)
                        {
                            break;
                        }
                        kew.push_back(t);
                        if (!last.contains(t))
                        {
                            break;
                        }
                        t = last.at(t);
                    }
                }
            }

            std::reverse(p_res.begin(), p_res.end());

            return p_res.size();
        }

        static const interval& root()
        {
            static const interval R = {0,0};
            return R;
        }

        void make_smaller(std::vector<interval>& p_intervals)
        {
            auto itr = p_intervals.begin();
            auto jtr = itr;
            while (itr != p_intervals.end())
            {
                auto ktr = itr;
                while (itr != p_intervals.end() && itr->first == ktr->first)
                {
                    ++itr;
                }
                auto ltr = itr - 1;
                while (ktr != ltr)
                {
                    smaller[*ltr].push_back(*ktr);
                    ++ktr;
                }
                if (jtr != ltr)
                {
                    *jtr = *ltr;
                }
                ++jtr;
            }
            p_intervals.erase(jtr, p_intervals.end());
        }

        const position Q;
        std::unordered_map<interval,std::vector<interval>> smaller;
        std::unordered_map<position,interval> start;
        std::unordered_map<interval,interval> parent;
        std::unordered_map<interval,interval> last;
        std::unordered_map<interval,interval> left;
    };
}
// namespace scindo

#endif // SCINDO_STABBY_SCHMIDT_HPP
