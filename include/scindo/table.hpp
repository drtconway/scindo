#ifndef SCINDO_TABLE_HPP
#define SCINDO_TABLE_HPP

#include <iostream>
#include <fstream>
#include <sstream>

namespace scindo
{
    namespace detail
    {
        template <size_t I, typename... Columns>
        std::enable_if<I == sizeof...(Columns), int>::type
        compare(const int& p_keyNum, const std::tuple<Columns...>& p_lhs, const std::tuple<Columns...>& p_rhs)
        {
            std::ostringstream s;
            s << "key number " << p_keyNum << " (" << (p_keyNum + 1) << ") out of " << sizeof...(Columns);
            
            throw std::runtime_error(s.str());
        }

        template <size_t I, typename... Columns>
        std::enable_if<I < sizeof...(Columns), int>::type
        compare(const int& p_keyNum, const std::tuple<Columns...>& p_lhs, const std::tuple<Columns...>& p_rhs)
        {
            if (p_keyNum == I)
            {
                const auto& lhs = std::get<I>(p_lhs);
                const auto& rhs = std::get<I>(p_rhs);
                if (lhs < rhs)
                {
                    return -1;
                }
                if (lhs > rhs)
                {
                    return 1;
                }
                return 0;
            }
            else
            {
                return compare<I+1>(p_keyNum, p_lhs, p_rhs);
            }
        }

        template <size_t I, typename... Columns>
        std::enable_if<I == sizeof...(Columns), std::ostream&>::type
        write_row(std::ostream& p_out, const std::tuple<Columns...>& p_row)
        {
            return p_out;
        }

        template <size_t I, typename... Columns>
        std::enable_if<I < sizeof...(Columns), std::ostream&>::type
        write_row(std::ostream& p_out, const std::tuple<Columns...>& p_row)
        {
            if (I > 0)
            {
                p_out << '\t';
            }
            p_out << std::get<I>(p_row);
            return write_row<I+1>(p_out, p_row);
        }
    }
    // namespace detail

    template <typename... Columns>
    struct table : std::vector<std::tuple<Columns...>>
    {
        using row_type = std::tuple<Columns...>;

        const std::vector<std::string> hdrs;
        std::unordered_map<std::string,size_t> idx;

        table(std::initializer_list<std::string> p_hdrs)
            : hdrs(p_hdrs)
        {
            if (p_hdrs.size() != sizeof...(Columns))
            {
                throw std::runtime_error("header and type are not the same width");
            }
            for (size_t i = 0; i < hdrs.size(); ++i)
            {
                idx[hdrs[i]] = i;
            }
        }

        void write(const std::string& p_filename) const
        {
            std::ofstream out(p_filename);
            write(out);
        }

        void write(std::ostream& p_out) const
        {
            for (size_t i = 0; i < hdrs.size(); ++i)
            {
                if (i > 0)
                {
                    p_out << '\t';
                }
                p_out << hdrs[i];
            }
            p_out << std::endl;

            for (auto itr = this->begin(); itr != this->end(); ++itr)
            {
                detail::write_row<0>(p_out, *itr);
                p_out << std::endl;
            }
        }

        void sort(const std::vector<std::string>& p_keys)
        {
            std::vector<int> keys;
            for (auto itr = p_keys.begin(); itr != p_keys.end(); ++itr)
            {
                const auto& k0 = *itr;
                if (k0.starts_with('-'))
                {
                    std::string k(k0.begin() + 1, k0.end());
                    if (!idx.contains(k))
                    {
                        throw std::runtime_error("sort key does not exist");
                    }
                    keys.push_back(-(idx[k] + 1));
                }
                else
                {
                    if (!idx.contains(k0))
                    {
                        throw std::runtime_error("sort key does not exist");
                    }
                    keys.push_back(idx[k0] + 1);
                }
            }
            sort(keys);
        }

        void sort(const std::vector<int>& p_keys)
        {
            std::sort(this->begin(), this->end(), [&](const auto& p_lhs, const auto& p_rhs) {
                for (auto itr = p_keys.begin(); itr != p_keys.end(); ++itr)
                {
                    int kn = *itr;
                    int n = (kn < 0 ? -kn : kn) - 1;
                    int c = detail::compare<0>(n, p_lhs, p_rhs);
                    //std::cerr << "cmp: " << nlohmann::json(p_lhs) << '\t' << nlohmann::json(p_rhs) << '\t' << c << std::endl;
                    if (kn > 0)
                    {
                        if (c != 0)
                        {
                            return c < 0;
                        }
                    }
                    else
                    {
                        if (c != 0)
                        {
                            return c > 0;
                        }
                    }
                }
                return false;
            });
        }

        template <typename X>
        void filter(X p_predicate)
        {
            static_assert(std::is_convertible<X, std::function<bool(const row_type&)>>::value);
            auto itr = this->begin();
            auto jtr = this->begin();
            while (itr != this->end())
            {
                if (p_predicate(*itr))
                {
                    if (jtr != itr)
                    {
                        *jtr = *itr;
                    }
                    ++jtr;
                }
                ++itr;
            }
            this->erase(jtr, this->end());
        }
    };
}
// namespace scindo

#endif // SCINDO_TABLE_HPP
