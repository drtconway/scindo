#ifndef SCINDO_TABLE_COLUMN_HPP
#define SCINDO_TABLE_COLUMN_HPP

#include <charconv>
#include <functional>
#include <string>
#include <string_view>
#include <tuple>

namespace scindo
{
    namespace column
    {
        template <typename X>
        std::enable_if<std::is_convertible<X, std::function<bool(const std::string_view&)>>::value, bool>::type
        make_list(const std::string_view& p_str, const char p_sep, X p_acceptor)
        {
            auto b = p_str.begin();
            while (b != p_str.end())
            {
                auto p = b;
                while (b != p_str.end() && *b != p_sep)
                {
                    ++b;
                }
                if (!p_acceptor(std::string_view(p, b)))
                {
                    return false;
                }
                if (b != p_str.end())
                {
                    // b == p_sep
                    ++b;
                }
            }
            return true;
        }

        template <char Sep, typename X>
        std::enable_if<std::is_convertible<X, std::function<bool(const std::string_view&)>>::value, bool>::type
        make_list(const std::string_view& p_str, X p_acceptor)
        {
            return make_list<>(p_str, Sep, p_acceptor);
        }

        template <typename X>
        std::enable_if<std::is_convertible<X, std::function<bool(const std::string_view&,
                                                                 const std::string_view&)>>::value, bool>::type
        make_pair(const std::string_view& p_str, const char p_sep, X p_acceptor)
        {
            size_t pos = p_str.find(p_sep);
            return pos != p_str.npos
                && p_acceptor(std::string_view(p_str.begin(), p_str.begin() + pos),
                              std::string_view(p_str.begin() + pos + 1, p_str.end()));
        }

        template <char Sep, typename X>
        std::enable_if<std::is_convertible<X, std::function<bool(const std::string_view&,
                                                                 const std::string_view&)>>::value, bool>::type
        make_pair(const std::string_view& p_str, X p_acceptor)
        {
            return make_pair<>(p_str, Sep, p_acceptor);
        }

        bool make(const std::string_view& p_str, std::string_view& p_res)
        {
            p_res = p_str;
            return true;
        }

        bool make(const std::string_view& p_str, std::string& p_res)
        {
            p_res = std::string(p_str.begin(), p_str.end());
            return true;
        }

        template <typename Int, std::enable_if_t<std::is_integral<Int>::value, bool> = true>
        bool make(const std::string_view& p_str, Int& p_res)
        {
            const char* b = std::data(p_str);
            const char* e = b + std::size(p_str);
            if(auto [p, ec] = std::from_chars(b, e, p_res);
               ec == std::errc() && p == e)
            {
                return true;
            }
            return false;
        }

        template <typename Flt, std::enable_if_t<std::is_floating_point<Flt>::value, bool> = true>
        bool make(const std::string_view& p_str, Flt& p_res)
        {
            const char* b = std::data(p_str);
            const char* e = b + std::size(p_str);
            //if(auto [p, ec] = std::from_chars(b, e, p_res);
            //   ec == std::errc() && p == e)
            //{
            //    return true;
            //}
            //return false;
            char* e0 = NULL;
            p_res = std::strtod(b, &e0);
            return e0 == e;
        }

        template<typename T>
        bool make(const std::string_view& p_str, std::vector<T>& p_res, const char p_sep = '\t')
        {
            p_res.clear();
            return make_list<>(p_str, p_sep, [&](const std::string_view& p_item) {
                p_res.push_back(T());
                return make(p_item, p_res.back());
            });
        }

        template <typename T, typename U>
        bool make(const std::string_view& p_str, std::pair<T,U>& p_res, const char p_sep = ',')
        {
            return make_pair<>(p_str, p_sep, [&](const std::string_view& p_first, const std::string_view& p_second) {
                return make(p_first, p_res.first) && make(p_second, p_res.second);
            });
        }

        template <size_t I, typename... Columns>
        std::enable_if<I == sizeof...(Columns), bool>::type
        make(std::vector<std::string_view>& p_strs, std::tuple<Columns...>& p_res)
        {
            return p_strs.size() == I;
        }

        template <size_t I, typename... Columns>
        std::enable_if<I < sizeof...(Columns), bool>::type
        make(std::vector<std::string_view>& p_strs, std::tuple<Columns...>& p_res)
        {
            if (p_strs.size() <= I)
            {
                return false;
            }
            if (!make(p_strs[I], std::get<I>(p_res)))
            {
                return false;
            }
            return make<I+1,Columns...>(p_strs, p_res);
        }

        template <typename... Columns>
        bool make(const std::string_view& p_strs, std::tuple<Columns...>& p_res, const char p_sep = '\t')
        {
            std::vector<std::string_view> parts;
            return make(p_strs, parts, p_sep) && make<0>(parts, p_res);
        }
    }
    // namespace column
}
// namespace scindo

#endif // SCINDO_TABLE_COLUMN_HPP
