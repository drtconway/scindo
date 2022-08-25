#ifndef SCINDO_ATOM_HPP
#define SCINDO_ATOM_HPP

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>

namespace scindo
{
    struct atom : std::string_view
    {
        atom(const char* p_str)
            : std::string_view(find_or_make(p_str))
        {
        }

        atom(const std::string& p_str)
            : std::string_view(find_or_make(p_str))
        {
        }

        static const std::string_view& find_or_make(const char* p_str)
        {
            std::string_view x = std::string_view(p_str);
            return find_or_make(x);
        }

        static const std::string_view& find_or_make(const std::string& p_str)
        {
            std::string_view x = std::string_view(p_str.begin(), p_str.end());
            return find_or_make(x);
        }

        static const std::string_view& find_or_make(const std::string_view& p_str)
        {
            std::unordered_map<std::string_view,std::unique_ptr<std::string>>& S = store();
            auto itr = S.find(p_str);
            if (itr != S.end())
            {
                std::unique_ptr<std::string> p = std::unique_ptr<std::string>(new std::string(p_str.begin(), p_str.end()));
                std::string_view k = std::string_view(p->begin(), p->end());
                S[k] = std::move(p);
                itr = S.find(p_str);
            }
            return itr->first;
        }

        static std::unordered_map<std::string_view,std::unique_ptr<std::string>>& store()
        {
            static std::unordered_map<std::string_view,std::unique_ptr<std::string>> S;
            return S;
        }
    };
}
// namespace scindo

#endif // SCINDO_ATOM_HPP
