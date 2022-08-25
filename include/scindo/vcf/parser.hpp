#ifndef SCINDO_VCF_PARSER_HPP
#define SCINDO_VCF_PARSER_HPP

#include <tao/pegtl.hpp>
#include <boost/log/trivial.hpp>

namespace scindo
{
    namespace vcf
    {
        struct number
        {
            enum kind { given, A, G, R, dot };
            size_t value;
        };

        enum type { Integer, Float, Character, String, Flag };

        struct info
        {
            number num;
            type typ;
            std::string description;
            std::string source;
            std::string version;
        };

        struct format
        {
            std::string id;
            number num;
            type typ;
            std::string description;
            std::string source;
            std::string version;
        };

        struct filter
        {
            std::string id;
            std::string description;
            std::string source;
            std::string version;
        };

        namespace detail
        {
            struct state
            {
                std::string name;
                std::string value;
                std::string key;
                std::vector<std::pair<std::string,std::string>> values;
                std::vector<std::string> headers;
                std::vector<std::string> fields;

                std::string str;

                std::unique_ptr<tao::pegtl::position> pos_ptr;

                std::function<void(const std::string& p_name, const std::string& p_value)> unstructured;
                std::function<void(const std::string& p_name, const std::vector<std::pair<std::string,std::string>>& p_value)> structured;
                std::function<void(const std::vector<std::string>& p_headers)> header;
                std::function<void(const std::vector<std::string>& p_headers)> variant;
            };

            template<typename Rule>
            struct action
                : tao::pegtl::nothing<Rule>
            {};

            struct vcf_eol : tao::pegtl::ascii::eol {};

            template <>
            struct action<vcf_eol>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    if (!s.pos_ptr)
                    {
                        s.pos_ptr = std::unique_ptr<tao::pegtl::position>(new tao::pegtl::position(in.position()));
                    }
                    else
                    {
                        *s.pos_ptr = in.position();
                    }
                }
            };

            struct meta_name
                : tao::pegtl::ascii::identifier
            {};

            template <>
            struct action<meta_name>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.name = in.string();
                }
            };

            struct meta_unstructured_value
                : tao::pegtl::plus<tao::pegtl::ascii::not_one<'\r', '\n'>>
            {};

            template <>
            struct action<meta_unstructured_value>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.value = in.string();
                    if (s.unstructured)
                    {
                        s.unstructured(s.name, s.value);
                    }
                }
            };

            struct meta_key : tao::pegtl::ascii::identifier {};

            template <>
            struct action<meta_key>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.key = in.string();
                }
            };

            struct meta_key_value_name : tao::pegtl::ascii::identifier {};

            template <>
            struct action<meta_key_value_name>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.values.push_back(std::make_pair(s.key, in.string()));
                }
            };

            struct meta_key_value_number
                : tao::pegtl::sor<
                    tao::pegtl::ascii::one<'.'>,
                    tao::pegtl::plus<tao::pegtl::ascii::digit>
                >
            {};

            template <>
            struct action<meta_key_value_number>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.values.push_back(std::make_pair(s.key, in.string()));
                }
            };

            struct meta_plain_char : tao::pegtl::ascii::not_one<'\\', '"'> {};

            template <>
            struct action<meta_plain_char>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.str += in.string();
                }
            };

            struct meta_escaped_char : tao::pegtl::ascii::one<'\\', '"'> {};

            template <>
            struct action<meta_escaped_char>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.str += in.string();
                }
            };

            struct meta_escape_sequence
                : tao::pegtl::seq<
                    tao::pegtl::ascii::one<'\\'>,
                    meta_escaped_char>
            {};

            struct meta_quoted_string
                : tao::pegtl::seq<
                    tao::pegtl::ascii::one<'"'>,
                    tao::pegtl::star<
                        tao::pegtl::sor<
                            meta_plain_char,
                            meta_escape_sequence
                        >
                    >,
                    tao::pegtl::ascii::one<'"'>
                >
            {};

            struct meta_key_value_string : meta_quoted_string {};

            template <>
            struct action<meta_key_value_string>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.values.push_back(std::make_pair(s.key, s.str));
                    s.str.clear();
                }
            };

            struct meta_structured_value
                : tao::pegtl::seq<
                    tao::pegtl::one<'<'>,
                    tao::pegtl::list<
                        tao::pegtl::seq<
                            meta_key,
                            tao::pegtl::one<'='>,
                            tao::pegtl::sor<
                                meta_key_value_name,
                                meta_key_value_number,
                                meta_key_value_string
                            >
                        >,
                        tao::pegtl::one<','>
                    >,
                    tao::pegtl::one<'>'>
                >
            {};

            template <>
            struct action<meta_structured_value>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    if (s.structured)
                    {
                        s.structured(s.name, s.values);
                    }
                    s.values.clear();
                }
            };

            struct meta_line
                : tao::pegtl::seq<
                    tao::pegtl::ascii::string<'#', '#'>,
                    meta_name,
                    tao::pegtl::ascii::one<'='>,
                    tao::pegtl::sor<
                        meta_structured_value,
                        meta_unstructured_value
                    >,
                    vcf_eol
                >
            {};

            struct header_name
                : tao::pegtl::plus<tao::pegtl::ascii::not_one<'\t', '\r', '\n'>>
            {};

            template <>
            struct action<header_name>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.headers.push_back(in.string());
                }
            };

            struct header_line
                : tao::pegtl::seq<
                    tao::pegtl::ascii::one<'#'>,
                    tao::pegtl::list<
                        header_name,
                        tao::pegtl::ascii::one<'\t'>
                    >,
                    vcf_eol
                >
            {};

            template <>
            struct action<header_line>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    if (s.header)
                    {
                        s.header(s.headers);
                    }
                    s.headers.clear();
                }
            };

            struct headers
                : tao::pegtl::seq<
                    tao::pegtl::star<meta_line>,
                    header_line
                >
            {};

            struct data_line_field
                : tao::pegtl::plus<tao::pegtl::ascii::not_one<'\t', '\r', '\n'>>
            {};

            template <>
            struct action<data_line_field>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.fields.push_back(in.string());
                }
            };

            struct data_line
                : tao::pegtl::seq<
                    tao::pegtl::list<
                        data_line_field,
                        tao::pegtl::ascii::one<'\t'>
                    >,
                    vcf_eol
                >
            {};

            template <>
            struct action<data_line>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    if (s.variant)
                    {
                        s.variant(s.fields);
                    }
                    s.fields.clear();
                }
            };

            struct vcf_file
                : tao::pegtl::seq<
                    headers,
                    tao::pegtl::star<data_line>
                >
            {};
        }
        // namespace detail

        struct vcf_file
        {
            using position = tao::pegtl::position;

            std::istream& in;
            std::unordered_map<std::string, scindo::vcf::info> info;
            std::unordered_map<std::string, scindo::vcf::filter> filter;
            std::unordered_map<std::string, scindo::vcf::format> format;

            vcf_file(std::istream& p_in)
                : in(p_in)
            {
            }

        private:
            void structured_header(const std::string& p_name,
                                   const std::vector<std::pair<std::string, std::string>>& p_items)
            {
                std::unordered_map<std::string,std::string> items;
                std::unordered_map<std::string,size_t> idx;
                size_t n = 0;
                for (auto itr = p_items.begin(); itr != p_items.end(); ++itr, ++n)
                {
                    if (items.contains(itr->first))
                    {
                        BOOST_LOG_TRIVIAL(warning)
                            << pos().source << ":" << pos().line << ": "
                            << "for metadata item `" << p_name << "`, "
                            << "property `" << itr->first << "` occurred more than once.";
                    }
                    else
                    {
                        items[itr->first] = itr->second;
                        idx[itr->first] = n;
                    }
                }

                if (p_name == "INFO")
                {
                    if (check(idx, {"ID", "Number", "Type", "Description"}, {"Source", "Version"}))
                    {
                        const std::string& id = items.at("ID");
                        if (info.contains(id))
                        {
                            BOOST_LOG_TRIVIAL(warning)
                                << pos().source << ":" << pos().line << ": "
                                << "duplicate info `" << id << "`.";
                        }
                    }
                }
                if (p_name == "FORMAT")
                {
                    if (check(idx, {"ID", "Number", "Type", "Description"}, {"Source", "Version"}))
                    {
                        const std::string& id = items.at("ID");
                        if (info.contains(id))
                        {
                            BOOST_LOG_TRIVIAL(warning)
                                << pos().source << ":" << pos().line << ": "
                                << "duplicate info `" << id << "`.";
                        }
                    }
                }
            }

            bool check(const std::unordered_map<std::string,size_t>& p_idx,
                      std::initializer_list<const char*> p_required,
                      std::initializer_list<const char*> p_optional) const
            {
                bool ok = true;
                size_t i = 0;
                for (auto itr = p_required.begin(); itr != p_required.end(); ++itr, ++i)
                {
                    if (!p_idx.contains(*itr) || p_idx.at(*itr) != i)
                    {
                        if (ok)
                        {
                            BOOST_LOG_TRIVIAL(warning)
                                << pos().source << ":" << pos().line << ": "
                                << "required property `" << *itr << "` not found at position " << (i + 1) << ".";
                            ok = false;
                        }
                    }
                }
                return ok;
            }

            const position& pos() const
            {
                return *S.pos_ptr;
            }

            scindo::vcf::detail::state S;
        };
    }
    // namespace vcf
}
// namespace scindo

#endif // SCINDO_VCF_PARSER_HPP
