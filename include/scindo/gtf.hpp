#ifndef SCINDO_GTF_HPP
#define SCINDO_GTF_HPP

#include <istream>
#include <variant>
#include <vector>
#include <unordered_map>

#include <tao/pegtl.hpp>
#include "scindo/files.hpp"

namespace scindo
{
    namespace gtf
    {
        struct missing {};

        using label = std::string;
        using value = std::variant<std::string,int64_t,double,missing>;
        using attribute = std::pair<label,value>;

        using handler = std::function<void(const std::string& p_seqname,
                                           const std::string& p_source,
                                           const std::string& p_feature, 
                                           const int64_t& p_start,
                                           const int64_t& p_end,
                                           const std::optional<double>& p_score,
                                           const char& p_strand,
                                           const std::optional<int>& p_frame,
                                           const std::vector<scindo::gtf::attribute>& p_attrs)>;

        namespace detail
        {
            struct tab : tao::pegtl::ascii::one<'\t'> {};
            struct eol : tao::pegtl::ascii::eol {};
            
            struct word
                : tao::pegtl::plus<
                    tao::pegtl::ascii::not_one<'\t', '\r', '\n'>
                >
            {};

            struct integer
                : tao::pegtl::plus<tao::pegtl::ascii::digit>
            {};

            struct real
                : tao::pegtl::sor<
                    tao::pegtl::seq<
                        tao::pegtl::plus<tao::pegtl::ascii::digit>,
                        tao::pegtl::opt<
                            tao::pegtl::ascii::one<'.'>,
                            tao::pegtl::star<tao::pegtl::ascii::digit>
                        >
                    >,
                    tao::pegtl::seq<
                        tao::pegtl::star<tao::pegtl::ascii::digit>,
                        tao::pegtl::ascii::one<'.'>,
                        tao::pegtl::plus<tao::pegtl::ascii::digit>
                    >
                >
            {};

            struct missing : tao::pegtl::ascii::one<'.'> {};

            struct seqname : word {};
            struct source : word {};
            struct feature : word {};
            struct start : tao::pegtl::seq<integer> {};
            struct end : tao::pegtl::seq<integer> {};
            struct score : tao::pegtl::sor<real,missing> {};
            struct strand : tao::pegtl::ascii::one<'+','-'> {};
            struct frame : tao::pegtl::ascii::one<'0','1','2'> {};

            struct label : tao::pegtl::ascii::identifier {};

            struct quoted_string_content
                : tao::pegtl::star<tao::pegtl::ascii::not_one<'"'>>
            {};

            struct quoted_string
                : tao::pegtl::seq<
                    tao::pegtl::ascii::one<'"'>,
                    quoted_string_content,
                    tao::pegtl::ascii::one<'"'>
                >
            {};

            struct label_value : tao::pegtl::ascii::identifier {};

            struct attr_value
                : tao::pegtl::sor<
                    quoted_string,
                    integer,
                    real,
                    label_value
                >
            {};

            struct attr
                : tao::pegtl::seq<
                    label,
                    tao::pegtl::plus<tao::pegtl::ascii::one<' '>>,
                    attr_value
                >
            {};

            struct attr_separator
                : tao::pegtl::seq<
                    tao::pegtl::ascii::one<';'>,
                    tao::pegtl::star<tao::pegtl::ascii::one<' '>>
                >
            {};

            struct attrs
                : tao::pegtl::seq<
                    tao::pegtl::list<
                        attr,
                        attr_separator
                    >,
                    tao::pegtl::opt<attr_separator>
                >
            {};

            struct row
                : tao::pegtl::seq<
                    seqname,    tab,
                    source,     tab,
                    feature,    tab,
                    start,      tab,
                    end,        tab,
                    tao::pegtl::sor<score,missing>,      tab,
                    tao::pegtl::sor<strand,missing>,     tab,
                    tao::pegtl::sor<frame,missing>,      tab,
                    attrs
                >
            {};

            struct comment
                : tao::pegtl::seq<
                    tao::pegtl::ascii::one<'#'>,
                    tao::pegtl::star<
                        tao::pegtl::ascii::not_one<'\r', '\n'>
                    >,
                    eol
                >
            {};

            struct gtf
                : tao::pegtl::star<
                    tao::pegtl::sor<
                        comment,
                        tao::pegtl::seq<row, eol>
                    >
                >
            {};

            struct state
            {
                enum field_name { seqname, source, feature, start, end, score, strand, frame };

                std::unordered_map<field_name,scindo::gtf::value> fields;
                std::vector<scindo::gtf::attribute> attrs;

                std::vector<std::string> words;
                scindo::gtf::label label;
                scindo::gtf::value value;

                scindo::gtf::handler handler;
            };

            template<typename Rule>
            struct action
                : tao::pegtl::nothing<Rule>
            {};

            template <>
            struct action<missing>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.value = scindo::gtf::missing{};
                }
            };

            template <>
            struct action<integer>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.value = std::stoll(t);
                }
            };

            template <>
            struct action<real>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.value = std::stod(t);
                }
            };

            template <>
            struct action<seqname>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.fields[scindo::gtf::detail::state::seqname] = t;
                }
            };

            template <>
            struct action<source>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.fields[scindo::gtf::detail::state::source] = t;
                }
            };

            template <>
            struct action<feature>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.fields[scindo::gtf::detail::state::feature] = t;
                }
            };

            template <>
            struct action<start>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.fields[scindo::gtf::detail::state::start] = s.value;
                }
            };

            template <>
            struct action<end>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.fields[scindo::gtf::detail::state::end] = s.value;
                }
            };

            template <>
            struct action<score>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    if (in.string() != ".")
                    {
                        s.fields[scindo::gtf::detail::state::score] = s.value;
                    }
                }
            };

            template <>
            struct action<strand>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.fields[scindo::gtf::detail::state::strand] = t;
                }
            };

            template <>
            struct action<frame>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.fields[scindo::gtf::detail::state::frame] = t;
                }
            };

            template <>
            struct action<quoted_string_content>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.value = t;
                }
            };

            template <>
            struct action<label_value>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.value = t;
                }
            };

            template <>
            struct action<label>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const std::string& t = in.string();
                    s.label = t;
                }
            };

            template <>
            struct action<attr>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    s.attrs.push_back(std::make_pair(s.label, s.value));
                    s.label.clear();
                    s.value = scindo::gtf::missing{};
                }
            };

            template <>
            struct action<row>
            {
                template<typename ActionInput>
                static void apply(const ActionInput& in, state& s)
                {
                    const scindo::gtf::value& sn = s.fields.at(scindo::gtf::detail::state::seqname);
                    const scindo::gtf::value& src = s.fields.at(scindo::gtf::detail::state::source);
                    const scindo::gtf::value& ft = s.fields.at(scindo::gtf::detail::state::feature);
                    int64_t st = std::get<int64_t>(s.fields.at(scindo::gtf::detail::state::start));
                    int64_t en = std::get<int64_t>(s.fields.at(scindo::gtf::detail::state::end));
                    std::optional<double> sc;
                    if (s.fields.contains(scindo::gtf::detail::state::score))
                    {
                        sc = std::get<double>(s.fields.at(scindo::gtf::detail::state::score));
                    }
                    char sd = '.';
                    if (s.fields.contains(scindo::gtf::detail::state::strand))
                    {
                        sd = std::get<std::string>(s.fields.at(scindo::gtf::detail::state::strand))[0];
                    }
                    std::optional<int> fr;
                    if (s.fields.contains(scindo::gtf::detail::state::frame))
                    {
                        char c = std::get<std::string>(s.fields.at(scindo::gtf::detail::state::frame))[0];
                        fr = c - '0';
                    }

                    if (s.handler)
                    {
                        s.handler(std::get<std::string>(sn),
                                  std::get<std::string>(src),
                                  std::get<std::string>(ft),
                                  st, en, sc, sd, fr, s.attrs);
                    }
                    s.fields.clear();
                    s.attrs.clear();
                }
            };
        }
        // namespace detail

        struct gtf_file
        {
            const std::string filename;

            gtf_file(const std::string& p_filename)
                : filename(p_filename)
            {
            }

            template<typename X>
            void parse(X p_handler)
            {
                input_file_holder_ptr in_ptr = files::in(filename);
                return parse(**in_ptr, p_handler);
            }

            template<typename X>
            void parse(std::istream& p_input, X p_handler)
            {
                scindo::gtf::detail::state s;
                s.handler = p_handler;

                size_t byte_num = 0;
                size_t line_num = 0;
                std::string line;
                while (std::getline(p_input, line))
                {
                    const size_t line_len = line.size();
                    if (starts_with(line, '#'))
                    {
                        ++line_num;
                        byte_num += line_len;
                        continue;
                    }
                    ++line_num;
                    tao::pegtl::memory_input in(line, filename);
                    if (!tao::pegtl::parse<scindo::gtf::detail::row, scindo::gtf::detail::action>(in, s))
                    {
                        throw std::runtime_error("error parsing gtf");
                    }
                    byte_num += line_len;
                    
                }
            }

            static bool starts_with(const std::string& p_str, char p_ch)
            {
                return p_str.size() > 0 && p_str.front() == p_ch;
            }

        };
    }
    // namespace gtf
}
// namespace scindo

#endif // SCINDO_GTF_HPP
