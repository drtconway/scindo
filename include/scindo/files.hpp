#ifndef SCINDO_FILES_HPP
#define SCINDO_FILES_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <sys/stat.h>

namespace scindo
{
    class input_file_holder
    {
    public:
        virtual std::istream& operator*() = 0;

        virtual ~input_file_holder() {}
    };
    typedef std::shared_ptr<input_file_holder> input_file_holder_ptr;

    class output_file_holder
    {
    public:
        virtual std::ostream& operator*() = 0;

        virtual ~output_file_holder() {}
    };
    typedef std::shared_ptr<output_file_holder> output_file_holder_ptr;

    class files_impl
    {
    public:
        virtual input_file_holder_ptr in(const std::string& p_name) = 0;

        virtual output_file_holder_ptr out(const std::string& p_name) = 0;

        virtual bool exists(const std::string& p_name) = 0;

        virtual void remove(const std::string& p_name) = 0;

        static bool ends_with(const std::string& p_str, const std::string& p_suffix)
        {
            auto s = p_str.rbegin();
            auto t = p_suffix.rbegin();
            while (s != p_str.rend() && t != p_suffix.rend())
            {
                if (*s != *t)
                {
                    return false;
                }
                ++s;
                ++t;
            }
            return t == p_suffix.rend();
        }
    };

    class tmp_file
    {
    public:
        tmp_file(files_impl& p_impl, const std::string& p_name)
            : m_impl(p_impl), m_name(p_name)
        {
        }

        input_file_holder_ptr in()
        {
            return m_impl.in(m_name);
        }

        output_file_holder_ptr out()
        {
            return m_impl.out(m_name);
        }

        void remove()
        {
            m_impl.remove(m_name);
        }

        ~tmp_file()
        {
            m_impl.remove(m_name);
        }

    private:
        files_impl& m_impl;
        const std::string m_name;
    };
    typedef std::shared_ptr<tmp_file> tmp_file_ptr;

    namespace detail
    {
        class cin_file_holder : public input_file_holder
        {
        public:
            std::istream& operator*()
            {
                return std::cin;
            };
        };

        class cout_file_holder : public output_file_holder
        {
        public:
            std::ostream& operator*()
            {
                return std::cout;
            };
        };

        class plain_input_file_holder : public input_file_holder
        {
        public:
            plain_input_file_holder(const std::string& p_name)
                : m_file(p_name)
            {
            }

            std::istream& operator*()
            {
                return m_file;
            }

        private:
            std::ifstream m_file;
        };

        class plain_output_file_holder : public output_file_holder
        {
        public:
            plain_output_file_holder(const std::string& p_name)
                : m_file(p_name)
            {
            }

            std::ostream& operator*()
            {
                return m_file;
            }

        private:
            std::ofstream m_file;
        };

        class gzip_input_file_holder : public input_file_holder
        {
        public:
            gzip_input_file_holder(const std::string& p_name)
                : m_file(p_name, std::ios_base::in | std::ios_base::binary)
            {
                m_filter.push(boost::iostreams::gzip_decompressor(boost::iostreams::gzip::default_window_bits, 1024*1024));
                m_filter.push(m_file);
            }

            std::istream& operator*()
            {
                return m_filter;
            }

        private:
            std::ifstream m_file;
            boost::iostreams::filtering_istream m_filter;
        };

        class gzip_output_file_holder : public output_file_holder
        {
        public:
            gzip_output_file_holder(const std::string& p_name)
                : m_file(p_name, std::ios_base::out | std::ios_base::binary)
            {
                m_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip::default_compression, 1024*1024));
                m_filter.push(m_file);
            }

            std::ostream& operator*()
            {
                return m_filter;
            }

            ~gzip_output_file_holder()
            {
                m_filter.flush();
                m_filter.pop();
            }

        private:
            std::ofstream m_file;
            boost::iostreams::filtering_ostream m_filter;
        };

        class plain_input_string_file_holder : public input_file_holder
        {
        public:
            plain_input_string_file_holder(const std::string& p_content)
                : m_file(p_content)
            {
            }

            std::istream& operator*()
            {
                return m_file;
            }

        private:
            std::istringstream m_file;
        };

        class plain_output_string_file_holder : public output_file_holder
        {
        public:
            plain_output_string_file_holder(std::map<std::string,std::string>& p_files, const std::string& p_name)
                : m_files(p_files), m_name(p_name)
            {
            }

            std::ostream& operator*()
            {
                return m_file;
            }

            ~plain_output_string_file_holder()
            {
                m_files[m_name] = m_file.str();
            }

        private:
            std::map<std::string,std::string>& m_files;
            const std::string m_name;
            std::ostringstream m_file;
        };

        class gzip_input_string_file_holder : public input_file_holder
        {
        public:
            gzip_input_string_file_holder(const std::string& p_content)
                : m_file(p_content, std::ios_base::in | std::ios_base::binary)
            {
                m_filter.push(boost::iostreams::gzip_decompressor());
                m_filter.push(m_file);
            }

            std::istream& operator*()
            {
                return m_filter;
            }

        private:
            std::istringstream m_file;
            boost::iostreams::filtering_istream m_filter;
        };

        class gzip_output_string_file_holder : public output_file_holder
        {
        public:
            gzip_output_string_file_holder(std::map<std::string,std::string>& p_files, const std::string& p_name)
                : m_files(p_files), m_name(p_name), m_file(std::ios_base::out | std::ios_base::binary)
            {
                m_filter.push(boost::iostreams::gzip_compressor());
                m_filter.push(m_file);
            }

            std::ostream& operator*()
            {
                return m_filter;
            }

            ~gzip_output_string_file_holder()
            {
                m_filter.flush();
                m_filter.pop();
                m_files[m_name] = m_file.str();
            }

        private:
            std::map<std::string,std::string>& m_files;
            const std::string m_name;
            std::ostringstream m_file;
            boost::iostreams::filtering_ostream m_filter;
        };

        class basic_files_impl : public files_impl
        {
        public:
            virtual input_file_holder_ptr in(const std::string& p_name)
            {
                if (p_name == "-")
                {
                    return input_file_holder_ptr(new detail::cin_file_holder);
                }
                if (ends_with(p_name, ".gz"))
                {
                    return input_file_holder_ptr(new detail::gzip_input_file_holder(p_name));
                }
                return input_file_holder_ptr(new detail::plain_input_file_holder(p_name));
            }

            virtual output_file_holder_ptr out(const std::string& p_name)
            {
                if (p_name == "-")
                {
                    return output_file_holder_ptr(new detail::cout_file_holder);
                }
                if (ends_with(p_name, ".gz"))
                {
                    return output_file_holder_ptr(new detail::gzip_output_file_holder(p_name));
                }
                return output_file_holder_ptr(new detail::plain_output_file_holder(p_name));
            }

            virtual bool exists(const std::string& p_name)
            {
                struct stat buf;
                return (stat(p_name.c_str(), &buf) == 0);
            }

            virtual void remove(const std::string& p_name)
            {
                boost::filesystem::remove(p_name);
            }
        };

        class string_files_impl : public files_impl
        {
        public:
            string_files_impl(std::map<std::string,std::string>& p_files)
                : m_files(p_files)
            {
            }

            virtual input_file_holder_ptr in(const std::string& p_name)
            {
                if (!exists(p_name))
                {
                    throw std::runtime_error("no such string file: " + p_name);
                }
                const std::string& content = m_files.find(p_name)->second;
                if (ends_with(p_name, ".gz"))
                {
                    return input_file_holder_ptr(new detail::gzip_input_string_file_holder(content));
                }
                return input_file_holder_ptr(new detail::plain_input_string_file_holder(content));
            }

            virtual output_file_holder_ptr out(const std::string& p_name)
            {
                if (ends_with(p_name, ".gz"))
                {
                    return output_file_holder_ptr(new detail::gzip_output_string_file_holder(m_files, p_name));
                }
                return output_file_holder_ptr(new detail::plain_output_string_file_holder(m_files, p_name));
            }

            virtual bool exists(const std::string& p_name)
            {
                return m_files.find(p_name) != m_files.end();
            }

            virtual void remove(const std::string& p_name)
            {
                m_files.erase(p_name);
            }

        private:
            std::map<std::string,std::string>& m_files;
        };
    }
    // namespace detail

    struct files
    {
        static input_file_holder_ptr in(const std::string& p_name)
        {
            return impl().in(p_name);
        }

        static output_file_holder_ptr out(const std::string& p_name)
        {
            return impl().out(p_name);
        }

        static tmp_file_ptr tmp(std::string p_suffix = "")
        {
            std::string name = tmp_name();
            return tmp_file_ptr(new tmp_file(impl(), name));
        }

        static bool exists(const std::string& p_name)
        {
            return impl().exists(p_name);
        }

        static void remove(const std::string& p_name)
        {
            return impl().remove(p_name);
        }

        static void file_impl()
        {
                auto p = std::unique_ptr<files_impl>(new detail::basic_files_impl);
                impl(std::move(p));
        }

        static void string_impl(std::map<std::string,std::string>& p_files)
        {
                auto p = std::unique_ptr<files_impl>(new detail::string_files_impl(p_files));
                impl(std::move(p));
        }

        static files_impl& impl(std::unique_ptr<files_impl> p_impl = std::unique_ptr<files_impl>())
        {
            static std::unique_ptr<files_impl> impl_ptr;
            if (p_impl)
            {
                impl_ptr = std::move(p_impl);
            }
            if (!impl_ptr)
            {
                impl_ptr = std::unique_ptr<files_impl>(new detail::basic_files_impl);
            }
            return *impl_ptr;
        }

        static std::string tmp_name(std::string p_suffix = "")
        {
            boost::uuids::random_generator gen;
            boost::uuids::uuid id = gen();
            std::ostringstream out;
            out << "/tmp/" << id << p_suffix;
            return out.str();
        }
    };
}
// namespace scindo

#endif // SCINDO_FILES_HPP
