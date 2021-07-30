#ifndef SCINDO_GZIP_HPP
#define SCINDO_GZIP_HPP

#include <cstdio>
#include <zlib.h>
#include <nlohmann/json.hpp>

namespace scindo
{
    namespace detail
    {
        struct buffer
        {
            template <typename X>
            static
            typename std::enable_if<std::is_convertible<X,std::function<void(std::string::const_iterator,std::string::const_iterator,bool&)>>::value,bool>::type
            lines(std::string& buf, std::string& tmp, X p_acceptor)
            {
                size_t pos = 0;
                while (pos < buf.size())
                {
                    auto nxt = buf.find('\n', pos);
                    if (nxt == std::string::npos)
                    {
                        break;
                    }
                    bool stop = false;
                    p_acceptor(buf.begin() + pos, buf.begin() + nxt, stop);
                    if (stop)
                    {
                        return true;
                    }
                    pos = nxt + 1;
                }
                if (pos > 0)
                {
                    tmp.insert(tmp.end(), buf.begin() + pos, buf.end());
                    buf.swap(tmp);
                    tmp.clear();
                }
                return false;
            }

            template <typename X>
            static
            typename std::enable_if<std::is_convertible<X,std::function<void(const std::pair<std::string::const_iterator,std::string::const_iterator>&)>>::value,std::string::iterator>::type
            lines(std::string::iterator p_begin, std::string::iterator p_end, X p_acceptor)
            {
                size_t len = p_end - p_begin;
                auto itr = p_begin;
                while (itr != p_end)
                {
                    auto jtr = find(itr, p_end, '\n');
                    if (jtr == p_end)
                    {
                        break;
                    }
                    p_acceptor(std::pair<std::string::const_iterator,std::string::const_iterator>(itr, jtr));
                    itr = jtr + 1;
                }
                return itr;
            }
        };

    }
    // namespace detail

    struct gzip
    {
        struct reader
        {
            using span = std::pair<std::string::const_iterator,std::string::const_iterator>;

            const std::string name;
            gzFile fp;

            reader(const std::string& p_name, size_t p_in_buffer_size = 65536)
                : name(p_name), fp(gzopen(name.c_str(), "rb"))
            {
                if (!fp)
                {
                    throw std::runtime_error("unable to open file");
                }
                gzbuffer(fp, p_in_buffer_size);
            }

            ~reader()
            {
                int err = gzclose(fp);
                if (err)
                {
                    std::cerr << "gzclose returned error code: " << err << std::endl;
                }
            }

            size_t next(std::string& buf, size_t off = 0)
            {
                int n = gzread(fp, reinterpret_cast<void*>(buf.data() + off), buf.size() - off);
                if (n < 0)
                {
                    throw std::runtime_error("gzread failed");
                }
                return n;
            }

            std::string::iterator next(std::string::iterator p_begin, std::string::iterator p_end)
            {
                size_t l = p_end - p_begin;
                int n = gzread(fp, reinterpret_cast<void*>(&*p_begin), l);
                if (n < 0)
                {
                    throw std::runtime_error("gzread failed");
                }
                return p_begin + n;
            }
        };

        struct writer
        {
            const std::string name;
            gzFile fp;

            writer(const std::string& p_name, size_t p_out_buffer_size = 65536)
                : name(p_name), fp(gzopen(name.c_str(), "wb"))
            {
                if (!fp)
                {
                    throw std::runtime_error("unable to open file");
                }
                gzbuffer(fp, p_out_buffer_size);
            }

            ~writer()
            {
                int err = gzclose(fp);
                if (err)
                {
                    std::cerr << "gzclose returned error code: " << err << std::endl;
                }
            }

            template <typename Itr>
            size_t write(Itr p_begin, Itr p_end)
            {
                constexpr size_t elem_bytes = sizeof(typename std::iterator_traits<Itr>::value_type);
                size_t num_elements = p_end - p_begin;
                size_t num_bytes = elem_bytes * num_elements;
                size_t n = gzwrite(fp, reinterpret_cast<const void*>(&*p_begin), num_bytes);
                if (n != num_bytes)
                {
                    throw std::runtime_error("gzwrite failed");
                }
                return n;
            }
        };

        struct lines
        {
            reader rdr;
            std::string buf;
            std::deque<reader::span> spans;
            std::string::iterator itr;
            std::string::iterator jtr;

            lines(const std::string& p_name, size_t p_in_buffer_size = 65536, size_t p_out_buffer_size = 65536)
                : rdr(p_name, p_in_buffer_size)
            {
                buf.resize(p_out_buffer_size);
                itr = buf.begin();
                jtr = buf.end();
                get_more();
            }

            bool more() const
            {
                return spans.size() > 0;
            }

            const reader::span& operator*() const
            {
                return spans.front();
            }

            lines& operator++()
            {
                spans.pop_front();
                if (spans.size() == 0)
                {
                    get_more();
                }
                return *this;
            }

            void get_more()
            {
                if (jtr != buf.end())
                {
                    return;
                }
                if (itr != buf.begin())
                {
                    itr = std::copy(itr, buf.end(), buf.begin());
                }
                jtr = rdr.next(itr, buf.end());
                itr = detail::buffer::lines(buf.begin(), jtr, [&](const reader::span& p_span) {
                    spans.push_back(p_span);
                });
            }
        };

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(const reader::span&)>>::value,void>::type
        with_lines(const std::string& p_name, X p_acceptor)
        {
            reader rdr(p_name);
            std::string buf;
            buf.resize(1024*1024);
            auto itr = buf.begin();
            auto jtr = buf.end();
            while (jtr == buf.end())
            {
                if (itr != buf.begin())
                {
                    itr = std::copy(itr, jtr, buf.begin());
                }
                jtr = rdr.next(itr, buf.end());
                itr = detail::buffer::lines(buf.begin(), jtr, p_acceptor);
            }
            if (itr != buf.begin())
            {
                itr = std::copy(itr, jtr, buf.begin());
            }
            if (itr != buf.begin())
            {
                p_acceptor(reader::span(buf.begin(), itr));
            }
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(const reader::span&,const reader::span&,bool&)>>::value,void>::type
        with_lines(const std::string& p_lhs_name, const std::string& p_rhs_name, X p_acceptor,
                   size_t p_in_buffer_size = 1024*1024, size_t p_out_buffer_size = 1024*1024)
        {
            lines lhs_rdr (p_lhs_name, p_in_buffer_size >> 1, p_out_buffer_size >> 1);
            lines rhs_rdr (p_rhs_name, p_in_buffer_size >> 1, p_out_buffer_size >> 1);
            bool stop = false;
            while (lhs_rdr.more() && rhs_rdr.more())
            {
                p_acceptor(*lhs_rdr, *rhs_rdr, stop);
                if (stop)
                {
                    return;
                }
                ++lhs_rdr;
                ++rhs_rdr;
            }
            if (lhs_rdr.more() != rhs_rdr.more())
            {
                throw std::runtime_error("files had unequal numbers of lines.");
            }
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(std::string::const_iterator,std::string::const_iterator)>>::value,void>::type
        with_lines(const std::string& p_name, X p_acceptor)
        {
            with_lines(p_name, [&](auto beg, auto end, bool& stop) { p_acceptor(beg, end); });
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(std::string::const_iterator,std::string::const_iterator,bool&)>>::value,void>::type
        with_lines(const std::string& p_name, X p_acceptor)
        {
            std::string buf;
            std::string tmp;
            bool stopped = false;
            with_blocks(p_name, [&](auto beg, auto end, bool& stop) {
                buf.insert(buf.end(), beg, end);
                stop = detail::buffer::lines(buf, tmp, p_acceptor);
                stopped = stop;
            });
            if (!stopped && buf.size() > 0)
            {
                bool stop = false;
                p_acceptor(buf.begin(), buf.end(), stop);
            }
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(std::string::const_iterator,std::string::const_iterator)>>::value,void>::type
        with_blocks(const std::string& p_name, X p_acceptor)
        {
            with_blocks(p_name, [&](auto beg, auto end, bool& stop) { p_acceptor(beg, end); });
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(std::string::const_iterator,std::string::const_iterator,bool&)>>::value,void>::type
        with_blocks(const std::string& p_name, X p_acceptor)
        {
            static const size_t out_buf_size = 1ULL << 20;

            std::string out;
            out.resize(out_buf_size);

            int errn = 0;
            gzFile f = gzopen(p_name.c_str(), "rb");
            if (!f)
            {
                throw std::runtime_error("unable to open file");
            }
            gzbuffer(f, 256*1024*1024);

            while (true)
            {
                int n = gzread(f, reinterpret_cast<void*>(out.data()), out.size());
                if (n < 0)
                {
                    std::cerr << "gzread returned error code: " << n << std::endl;
                    std::cerr << gzerror(f, &errn) << std::endl;
                    break;
                }
                if (n > 0)
                {
                    bool stop = false;
                    p_acceptor(out.begin(), out.begin() + n, stop);
                    if (stop)
                    {
                        break;
                    }
                }
                if (n < out.size())
                {
                    break;
                }
            }
            int err = gzclose(f);
            if (err < 0)
            {
                std::cerr << "gzclose returned error code: " << err << std::endl;
                std::cerr << gzerror(f, &errn) << std::endl;
            }
        }
    };
}
// namespace scindo

#endif // SCINDO_GZIP_HPP
