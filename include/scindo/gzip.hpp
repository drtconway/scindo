#ifndef SCINDO_GZIP_HPP
#define SCINDO_GZIP_HPP

#include <cstdio>
#include <zlib.h>
#include <nlohmann/json.hpp>

namespace scindo
{
    struct gzip
    {
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
            size_t pos = 0;
            bool stopped = false;
            with_blocks(p_name, [&](auto beg, auto end, bool& stop) {
                buf.insert(buf.end(), beg, end);
                while (pos < buf.size())
                {
                    auto nxt = buf.find('\n', pos);
                    if (nxt == std::string::npos)
                    {
                        break;
                    }
                    p_acceptor(buf.begin() + pos, buf.begin() + nxt, stop);
                    if (stop)
                    {
                        stopped = false;
                        return;
                    }
                    pos = nxt + 1;
                }
                if (pos > 0)
                {
                    tmp.insert(tmp.end(), buf.begin() + pos, buf.end());
                    buf.swap(tmp);
                    tmp.clear();
                    pos = 0;
                }
            });
            if (!stopped && buf.size() > 0)
            {
                bool stop = false;
                p_acceptor(buf.begin() + pos, buf.end(), stop);
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
            gzbuffer(f, 1024*1024);

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
