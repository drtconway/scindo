#ifndef SCINDO_LINES_OF_FILE_HPP
#define SCINDO_LINES_OF_FILE_HPP

#include <istream>
#include <thread>

#include <iostream>

#include "scindo/concurrent_deque.hpp"
#include "scindo/gzip.hpp"

namespace scindo
{
    namespace detail
    {
        struct block_of_lines
        {
            static void make(std::istream& p_in, size_t p_n, std::deque<std::string>& p_res)
            {
                std::string line;
                while (p_res.size() < p_n)
                {
                    if (!std::getline(p_in, line))
                    {
                        return;
                    }
                    p_res.push_back(std::move(line));
                    line.clear();
                }
            }
        };
    }
    // namespace detail

    template <bool BG = false>
    class lines_of_file;

    template <>
    class lines_of_file<false>
    {
    public:
        lines_of_file(std::istream& p_in)
            : m_in(p_in)
        {
        }

        bool next(std::string& p_line)
        {
            if (std::getline(m_in, p_line))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

    private:
        std::istream& m_in;
    };

    template <>
    class lines_of_file<true>
    {
    public:
        lines_of_file(std::istream& p_in, size_t p_bufferSize = 100)
            : m_deque(p_bufferSize), m_thread([&]() { read_lines(p_in); })
        {
        }

        lines_of_file(std::string& p_filename, size_t p_bufferSize = 100)
            : m_deque(p_bufferSize), m_thread([&]() { read_lines(p_filename); })
        {
        }

        bool next(std::string& p_line)
        {
            while (m_block.size() == 0)
            {
                if (!m_deque.pop_front(m_block))
                {
                    return false;
                }
            }
            std::swap(p_line, m_block.front());
            m_block.pop_front();
            return true;
        }

        ~lines_of_file()
        {
            m_deque.end();
            m_thread.join();
        }

    private:
        void read_lines(const std::string& p_filename)
        {
            if (ends_with(p_filename, ".gz"))
            {
                const size_t N = 1000;
                std::deque<std::string> blk;
                size_t bn = 0;
                gzip::with_lines(p_filename, [&](auto beg, auto end, bool& stop) {
                    blk.push_back(std::string(beg, end));
                    if (blk.size() == N)
                    {
                        bn += 1;
                        if (!m_deque.push_back(std::move(blk)))
                        {
                            blk.clear();
                            stop = true;
                            return;
                        }
                        blk.clear();
                    }
                });
                if (blk.size() > 0)
                {
                    m_deque.push_back(std::move(blk));
                }
                m_deque.end();
            }
            else
            {
                input_file_holder_ptr inp = files::in(p_filename);
                read_lines(**inp);
            }
        }

        void read_lines(std::istream& p_in)
        {
            const size_t N = 1000;
            while (true)
            {
                std::deque<std::string> blk;
                detail::block_of_lines::make(p_in, N, blk);
                size_t z = blk.size();
                if (z > 0)
                {
                    if (!m_deque.push_back(std::move(blk)))
                    {
                        break;
                    }
                }
                if (z < N)
                {
                    break;
                }
            }
            m_deque.end();
        }

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
        scindo::concurrent_deque<std::deque<std::string>> m_deque;
        std::thread m_thread;
        std::deque<std::string> m_block;
    };
}
// namespace scindo

#endif // SCINDO_LINES_OF_FILE_HPP
