#ifndef SCINDO_LINES_OF_FILE_HPP
#define SCINDO_LINES_OF_FILE_HPP

#include <istream>
#include <thread>

#include <iostream>

#include "scindo/concurrent_deque.hpp"

namespace scindo
{
    namespace detail
    {
        struct block_of_lines : public std::deque<std::string>
        {
            block_of_lines() {}

            block_of_lines(std::istream& p_in, size_t p_n)
            {
                std::string line;
                while (size() < p_n)
                {
                    if (!std::getline(p_in, line))
                    {
                        return;
                    }
                    push_back(line);
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
            : m_in(p_in), m_deque(p_bufferSize), m_thread([this]() { read_lines(); })
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
            m_thread.join();
        }

    private:
        void read_lines()
        {
            const size_t N = 1000;
            while (true)
            {
                detail::block_of_lines blk(m_in, N);
                size_t z = blk.size();
                if (z > 0)
                {
                    m_deque.push_back(std::move(blk));
                }
                if (z < N)
                {
                    break;
                }
            }
            m_deque.end();
        }

        std::istream& m_in;
        scindo::concurrent_deque<detail::block_of_lines> m_deque;
        std::thread m_thread;
        detail::block_of_lines m_block;
    };
}
// namespace scindo

#endif // SCINDO_LINES_OF_FILE_HPP
