#ifndef SCINDO_FASTQ_HPP
#define SCINDO_FASTQ_HPP

#include <istream>
#include <string>
#include <tuple>
#include <boost/algorithm/string/trim.hpp>

#include "scindo/lines_of_file.hpp"

namespace scindo
{
    using fastq_read = std::tuple<std::string,std::string,std::string,std::string>;
    using fastq_text = std::pair<std::string::const_iterator,std::string::const_iterator>;
    using fastq_tuple = std::tuple<fastq_text,fastq_text,fastq_text,fastq_text>;

    namespace detail
    {
        struct text_helpers {
            static bool isspace(char c) {
                switch (c) {
                    case ' ':
                    case '\t':
                    case '\r':
                    case '\n':
                    case '\f':
                    case '\v': {
                        return true;
                    }
                    default: {
                        return false;
                    }
                }
            }

            static void trim(std::string& p_str)
            {
                if (p_str.size() > 0 && isspace(p_str[0])) {
                    boost::algorithm::trim_if(p_str, isspace);
                    return;
                }
                
                // Use a simplified algorithm when there is no space at the front.
                //
                while (p_str.size() > 0 && isspace(p_str.back())) {
                    p_str.pop_back();
                }
            }
        };
    }
    // namespace detail

    class fastq_reader
    {
    public:
        typedef fastq_read item_type;
        typedef item_type const& item_result_type;

        fastq_reader(std::string& p_name, size_t p_bufferSize = 100)
            : m_lines(p_name, p_bufferSize), m_more(true)
        {
            next();
        }

        fastq_reader(std::istream& p_in, size_t p_bufferSize = 100)
            : m_lines(p_in, p_bufferSize), m_more(true)
        {
            next();
        }

        bool more() const
        {
            return m_more;
        }

        item_result_type operator*() const
        {
            return m_curr;
        }

        fastq_reader& operator++()
        {
            next();
            return *this;
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(const fastq_tuple&)>>::value,void>::type
        with(const std::string& p_name, X p_acceptor)
        {
            std::string id;
            std::string seq;
            std::string id2;
            std::string qal;
            std::vector<std::string*> r0(4, NULL);
            r0[0] = &id;
            r0[1] = &seq;
            r0[2] = &id2;
            r0[3] = &qal;

            size_t n = 0;
            gzip::with_lines(p_name, [&](auto beg, auto end) {
                std::string& s = r0[n];
                s.clear();
                s.insert(s.end(), beg, end);
                n += 1;

                if (n == 4)
                {
                    for (size_t i = 0; i < 4; ++i)
                    {
                        boost::algorithm::trim(*r0[i]);
                    }

                    if (!starts_with(id, '@'))
                    {
                        std::cerr << id << std::endl;
                        throw std::domain_error("malformed read identifier line.");
                    }

                    if (!starts_with(id2, '+'))
                    {
                        std::cerr << id2 << std::endl;
                        throw std::domain_error("malformed second read identifier line.");
                    }

                    fastq_tuple tup(fastq_text(id.begin() + 1, id.end()),
                                    fastq_text(seq.begin(), seq.end()),
                                    fastq_text(id2.begin() + 1, id2.end()),
                                    fastq_text(qal.begin(), qal.end()));
                    p_acceptor(tup);

                    n = 0;
                }
            });
            if (n != 0)
            {
                throw std::domain_error("premature end of file");
            }
        }

        template <typename X>
        static
        typename std::enable_if<std::is_convertible<X,std::function<void(const fastq_tuple&,const fastq_tuple&,bool&)>>::value,void>::type
        with(const std::string& p_lhs_name, const std::string& p_rhs_name, X p_acceptor)
        {
            std::string lhs_id;
            std::string lhs_seq;
            std::string lhs_id2;
            std::string lhs_qal;
            std::vector<std::string*> lhs_r0(4, NULL);
            lhs_r0[0] = &lhs_id;
            lhs_r0[1] = &lhs_seq;
            lhs_r0[2] = &lhs_id2;
            lhs_r0[3] = &lhs_qal;

            std::string rhs_id;
            std::string rhs_seq;
            std::string rhs_id2;
            std::string rhs_qal;
            std::vector<std::string*> rhs_r0(4, NULL);
            rhs_r0[0] = &rhs_id;
            rhs_r0[1] = &rhs_seq;
            rhs_r0[2] = &rhs_id2;
            rhs_r0[3] = &rhs_qal;

            size_t n = 0;
            gzip::with_lines(p_lhs_name, p_rhs_name, [&](auto lhs_span, auto rhs_span, bool& stop) {
                std::string& lhs = *lhs_r0[n];
                lhs.clear();
                lhs.insert(lhs.end(), lhs_span.first, lhs_span.second);

                std::string& rhs = *rhs_r0[n];
                rhs.clear();
                rhs.insert(rhs.end(), rhs_span.first, rhs_span.second);

                n += 1;

                if (n == 4)
                {
                    n = 0;
                    for (size_t i = 0; i < 4; ++i)
                    {
                        detail::text_helpers::trim(*lhs_r0[i]);
                        detail::text_helpers::trim(*rhs_r0[i]);
                    }

                    if (!starts_with(lhs_id, '@'))
                    {
                        std::cerr << lhs_id << std::endl;
                        throw std::domain_error("malformed read identifier line.");
                    }
                    if (!starts_with(rhs_id, '@'))
                    {
                        std::cerr << rhs_id << std::endl;
                        throw std::domain_error("malformed read identifier line.");
                    }

                    if (!starts_with(lhs_id2, '+'))
                    {
                        std::cerr << lhs_id2 << std::endl;
                        throw std::domain_error("malformed second read identifier line.");
                    }
                    if (!starts_with(rhs_id2, '+'))
                    {
                        std::cerr << rhs_id2 << std::endl;
                        throw std::domain_error("malformed second read identifier line.");
                    }

                    fastq_tuple lhs_tup(
                        fastq_text(lhs_id.begin() + 1, lhs_id.end()),
                        fastq_text(lhs_seq.begin(), lhs_seq.end()),
                        fastq_text(lhs_id2.begin() + 1, lhs_id2.end()),
                        fastq_text(lhs_qal.begin(), lhs_qal.end()));

                    fastq_tuple rhs_tup(
                        fastq_text(rhs_id.begin() + 1, rhs_id.end()),
                        fastq_text(rhs_seq.begin(), rhs_seq.end()),
                        fastq_text(rhs_id2.begin() + 1, rhs_id2.end()),
                        fastq_text(rhs_qal.begin(), rhs_qal.end()));

                    p_acceptor(lhs_tup, rhs_tup, stop);
                    if (stop)
                    {
                        return;
                    }
                }
            }, 128*1024*1024, 128*1024*1024);
            if (n != 0)
            {
                throw std::domain_error("premature end of file");
            }
        }

    private:
        void next()
        {
            if (!m_lines.next(m_next))
            {
                m_more = false;
                return;
            }
            if (!starts_with(m_next, '@'))
            {
                std::cerr << m_next << std::endl;
                throw std::domain_error("missing read identifier line.");
            }

            detail::text_helpers::trim(m_next);
            std::string& id1 = std::get<0>(m_curr);
            id1.clear();
            id1.insert(id1.end(), m_next.begin() + 1, m_next.end()); // drop the @

            if (!m_lines.next(m_next))
            {
                throw std::domain_error("missing read sequence line.");
            }

            detail::text_helpers::trim(m_next);
            std::string& seq = std::get<1>(m_curr);
            seq.clear();
            seq.insert(seq.end(), m_next.begin(), m_next.end());

            if (!m_lines.next(m_next) || !starts_with(m_next, '+'))
            {
                throw std::domain_error("missing read spacer line.");
            }

            detail::text_helpers::trim(m_next);
            std::string& id2 = std::get<2>(m_curr);
            id2.clear();
            id2.insert(id2.end(), m_next.begin() + 1, m_next.end()); // drop the +

            if (!m_lines.next(m_next))
            {
                throw std::domain_error("missing read quality score line.");
            }

            detail::text_helpers::trim(m_next);
            std::string& qual = std::get<3>(m_curr);
            qual.clear();
            qual.insert(qual.end(), m_next.begin(), m_next.end());
        }

        static bool starts_with(const std::string& p_str, char p_ch)
        {
            return p_str.size() > 0 && p_str.front() == p_ch;
        }

        scindo::lines_of_file<true> m_lines;
        bool m_more;
        std::string m_next;
        fastq_read m_curr;
    };

    class fastq_writer
    {
    private:
        struct output_buffer
        {
            std::string buf;
            std::string::iterator cur;

            output_buffer(size_t p_buffer_size)
            {
                buf.resize(p_buffer_size);
                cur = buf.begin();
            }

            std::string::const_iterator append(std::string::const_iterator p_begin, std::string::const_iterator p_end)
            {
                size_t src_len = p_end - p_begin;
                size_t free_space = buf.end() - cur;
                size_t amount_to_copy = std::min(src_len, free_space);
                cur = std::copy(p_begin, p_begin + amount_to_copy, cur);
                return p_begin + amount_to_copy;
            }

            template <typename X>
            void flush(X p_acceptor)
            {
                p_acceptor(buf.begin(), cur);
                cur = buf.begin();
            }
        };
        using output_buffer_ptr = std::shared_ptr<output_buffer>;

    public:
        fastq_writer(std::ostream& p_out, size_t p_queue_size = 5, size_t p_in_buffer_size = 4*1024*1024)
            : m_spare_buffers(p_queue_size), m_full_buffers(p_queue_size), m_thread([&]() { write_strings(p_out); })
        {
            size_t buf_size = p_in_buffer_size / p_queue_size;
            for (size_t i = 0; i < p_queue_size; ++i)
            {
                m_spare_buffers.push_back(output_buffer_ptr(new output_buffer(buf_size)));
            }
        }

        fastq_writer(const std::string& p_out_name, size_t p_queue_size = 4, size_t p_in_buffer_size = 64*1024*1024, size_t p_out_buffer_size = 8*1024*1024)
            : m_spare_buffers(p_queue_size), m_full_buffers(p_queue_size), m_thread([&]() { write_strings(p_out_name, p_out_buffer_size); })
        {
            size_t buf_size = p_in_buffer_size / p_queue_size;
            for (size_t i = 0; i < p_queue_size; ++i)
            {
                m_spare_buffers.push_back(output_buffer_ptr(new output_buffer(buf_size)));
            }
        }

        ~fastq_writer()
        {
            m_spare_buffers.end();
            m_full_buffers.end();
            m_thread.join();
        }

        void write(const fastq_tuple& p_read)
        {
            static const std::string endl("\n");
            auto r0 = std::get<0>(p_read);
            write(r0.first, r0.second);
            write(endl.begin(), endl.end());
            auto r1 = std::get<1>(p_read);
            write(r1.first, r1.second);
            write(endl.begin(), endl.end());
            auto r2 = std::get<2>(p_read);
            write(r2.first, r2.second);
            write(endl.begin(), endl.end());
            auto r3 = std::get<3>(p_read);
            write(r3.first, r3.second);
            write(endl.begin(), endl.end());
        }

        void write(const fastq_read& p_read)
        {
            static const std::string endl("\n");
            auto r0 = std::get<0>(p_read);
            write(r0.begin(), r0.end());
            write(endl.begin(), endl.end());
            auto r1 = std::get<1>(p_read);
            write(r1.begin(), r1.end());
            write(endl.begin(), endl.end());
            auto r2 = std::get<2>(p_read);
            write(r2.begin(), r2.end());
            write(endl.begin(), endl.end());
            auto r3 = std::get<3>(p_read);
            write(r3.begin(), r3.end());
            write(endl.begin(), endl.end());
        }

        void write(std::string::const_iterator p_begin, std::string::const_iterator p_end)
        {
            if (!cur_buf)
            {
                if (!m_spare_buffers.pop_front(cur_buf))
                {
                    throw std::runtime_error("couldn't get a spare buffer");
                }
            }

            auto itr = cur_buf->append(p_begin, p_end);
            while (itr != p_end)
            {
                m_full_buffers.push_back(cur_buf);
                if (!m_spare_buffers.pop_front(cur_buf))
                {
                    throw std::runtime_error("couldn't get a spare buffer");
                }
                itr = cur_buf->append(itr, p_end);
            }
        }

        static void write(std::ostream& p_out, const fastq_read& p_read)
        {
            p_out << '@' << std::get<0>(p_read) << std::endl;
            p_out << std::get<1>(p_read) << std::endl;
            p_out << '+' << std::get<2>(p_read) << std::endl;
            p_out << std::get<3>(p_read) << std::endl;
        }

    private:
        void write_strings(std::ostream& p_out)
        {
            output_buffer_ptr buf_ptr;
            while (m_full_buffers.pop_front(buf_ptr))
            {
                buf_ptr->flush([&](auto p_begin, auto p_end) {
                    size_t len = p_end - p_begin;
                    p_out.write(&*p_begin, len);
                });
                m_spare_buffers.push_back(buf_ptr);
            }
        }

        void write_strings(const std::string& p_out_name, size_t p_out_buffer_size)
        {
            gzip::writer out(p_out_name, p_out_buffer_size);
            output_buffer_ptr buf_ptr;
            while (m_full_buffers.pop_front(buf_ptr))
            {
                buf_ptr->flush([&](auto p_begin, auto p_end) {
                    out.write(p_begin, p_end);
                });
                m_spare_buffers.push_back(buf_ptr);
            }
        }

        scindo::concurrent_deque<output_buffer_ptr> m_spare_buffers;
        scindo::concurrent_deque<output_buffer_ptr> m_full_buffers;
        std::thread m_thread;
        output_buffer_ptr cur_buf;
    };
}
// namespace scindo

#endif // SCINDO_FASTQ_HPP

