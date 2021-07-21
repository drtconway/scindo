#ifndef SCINDO_FASTQ_HPP
#define SCINDO_FASTQ_HPP

#include <istream>
#include <string>
#include <tuple>
#include <boost/algorithm/string/trim.hpp>

namespace scindo
{
    typedef std::tuple<std::string,std::string,std::string,std::string> fastq_read;

    class fastq_reader
    {
    public:
        typedef fastq_read item_type;
        typedef item_type const& item_result_type;

        fastq_reader(std::istream& p_in)
            : m_in(p_in), m_more(true)
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

    private:
        void next()
        {
            if (!std::getline(m_in, m_next))
            {
                m_more = false;
                return;
            }
            if (!starts_with(m_next, '@'))
            {
                std::cerr << m_next << std::endl;
                throw std::domain_error("missing read identifier line.");
            }

            boost::algorithm::trim(m_next);
            std::string& id1 = std::get<0>(m_curr);
            id1.clear();
            id1.insert(id1.end(), m_next.begin() + 1, m_next.end()); // drop the @

            if (!std::getline(m_in, m_next))
            {
                throw std::domain_error("missing read sequence line.");
            }

            boost::algorithm::trim(m_next);
            std::string& seq = std::get<1>(m_curr);
            seq.clear();
            seq.insert(seq.end(), m_next.begin(), m_next.end());

            if (!std::getline(m_in, m_next) || !starts_with(m_next, '+'))
            {
                throw std::domain_error("missing read spacer line.");
            }

            boost::algorithm::trim(m_next);
            std::string& id2 = std::get<2>(m_curr);
            id2.clear();
            id2.insert(id2.end(), m_next.begin() + 1, m_next.end()); // drop the +

            if (!std::getline(m_in, m_next))
            {
                throw std::domain_error("missing read quality score line.");
            }

            boost::algorithm::trim(m_next);
            std::string& qual = std::get<3>(m_curr);
            qual.clear();
            qual.insert(qual.end(), m_next.begin(), m_next.end());
        }

        static bool starts_with(const std::string& p_str, char p_ch)
        {
            return p_str.size() > 0 && p_str.front() == p_ch;
        }

        std::istream& m_in;
        bool m_more;
        std::string m_next;
        fastq_read m_curr;
    };

    class fastq_writer
    {
    public:
        static void write(std::ostream& p_out, const fastq_read& p_read)
        {
            p_out << '@' << std::get<0>(p_read) << std::endl;
            p_out << std::get<1>(p_read) << std::endl;
            p_out << '+' << std::get<2>(p_read) << std::endl;
            p_out << std::get<3>(p_read) << std::endl;
        }
    };
}
// namespace scindo

#endif // SCINDO_FASTQ_HPP
