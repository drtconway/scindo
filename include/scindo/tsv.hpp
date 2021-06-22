#ifndef SCINDO_TSV_HPP
#define SCINDO_TSV_HPP

#include <istream>

namespace scindo
{
    namespace detail
    {
        struct tsv_reader
        {
            std::istream& in;
            std::string line;

            tsv_reader(std::istream& p_in)
                : in(p_in)
            {
            }

            bool next(std::vector<std::string>& p_row)
            {
                if (!std::getline(in, line))
                {
                    return false;
                }
                p_row.clear();
                size_t i = 0;
                size_t n = line.find('\t');
                while (n != std::string::npos)
                {
                    p_row.push_back(std::string(line.begin() + i, line.begin() + n));
                    i = n + 1;
                    n = line.find('\t', i);
                }
                p_row.push_back(std::string(line.begin() + i, line.end()));
                return true;
            }
        };
    }
    // namespace detail

    struct tsv : std::vector<std::vector<std::string>>
    {
        std::vector<std::string> hdr;
        std::unordered_map<std::string,size_t> idx;

        tsv(std::istream& p_in)
        {
            detail::tsv_reader r(p_in);
            if (!r.next(hdr))
            {
                throw std::runtime_error("no header line");
            }
            for (size_t i = 0; i < hdr.size(); ++i)
            {
                idx[hdr[i]] = i;
            }

            std::vector<std::string> x;
            while (r.next(x))
            {
                this->push_back(std::vector<std::string>());
                std::swap(this->back(), x);
            }
        }
    };
}
// namespace scindo

#endif // SCINDO_TSV_HPP
