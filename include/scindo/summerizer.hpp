#ifndef SCINDO_SUMMERIZER_HPP
#define SCINDO_SUMMERIZER_HPP

namespace scindo
{
    struct summerizer
    {
        uint64_t n;
        double s;
        double s2;

        summerizer()
            : n(0), s(0), s2(0)
        {
        }

        void push_back(double p_x)
        {
            n += 1;
            s += p_x;
            s2 += p_x*p_x;
        }

        uint64_t size() const
        {
            return n;
        }

        double sum() const
        {
            return s;
        }

        double mean() const
        {
            return s/n;
        }

        double var() const
        {
            double m = mean();
            return s2/n - m*m;
        }

        double sd() const
        {
            return std::sqrt(var());
        }

        void clear()
        {
            n = 0;
            s = 0;
            s2 = 0;
        }
    };
}
// namespace scindo

#endif // SCINDO_SUMMERIZER_HPP
