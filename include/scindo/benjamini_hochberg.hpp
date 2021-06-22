#ifndef SCINDO_BENJAMINI_HOCHBERG_HPP
#define SCINDO_BENJAMINI_HOCHBERG_HPP

namespace scindo
{
    struct benjamini_hochberg
    {
        const double Q;

        benjamini_hochberg(const double& p_Q)
            : Q(p_Q)
        {
        }

        void operator()(const std::vector<double>& p_pValues, std::vector<double>& p_qValue) const
        {
            std::vector<size_t> perm;
            perm.reserve(p_pValues.size());
            for (size_t i = 0; i < p_pValues.size(); ++i)
            {
                perm.push_back(i);
            }

            std::sort(perm.begin(), perm.end(), [&](const auto& p_lhs, const auto& p_rhs) {
                return p_pValues[p_lhs] < p_pValues[p_rhs];
            });

            const size_t N = perm.size();
            p_qValue.resize(N, -1);
            for (size_t i = 0; i < N; ++i)
            {
                double q = double(i)/double(N) * Q;
                p_qValue[perm[i]] = q;
            }
        }
    };
}
// namespace scindo

#endif // SCINDO_BENJAMINI_HOCHBERG_HPP
