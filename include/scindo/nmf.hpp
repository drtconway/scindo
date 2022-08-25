#ifndef SCINDO_NMF_HPP
#define SCINDO_NMF_HPP

namespace scindo
{
    template <typename T>
    struct is_matrix
    {
        static constexpr bool value = false;
    };

    template <typename T, typename U>
    struct most_general_type
    {
        using type = decltype(T() + U());
    };

    template <typename ValueType>
    struct basic_matrix
    {
        using value_type = ValueType;

        size_t theHeight;
        size_t theWidth;
        std::vector<std::vector<value_type>> data;

        basic_matrix(const size_t& p_height, const size_t& p_width)
            : theHeight(p_height), theWidth(p_width)
        {
            data.resize(theHeight);
            for (size_t i = 0; i < theHeight; ++i)
            {
                data[i].resize(theWidth);
            }
        }

        size_t width() const
        {
            return theWidth;
        }

        size_t height() const
        {
            return theHeight;
        }


        value_type operator()(const size_t p_r, const size_t& p_c) const
        {
            return data[p_r][p_c];
        }

        value_type& operator()(const size_t p_r, const size_t& p_c)
        {
            return data[p_r][p_c];
        }
    };

    template <typename ValueType>
    struct is_matrix<basic_matrix<ValueType>>
    {
        static constexpr bool value = true;
    };

    template <typename MatrixType>
    struct transpose
    {
        using value_type = typename MatrixType::value_type;

        const MatrixType& M;

        transpose(const MatrixType& p_M)
            : M(p_M)
        {}

        size_t width() const
        {
            return M.height();
        }

        size_t height() const
        {
            return M.width();
        }

        typename MatrixType::value_type operator()(size_t p_r, size_t p_c) const
        {
            return M(p_c, p_r);
        }
    };

    template <typename MatrixType>
    struct is_matrix<transpose<MatrixType>>
    {
        static constexpr bool value = true;
    };

    template <typename MatrixType>
    transpose<MatrixType> t(const MatrixType& M)
    {
        return transpose<MatrixType>(M);
    }

    template <typename AType, typename BType, typename ValueType = most_general_type<typename AType::value_type, typename BType::value_type>::type>
    struct mult
    {
        using value_type = ValueType;

        const AType& A;
        const BType& B;

        mult(const AType& p_A, const BType& p_B)
            : A(p_A), B(p_B)
        {}

        size_t height() const
        {
            return A.height();
        }

        size_t width() const
        {
            return B.width();
        }

        value_type operator()(size_t p_r, size_t p_c) const
        {
            value_type v = 0;
            for (size_t p = 0; p < A.width(); ++p)
            {
                v += A(p_r, p) * B(p, p_c);
            }
            return v;
        }
    };

    template <typename AType, typename BType, typename ValueType>
    struct is_matrix<mult<ValueType,AType,BType>>
    {
        static constexpr bool value = true;
    };

    template <typename AType, typename BType>
    typename std::enable_if<is_matrix<AType>::value && is_matrix<BType>::value, mult<AType,BType>>::type operator*(const AType& a, const BType& b)
    {
        return mult<AType,BType>(a, b);
    }

    template <typename VType, typename WType, typename HType>
    struct nmf
    {
        static_assert(std::is_integral<typename VType::value_type>::value || std::is_floating_point<typename VType::value_type>::value);
        static_assert(std::is_integral<typename WType::value_type>::value || std::is_floating_point<typename WType::value_type>::value);
        static_assert(std::is_integral<typename HType::value_type>::value || std::is_floating_point<typename HType::value_type>::value);
        static_assert(std::is_same<typename most_general_type<bool, float>::type, float>::value);
        static_assert(!std::is_same<typename most_general_type<typename VType::value_type, typename WType::value_type>::type, void>::value);
        static_assert(!std::is_same<typename most_general_type<typename VType::value_type, typename HType::value_type>::type, void>::value);
        const VType& V;
        WType W;
        HType H;

        nmf(const VType& p_V, size_t p_p)
            : V(p_V), W(V.height(), p_p), H(p_p, V.width())
        {}

        template <typename ValueSource>
        void init(ValueSource p_values)
        {
            // Fill W with values.
            //
            for (size_t r = 0; r < W.height(); ++r)
            {
                for (size_t c = 0; c < W.width(); ++c)
                {
                    W(r,c) = p_values();
                }
            }

            // Fill W with values.
            //
            for (size_t r = 0; r < H.height(); ++r)
            {
                for (size_t c = 0; c < H.width(); ++c)
                {
                    W(r,c) = p_values();
                }
            }
        }

        double iterate()
        {
            double delta = 0;
            {
                HType H1(H.height(), H.width());
                for (size_t r = 0; r < H.height(); ++r)
                {
                    for (size_t c = 0; c < H.width(); ++c)
                    {
                        H1(r,c) = H(r,c) * ((t(W)*V)(r,c)/(t(W)*W*H)(r,c));
                        double drc = H1(r,c) - H(r,c);
                        delta += drc*drc;
                    }
                }
                std::swap(H.data, H1.data);
            }
            {
                WType W1(W.height(), W.width());
                for (size_t r = 0; r < W.height(); ++r)
                {
                    for (size_t c = 0; c < W.width(); ++c)
                    {
                        W1(r,c) = W(r,c) * ((V*t(H))(r,c) / (W*H*t(H))(r,c));
                        double drc = W1(r,c) - W(r,c);
                        delta += drc*drc;
                    }
                }
                std::swap(W.data, W1.data);
            }
            return delta;
        }
    };
}
// namespace scindo

#endif // SCINDO_NMF_HPP
