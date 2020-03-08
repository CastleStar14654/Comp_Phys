#ifndef MISC_LOW_TRI_MATRIX
#define MISC_LOW_TRI_MATRIX

#include <iterator>

#include "Base_Matrix.h"
#include "Diag_Matrix.h"
#include "Base_Tri_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Lower Triangular Matrix; only half of the items are stored
template <typename T, size_t N>
class Low_Tri_Matrix : public Base_Tri_Matrix<T, N>
{
private:
    // using Base_Tri_Matrix<T, N>::rs;
    // using Base_Tri_Matrix<T, N>::cs;
    using Base_Tri_Matrix<T, N>::elem;
    using Base_Tri_Matrix<T, N>::data_ln;

public:
    using typename Base_Tri_Matrix<T, N>::size_type;
    using Base_Tri_Matrix<T, N>::Base_Tri_Matrix;
    using Base_Tri_Matrix<T, N>::operator=;
    using Base_Tri_Matrix<T, N>::operator();

    Low_Tri_Matrix(const Low_Tri_Matrix &mat) = default;
    Low_Tri_Matrix(Low_Tri_Matrix &&mat) = default;

    Low_Tri_Matrix &operator=(const Low_Tri_Matrix &mat) = default;
    Low_Tri_Matrix &operator=(Low_Tri_Matrix &&mat) = default;

    const T &operator()(size_type row, size_type col) const override;
};

// -------------------------------------------------------------------------

template <typename T, size_t N>
Low_Tri_Matrix<T, N> operator*(const Low_Tri_Matrix<T, N> &a, const Low_Tri_Matrix<T, N> &b)
{
    Low_Tri_Matrix<T, N> res{};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t k = 0; k <= i; k++)
        {
            T temp {a(i, k)};
            for (std::size_t j = 0; j <= k; j++)
            {
                res(i, j) += temp * b(k, j);
            }
        }
    return res;
}

// ---------------------------------------------------------------------------

template <typename T, size_t N>
Low_Tri_Matrix<T, N> operator*(const Low_Tri_Matrix<T, N> &a, const Diag_Matrix<T, N> &b)
{
    Low_Tri_Matrix<T, N> res{};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t k = 0; k <= i; k++)
        {
            res(i, k) = a(i, k) * b(k);
        }
    return res;
}

template <typename T, size_t N>
Low_Tri_Matrix<T, N> operator*(const Diag_Matrix<T, N> &a, const Low_Tri_Matrix<T, N> &b)
{
    Low_Tri_Matrix<T, N> res{};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = 0; j <= i; j++)
        {
            res(i, j) += a(i) * b(i, j);
        }
    return res;
}

// ========================== Low_Tri_Matrix =================================

// -------------------- Low_Tri_Matrix: row & column ----------------------------

template <typename T, size_t N>
const T &Low_Tri_Matrix<T, N>::operator()(size_type row, size_type col) const
{
    if (row < col)
    {
        return Base_Matrix<T, N, N>::zero;
    }
    else
    {
        return Base_Tri_Matrix<T, N>::operator()(row, col);
    }
}

} // namespace Misc

#endif // MISC_LOW_TRI_MATRIX
