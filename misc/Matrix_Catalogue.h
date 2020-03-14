#ifndef MISC_MATRIX_CATALOGUE
#define MISC_MATRIX_CATALOGUE

#include "Base_Matrix.h"

#include "Band_Matrix.h"
#include "Base_Tri_Matrix.h"
#include "Diag_Matrix.h"
#include "Matrix.h"

#include "Base_Half_Band_Matrix.h"
#include "Low_Tri_Matrix.h"
#include "Symm_Matrix.h"
#include "Up_Tri_Matrix.h"

#include "Low_Band_Matrix.h"
#include "Symm_Band_Matrix.h"
#include "Up_Band_Matrix.h"

#include "Sparse_Matrix.h"

namespace Misc
{

template <typename T, size_t N>
T operator*(const std::array<T, N> &a, const std::array<T, N> &b)
{
    T res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

template <typename T, size_t N, size_t C>
std::array<T, C> operator*(const std::array<T, N> &a, const Base_Matrix<T, N, C> &b)
{
    std::array<T, C> res{};
    for (std::size_t j = 0; j < N; j++)
    {
        T temp {a[j]};
        for (std::size_t i = 0; i < C; i++)
        {
            res[i] += temp*b(j, i);
        }
    }
    return res;
}

template <typename T, size_t R, size_t N>
std::array<T, R> operator*(const Base_Matrix<T, R, N> &a, const std::array<T, N> &b)
{
    std::array<T, R> res{};
    for (std::size_t i = 0; i < R; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            res[i] += a(i, j)*b[j];
        }
    return res;
}

template <typename T, size_t N, size_t C>
std::array<T, C> operator*(const std::array<T, N> &a, const Sparse_Matrix<T, N, C> &b)
{
    std::array<T, C> res{};
    for (std::size_t j = 0; j < N; j++)
    {
        T temp {a[j]};
        for (const auto& p: b[j])
        {
            res[p.first] += temp*p.second;
        }
    }
    return res;
}

template <typename T, size_t R, size_t N>
std::array<T, R> operator*(const Sparse_Matrix<T, R, N> &a, const std::array<T, N> &b)
{
    std::array<T, R> res{};
    for (std::size_t i = 0; i < R; i++)
        for (const auto& p: a[i])
        {
            res[i] += p.second*b[p.first];
        }
    return res;
}

} // namespace Misc

#endif // MISC_MATRIX_CATALOGUE
