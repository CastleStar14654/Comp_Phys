#ifndef MISC_MATRIXLIB
#define MISC_MATRIXLIB

#include <complex>

#include "matrixlib/Base_Matrix.h"

#include "matrixlib/Band_Matrix.h"
#include "matrixlib/Base_Tri_Matrix.h"
#include "matrixlib/Diag_Matrix.h"
#include "matrixlib/Matrix.h"

#include "matrixlib/Base_Half_Band_Matrix.h"
#include "matrixlib/Hermite_Matrix.h"
#include "matrixlib/Low_Tri_Matrix.h"
#include "matrixlib/Symm_Matrix.h"
#include "matrixlib/Up_Tri_Matrix.h"

#include "matrixlib/Low_Band_Matrix.h"
#include "matrixlib/Symm_Band_Matrix.h"
#include "matrixlib/Up_Band_Matrix.h"

#include "matrixlib/Sparse_Matrix.h"

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

template <typename T, size_t N>
std::complex<T> operator*(const std::array<std::complex<T>, N> &a, const std::array<std::complex<T>, N> &b)
{
    std::complex<T> res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += std::conj(a[i]) * b[i];
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

template <typename T, size_t N, size_t C>
std::array<std::complex<T>, C> operator*(const std::array<std::complex<T>, N> &a, const Base_Matrix<std::complex<T>, N, C> &b)
{
    std::array<std::complex<T>, C> res{};
    for (std::size_t j = 0; j < N; j++)
    {
        std::complex<T> temp {std::conj(a[j])};
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

#endif // MISC_MATRIXLIB
