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
inline T operator*(const std::array<T, N> &a, const std::array<T, N> &b)
{
    T res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}
template <typename T, size_t N, size_t _C>
inline T operator*(const std::array<T, N> &a, const Column<T, N, _C> &b)
{
    T res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}
template <typename T, size_t N, size_t _R>
inline T operator*(const Row<T, N, _R> &a, const std::array<T, N> &b)
{
    T res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

template <typename T, size_t N>
inline std::complex<T> operator*(const std::array<std::complex<T>, N> &a, const std::array<std::complex<T>, N> &b)
{
    std::complex<T> res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += std::conj(a[i]) * b[i];
    }
    return res;
}

template <typename T, size_t N, size_t C>
inline std::array<T, C> operator*(const std::array<T, N> &a, const Base_Matrix<T, N, C> &b)
{
    std::array<T, C> res;
    res.fill(T{});
    for (std::size_t j = 0; j < N; j++)
    {
        T temp{a[j]};
        for (std::size_t i = 0; i < C; i++)
        {
            res[i] += temp * b(j, i);
        }
    }
    return res;
}

template <typename T, size_t N, size_t C>
inline std::array<std::complex<T>, C> operator*(const std::array<std::complex<T>, N> &a, const Base_Matrix<std::complex<T>, N, C> &b)
{
    std::array<std::complex<T>, C> res;
    res.fill(std::complex<T>{});
    for (std::size_t j = 0; j < N; j++)
    {
        std::complex<T> temp{std::conj(a[j])};
        for (std::size_t i = 0; i < C; i++)
        {
            res[i] += temp * b(j, i);
        }
    }
    return res;
}

template <typename T, size_t R, size_t N>
inline std::array<T, R> operator*(const Base_Matrix<T, R, N> &a, const std::array<T, N> &b)
{
    std::array<T, R> res;
    res.fill(T{});
    for (std::size_t i = 0; i < R; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            res[i] += a(i, j) * b[j];
        }
    return res;
}

} // namespace Misc

#endif // MISC_MATRIXLIB
