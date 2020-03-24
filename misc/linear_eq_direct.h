#ifndef MISC_LINEAR_EQ_DIRECT
#define MISC_LINEAR_EQ_DIRECT

#include <stdexcept>
#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <iterator>
#include <iostream>

#include "Matrix_Catalogue.h"

namespace Misc
{

// ================== DECLEARATIONS ==================
// ------------------ implementations ----------------

template <size_t B = 0, typename T, size_t N>
inline void _up_back_sub(const Base_Matrix<T, N, N> &in_mat,
                         const std::array<T, N> &in_b,
                         std::array<T, N> &out_x)
{
    constexpr size_t M{B ? B : N - 1};
    std::copy(in_b.begin(), in_b.end(), out_x.begin());
    size_t up_bound{N - 1};
    for (size_t i = N - 1; i != -1; i--)
    {
        if (in_mat(i, i) == 0)
        {
            throw std::runtime_error("back_sub(): singular matrix");
        }
        up_bound = N - 1 > i + M ? i + M : N - 1;
        for (size_t j = i + 1; j <= up_bound; j++)
        {
            out_x[i] -= in_mat(i, j) * out_x[j];
        }
        out_x[i] /= in_mat(i, i);
    }
}

template <size_t B = 0, typename T, size_t N>
inline void _low_back_sub(const Base_Matrix<T, N, N> &in_mat,
                          const std::array<T, N> &in_b,
                          std::array<T, N> &out_x)
{
    constexpr size_t M{B ? B : N - 1};
    std::copy(in_b.begin(), in_b.end(), out_x.begin());
    size_t low_bound{0};
    for (size_t i = 0; i < N; i++)
    {
        if (in_mat(i, i) == 0)
        {
            throw std::runtime_error("back_sub(): singular matrix");
        }
        low_bound = i > M ? i - M : 0;
        for (size_t j = low_bound; j < i; j++)
        {
            out_x[i] -= in_mat(i, j) * out_x[j];
        }
        out_x[i] /= in_mat(i, i);
    }
}

template <typename T, size_t N>
inline void _l_u_decomposition(const Base_Matrix<T, N, N> &in_mat,
                               Base_Matrix<T, N, N> &out_l,
                               Base_Matrix<T, N, N> &out_u,
                               std::array<size_t, N> *pivot = nullptr);

template <typename T, size_t N>
inline void _tridiagonal_decomposition(const Band_Matrix<T, N, 1> &in_mat,
                                       Base_Matrix<T, N, N> &out_l,
                                       Base_Matrix<T, N, N> &out_u);

// this function would sort the original container!
template <typename It>
inline std::size_t _inversion_number(It b, It e);

// IMPLEMENTATION, would change the input into an identity
// return the inverse matrix of an Upper Triangular Matrix
// using Gauss-Jordan method
template <typename T, size_t N>
inline Up_Tri_Matrix<T, N> _inv(Up_Tri_Matrix<T, N> &in_mat);

// IMPLEMENTATION, would change the input into an identity
// return the inverse matrix of a Symmetric Matrix
// using Gauss-Jordan method
// template <typename T, size_t N>
// Symm_Matrix<T, N> _inv(Symm_Matrix<T, N>& in_mat);

// ------------------ APIs ----------------------------

// back-substitution
template <typename T, size_t N>
inline void back_sub(const Diag_Matrix<T, N> &in_mat,
                     std::array<T, N> &in_b)
{
    for (size_t i = 0; i < N; i++)
    {
        in_b[i] /= in_mat(i);
    }
}

template <typename T, size_t N>
inline void back_sub(const Diag_Matrix<T, N> &in_mat,
                     const std::array<T, N> &in_b,
                     std::array<T, N> &out_x)
{
    for (size_t i = 0; i < N; i++)
    {
        out_x[i] = in_b[i] / in_mat(i);
    }
}

template <typename T, size_t N>
inline void back_sub(const std::array<T, N> &in_mat,
                     std::array<T, N> &in_b)
{
    for (size_t i = 0; i < N; i++)
    {
        in_b[i] /= in_mat[i];
    }
}

template <typename T, size_t N>
inline void back_sub(const std::array<T, N> &in_mat,
                     const std::array<T, N> &in_b,
                     std::array<T, N> &out_x)
{
    for (size_t i = 0; i < N; i++)
    {
        out_x[i] = in_b[i] / in_mat[i];
    }
}

template <typename T, size_t N>
inline void back_sub(const Up_Tri_Matrix<T, N> &in_mat,
                     std::array<T, N> &in_b,
                     bool reverse = false)
{
    if (reverse)
    {
        _low_back_sub(in_mat.trans(), in_b, in_b);
    }
    else
    {
        _up_back_sub(in_mat, in_b, in_b);
    }
}

template <typename T, size_t N>
inline void back_sub(const Up_Tri_Matrix<T, N> &in_mat,
                     const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool reverse = false)
{
    if (reverse)
    {
        _low_back_sub(in_mat.trans(), in_b, out_x);
    }
    else
    {
        _up_back_sub(in_mat, in_b, out_x);
    }
}

template <typename T, size_t N, size_t M>
inline void back_sub(const Up_Band_Matrix<T, N, M> &in_mat,
                     std::array<T, N> &in_b,
                     bool reverse = false)
{
    if (reverse)
    {
        _low_back_sub<M>(in_mat.trans(), in_b, in_b);
    }
    else
    {
        _up_back_sub<M>(in_mat, in_b, in_b);
    }
}

template <typename T, size_t N, size_t M>
inline void back_sub(const Up_Band_Matrix<T, N, M> &in_mat,
                     const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool reverse = false)
{
    if (reverse)
    {
        _low_back_sub<M>(in_mat.trans(), in_b, out_x);
    }
    else
    {
        _up_back_sub<M>(in_mat, in_b, out_x);
    }
}

template <typename T, size_t N>
inline void back_sub(const Low_Tri_Matrix<T, N> &in_mat,
                     std::array<T, N> &in_b,
                     bool reverse = false)
{
    if (reverse)
    {
        _up_back_sub(in_mat.trans(), in_b, in_b);
    }
    else
    {
        _low_back_sub(in_mat, in_b, in_b);
    }
}

template <typename T, size_t N>
inline void back_sub(const Low_Tri_Matrix<T, N> &in_mat,
                     const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool reverse = false)
{
    if (reverse)
    {
        _up_back_sub(in_mat.trans(), in_b, out_x);
    }
    else
    {
        _low_back_sub(in_mat, in_b, out_x);
    }
}

template <typename T, size_t N, size_t M>
inline void back_sub(const Low_Band_Matrix<T, N, M> &in_mat,
                     std::array<T, N> &in_b,
                     bool reverse = false)
{
    if (reverse)
    {
        _up_back_sub<M>(in_mat.trans(), in_b, in_b);
    }
    else
    {
        _low_back_sub<M>(in_mat, in_b, in_b);
    }
}

template <typename T, size_t N, size_t M>
inline void back_sub(const Low_Band_Matrix<T, N, M> &in_mat,
                     const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool reverse = false)
{
    if (reverse)
    {
        _up_back_sub<M>(in_mat.trans(), in_b, out_x);
    }
    else
    {
        _low_back_sub<M>(in_mat, in_b, out_x);
    }
}

// LU decomposition
template <typename T, size_t N>
inline void lu_factor(const Base_Matrix<T, N, N> &in_mat,
                      Low_Tri_Matrix<T, N> &out_l,
                      Up_Tri_Matrix<T, N> &out_u)
{
    _l_u_decomposition(in_mat, out_l, out_u);
}

template <typename T, size_t N>
inline void lu_factor(Base_Matrix<T, N, N> &in_mat)
{
    _l_u_decomposition(in_mat, in_mat, in_mat);
}

// Pivoted LU decomposition
template <typename T, size_t N>
inline void lu_factor(const Base_Matrix<T, N, N> &in_mat,
                      Low_Tri_Matrix<T, N> &out_l,
                      Up_Tri_Matrix<T, N> &out_u,
                      std::array<size_t, N> &pivot)
{
    _l_u_decomposition(in_mat, out_l, out_u, &pivot);
}

template <typename T, size_t N>
inline void lu_factor(Base_Matrix<T, N, N> &in_mat,
                      std::array<size_t, N> &pivot)
{
    _l_u_decomposition(in_mat, in_mat, in_mat, &pivot);
}

// determinant
template <typename T, size_t N>
inline T det(const Base_Matrix<T, N, N> &in_mat);

template <typename T, size_t N>
inline T det(const Symm_Matrix<T, N> &in_mat)
{
    Low_Tri_Matrix<T, N> out_l{};
    std::array<T, N> out_d{};
    ldl_factor(in_mat, out_l, out_d);
    T res{1};
    for (size_t i = 0; i < N; i++)
    {
        res *= out_d[i];
    }
    return res;
}

template <typename T, size_t N, size_t M>
inline T det(const Symm_Band_Matrix<T, N, M> &in_mat)
{
    Low_Band_Matrix<T, N, M> out_l{};
    std::array<T, N> out_d{};
    ldl_factor(in_mat, out_l, out_d);
    T res{1};
    for (size_t i = 0; i < N; i++)
    {
        res *= out_d[i];
    }
    return res;
}

// inversion number of a container
template <typename Container>
inline std::size_t inv_num(const Container &c)
{
    std::vector<typename Container::value_type> v{c.begin(), c.end()};
    return _inversion_number(v.begin(), v.end());
}

// Thomas tridiagonal decomposition
template <typename T, size_t N>
inline void tri_factor(const Band_Matrix<T, N, 1> &in_mat,
                       Low_Band_Matrix<T, N, 1> &out_l,
                       Up_Band_Matrix<T, N, 1> &out_u)
{
    _tridiagonal_decomposition(in_mat, out_l, out_u);
}

template <typename T, size_t N>
inline void tri_factor(const Band_Matrix<T, N, 1> &in_mat,
                       Band_Matrix<T, N, 1> &out_mat)
{
    _tridiagonal_decomposition(in_mat, out_mat, out_mat);
}

template <typename T, size_t N>
inline void tri_factor(Band_Matrix<T, N, 1> &in_mat)
{
    _tridiagonal_decomposition(in_mat, in_mat, in_mat);
}

// return the inverse matrix
// using Gauss-Jordan method
template <typename T, size_t N>
inline Matrix<T, N, N> inv(const Base_Matrix<T, N, N> &in_mat);

// return the inverse matrix of an Upper Triangular Matrix
// using Gauss-Jordan method
template <typename T, size_t N>
inline Up_Tri_Matrix<T, N> inv(const Up_Tri_Matrix<T, N> &in_mat)
{
    Up_Tri_Matrix<T, N> temp{in_mat};
    return _inv(temp);
}

// return the inverse matrix of a Lower Triangular Matrix
// using Gauss-Jordan method
template <typename T, size_t N>
inline Low_Tri_Matrix<T, N> inv(const Low_Tri_Matrix<T, N> &in_mat);

template <typename T, size_t N>
inline Diag_Matrix<T, N> inv(const Diag_Matrix<T, N> &in_mat)
{
    Diag_Matrix<T, N> res{};
    for (size_t i = 0; i < N; i++)
    {
        res(i) = 1. / in_mat(i);
    }
    return res;
}

// // return the inverse matrix of a Symmetric Matrix
// // using Gauss-Jordan method
// template <typename T, size_t N>
// Symm_Matrix<T, N> inv(const Symm_Matrix<T, N>& in_mat)
// {
//     Symm_Matrix<T, N> temp {in_mat};
//     return _inv(temp);
// }

// return the inverse matrix of an Upper Band Matrix
// using Gauss-Jordan method
template <typename T, size_t N, size_t M>
inline Up_Tri_Matrix<T, N> inv(const Up_Band_Matrix<T, N, M> &in_mat)
{
    Up_Tri_Matrix<T, N> temp{in_mat};
    return _inv(temp);
}

// return the inverse matrix of a Lower Band Matrix
// using Gauss-Jordan method
template <typename T, size_t N, size_t M>
inline Low_Tri_Matrix<T, N> inv(const Low_Band_Matrix<T, N, M> &in_mat)
{
    Low_Tri_Matrix<T, N> temp{in_mat};
    return inv(temp);
}

// // return the inverse matrix of a Symmetric Band Matrix
// // using Gauss-Jordan method
// template <typename T, size_t N, size_t M>
// Symm_Matrix<T, N> inv(const Symm_Band_Matrix<T, N, M>& in_mat)
// {
//     Symm_Matrix<T, N> temp {in_mat};
//     return _inv(temp);
// }

// LDL Decomposition
template <typename T, size_t N>
inline void ldl_factor(const Symm_Matrix<T, N> &in_mat,
                       Low_Tri_Matrix<T, N> &out_l,
                       std::array<T, N> &out_d);

template <typename T, size_t N, size_t M>
inline void ldl_factor(const Symm_Band_Matrix<T, N, M> &in_mat,
                       Low_Band_Matrix<T, N, M> &out_l,
                       std::array<T, N> &out_d);

// Cholesky Decomposition
template <typename T, size_t N>
inline void cholesky(const Symm_Matrix<T, N> &in_mat,
                     Low_Tri_Matrix<T, N> &out_l)
{
    std::array<T, N> out_d{};
    ldl_factor(in_mat, out_l, out_d);
    size_t count;
    for (count = 0; count < N && out_d[count] >= 0; count++)
    {
        out_d[count] = std::sqrt(out_d[count]);
    }
    if (count != N)
    {
        throw std::runtime_error("cholesky(): non positive-semidefinite matrix");
    }

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j <= i; j++)
        {
            out_l(i, j) *= out_d[j];
        }
}

template <typename T, size_t N, size_t M>
inline void cholesky(const Symm_Band_Matrix<T, N, M> &in_mat,
                     Low_Band_Matrix<T, N, M> &out_l)
{
    std::array<T, N> out_d{};
    ldl_factor(in_mat, out_l, out_d);
    size_t count;
    for (count = 0; count < N && out_d[count] >= 0; count++)
    {
        out_d[count] = std::sqrt(out_d[count]);
    }
    if (count != N)
    {
        throw std::runtime_error("cholesky(): non positive-semidefinite matrix");
    }

    for (size_t i = 0; i < N; i++)
        for (size_t j = (i > M) ? (i - M) : 0; j <= i; j++)
        {
            out_l(i, j) *= out_d[j];
        }
}

// ================== DEFINITIONS =====================

template <typename T, size_t N>
inline void _l_u_decomposition(const Base_Matrix<T, N, N> &in_mat,
                               Base_Matrix<T, N, N> &out_l,
                               Base_Matrix<T, N, N> &out_u,
                               std::array<size_t, N> *p_pivot)
{
    // prepare the out_l
    if (&in_mat != &out_l && &out_l != &out_u)
        for (size_t i = 0; i < N; i++)
        {
            out_l(i, i) = 1.;
        }
    // prepare the pivot
    if (p_pivot)
        for (size_t i = 0; i < N; i++)
        {
            (*p_pivot)[i] = i;
        }
    // whether the in_mat should be referred via p_pivot
    bool via{p_pivot && (&in_mat != &out_l)};

    // calculation
    for (size_t i = 0; i < N; i++)
    {
        // calc the current pivot
        out_u(i, i) = in_mat(via ? (*p_pivot)[i] : i, i);
        for (size_t k = 0; k < i; k++)
        {
            out_u(i, i) -= out_l(i, k) * out_u(k, i);
        }

        // calc l (but the division)
        for (size_t j = i + 1; j < N; j++)
        {
            out_l(j, i) = in_mat(via ? (*p_pivot)[j] : j, i);
            for (size_t k = 0; k < i; k++)
            {
                out_l(j, i) -= out_l(j, k) * out_u(k, i);
            }
        }

        // find the index of the pivot line
        if (p_pivot)
        {
            static T pivot_val;
            static size_t pivot_index;
            pivot_val = std::abs(out_u(i, i));
            pivot_index = i;
            for (size_t p = i + 1; p < N; p++)
            {
                if (std::abs(out_l(p, i)) > pivot_val)
                {
                    pivot_val = std::abs(out_l(p, i));
                    pivot_index = p;
                }
            }
            // swap two lines if necessary
            if (pivot_index != i)
            {
                std::swap((*p_pivot)[pivot_index], (*p_pivot)[i]);
                if (via)
                {
                    std::swap_ranges(&out_l(i, 0), &out_l(i, i), &out_l(pivot_index, 0));
                    std::swap(out_u(i, i), out_l(pivot_index, i));
                }
                else
                {
                    std::swap_ranges(&out_l(i, 0), &out_l(i, N), &out_l(pivot_index, 0));
                }
            }
        }

        if (out_u(i, i) == 0.)
        {
            throw std::runtime_error("_l_u_decomposition(): divided by zero");
        }

        // calc u
        for (size_t j = i + 1; j < N; j++)
        {
            // the division of l(j, i) here
            out_l(j, i) /= out_u(i, i);
            out_u(i, j) = in_mat(via ? (*p_pivot)[i] : i, j);
            for (size_t k = 0; k < i; k++)
            {
                out_u(i, j) -= out_l(i, k) * out_u(k, j);
            }
        }
    }
}

template <typename T, size_t N>
inline void _tridiagonal_decomposition(const Band_Matrix<T, N, 1> &in_mat,
                                       Band_Matrix<T, N, 1> &out_l,
                                       Band_Matrix<T, N, 1> &out_u)
{
    // prepare the output matrix
    if (&out_l != &out_u)
        for (size_t i = 0; i < N; i++)
        {
            out_l(i, i) = 1;
        }
    if (&in_mat != &out_u)
        for (size_t i = 0; i < N - 1; i++)
        {
            out_u(i, i + 1) = in_mat(i, i + 1);
        }

    if (&in_mat == &out_l && &in_mat == &out_u)
    {
        for (size_t i = 1; i < N; i++)
        {
            out_l(i, i - 1) /= out_u(i - 1, i - 1);
            out_u(i, i) -= out_l(i, i - 1) * out_u(i - 1, i);
        }
    }
    else
    {
        out_u(0, 0) = in_mat(0, 0);
        for (size_t i = 1; i < N; i++)
        {
            out_l(i, i - 1) = in_mat(i, i - 1) / out_u(i - 1, i - 1);
            out_u(i, i) = in_mat(i, i) - out_l(i, i - 1) * out_u(i - 1, i);
        }
    }
}

// determinant
template <typename T, size_t N>
inline T det(const Base_Matrix<T, N, N> &in_mat)
{
    // lu_factor in situ
    Matrix<T, N, N> temp{in_mat};
    // using the pivoted one for higher precision
    std::array<size_t, N> pivot;

    lu_factor(temp, pivot);
    T res{1};
    for (size_t i = 0; i < N; i++)
    {
        res *= temp(i, i);
    }
    // multiply the inversion number of the pivot matrix
    return res * (1 - 2 * int(inv_num(pivot) % 2));
}

// this function would sort the original container!
// using MERGE_SORT
template <typename It>
inline std::size_t _inversion_number(It b, It e)
{
    // trivial case
    if (next(b) == e)
    {
        return 0;
    }
    // half-partition
    It m{b + (e - b) / 2};
    // merge_sort two sub-array
    std::size_t res{_inversion_number(b, m)};
    res += _inversion_number(m, e);

    // copy two sub-array out
    auto value{*b};
    std::vector<decltype(value)> left(b, m);
    std::vector<decltype(value)> right(m, e);
    auto l{left.begin()};
    auto r{right.begin()};
    // merge
    for (It current = b; current != e; ++current)
    {
        if (l != left.end() && (r == right.end() || *l <= *r))
        {
            *current = *l;
            ++l;
        }
        else
        {
            *current = *r;
            ++r;
            res += left.end() - l;
        }
    }
    return res;
}

// inverse matrix
template <typename T, size_t N>
inline Matrix<T, N, N> inv(const Base_Matrix<T, N, N> &in_mat)
{
    // temp will be transformed into an identity
    Matrix<T, N, N> temp{in_mat};
    // res will be returned
    Matrix<T, N, N> res{};
    for (size_t i = 0; i < N; i++)
    {
        res(i, i) = 1.;
    }

    for (size_t i = 0; i < N; i++)
    {
        // find the pivot;
        static size_t pivot_index;
        // static T pivot_value;
        pivot_index = i;
        auto pivot_value{std::abs(temp(i, i))};
        for (size_t p = i + 1; p < N; p++)
            if (std::abs(temp(p, i)) > pivot_value)
            {
                pivot_value = std::abs(temp(p, i));
                pivot_index = p;
            }
        if (pivot_index != i)
        {
            std::swap(temp.data()[i], temp.data()[pivot_index]);
            std::swap(res.data()[i], res.data()[pivot_index]);
        }

        // singular
        static T t_ii;
        t_ii = temp(i, i);
        if (t_ii == 0.)
        {
            throw std::runtime_error("inv(): divided by zero");
        }

        // refresh two matrices
        // unnecessary calculations are avoided
        for (size_t j = i + 1; j < N; j++)
        {
            temp(i, j) /= t_ii;
        }
        for (size_t j = 0; j < N; j++)
        {
            res(i, j) /= t_ii;
        }

        for (size_t k = 0; k < i; k++)
        {
            static T t_ki;
            t_ki = temp(k, i);
            for (size_t j = 0; j < N; j++)
            {
                res(k, j) -= t_ki * res(i, j);
            }
            for (size_t j = i + 1; j < N; j++)
            {
                temp(k, j) -= t_ki * temp(i, j);
            }
        }

        for (size_t k = i + 1; k < N; k++)
        {
            static T t_ki;
            t_ki = temp(k, i);
            for (size_t j = 0; j < N; j++)
            {
                res(k, j) -= t_ki * res(i, j);
            }
            for (size_t j = i + 1; j < N; j++)
            {
                temp(k, j) -= t_ki * temp(i, j);
            }
        }
    }
    return res;
}

// IMPLEMENTATION, would change the input into an identity
// return the inverse matrix of an Upper Triangular Matrix
// using Gauss-Jordan method
template <typename T, size_t N>
inline Up_Tri_Matrix<T, N> _inv(Up_Tri_Matrix<T, N> &in_mat)
{
    // res will be returned
    Up_Tri_Matrix<T, N> res{};
    for (size_t i = 0; i < N; i++)
    {
        res(i, i) = 1.;
    }

    for (size_t i = 0; i < N; i++)
    {
        // singular
        static T t_ii;
        t_ii = in_mat(i, i);
        if (t_ii == 0.)
        {
            throw std::runtime_error("inv(): divided by zero");
        }

        // refresh two matrices
        // unnecessary calculations are avoided
        for (size_t j = i + 1; j < N; j++)
        {
            in_mat(i, j) /= t_ii;
        }
        for (size_t j = i; j < N; j++)
        {
            res(i, j) /= t_ii;
        }

        for (size_t k = 0; k < i; k++)
        {
            static T t_ki;
            t_ki = in_mat(k, i);
            for (size_t j = i; j < N; j++)
            {
                res(k, j) -= t_ki * res(i, j);
            }
            for (size_t j = i + 1; j < N; j++)
            {
                in_mat(k, j) -= t_ki * in_mat(i, j);
            }
        }
    }
    return res;
}

// IMPLEMENTATION, would change the input into an identity
// return the inverse matrix of a Lower Triangular Matrix
// using Gauss-Jordan method
template <typename T, size_t N>
inline Low_Tri_Matrix<T, N> inv(const Low_Tri_Matrix<T, N> &in_mat)
{
    // res will be returned
    Low_Tri_Matrix<T, N> res{};
    for (size_t i = 0; i < N; i++)
    {
        res(i, i) = 1.;
    }

    for (size_t i = 0; i < N; i++)
    {
        // singular
        static T t_ii;
        t_ii = in_mat(i, i);
        if (t_ii == 0.)
        {
            throw std::runtime_error("inv(): divided by zero");
        }

        // refresh two matrices
        // unnecessary calculations are avoided
        for (size_t j = 0; j <= i; j++)
        {
            res(i, j) /= t_ii;
        }

        for (size_t k = i + 1; k < N; k++)
        {
            static T t_ki;
            t_ki = in_mat(k, i);
            for (size_t j = 0; j <= i; j++)
            {
                res(k, j) -= t_ki * res(i, j);
            }
        }
    }
    return res;
}

// IMPLEMENTATION, would change the input into an identity
// return the inverse matrix of a Symmetric Matrix
// using Gauss-Jordan method
// template <typename T, size_t N>
// Symm_Matrix<T, N> _inv(Symm_Matrix<T, N>& in_mat)
// {

// }

template <typename T, size_t N>
inline void ldl_factor(const Symm_Matrix<T, N> &in_mat,
                       Low_Tri_Matrix<T, N> &out_l,
                       std::array<T, N> &out_d)
{
    Up_Tri_Matrix<T, N> temp_t{};

    for (size_t i = 0; i < N; i++)
    {
        out_l(i, i) = 1.;
        for (size_t j = 0; j < i; j++)
        {
            static T t_jj;

            t_jj = temp_t(j, j);
            if (t_jj == 0.)
            {
                temp_t(j, i) = 0.;
            }
            else
            {
                temp_t(j, i) = in_mat(i, j);
                for (size_t k = 0; k < j; k++)
                {
                    temp_t(j, i) -= out_l(i, k) * temp_t(k, j);
                }
            }

            out_l(i, j) = (t_jj == 0.) ? 0. : (temp_t(j, i) / temp_t(j, j));
        }

        temp_t(i, i) = in_mat(i, i);
        for (size_t k = 0; k < i; k++)
        {
            temp_t(i, i) -= out_l(i, k) * temp_t(k, i);
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        out_d[i] = temp_t(i, i);
    }
}

/*beg:ldl_factor*/
template <typename T, size_t N, size_t M>
inline void ldl_factor(const Symm_Band_Matrix<T, N, M> &in_mat,
                       Low_Band_Matrix<T, N, M> &out_l,
                       std::array<T, N> &out_d)
{
    Up_Band_Matrix<T, N, M> temp_t{};

    for (size_t i = 0; i < N; i++)
    {
        out_l(i, i) = 1;
        static size_t limit;
        limit = (i > M) ? (i - M) : 0;
        for (size_t j = limit; j < i; j++)
        {
            static T t_jj;

            t_jj = temp_t(j, j);
            if (t_jj == 0.)
            {
                temp_t(j, i) = 0.;
            }
            else
            {
                temp_t(j, i) = in_mat(i, j);
                for (size_t k = limit; k < j; k++)
                {
                    temp_t(j, i) -= out_l(i, k) * temp_t(k, j);
                }
            }

            out_l(i, j) = (t_jj == 0.) ? 0. : (temp_t(j, i) / temp_t(j, j));
        }

        temp_t(i, i) = in_mat(i, i);
        for (size_t k = limit; k < i; k++)
        {
            temp_t(i, i) -= out_l(i, k) * temp_t(k, i);
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        out_d[i] = temp_t(i, i);
    }
}
/*end:ldl_factor*/

} // namespace Misc

#endif
