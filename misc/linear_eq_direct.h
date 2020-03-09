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

template <typename T, size_t N>
void _l_u_decomposition(const Base_Matrix<T, N, N> &in_mat,
                        Base_Matrix<T, N, N> &out_l,
                        Base_Matrix<T, N, N> &out_u,
                        std::array<size_t, N> *pivot = nullptr);

template <typename T, size_t N>
void _tridiagonal_decomposition(const Band_Matrix<T, N, 1> &in_mat,
                                Base_Matrix<T, N, N> &out_l,
                                Base_Matrix<T, N, N> &out_u);

// this function would sort the original container!
template <typename It>
std::size_t _inversion_number(It b, It e);

// ------------------ APIs ----------------------------

// LU decomposition
template <typename T, size_t N>
void l_u_factor(const Base_Matrix<T, N, N> &in_mat,
                Low_Tri_Matrix<T, N> &out_l,
                Up_Tri_Matrix<T, N> &out_u)
{
    _l_u_decomposition(in_mat, out_l, out_u);
}

template <typename T, size_t N>
void l_u_factor(Base_Matrix<T, N, N> &in_mat)
{
    _l_u_decomposition(in_mat, in_mat, in_mat);
}

// determinant
template <typename T, size_t N>
T det(const Base_Matrix<T, N, N> &in_mat);

// inversion number of a container
template <typename Container>
std::size_t inv_num(const Container& c)
{
    std::vector<typename Container::value_type> v {c.begin(), c.end()};
    return _inversion_number(v.begin(), v.end());
}

// Pivoted LU decomposition
template <typename T, size_t N>
void l_u_factor(const Base_Matrix<T, N, N> &in_mat,
                Low_Tri_Matrix<T, N> &out_l,
                Up_Tri_Matrix<T, N> &out_u,
                std::array<size_t, N> &pivot)
{
    _l_u_decomposition(in_mat, out_l, out_u, &pivot);
}

template <typename T, size_t N>
void l_u_factor(Base_Matrix<T, N, N> &in_mat,
                std::array<size_t, N> &pivot)
{
    _l_u_decomposition(in_mat, in_mat, in_mat, &pivot);
}

// Thomas tridiagonal decomposition
template <typename T, size_t N>
void tri_factor(const Band_Matrix<T, N, 1> &in_mat,
                Low_Band_Matrix<T, N, 1> &out_l,
                Up_Band_Matrix<T, N, 1> &out_u)
{
    _tridiagonal_decomposition(in_mat, out_l, out_u);
}

template <typename T, size_t N>
void tri_factor(const Band_Matrix<T, N, 1> &in_mat,
                Band_Matrix<T, N, 1> &out_mat)
{
    _tridiagonal_decomposition(in_mat, out_mat, out_mat);
}

template <typename T, size_t N>
void tri_factor(Band_Matrix<T, N, 1> &in_mat)
{
    _tridiagonal_decomposition(in_mat, in_mat, in_mat);
}

// ================== DEFINITIONS =====================

template <typename T, size_t N>
void _l_u_decomposition(const Base_Matrix<T, N, N> &in_mat,
                        Base_Matrix<T, N, N> &out_l,
                        Base_Matrix<T, N, N> &out_u,
                        std::array<size_t, N> *p_pivot)
{
    // prepare the out_l
    if (&in_mat != &out_l && &out_l != &out_u)
        for (size_t i = 0; i < N; i++)
        {
            out_l(i, i) = 1;
        }
    // prepare the pivot
    if (p_pivot)
        for (size_t i = 0; i < N; i++)
        {
            (*p_pivot)[i] = i;
        }
    // whether the in_mat should be referred via p_pivot
    bool via {p_pivot && (&in_mat != &out_l)};

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

        if (out_u(i, i) == 0)
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
void _tridiagonal_decomposition(const Band_Matrix<T, N, 1> &in_mat,
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
T det(const Base_Matrix<T, N, N> &in_mat)
{
    // l_u_factor in situ
    Matrix<T, N, N> temp {in_mat};
    // using the pivoted one for higher precision
    std::array<size_t, N> pivot;

    l_u_factor(temp, pivot);
    T res {1};
    for (size_t i = 0; i < N; i++)
    {
        res *= temp(i, i);
    }
    // multiply the inversion number of the pivot matrix
    return res*(1 - 2*int(inv_num(pivot)%2));
}

// this function would sort the original container!
// using MERGE_SORT
template <typename It>
std::size_t _inversion_number(It b, It e)
{
    // trivial case
    if (next(b) == e)
    {
        return 0;
    }
    // half-partition
    It m {b + (e-b)/2};
    // merge_sort two sub-array
    std::size_t res{_inversion_number(b, m)};
    res += _inversion_number(m, e);

    // copy two sub-array out
    auto value {*b};
    std::vector<decltype(value)> left(b, m);
    std::vector<decltype(value)> right(m, e);
    auto l {left.begin()};
    auto r {right.begin()};
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

} // namespace Misc

#endif
