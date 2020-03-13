#ifndef MISC_LINEAR_EQ_ITERATIVE
#define MISC_LINEAR_EQ_ITERATIVE

#include <stdexcept>
#include <array>
#include <map>
#include <vector>
#include <cmath>
#include <utility>
#include <iterator>
#include <iostream>

#include "Matrix_Catalogue.h"

namespace Misc
{

// ================== DECLEARATIONS ==================

template <typename T, size_t N>
T norm_1(const std::array<T, N> &x)
{
    T res{};
    for (auto &i : x)
    {
        res += std::abs(i);
    }
    return res;
}

template <typename T, size_t N>
T norm_1(const std::array<T, N> &x1, const std::array<T, N> &x2)
{
    T res{};
    for (size_t i = 0; i < N; i++)
    {
        res += std::abs(x1[i] - x2[i]);
    }
    return res;
}

// Jacobi iterative method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
//      2: rho(A) > 1 (i.e. not converge)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
template <typename T, size_t N>
int jacobi(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
           std::array<T, N> &out_x,
           bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-10);

// Gauss-Seidel iterative method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
//      2: rho(A) > 1 (i.e. not converge)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
template <typename T, size_t N>
int gauss_seidel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                 std::array<T, N> &out_x,
                 bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-10);

// successive over relaxation iterative method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
//      2: rho(A) > 1 (i.e. not converge)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
template <typename T, size_t N>
int suc_over_rel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                 std::array<T, N> &out_x, T omega,
                 bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-10);

// ================== DEFINITIONS ==================

template <typename T, size_t N>
int jacobi(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
           std::array<T, N> &out_x,
           bool sparse, size_t max_times, double rel_epsilon)
{
    std::array<T, N> prev_x{std::move(out_x)};
    T delta_norm;
    T new_norm;
    T prev_norm{norm_1(prev_x)};
    static double not_conv_cri{2};
    if (sparse)
    {
        std::array<std::map<size_t, T>, N> sparse_mat;
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++)
                if (in_mat(i, j) != T{})
                {
                    sparse_mat[i][j] = in_mat(i, j);
                }
        for (size_t count = 0; count < max_times; count++)
        {
            out_x = in_b;
            for (size_t i = 0; i < N; i++)
            {
                T &obj{out_x[i]};
                for (auto &p : sparse_mat[i])
                {
                    if (p.first != i)
                    {
                        obj -= p.second * prev_x[p.first];
                    }
                }
                obj /= in_mat(i, i);
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
    }
    else
    {
        for (size_t count = 0; count < max_times; count++)
        {
            out_x = in_b;
            for (size_t i = 0; i < N; i++)
            {
                T &obj{out_x[i]};
                for (size_t j = 0; j < i; j++)
                {
                    obj -= in_mat(i, j) * prev_x[j];
                }
                for (size_t j = i + 1; j < N; j++)
                {
                    obj -= in_mat(i, j) * prev_x[j];
                }
                obj /= in_mat(i, i);
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);
for (auto i : out_x)
{
    std::cout << i << '\t';
}
std::cout << std::endl;

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
    }
    return 1;
}

template <typename T, size_t N>
int gauss_seidel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                 std::array<T, N> &out_x,
                 bool sparse, size_t max_times, double rel_epsilon)
{
    std::array<T, N> prev_x{std::move(out_x)};
    T delta_norm;
    T new_norm;
    T prev_norm{norm_1(prev_x)};
    static double not_conv_cri{2};
    if (sparse)
    {
        std::array<std::map<size_t, T>, N> sparse_mat;
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++)
                if (in_mat(i, j) != T{})
                {
                    sparse_mat[i][j] = in_mat(i, j);
                }
        for (size_t count = 0; count < max_times; count++)
        {
            out_x = in_b;
            for (size_t i = 0; i < N; i++)
            {
                T &obj{out_x[i]};
                for (auto &p : sparse_mat[i])
                {
                    if (p.first < i)
                    {
                        obj -= p.second * out_x[p.first];
                    }
                    else if (p.first > i)
                    {
                        obj -= p.second * prev_x[p.first];
                    }
                }
                obj /= in_mat(i, i);
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
    }
    else
    {
        for (size_t count = 0; count < max_times; count++)
        {
            out_x = in_b;
            for (size_t i = 0; i < N; i++)
            {
                T &obj{out_x[i]};
                for (size_t j = 0; j < i; j++)
                {
                    obj -= in_mat(i, j) * out_x[j];
                }
                for (size_t j = i + 1; j < N; j++)
                {
                    obj -= in_mat(i, j) * prev_x[j];
                }
                obj /= in_mat(i, i);
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);
for (auto i : out_x)
{
    std::cout << i << '\t';
}
std::cout << std::endl;

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
    }
    return 1;
}

template <typename T, size_t N>
int suc_over_rel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                 std::array<T, N> &out_x, T omega,
                 bool sparse, size_t max_times, double rel_epsilon)
{
    if (omega <= 0 || omega >= 2)
    {
        throw std::runtime_error("suc_over_rel(): invalid omega");
    }

    std::array<T, N> prev_x{std::move(out_x)};
    T delta_norm;
    T new_norm;
    T prev_norm{norm_1(prev_x)};
    static double not_conv_cri{2};
    if (sparse)
    {
        std::array<std::map<size_t, T>, N> sparse_mat;
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++)
                if (in_mat(i, j) != T{})
                {
                    sparse_mat[i][j] = in_mat(i, j);
                }
        for (size_t count = 0; count < max_times; count++)
        {
            out_x = in_b;
            for (size_t i = 0; i < N; i++)
            {
                T &obj{out_x[i]};
                for (auto &p : sparse_mat[i])
                {
                    if (p.first < i)
                    {
                        obj -= p.second * out_x[p.first];
                    }
                    else if (p.first > i)
                    {
                        obj -= p.second * prev_x[p.first];
                    }
                }
                obj *= omega/in_mat(i, i);
                obj += (1-omega)*prev_x[i];
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
    }
    else
    {
        for (size_t count = 0; count < max_times; count++)
        {
            out_x = in_b;
            for (size_t i = 0; i < N; i++)
            {
                T &obj{out_x[i]};
                for (size_t j = 0; j < i; j++)
                {
                    obj -= in_mat(i, j) * out_x[j];
                }
                for (size_t j = i + 1; j < N; j++)
                {
                    obj -= in_mat(i, j) * prev_x[j];
                }
                obj *= omega/in_mat(i, i);
                obj += (1-omega)*prev_x[i];
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);
for (auto i : out_x)
{
    std::cout << i << '\t';
}
std::cout << std::endl;

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
    }
    return 1;
}

} // namespace Misc

#endif // MISC_LINEAR_EQ_ITERATIVE
