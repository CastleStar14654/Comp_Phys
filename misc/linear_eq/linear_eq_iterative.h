#ifndef MISC_LINEAR_EQ_ITERATIVE
#define MISC_LINEAR_EQ_ITERATIVE

#include <stdexcept>
#include <array>
#include <vector>
#include <complex>
#include <cmath>
#include <utility>
#include <iterator>
#include <iostream>
#include <memory>

#include "../matrixlib.h"

namespace Misc
{

// ================== DECLEARATIONS ==================

template <typename T, size_t N>
inline T norm_1(const std::array<T, N> &x)
{
    T res{};
    for (auto &i : x)
    {
        res += std::abs(i);
    }
    return res;
}

template <typename T, size_t N>
inline T norm_1(const std::array<T, N> &x1, const std::array<T, N> &x2)
{
    T res{};
    for (size_t i = 0; i < N; i++)
    {
        res += std::abs(x1[i] - x2[i]);
    }
    return res;
}

template <typename T, size_t N>
inline T norm_1(const std::array<std::complex<T>, N> &x)
{
    T res{};
    for (auto &i : x)
    {
        res += std::abs(i);
    }
    return res;
}

template <typename T, size_t N>
inline T norm_1(const std::array<std::complex<T>, N> &x1, const std::array<std::complex<T>, N> &x2)
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
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
template <typename T, size_t N>
inline int jacobi(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                  std::array<T, N> &out_x,
                  bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-15);

// Gauss-Seidel iterative method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
//      2: rho(A) > 1 (i.e. not converge)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
template <typename T, size_t N>
inline int gauss_seidel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                        std::array<T, N> &out_x,
                        bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-15);

// successive over relaxation iterative method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
//      2: rho(A) > 1 (i.e. not converge)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// omega:       parameter;
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
template <typename T, size_t N>
inline int suc_over_rel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                        std::array<T, N> &out_x, T omega,
                        bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-15);

// gradient descent iterative method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
template <typename T, size_t N>
inline int grad_des(const Symm_Matrix<T, N> &in_mat, const std::array<T, N> &in_b,
                    std::array<T, N> &out_x,
                    bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-15);

/*beg:conj_grad_dec*/
// conjugate gradient method method to solve eqs Ax = b
// ---------- return ------------
// int
//      0: end normally
//      1: meet max_times (i.e. converge too slowly)
// ---------- parameters -------------
// in_mat:      A;          in_b:   b
// out_x:       x for output, the original value is used as the initial value
// sparse:      whether A is sparse; some optimization of speed will be applied; will consume more RAM
// max_times:   maximum iteration times;
// rel_epsilon: when the relative error is less than this, iteration will be stopped
template <typename T, size_t N>
inline int conj_grad(const Symm_Matrix<T, N> &in_mat, const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-15);
/*end:conj_grad_dec*/

// ================== DEFINITIONS ==================

template <typename T, size_t N>
inline int jacobi(const Sparse_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                  std::array<T, N> &out_x,
                  bool sparse = true, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    std::array<T, N> prev_x{std::move(out_x)};
    T delta_norm;
    T new_norm;
    T prev_norm{norm_1(prev_x)};
    static double not_conv_cri{2};

    for (size_t count = 0; count < max_times; count++)
    {
        out_x = in_b;
        for (size_t i = 0; i < N; i++)
        {
            T &obj{out_x[i]};
            for (auto &p : in_mat[i])
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
            std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: new_norm / prev_norm=" << new_norm / prev_norm << std::endl;
            return 2;
        }
        else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
        {
            return 0;
        }
        prev_x = std::move(out_x);
        prev_norm = new_norm;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: delta_norm / prev_norm=" << delta_norm / prev_norm << std::endl;
    return 1;
}

template <typename T, size_t N>
inline int jacobi(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                  std::array<T, N> &out_x,
                  bool sparse, size_t max_times, double rel_epsilon)
{
    if (sparse)
    {
        return jacobi(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x, sparse, max_times, rel_epsilon);
    }
    else
    {
        std::array<T, N> prev_x{std::move(out_x)};
        T delta_norm;
        T new_norm;
        T prev_norm{norm_1(prev_x)};
        static double not_conv_cri{2};

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

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: new_norm / prev_norm=" << new_norm / prev_norm << std::endl;
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
        std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: delta_norm / prev_norm=" << delta_norm / prev_norm << std::endl;
        return 1;
    }
}

template <typename T, size_t N>
inline int gauss_seidel(const Sparse_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                        std::array<T, N> &out_x,
                        bool sparse = true, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    std::array<T, N> prev_x{std::move(out_x)};
    T delta_norm;
    T new_norm;
    T prev_norm{norm_1(prev_x)};
    static double not_conv_cri{2};

    for (size_t count = 0; count < max_times; count++)
    {
        out_x = in_b;

        for (size_t i = 0; i < N; i++)
        {
            T &obj{out_x[i]};
            for (auto &p : in_mat[i])
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
            std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: new_norm / prev_norm=" << new_norm / prev_norm << std::endl;
            return 2;
        }
        else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
        {
            return 0;
        }
        prev_x = std::move(out_x);
        prev_norm = new_norm;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: delta_norm / prev_norm=" << delta_norm / prev_norm << std::endl;
    return 1;
}

template <typename T, size_t N>
inline int gauss_seidel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                        std::array<T, N> &out_x,
                        bool sparse, size_t max_times, double rel_epsilon)
{
    if (sparse)
    {
        return gauss_seidel(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x, sparse, max_times, rel_epsilon);
    }
    else
    {
        std::array<T, N> prev_x{std::move(out_x)};
        T delta_norm;
        T new_norm;
        T prev_norm{norm_1(prev_x)};
        static double not_conv_cri{2};

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

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: new_norm / prev_norm=" << new_norm / prev_norm << std::endl;
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
        std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: delta_norm / prev_norm=" << delta_norm / prev_norm << std::endl;
        return 1;
    }
}

template <typename T, size_t N>
inline int suc_over_rel(const Sparse_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                        std::array<T, N> &out_x, T omega,
                        bool sparse = true, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    if (omega <= 0 || omega >= 2)
    {
        throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": invalid omega");
    }

    std::array<T, N> prev_x{std::move(out_x)};
    T delta_norm;
    T new_norm;
    T prev_norm{norm_1(prev_x)};
    static double not_conv_cri{2};

    for (size_t count = 0; count < max_times; count++)
    {
        out_x = in_b;
        for (size_t i = 0; i < N; i++)
        {
            T &obj{out_x[i]};
            for (auto &p : in_mat[i])
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
            obj *= omega / in_mat(i, i);
            obj += (1 - omega) * prev_x[i];
        }

        delta_norm = norm_1(out_x, prev_x);
        new_norm = norm_1(out_x);

        if (prev_norm && new_norm / prev_norm > not_conv_cri)
        {
            std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: new_norm / prev_norm=" << new_norm / prev_norm << std::endl;
            return 2;
        }
        else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
        {
            return 0;
        }
        prev_x = std::move(out_x);
        prev_norm = new_norm;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: delta_norm / prev_norm=" << delta_norm / prev_norm << std::endl;
    return 1;
}

template <typename T, size_t N>
inline int suc_over_rel(const Base_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                        std::array<T, N> &out_x, T omega,
                        bool sparse, size_t max_times, double rel_epsilon)
{
    if (omega <= 0 || omega >= 2)
    {
        throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": invalid omega");
    }
    if (sparse)
    {
        return suc_over_rel(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x, omega, sparse, max_times, rel_epsilon);
    }
    else
    {
        std::array<T, N> prev_x{std::move(out_x)};
        T delta_norm;
        T new_norm;
        T prev_norm{norm_1(prev_x)};
        static double not_conv_cri{2};

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
                obj *= omega / in_mat(i, i);
                obj += (1 - omega) * prev_x[i];
            }
            delta_norm = norm_1(out_x, prev_x);
            new_norm = norm_1(out_x);

            if (prev_norm && new_norm / prev_norm > not_conv_cri)
            {
                std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: new_norm / prev_norm=" << new_norm / prev_norm << std::endl;
                return 2;
            }
            else if (prev_norm && delta_norm / prev_norm < rel_epsilon)
            {
                return 0;
            }
            prev_x = std::move(out_x);
            prev_norm = new_norm;
        }
        std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: delta_norm / prev_norm=" << delta_norm / prev_norm << std::endl;
        return 1;
    }
}

template <typename T, size_t N>
inline int _grad_des_sparse(const Sparse_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                            std::array<T, N> &out_x,
                            bool, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    T alpha;
    // calculate residue
    std::array<T, N> residue{in_b};
    // !!! this Ar is actually Ax
    std::array<T, N> Ar{in_mat * out_x};
    for (size_t i = 0; i < N; i++)
    {
        residue[i] -= Ar[i];
    }
    // finish residue
    // this is the true Ar
    Ar = in_mat * residue;

    const T res_b{norm_1(in_b)};
    T res_norm;

    for (size_t count = 0; count < max_times; count++)
    {
        alpha = (residue * residue) / (residue * Ar);

        for (size_t i = 0; i < N; i++)
        {
            out_x[i] += alpha * residue[i];
            residue[i] -= alpha * Ar[i];
        }

        res_norm = norm_1(residue);

        if (res_norm / res_b < rel_epsilon)
        {
            return 0;
        }

        Ar = in_mat * residue;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: res_norm / res_b=" << res_norm / res_b << std::endl;
    return 1;
}

template <typename T, size_t N>
inline int grad_des(const Symm_Matrix<T, N> &in_mat, const std::array<T, N> &in_b,
                    std::array<T, N> &out_x,
                    bool sparse, size_t max_times, double rel_epsilon)
{
    if (sparse)
    {
        return _grad_des_sparse(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x, sparse, max_times, rel_epsilon);
    }
    else
    {
        T alpha;
        // calculate residue
        std::array<T, N> residue{in_b};
        // !!! this Ar is actually Ax
        std::array<T, N> Ar{in_mat * out_x};
        for (size_t i = 0; i < N; i++)
        {
            residue[i] -= Ar[i];
        }
        // finish residue
        // this is the true Ar
        Ar = in_mat * residue;

        const T res_b{norm_1(in_b)};
        T res_norm;

        for (size_t count = 0; count < max_times; count++)
        {
            alpha = (residue * residue) / (residue * Ar);

            for (size_t i = 0; i < N; i++)
            {
                out_x[i] += alpha * residue[i];
                residue[i] -= alpha * Ar[i];
            }

            res_norm = norm_1(residue);

            if (res_norm / res_b < rel_epsilon)
            {
                return 0;
            }

            Ar = in_mat * residue;
        }
        std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: res_norm / res_b=" << res_norm / res_b << std::endl;
        return 1;
    }
}

template <typename T, size_t N>
inline int _conj_grad_sparse(const Sparse_Matrix<T, N, N> &in_mat, const std::array<T, N> &in_b,
                             std::array<T, N> &out_x,
                             bool, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    T alpha;
    // calculate residue
    std::array<T, N> residue{in_b};
    // !!! this Ap is actually Ax
    std::array<T, N> Ap{in_mat * out_x};
    for (size_t i = 0; i < N; i++)
    {
        residue[i] -= Ap[i];
    }
    // finish residue
    std::array<T, N> search_p{residue};
    // this is the true Ap
    Ap = in_mat * search_p;

    const T res_b{norm_1(in_b)};
    T res_norm;
    T rr{residue * residue};
    T prev_rr;
    T beta;

    for (size_t count = 0; count < max_times; count++)
    {
        alpha = rr / (residue * Ap);

        for (size_t i = 0; i < N; i++)
        {
            out_x[i] += alpha * search_p[i];
            residue[i] -= alpha * Ap[i];
        }

        res_norm = norm_1(residue);

        if (res_norm / res_b < rel_epsilon)
        {
            return 0;
        }

        prev_rr = rr;
        rr = residue * residue;

        beta = rr / prev_rr;

        for (size_t i = 0; i < N; i++)
        {
            search_p[i] = beta * search_p[i] + residue[i];
        }

        Ap = in_mat * search_p;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: res_norm / res_b=" << res_norm / res_b << std::endl;
    return 1;
}

/*beg:conj_grad_imp*/
template <typename T, size_t N>
inline int conj_grad(const Symm_Matrix<T, N> &in_mat, const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool sparse, size_t max_times, double rel_epsilon)
{
    if (sparse)
    {
        return _grad_des_sparse(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x, sparse, max_times, rel_epsilon);
    }
    else
    {
        T alpha;
        // calculate residue
        std::array<T, N> residue{in_b};
        // !!! this Ap is actually Ax
        std::array<T, N> Ap{in_mat * out_x};
        for (size_t i = 0; i < N; i++)
        {
            residue[i] -= Ap[i];
        }
        // finish residue
        std::array<T, N> search_p{residue};
        // this is the true Ap
        Ap = in_mat * search_p;

        const T res_b{norm_1(in_b)};
        T res_norm;
        T rr{residue * residue};
        T prev_rr;
        T beta;

        for (size_t count = 0; count < max_times; count++)
        {
            alpha = rr / (residue * Ap);

            for (size_t i = 0; i < N; i++)
            {
                out_x[i] += alpha * search_p[i];
                residue[i] -= alpha * Ap[i];
            }

            res_norm = norm_1(residue);

            if (res_norm / res_b < rel_epsilon)
            {
                return 0;
            }

            prev_rr = rr;
            rr = residue * residue;

            beta = rr / prev_rr;

            for (size_t i = 0; i < N; i++)
            {
                search_p[i] = beta * search_p[i] + residue[i];
            }

            Ap = in_mat * search_p;
        }
        std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: res_norm / res_b=" << res_norm / res_b << std::endl;
        return 1;
    }
}
/*end:conj_grad_imp*/

template <typename T, size_t N>
inline int conj_grad(const Hermite_Matrix<std::complex<T>, N> &in_mat, const std::array<std::complex<T>, N> &in_b,
                     std::array<std::complex<T>, N> &out_x,
                     bool sparse = false, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    std::unique_ptr<Sparse_Matrix<std::complex<T>, N, N>> p_sparse_mat{};
    if (sparse)
    {
        p_sparse_mat = std::make_unique<Sparse_Matrix<std::complex<T>, N, N>>(in_mat);
    }

    T alpha;
    // calculate residue
    std::array<std::complex<T>, N> residue{in_b};
    // !!! this Ap is actually Ax
    std::array<std::complex<T>, N> Ap{sparse ? (*p_sparse_mat) * out_x : in_mat * out_x};
    for (size_t i = 0; i < N; i++)
    {
        residue[i] -= Ap[i];
    }
    // finish residue
    std::array<std::complex<T>, N> search_p{residue};
    // this is the true Ap
    Ap = sparse ? (*p_sparse_mat) * search_p : in_mat * search_p;

    const T res_b{norm_1(in_b)};
    T res_norm;
    T rr{std::abs(residue * residue)};
    T prev_rr;
    std::complex<T> beta;

    for (size_t count = 0; count < max_times; count++)
    {
        alpha = rr / std::abs(residue * Ap);

        for (size_t i = 0; i < N; i++)
        {
            out_x[i] += alpha * search_p[i];
            residue[i] -= alpha * Ap[i];
        }

        res_norm = norm_1(residue);

        if (res_norm / res_b < rel_epsilon)
        {
            return 0;
        }

        prev_rr = rr;
        rr = std::abs(residue * residue);

        beta = rr / prev_rr;

        for (size_t i = 0; i < N; i++)
        {
            search_p[i] = beta * search_p[i] + residue[i];
        }

        Ap = sparse ? (*p_sparse_mat) * search_p : in_mat * search_p;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": iterative method not converged: res_norm / res_b=" << res_norm / res_b << std::endl;
    return 1;
}

template <typename T, size_t N, size_t M>
inline int grad_des(const Symm_Band_Matrix<T, N, M> &in_mat, const std::array<T, N> &in_b,
                    std::array<T, N> &out_x,
                    bool sparse = true, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    return _grad_des_sparse(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x,
                            true, max_times, rel_epsilon);
}

template <typename T, size_t N, size_t M>
inline int conj_grad(const Symm_Band_Matrix<T, N, M> &in_mat, const std::array<T, N> &in_b,
                     std::array<T, N> &out_x,
                     bool sparse = true, size_t max_times = 1000, double rel_epsilon = 1e-15)
{
    return _conj_grad_sparse(Sparse_Matrix<T, N, N>{in_mat}, in_b, out_x,
                             true, max_times, rel_epsilon);
}

} // namespace Misc

#endif // MISC_LINEAR_EQ_ITERATIVE
