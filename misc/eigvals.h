#ifndef MISC_EIGVALS
#define MISC_EIGVALS

#include <array>
#include <utility>
#include <vector>

#include "matrixlib.h"

namespace Misc
{

// ======================= APIs =============================

/* find the Hessenberg matrix of `in_mat`. `in_mat` == `q_mat` H `q_mat`^T, where `q_mat` is orthogonal
 * RETURN:
 *      H when `in_situ` is not passed, or else `in_mat` is modified;
 *      `q_mat` is modified.
 * PARAMETERS:
 *      in_mat: the matrix to be calculated. if `in_situ` is passed, then it will be modified
 *      q_mat:  optional; will be modified in situ
 *      in_situ:    its value is not used
 */
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat);
template <typename T, size_t N>
inline void hessenberg(Matrix<T, N, N> &in_mat, bool in_situ);
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat, Matrix<T, N, N> &q_mat);
template <typename T, size_t N>
inline void hessenberg(Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat, bool in_situ);

template <typename T, size_t N>
inline void _hessenberg_impl(Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr);

/* eigenvalues of a *Hessenberg matrix* H (`in_mat`), using shifted implicit QR
 * H = Q R Q^T, R is a block diagonal matrix
 * intended to be inner implementation
 * ASSUMPTION: H = Q0^T A Q0, where Q0 is `*p_q_mat` (if `p_q_mat` != nullptr)
 * RETURN:
 *      `in_mat` & `*p_q_mat` are modified.
 *      H = Q R Q^T = Q0^T A Q0, so A = Q0 Q R Q^T Q0^T. `*p_q_mat` would be Q0 Q
 * parameters:
 *      in_mat: the matrix H to be calculated. if `in_situ` is passed, then it will be modified
 *      p_q_mat:pointer to q_mat, optional; will be modified in situ
 *      tol:    relative tolerance
 */
template <typename T, size_t N>
inline void _eig_hessenberg(Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr,
                            const T tol = std::sqrt(std::numeric_limits<T>::epsilon()));

/* eigenvalues of a matrix A (`in_mat`),
 * modifying A into a Hessenberg Matrix and then using shifted implicit QR
 * RETURN:
 *      array of eigenvalues;
 *      vector of index where eigenvalues are complex pairs (Re first and Im second)
 *          for example, if 2 is in vector, then array[2] +- i array[3] are eigenvalues
 * parameters:
 *      in_mat: the matrix A to be calculated.
 */
template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(Matrix<T, N, N> in_mat,
                                                                 const T tol = std::sqrt(std::numeric_limits<T>::epsilon()));

// ======================= implementations ========================

template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat)
{
    _hessenberg_impl(in_mat);
    return in_mat;
}
template <typename T, size_t N>
inline void hessenberg(Matrix<T, N, N> &in_mat, bool in_situ)
{
    _hessenberg_impl(in_mat);
}
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat, Matrix<T, N, N> &q_mat)
{
    _hessenberg_impl(in_mat, &q_mat);
    return in_mat;
}
template <typename T, size_t N>
inline void hessenberg(Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat, bool in_situ)
{
    _hessenberg_impl(in_mat, &q_mat);
}

template <typename T, size_t N>
inline void _hessenberg_impl(Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat)
{
    // initialize Q as Identity
    if (p_q_mat)
    {
        *p_q_mat = Matrix<T, N, N>{};
        for (size_t i = 0; i < N; i++)
        {
            (*p_q_mat)(i, i) = 1.;
        }
    }

    std::array<T, N> u; // in_mat(i+1:, i) + sigma * e_{i+1}

    for (size_t i = 0; i < N - 2; i++)
    {
        u[i] = 0.;
        for (size_t j = i + 1; j < N; j++)
        {
            u[j] = in_mat(j, i);
        }
        T sigma{std::sqrt(u * u)}; // sign(in_mat(i+1, i)) * |in_mat(i+1:, i)|
        sigma *= std::signbit(u[i + 1]) ? -1. : 1.;
        u[i + 1] += sigma;
        T beta{sigma * u[i + 1]}; // |u|^2 / 2
        if (beta == 0.)
        {
            continue;
        }
        // left multiply H
        in_mat(i + 1, i) = -sigma;
        for (size_t j = i + 2; j < N; j++)
        {
            in_mat(j, i) = 0.;
        }
        for (size_t k = i + 1; k < N; k++)
        {
            T alpha{};
            for (size_t j = i + 1; j < N; j++)
            {
                alpha += in_mat(j, k) * u[j];
            }
            alpha /= beta;
            for (size_t j = i + 1; j < N; j++)
            {
                in_mat(j, k) -= alpha * u[j];
            }
        }
        // right multiply H
        for (size_t k = 0; k < N; k++)
        {
            T alpha{};
            for (size_t j = i + 1; j < N; j++)
            {
                alpha += in_mat(k, j) * u[j];
            }
            alpha /= beta;
            for (size_t j = i + 1; j < N; j++)
            {
                in_mat(k, j) -= alpha * u[j];
            }
        }
        // get Q
        if (p_q_mat)
            for (size_t k = 1; k < N; k++)
            {
                T alpha{};
                for (size_t j = i + 1; j < N; j++)
                {
                    alpha += (*p_q_mat)(k, j) * u[j];
                }
                alpha /= beta;
                for (size_t j = i + 1; j < N; j++)
                {
                    (*p_q_mat)(k, j) -= alpha * u[j];
                }
            }
    }
}

template <typename T, size_t N>
inline void _eig_hessenberg(Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat, const T tol)
{
    std::array<T, 3> u;
    size_t mat_range{N - 1}; // the maximum index of current sub-matrix

    while (mat_range >= 1)
    {
        // sigma^2 = det(2*2 square at right below corner)
        T sigma2{in_mat(mat_range, mat_range) * in_mat(mat_range - 1, mat_range - 1) - in_mat(mat_range, mat_range - 1) * in_mat(mat_range - 1, mat_range)};
        // 2*Re(sigma) = tr(2*2 square at right below corner)
        T re_2sigma{in_mat(mat_range, mat_range) + in_mat(mat_range - 1, mat_range - 1)};
        // initialize
        if (mat_range != 1)
        {
            // shifted implicit QR
            u[0] = in_mat(0, 0) * (in_mat(0, 0) - re_2sigma) + in_mat(0, 1) * in_mat(1, 0) + sigma2;
            u[1] = in_mat(1, 0) * (in_mat(0, 0) - re_2sigma + in_mat(1, 1));
            u[2] = in_mat(2, 1) * in_mat(1, 0);
        }
        // complex eigenvalue pair
        else
        {
            T delta{re_2sigma * re_2sigma - 4. * sigma2};
            if (delta < 0.)
            {
                break;
            }
            // real eigenvalues
            else
            {
                // shifted QR; the shift is the eigenvalue with larger abs
                u[0] = in_mat(0, 0) - .5 * re_2sigma - (std::signbit(re_2sigma) ? .5 : -.5) * std::sqrt(delta);
                u[1] = in_mat(1, 0);
                u[2] = 0.;
            }
        }
        size_t loop_count{3};
        for (size_t i = -1; i != mat_range - 1; i++)
        {
            if (i + 3 == N)
            // i + 3 will be out of range, so loop_count is set to 2
            {
                loop_count = 2;
                u[2] = 0.;
            }
            if (i != -1)
                for (size_t j = 0; j < loop_count; j++)
                {
                    u[j] = in_mat(i + j + 1, i);
                }
            T sigma{std::sqrt(u * u)}; // sign(in_mat(i+1, i)) * |in_mat(i+1:, i)|
            sigma *= std::signbit(u[0]) ? -1. : 1.;
            u[0] += sigma;
            T beta{sigma * u[0]}; // |u|^2 / 2
            if (beta == 0.)
            {
                continue;
            }
            // left multiply H
            if (i != -1)
            {
                in_mat(i + 1, i) = -sigma;
                in_mat(i + 2, i) = 0.;
                if (i + 3 != N)
                {
                    in_mat(i + 3, i) = 0.; // in case out of range
                }
            }
            for (size_t k = i + 1; k < N; k++)
            {
                T alpha{};
                for (size_t j = 0; j < loop_count; j++)
                {
                    alpha += in_mat(i + j + 1, k) * u[j];
                }
                alpha /= beta;
                for (size_t j = 0; j < loop_count; j++)
                {
                    in_mat(i + j + 1, k) -= alpha * u[j];
                }
            }

            // right multiply H
            size_t bound{std::min(i + 4, N - 1)}; // in case out of range
            for (size_t k = 0; k <= bound; k++)
            {
                T alpha{};
                for (size_t j = 0; j < loop_count; j++)
                {
                    alpha += in_mat(k, i + j + 1) * u[j];
                }
                alpha /= beta;
                for (size_t j = 0; j < loop_count; j++)
                {
                    in_mat(k, i + j + 1) -= alpha * u[j];
                }
            }

            // get Q
            if (p_q_mat)
                for (size_t k = 0; k < N; k++)
                {
                    T alpha{};
                    for (size_t j = 0; j < loop_count; j++)
                    {
                        alpha += (*p_q_mat)(k, i + j + 1) * u[j];
                    }
                    alpha /= beta;
                    for (size_t j = 0; j < loop_count; j++)
                    {
                        (*p_q_mat)(k, i + j + 1) -= alpha * u[j];
                    }
                }
        }
        // check convergence
        T criteria{tol * std::max({1., std::abs(in_mat(mat_range, mat_range)), std::abs(in_mat(mat_range - 1, mat_range - 1))})};
        if (std::abs(in_mat(mat_range, mat_range - 1)) < criteria)
        {
            in_mat(mat_range, mat_range - 1) = 0.;
            mat_range--;
        }
        else if (mat_range >= 2 && std::abs(in_mat(mat_range - 1, mat_range - 2)) < criteria)
        {
            in_mat(mat_range - 1, mat_range - 2) = 0.;
            mat_range -= 2;
        }
    }
}

template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(Matrix<T, N, N> in_mat, const T tol)
{
    hessenberg(in_mat, true);
    _eig_hessenberg(in_mat, (Matrix<T, N, N> *)nullptr, tol);
    std::array<T, N> res_array;
    std::vector<size_t> res_vector{};

    for (size_t i = 0; i < N - 1; i++)
    {
        if (in_mat(i + 1, i) == 0.)
        {
            res_array[i] = in_mat(i, i);
        }
        else
        {
            T sigma2{in_mat(i, i) * in_mat(i + 1, i + 1) - in_mat(i, i + 1) * in_mat(i + 1, i)};
            T re_2sigma{in_mat(i, i) + in_mat(i + 1, i + 1)};
            T delta{re_2sigma * re_2sigma - 4. * sigma2};
            if (delta < 0.)
            {
                res_vector.push_back(i);
                res_array[i] = re_2sigma / 2.;
                res_array[i + 1] = std::sqrt(-delta) / 2.;
            }
            else
            {
                delta = std::sqrt(delta);
                res_array[i] = (re_2sigma - delta) / 2.;
                res_array[i + 1] = (re_2sigma + delta) / 2.;
            }
            i++;
        }
    }
    if (in_mat(N - 1, N - 2) == 0.)
    {
        res_array[N - 1] = in_mat(N - 1, N - 1);
    }
    return {res_array, res_vector};
}

} // namespace Misc

#endif // MISC_EIGVALS
