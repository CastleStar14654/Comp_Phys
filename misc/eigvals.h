#ifndef MISC_EIGVALS
#define MISC_EIGVALS

#include <array>
#include <utility>
#include <vector>
#include <random>
#include <thread>
#include <functional>
#include <stdexcept>

#include "matrixlib.h"
#include "rootfind.h"

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
 * for Sparse_Matrix, Arnold algorithms is used;
 * for Matrix or Base_Matrix, Householder transformation instead.
 */
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat);
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(const Base_Matrix<T, N, N> &in_mat);
template <typename T, size_t N>
inline void hessenberg(Matrix<T, N, N> &in_mat, bool in_situ);
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat, Matrix<T, N, N> &q_mat);
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(const Base_Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat);
template <typename T, size_t N>
inline void hessenberg(Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat, bool in_situ);

template <typename T, size_t N>
inline void _hessenberg_impl(Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr);

template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(const Sparse_Matrix<T, N, N> &in_mat);
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Sparse_Matrix<T, N, N> &in_mat, bool symmetrical);
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat);
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat, bool symmetrical);

template <typename T, size_t N>
inline Matrix<T, N, N> _hessenberg_impl(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr);

/* find the Hessenberg matrix of a symmetrical `in_mat`. `in_mat` == `q_mat` H `q_mat`^T, where `q_mat` is orthogonal
 * RETURN:
 *      a Symm_Band_Matrix H;
 *      `q_mat` is modified.
 * PARAMETERS:
 *      in_mat: the matrix to be calculated. if `in_situ` is passed, then it will be modified
 *      q_mat:  optional; will be modified in situ
 *      in_situ:    its value is not used
 * for Sparse_Matrix, Lanczos algorithm is used;
 * for Symm_Matrix or Symm_Band_Matrix, Householder transformation instead.
 */
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(Symm_Matrix<T, N> in_mat);
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(Symm_Matrix<T, N> in_mat, Matrix<T, N, N> &q_mat);

template <typename T, size_t N, size_t M>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Symm_Band_Matrix<T, N, M> &in_mat);
template <typename T, size_t N, size_t M>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Symm_Band_Matrix<T, N, M> &in_mat, Matrix<T, N, N> &q_mat);

template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> _hessenberg_impl(Symm_Matrix<T, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr);
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> _hessenberg_symm_impl(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr);

/* eigenvalues of a *Hessenberg matrix* H (`in_mat`),
 * using Francis shifted implicit QR.
 * H = Q R Q^T, R is a block diagonal matrix
 * intended to be inner implementation
 * ASSUMPTION: H = Q0^T A Q0, where Q0 is `*p_q_mat` (if `p_q_mat` != nullptr)
 * RETURN:
 *      `in_mat` & `*p_q_mat` are modified.
 *      H = Q R Q^T = Q0^T A Q0, so A = Q0 Q R Q^T Q0^T. `*p_q_mat` would be Q0 Q
 * parameters:
 *      in_mat: the matrix H to be calculated.
 *      p_q_mat:pointer to q_mat, optional; will be modified in situ
 *      tol:    relative tolerance
 */
template <typename T, size_t N>
inline void _eig_hessenberg(Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr,
                            const T tol = std::sqrt(std::numeric_limits<T>::epsilon()));
/* eigenvalues of a *tridiagonal Hessenberg matrix* H (`in_mat`),
 * using QR with implicit Wilkinson shift.
 * H = Q R Q^T, R is a diagonal matrix
 * intended to be inner implementation
 * ASSUMPTION: H = Q0^T A Q0, where Q0 is `*p_q_mat` (if `p_q_mat` != nullptr)
 * RETURN:
 *      `in_mat` & `*p_q_mat` are modified.
 *      H = Q R Q^T = Q0^T A Q0, so A = Q0 Q R Q^T Q0^T. `*p_q_mat` would be Q0 Q
 * parameters:
 *      in_mat: the matrix H to be calculated.
 *      p_q_mat:pointer to q_mat, optional; will be modified in situ
 *      tol:    relative tolerance
 */
template <typename T, size_t N>
inline void _eig_hessenberg_symm_qr(Symm_Band_Matrix<T, N, 1> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr,
                                    const T tol = std::sqrt(std::numeric_limits<T>::epsilon()));
/* eigenvalues of a *tridiagonal Hessenberg matrix* H (`in_mat`),
 * using bisection method.
 * intended to be inner implementation
 * RETURN:
 *      array of eigenvalues
 * parameters:
 *      in_mat: the matrix H to be calculated.
 *      tol:    relative tolerance
 */
template <typename T, size_t N>
inline std::array<T, N> _eig_hessenberg_symm_bisection(const Symm_Band_Matrix<T, N, 1> &in_mat,
                                                       const T tol = std::numeric_limits<T>::epsilon());
/* eigenvalues of a *tridiagonal Hessenberg matrix* H (`in_mat`),
 * using QR with implicit Wilkinson shift.
 * H = Q R Q^T, R is a diagonal matrix
 * intended to be inner implementation
 * ASSUMPTION: H = Q0^T A Q0, where Q0 is `*p_q_mat` (if `p_q_mat` != nullptr)
 * RETURN:
 *      `in_mat` & `*p_q_mat` are modified.
 *      H = Q R Q^T = Q0^T A Q0, so A = Q0 Q R Q^T Q0^T. `*p_q_mat` would be Q0 Q
 * parameters:
 *      in_mat: the matrix H to be calculated.
 *      p_q_mat:pointer to q_mat, optional; will be modified in situ
 *      tol:    relative tolerance
 */
template <typename T, size_t N>
inline std::array<T, N> _eig_hessenberg_symm_div_conq(Symm_Band_Matrix<T, N, 1> &in_mat, Matrix<T, N, N> *p_q_mat = nullptr,
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
template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(const Base_Matrix<T, N, N> &in_mat,
                                                                 const T tol = std::sqrt(std::numeric_limits<T>::epsilon()));
template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(const Sparse_Matrix<T, N, N> &in_mat,
                                                                 const T tol = std::sqrt(std::numeric_limits<T>::epsilon()));

template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> _eig_vals_impl(Matrix<T, N, N> &hessenberg_mat,
                                                                       const T tol);

enum class Eig_Symm_Method
{
    QR,
    bisection,
    divide_conquer
};

/* eigenvalues of a symmetrical matrix A (`in_mat`),
 * modifying A into a Hessenberg Matrix and then using shifted implicit QR
 * RETURN:
 *      array of eigenvalues;
 *      vector of index where eigenvalues are complex pairs (Re first and Im second)
 *          for example, if 2 is in vector, then array[2] +- i array[3] are eigenvalues
 * parameters:
 *      in_mat: the matrix A to be calculated.
 */
template <typename T, size_t N>
inline std::array<T, N> eig_vals_symm(const Sparse_Matrix<T, N, N> &in_mat,
                                      const Eig_Symm_Method method = Eig_Symm_Method::QR,
                                      const T tol = std::numeric_limits<T>::epsilon());
template <typename T, size_t N>
inline std::array<T, N> eig_vals_symm(const Symm_Matrix<T, N> &in_mat,
                                      const Eig_Symm_Method method = Eig_Symm_Method::QR,
                                      const T tol = std::numeric_limits<T>::epsilon());
template <typename T, size_t N, size_t M>
inline std::array<T, N> eig_vals_symm(const Symm_Band_Matrix<T, N, M> &in_mat,
                                      const Eig_Symm_Method method = Eig_Symm_Method::QR,
                                      const T tol = std::numeric_limits<T>::epsilon());

template <typename T, size_t N>
inline std::array<T, N> _eig_vals_symm_impl(Symm_Band_Matrix<T, N, 1> &hessenberg_mat,
                                            const Eig_Symm_Method method, const T tol);
// ======================= implementations ========================

template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(Matrix<T, N, N> in_mat)
{
    _hessenberg_impl(in_mat);
    return in_mat;
}
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(const Base_Matrix<T, N, N> &in_mat)
{
    Matrix<T, N, N> res{in_mat};
    _hessenberg_impl(res);
    return res;
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
inline Matrix<T, N, N> hessenberg(const Base_Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat)
{
    Matrix<T, N, N> res{in_mat};
    _hessenberg_impl(res, &q_mat);
    return res;
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
inline Matrix<T, N, N> hessenberg(const Sparse_Matrix<T, N, N> &in_mat)
{
    return _hessenberg_impl(in_mat);
}
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Sparse_Matrix<T, N, N> &in_mat, bool symmetrical)
{
    return _hessenberg_symm_impl(in_mat);
}
template <typename T, size_t N>
inline Matrix<T, N, N> hessenberg(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat)
{
    return _hessenberg_impl(in_mat, &q_mat);
}
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> &q_mat, bool symmetrical)
{
    return _hessenberg_symm_impl(in_mat, &q_mat);
}

template <typename T, size_t N>
inline Matrix<T, N, N> _hessenberg_impl(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat)
{
    // initialize out_mat
    Matrix<T, N, N> out_mat{};
    // initialize Q
    bool need_delete{false};
    static auto ran_eng{std::default_random_engine{}};
    static auto rand_double{std::uniform_real_distribution<T>{0., 1.}};
    if (!p_q_mat)
    {
        p_q_mat = new Matrix<T, N, N>{};
        need_delete = true;
    }
    std::array<T, N> r;
    r.fill(0.);
    r[0] = 1.;

    for (size_t i = 0; i < N; i++)
    {
        // add q
        (*p_q_mat).column(i) = r;
        r = in_mat * (*p_q_mat).column(i);
        for (size_t j = 0; j <= i; j++)
        {
            out_mat(j, i) = r * (*p_q_mat).column(j);
        }
        // get q_{i+1} from r
        for (size_t j = 0; j <= i; j++)
        {
            T inner_prod;
            inner_prod = j ? (r * (*p_q_mat).column(j)) : out_mat(0, i);
            for (size_t k = 0; k < N; k++)
            {
                r[k] -= inner_prod * (*p_q_mat)(k, j);
            }
        }
        if (i + 1 < N)
        {
            T q_abs{std::sqrt(r * r)};
            out_mat(i + 1, i) = q_abs;
            while (q_abs < 1. / (1 << 20))
            {
                for (auto &x : r)
                {
                    x = rand_double(ran_eng);
                }
                // Grand-Schmidt
                for (size_t j = 0; j <= i; j++)
                {
                    T inner_prod;
                    inner_prod = r * (*p_q_mat).column(j);
                    for (size_t k = 0; k < N; k++)
                    {
                        r[k] -= inner_prod * (*p_q_mat)(k, j);
                    }
                }
                q_abs = std::sqrt(r * r);
            }
            for (size_t k = 0; k < N; k++)
            {
                r[k] /= q_abs;
            }
        }
    }
    if (need_delete)
    {
        delete p_q_mat;
    }
    return out_mat;
}

// ------------------------------------------

template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(Symm_Matrix<T, N> in_mat)
{
    return _hessenberg_impl(in_mat);
}
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> hessenberg(Symm_Matrix<T, N> in_mat, Matrix<T, N, N> &q_mat)
{
    return _hessenberg_impl(in_mat, &q_mat);
}

template <typename T, size_t N, size_t M>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Symm_Band_Matrix<T, N, M> &in_mat)
{
    if (M == 1)
    {
        return Symm_Band_Matrix<T, N, 1>{in_mat};
    }
    else
    {
        Symm_Matrix<T, N> temp{in_mat};
        return _hessenberg_impl(temp);
    }
}
template <typename T, size_t N, size_t M>
inline Symm_Band_Matrix<T, N, 1> hessenberg(const Symm_Band_Matrix<T, N, M> &in_mat, Matrix<T, N, N> &q_mat)
{
    if (M == 1)
    {
        q_mat = Matrix<T, N, N>{};
        for (size_t i = 0; i < N; i++)
        {
            q_mat(i, i) = 1.;
        }
        return Symm_Band_Matrix<T, N, 1>{in_mat};
    }
    else
    {
        Symm_Matrix<T, N> temp{in_mat};
        return _hessenberg_impl(temp, &q_mat);
    }
}

template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> _hessenberg_impl(Symm_Matrix<T, N> &in_mat, Matrix<T, N, N> *p_q_mat)
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
    std::array<T, N> p;

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
        // left and right multiply H
        in_mat(i + 1, i) = -sigma;
        std::fill(&p[i], p.end(), 0.);
        for (size_t j = i + 1; j < N; j++)
        {
            for (size_t k = i + 1; k < N; k++)
            {
                p[j] += const_cast<const Symm_Matrix<T, N> &>(in_mat)(j, k) * u[k];
            }
        }
        T coef{u * p / (2. * beta)};
        for (size_t j = i + 1; j < N; j++)
        {
            p[j] -= coef * u[j];
            p[j] /= beta;
        }
        for (size_t j = i + 1; j < N; j++)
            for (size_t k = j; k < N; k++)
            {
                in_mat(k, j) -= p[k] * u[j] + p[j] * u[k];
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

    Symm_Band_Matrix<T, N, 1> res{};
    for (size_t i = 0; i < N; i++)
    {
        res(i, i) = in_mat(i, i);
    }
    for (size_t i = 0; i < N - 1; i++)
    {
        res(i + 1, i) = in_mat(i + 1, i);
    }
    return res;
}
template <typename T, size_t N>
inline Symm_Band_Matrix<T, N, 1> _hessenberg_symm_impl(const Sparse_Matrix<T, N, N> &in_mat, Matrix<T, N, N> *p_q_mat)
{
    // initialize out_mat
    Symm_Band_Matrix<T, N, 1> out_mat{};
    // initialize Q
    bool need_delete{false};
    static auto ran_eng{std::default_random_engine{}};
    static auto rand_double{std::uniform_real_distribution<T>{0., 1.}};
    if (!p_q_mat)
    {
        p_q_mat = new Matrix<T, N, N>{};
        need_delete = true;
    }
    std::array<T, N> r;
    r.fill(0.);
    r[0] = 1.;

    for (size_t i = 0; i < N; i++)
    {
        // add q
        (*p_q_mat).column(i) = r;
        r = in_mat * (*p_q_mat).column(i);
        out_mat(i, i) = r * (*p_q_mat).column(i);
        // get q_{i+1} from r
        for (size_t j = i ? i - 1 : i; j <= i; j++)
        {
            T inner_prod;
            inner_prod = (j == i) ? (r * (*p_q_mat).column(j)) : out_mat(i, i - 1);
            for (size_t k = 0; k < N; k++)
            {
                r[k] -= inner_prod * (*p_q_mat)(k, j);
            }
        }
        if (i + 1 < N)
        {
            T q_abs{std::sqrt(r * r)};
            out_mat(i + 1, i) = q_abs;
            while (q_abs < 1. / (1 << 20))
            {
                for (auto &x : r)
                {
                    x = rand_double(ran_eng);
                }
                // Grand-Schmidt
                for (size_t j = 0; j <= i; j++)
                {
                    T inner_prod;
                    inner_prod = r * (*p_q_mat).column(j);
                    for (size_t k = 0; k < N; k++)
                    {
                        r[k] -= inner_prod * (*p_q_mat)(k, j);
                    }
                }
                q_abs = std::sqrt(r * r);
            }
            for (size_t k = 0; k < N; k++)
            {
                r[k] /= q_abs;
            }
        }
    }
    if (need_delete)
    {
        delete p_q_mat;
    }
    return out_mat;
}

// ------------------------------------------

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
/*beg:symm_qr*/
template <typename T, size_t N>
inline void _eig_hessenberg_symm_qr_(Symm_Band_Matrix<T, N, 1> &in_mat, Matrix<T, N, N> *p_q_mat,
                                     const T tol, const size_t start = 0, const size_t end = N)
{
    std::vector<size_t> zero_break{start};
    for (size_t i = start + 1; i < end; i++)
        if (std::abs(in_mat(i, i - 1)) < tol * std::max({1.,
                                                         std::abs(in_mat(i, i)),
                                                         std::abs(in_mat(i - 1, i - 1))}))
        {
            in_mat(i, i - 1) = 0.;
            zero_break.push_back(i);
        }
    size_t mat_range{end - 1}; // the maximum index of current sub-matrix
    size_t lower_bound{zero_break.back()};
    zero_break.pop_back();

    while (mat_range >= start + 1)
    {
        if (mat_range <= lower_bound)
        {
            mat_range = lower_bound;
            lower_bound = zero_break.back();
            zero_break.pop_back();
        }
        // d = (a_{n-1, n-1} - a_{n, n}) / 2.
        T d{(in_mat(mat_range - 1, mat_range - 1) - in_mat(mat_range, mat_range)) / 2.};
        T shift{d
                    ? in_mat(mat_range, mat_range - 1) / d
                    : in_mat(mat_range, mat_range) - std::abs(in_mat(mat_range, mat_range - 1))};
        if (d)
        {
            shift *= shift;
            shift = in_mat(mat_range, mat_range) - d * shift / (1. + std::sqrt(1. + shift));
        }
        // initialize
        T bulge{in_mat(lower_bound + 1, lower_bound)};
        T offset;
        for (size_t i = lower_bound; i < mat_range; i++)
        {
            if (i == lower_bound)
            {
                offset = in_mat(lower_bound, lower_bound) - shift;
            }
            else
            {
                offset = in_mat(i, i - 1);
            }
            T sin_theta{std::sqrt(bulge * bulge + offset * offset)};
            T cos_theta{std::abs(offset) / sin_theta};
            sin_theta = bulge / sin_theta;
            if (std::signbit(offset))
            {
                sin_theta *= -1.;
            }

            // iterate
            if (i != lower_bound)
            {
                in_mat(i, i - 1) = cos_theta * offset + sin_theta * bulge;
            }
            T delta{in_mat(i + 1, i + 1) - in_mat(i, i)};
            T inc{sin_theta * (2 * in_mat(i + 1, i) * cos_theta + delta * sin_theta)};
            in_mat(i, i) += inc;
            in_mat(i + 1, i + 1) -= inc;
            in_mat(i + 1, i) *= (cos_theta - sin_theta) * (cos_theta + sin_theta);
            in_mat(i + 1, i) += delta * cos_theta * sin_theta;
            if (i + 1 != mat_range)
            {
                bulge = sin_theta * in_mat(i + 2, i + 1);
                in_mat(i + 2, i + 1) *= cos_theta;
            }

            // get Q
            if (p_q_mat)
                for (size_t k = start; k < end; k++)
                {
                    T temp{(*p_q_mat)(k, i)};
                    (*p_q_mat)(k, i) *= cos_theta;
                    (*p_q_mat)(k, i) += sin_theta * (*p_q_mat)(k, i + 1);
                    (*p_q_mat)(k, i + 1) *= cos_theta;
                    (*p_q_mat)(k, i + 1) -= sin_theta * temp;
                }
        }
        // check convergence
        T criteria{tol * std::max({1.,
                                   std::abs(in_mat(mat_range, mat_range)),
                                   std::abs(in_mat(mat_range - 1, mat_range - 1))})};
        if (std::abs(in_mat(mat_range, mat_range - 1)) < criteria)
        {
            in_mat(mat_range, mat_range - 1) = 0.;
            mat_range--;
        }
    }
}
/*end:symm_qr*/
template <typename T, size_t N>
inline void _eig_hessenberg_symm_qr(Symm_Band_Matrix<T, N, 1> &in_mat, Matrix<T, N, N> *p_q_mat, const T tol)
{
    _eig_hessenberg_symm_qr_(in_mat, p_q_mat, tol);
}
/*beg:symm_bis*/
template <typename T, size_t N>
inline std::array<T, N> _eig_hessenberg_symm_bisection(const Symm_Band_Matrix<T, N, 1> &in_mat, const T tol)
{
    static auto sturm_change_of_sign{
        [&](T x) -> size_t {
            T q{in_mat(0, 0) - x};
            size_t res{std::signbit(q) ? 1ULL : 0ULL};
            for (size_t i = 1; i < N; i++)
            {
                if (q == 0.)
                {
                    q = -std::numeric_limits<double>::infinity();
                }
                else
                {
                    q = std::isinf(q)
                            ? in_mat(i, i) - x
                            : in_mat(i, i) - x - in_mat(i, i - 1) * in_mat(i, i - 1) / q;
                }
                if (std::signbit(q))
                // q < +0.0
                {
                    res++;
                }
            }
            return res;
        }};
    std::array<T, N> res;
    // calc norm_inf
    T norm_inf{std::abs(in_mat(0, 0)) + std::abs(in_mat(0, 1))};
    for (size_t i = 1; i < N; i++)
    {
        T temp{std::abs(in_mat(i, i - 1))};
        temp += std::abs(in_mat(i, i));
        if (i != N - 1)
            temp += std::abs(in_mat(i, i + 1));
        if (temp > norm_inf)
        {
            norm_inf = temp;
        }
    }
    T left_bound{-norm_inf};
    for (size_t i = 1; i <= N; i++)
    {
        T right_bound{norm_inf};
        T del{(right_bound - left_bound) / 2.};
        while (del > tol * std::max(std::abs(right_bound), std::abs(right_bound)))
        {
            T mid{left_bound + del};
            if (sturm_change_of_sign(mid) < i)
            {
                left_bound = mid;
            }
            else
            {
                right_bound = mid;
            }
            del /= 2.;
        }
        left_bound += del;
        res[i - 1] = left_bound;
    }
    return res;
}
/*end:symm_bis*/
template <typename T, size_t N>
inline std::array<T, N> _eig_hessenberg_symm_div_conq(Symm_Band_Matrix<T, N, 1> &in_mat, Matrix<T, N, N> *p_q_mat, const T tol)
{
    static_assert(N >= 4);
    Matrix<T, N, N> q_mat{};
    Matrix<T, N, N> alter_q_mat{};
    for (size_t i = 0; i < N; i++)
    {
        q_mat(i, i) = 1.;
    }
    std::array<T, N> res;
    constexpr std::array<size_t, 4> mats_end{(N / 2) / 2, N / 2, N / 2 + (N - N / 2) / 2, N};
    constexpr T theta{1.};
    std::array<T, 3> rhos{
        in_mat(mats_end[0], mats_end[0] - 1) / theta,
        in_mat(mats_end[1], mats_end[1] - 1) / theta,
        in_mat(mats_end[2], mats_end[2] - 1) / theta,
    };

    std::array<T, N> z;
    z.fill(0.);
    z[mats_end[0] - 1] = 1.;
    z[mats_end[0]] = theta;
    z[mats_end[2] - 1] = 1.;
    z[mats_end[2]] = theta;
    // divide
    in_mat(mats_end[0] - 1, mats_end[0] - 1) -= rhos[0];
    in_mat(mats_end[0], mats_end[0] - 1) = 0.;
    std::thread t11(_eig_hessenberg_symm_qr_<T, N>, std::ref(in_mat), &q_mat, tol, 0, mats_end[0]);
    in_mat(mats_end[0], mats_end[0]) -= rhos[0] * (theta * theta);
    in_mat(mats_end[1] - 1, mats_end[1] - 1) -= rhos[1];
    in_mat(mats_end[1], mats_end[1] - 1) = 0.;
    std::thread t12(_eig_hessenberg_symm_qr_<T, N>, std::ref(in_mat), &q_mat, tol, mats_end[0], mats_end[1]);
    in_mat(mats_end[1], mats_end[1]) -= rhos[1] * (theta * theta);
    in_mat(mats_end[2] - 1, mats_end[2] - 1) -= rhos[2];
    in_mat(mats_end[2], mats_end[2] - 1) = 0.;
    std::thread t13(_eig_hessenberg_symm_qr_<T, N>, std::ref(in_mat), &q_mat, tol, mats_end[1], mats_end[2]);
    in_mat(mats_end[2], mats_end[2]) -= rhos[2] * (theta * theta);
    std::thread t14(_eig_hessenberg_symm_qr_<T, N>, std::ref(in_mat), &q_mat, tol, mats_end[2], mats_end[3]);
    t11.join();
    t12.join();
    t13.join();
    t14.join();
    // glue small
    z = z * q_mat;
    for (size_t i = 0; i < N; i++)
    {
        res[i] = in_mat(i, i);
    }
    std::sort(res.begin(), &res[mats_end[1]]);
    std::sort(&res[mats_end[1]], res.end());
    for (auto x : res)
    {
        std::cerr << x << ',';
    }
    std::cerr << '\n';

    std::cerr << alter_q_mat << std::endl;

    struct Glue
    {
        const Symm_Band_Matrix<T, N, 1> &in_mat;
        std::array<T, N> diag;
        std::array<T, N> &res;
        const std::array<T, N> &z;
        const Matrix<T, N, N> &q_mat;
        Matrix<T, N, N> &new_q_mat;
        std::array<T, N> &q_temp_array;
        T rho;
        size_t start;
        size_t end;
        size_t mid;
        Glue(const Symm_Band_Matrix<T, N, 1> &in_mat, std::array<T, N> &res, const std::array<T, N> &z,
             const Matrix<T, N, N> &q_mat, Matrix<T, N, N> &new_q_mat, std::array<T, N> &q_temp_array, T rho, size_t start, size_t end)
            : in_mat{in_mat}, diag{res}, res{res}, z{z}, q_mat{q_mat}, new_q_mat{new_q_mat}, q_temp_array{q_temp_array}, rho{rho}, start{start}, end{end}, mid{start + (end - start) / 2}
        {
        }
        T f(T lambda)
        {
            T ans{0.};
            for (size_t i = start; i < end; i++)
            {
                ans += z[i] * z[i] / (in_mat(i, i) - lambda);
            }
            ans *= rho;
            ans += 1.;
            return ans;
        }
        void less_half()
        {
            std::function<T(T)> f_lmd{
                [&](T lambda) {
                    return f(lambda);
                }};
            if (rho)
            {
                if (rho < 0.)
                {
                    T del = diag[start + 1] - diag[start];
                    T f_right{f(diag[start] - 3.0517578125e-05 * std::max(del, .5))};
                    // T left_bound{diag[start] - del};
                    // while (std::signbit(f(left_bound)) == std::signbit(f_right))
                    // {
                    //     left_bound -= del;
                    // }
                    res[start] = dekker_brent(f_lmd, diag[start] - z * z, diag[start] - 3.0517578125e-05 * del);
                }
                for (size_t i = rho < 0. ? start + 1 : start, j = start; i < mid; i++, j++)
                {
                    T del{diag[j + 1] - diag[j]};
                    if (del > 2. * std::numeric_limits<T>::epsilon() * std::max(std::abs(diag[j + 1]), std::abs(diag[j])))
                    {
                        del *= 3.0517578125e-05;
                        T del_left{std::max(del, std::numeric_limits<T>::epsilon() * std::abs(diag[j]))};
                        T del_right{std::max(del, std::numeric_limits<T>::epsilon() * std::abs(diag[j + 1]))};
                        res[i] = dekker_brent(f_lmd, diag[j] + del, diag[j + 1] - del);
                    }
                    else
                    {
                        res[i] = diag[j];
                    }
                }
            }
        }
        void greater_half()
        {
            std::function<T(T)> f_lmd{
                [&](T lambda) {
                    return f(lambda);
                }};
            if (rho)
            {
                if (rho > 0.)
                {
                    T del = diag[end - 1] - diag[end - 2];
                    T f_left{f(diag[end - 1] + 3.0517578125e-05 * std::max(del, .5))};
                    // T right_bound{diag[end - 1] + del};
                    // while (std::signbit(f(right_bound)) == std::signbit(f_left))
                    // {
                    //     right_bound += del;
                    // }
                    res[end - 1] = dekker_brent(f_lmd, diag[end - 1] + 3.0517578125e-05 * del, diag[end - 1] + z * z);
                }
                for (size_t i = mid, j = rho > 0. ? mid : mid - 1; j < end - 1; i++, j++)
                {
                    T del{diag[j + 1] - diag[j]};
                    if (del > 2. * std::numeric_limits<T>::epsilon() * std::max(std::abs(diag[j + 1]), std::abs(diag[j])))
                    {
                        del *= 3.0517578125e-05;
                        T del_left{std::max(del, std::numeric_limits<T>::epsilon() * std::abs(diag[j]))};
                        T del_right{std::max(del, std::numeric_limits<T>::epsilon() * std::abs(diag[j + 1]))};
                        res[i] = dekker_brent(f_lmd, diag[j] + del, diag[j + 1] - del);
                    }
                    else
                    {
                        res[i] = diag[j];
                    }
                }
            }
        }
        // refresh q_mat in less half or greater half
        void refresh_q_mat(bool less)
        {
            size_t ref_start{less ? start : mid};
            size_t ref_end{less ? mid : end};
            for (size_t j = ref_start; j < ref_end; j++)
            {
                T lambda{res[j]};
                T criteria{std::numeric_limits<T>::epsilon() * std::abs(lambda)};
                std::vector<size_t> i_s_equal{};
                for (size_t i = start; i < end; i++)
                {
                    if (std::abs(in_mat(i, i) - res[j]) > criteria || std::find(i_s_equal.begin(), i_s_equal.end(), i) != i_s_equal.end())
                    {
                        new_q_mat(i, j) = z[i] / (in_mat(i, i) - res[j]);
                    }
                    else
                    {
                        for (size_t i_ = start; i_ < end; i_++)
                        {
                            new_q_mat(i_, j) = 0.;
                        }
                        new_q_mat(i, j) = 1.;
                    }
                    // if (z[i])
                    // {
                    //     new_q_mat(i, j) = z[i] / (d - res[j]);
                    // }
                    // else
                    // {
                    //     new_q_mat(i, j) = i == j ? 1. : 0.;
                    // }
                }
            }
            for (size_t j = ref_start; j < ref_end; j++)
            {
                T norm_inv{0.};
                for (size_t i = start; i < end; i++)
                {
                    norm_inv += new_q_mat(i, j) * new_q_mat(i, j);
                }
                norm_inv = 1. / std::sqrt(norm_inv);
                for (size_t i = start; i < end; i++)
                {
                    new_q_mat(i, j) *= norm_inv;
                }
            }
            for (size_t j = ref_start; j < ref_end; j++)
            {
                std::fill(&q_temp_array[start], &q_temp_array[end], 0.);
                for (size_t k = start; k < end; k++)
                {
                    T q_k_j{new_q_mat(k, j)};
                    for (size_t i = start; i < end; i++)
                    {
                        q_temp_array[i] += q_mat(i, k) * q_k_j;
                    }
                }
                for (size_t i = start; i < end; i++)
                {
                    new_q_mat(i, j) = q_temp_array[i];
                }
            }
        }
    };

    std::array<T, N> q_temp_array; // used when calculating new q_mat
    {
        Glue up{in_mat, res, z, q_mat, alter_q_mat, q_temp_array, rhos[0], 0, mats_end[1]};
        Glue down{in_mat, res, z, q_mat, alter_q_mat, q_temp_array, rhos[2], mats_end[1], N};
        std::thread t21(up.less_half, &up);
        std::thread t22(up.greater_half, &up);
        std::thread t23(down.less_half, &down);
        std::thread t24(down.greater_half, &down);
        t21.join();
        t22.join();
        t23.join();
        t24.join();
        std::thread t25(up.refresh_q_mat, &up, true);
        std::thread t26(down.refresh_q_mat, &down, true);
        t25.join();
        t26.join();
        std::thread t27(up.refresh_q_mat, &up, false);
        std::thread t28(down.refresh_q_mat, &down, false);
        t27.join();
        t28.join();
    }
    // glue big
    for (size_t i = 0; i < N; i++)
    {
        in_mat(i, i) = res[i];
    }
    std::sort(res.begin(), res.end());
    z.fill(0.);
    z[mats_end[1] - 1] = 1.;
    z[mats_end[1]] = theta;
    for (auto x : res)
    {
        std::cerr << x << ',';
    }
    std::cerr << '\n';

    std::cerr << alter_q_mat << std::endl;
    z = z * alter_q_mat;

    Glue whole{in_mat, res, z, alter_q_mat, q_mat, q_temp_array, rhos[1], 0, N};
    std::thread t31(whole.less_half, &whole);
    std::thread t32(whole.greater_half, &whole);
    t31.join();
    t32.join();

    if (p_q_mat)
    {
        whole.refresh_q_mat(true);
        whole.refresh_q_mat(false);
        (*p_q_mat) = (*p_q_mat) * q_mat;
    }

    for (size_t i = 0; i < N; i++)
    {
        in_mat(i, i) = res[i];
    }
    return res;
}

// ------------------------------------------

template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(Matrix<T, N, N> in_mat, const T tol)
{
    hessenberg(in_mat, true);
    return _eig_vals_impl(in_mat, tol);
}
template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(const Base_Matrix<T, N, N> &in_mat, const T tol)
{
    Matrix<T, N, N> hessenberg_mat{hessenberg(in_mat)};
    return _eig_vals_impl(hessenberg_mat, tol);
}
template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> eig_vals(const Sparse_Matrix<T, N, N> &in_mat, const T tol)
{
    Matrix<T, N, N> hessenberg_mat{hessenberg(in_mat)};
    return _eig_vals_impl(hessenberg_mat, tol);
}

template <typename T, size_t N>
inline std::pair<std::array<T, N>, std::vector<size_t>> _eig_vals_impl(Matrix<T, N, N> &hessenberg_mat, const T tol)
{
    _eig_hessenberg(hessenberg_mat, static_cast<Matrix<T, N, N> *>(nullptr), tol);
    std::array<T, N> res_array;
    std::vector<size_t> res_vector{};

    for (size_t i = 0; i < N - 1; i++)
    {
        if (hessenberg_mat(i + 1, i) == 0.)
        {
            res_array[i] = hessenberg_mat(i, i);
        }
        else
        {
            T sigma2{hessenberg_mat(i, i) * hessenberg_mat(i + 1, i + 1) - hessenberg_mat(i, i + 1) * hessenberg_mat(i + 1, i)};
            T re_2sigma{hessenberg_mat(i, i) + hessenberg_mat(i + 1, i + 1)};
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
    if (hessenberg_mat(N - 1, N - 2) == 0.)
    {
        res_array[N - 1] = hessenberg_mat(N - 1, N - 1);
    }
    return {res_array, res_vector};
}

// ------------------------------------------

template <typename T, size_t N>
inline std::array<T, N> eig_vals_symm(const Sparse_Matrix<T, N, N> &in_mat,
                                      const Eig_Symm_Method method, const T tol)
{
    auto temp{hessenberg(in_mat, true)};
    return _eig_vals_symm_impl(temp, method, tol);
}
template <typename T, size_t N>
inline std::array<T, N> eig_vals_symm(const Symm_Matrix<T, N> &in_mat,
                                      const Eig_Symm_Method method, const T tol)
{
    auto temp{hessenberg(in_mat)};
    return _eig_vals_symm_impl(temp, method, tol);
}
template <typename T, size_t N, size_t M>
inline std::array<T, N> eig_vals_symm(const Symm_Band_Matrix<T, N, M> &in_mat,
                                      const Eig_Symm_Method method, const T tol)
{
    auto temp{hessenberg(in_mat)};
    return _eig_vals_symm_impl(temp, method, tol);
}
template <typename T, size_t N>
inline std::array<T, N> eig_vals_symm(const Symm_Band_Matrix<T, N, 1> &in_mat,
                                      const Eig_Symm_Method method, const T tol)
{
    return _eig_vals_symm_impl(in_mat, method, tol);
}

template <typename T, size_t N>
inline std::array<T, N> _eig_vals_symm_impl(Symm_Band_Matrix<T, N, 1> &hessenberg_mat,
                                            const Eig_Symm_Method method, const T tol)
{
    switch (method)
    {
    case Eig_Symm_Method::QR:
    {
        _eig_hessenberg_symm_qr(hessenberg_mat, static_cast<Matrix<T, N, N> *>(nullptr), std::max(tol, std::cbrt(std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon())));
        std::array<T, N> res;
        for (size_t i = 0; i < N; i++)
        {
            res[i] = hessenberg_mat(i, i);
        }
        std::sort(res.begin(), res.end());
        return res;
        break;
    }
    case Eig_Symm_Method::bisection:
        return _eig_hessenberg_symm_bisection(hessenberg_mat, tol);
        break;
    case Eig_Symm_Method::divide_conquer:
        return _eig_hessenberg_symm_div_conq(hessenberg_mat, static_cast<Matrix<T, N, N> *>(nullptr), std::max(tol, std::sqrt(std::numeric_limits<T>::epsilon())));
    default:
    {
        throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": wrong method");
    }
    }
}

} // namespace Misc

#endif // MISC_EIGVALS
