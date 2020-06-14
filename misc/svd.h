#ifndef MISC_SVD
#define MISC_SVD

#include <array>
#include <algorithm>
#include <utility>
#include <vector>
#include <random>
#include <thread>
#include <functional>
#include <stdexcept>

#include "eigvals.h"

namespace Misc
{
    /* return singular values of in_mat.
     * we would have `in_mat == u_mat * sv_diag * v_mat^T`.
     * assume R > C
     * using QR method.
     */
    template <typename T, size_t R, size_t C>
    inline std::array<T, C> svd(const Matrix<T, R, C> &in_mat,
                                Matrix<T, R, C> *p_u_mat = nullptr,
                                Matrix<T, C, C> *p_v_mat = nullptr,
                                const T tol = std::numeric_limits<T>::epsilon())
    {
        std::array<T, C> res;

        // get singular value
        // get AT A
        Symm_Matrix<T, C> ata{};
        for (std::size_t k = 0; k < R; k++)
            for (std::size_t i = 0; i < C; i++)
            {
                T r = in_mat(k, i);
                for (std::size_t j = 0; j <= i; j++)
                {
                    ata(i, j) += r * in_mat(k, j);
                }
            }
        Symm_Band_Matrix<T, C, 1> hessenberg_mat{_hessenberg_impl(ata, p_v_mat)};
        _eig_hessenberg_symm_qr_(hessenberg_mat, p_v_mat,
                                 std::max(tol, std::cbrt(std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon())));
        // sort
        std::array<size_t, C> order;
        for (size_t i = 0; i < C; i++)
            order[i] = i;
        std::stable_sort(order.begin(), order.end(),
                         [&hessenberg_mat](const size_t &a, const size_t &b) -> bool {
                             return hessenberg_mat(a, a) > hessenberg_mat(b, b);
                         });
        std::array<T, C> temp;
        if (p_v_mat)
            for (size_t i = 0; i < C; i++)
            {
                for (size_t j = 0; j < C; j++)
                    temp[j] = (*p_v_mat)(i, order[j]);
                std::copy(temp.begin(), temp.end(), &((*p_v_mat)(i, 0)));
            }
        // singular value
        for (size_t i = 0; i < C; i++)
        {
            res[i] = std::sqrt(hessenberg_mat(order[i], order[i]));
        }
        // get u_mat
        if (p_u_mat)
        {
            if (!p_v_mat)
                throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": cannot get u without v");
            // get `in_mat * v_mat`
            auto ud {in_mat * (*p_v_mat)};
            for (size_t j = 0; j < C; j++)
            {
                T coef{1. / res[j]};
                for (size_t i = 0; i < R; i++)
                {
                    (*p_u_mat)(i, j) = ud(i, j) * coef;
                }
            }
        }
        return res;
    }

} // namespace Misc

#endif // MISC_SVD
