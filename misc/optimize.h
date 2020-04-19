#ifndef MISC_OPTIMIZE
#define MISC_OPTIMIZE

#include <functional>
#include <cmath>
#include <array>
#include <map>
#include <iostream>
#include <iterator>
#include <utility>

#include "matrixlib.h"
#include "linear_eq.h"

namespace Misc
{

/* the Quadratic interpolation part
 */
template <typename T>
T minimize1d_quad(std::function<T(T)> func, T left, T mid, T right,
                  const T tol = std::numeric_limits<T>::epsilon())
{
    if (!(left < mid && mid < right))
    {
        T the_min = std::min(std::min(left, right), mid);
        right = std::max(std::max(left, right), mid);
        left = the_min;
        mid = left + (right - left) / 2.;
    }

    T diff{right - left};
    T f_left{func(left)};
    T f_mid{func(mid)};
    T f_right{func(right)};
    T del_left{left - mid};
    T del_right{right - mid};
    const T tol_x{std::sqrt(tol)};
    const T tol_y{tol};

    while (std::abs(right - left) > tol_x * (tol_x + std::abs(left) + std::abs(right)))
    {
        T F_left{(f_left - f_mid) / del_left};
        T F_right((f_right - f_mid) / del_right);
        if (std::abs(F_left) < tol_y * (tol_y + std::abs(f_mid / del_left)) && std::abs(F_right) < tol_y * (tol_y + std::abs(f_mid / del_right)))
        {
            break;
        }
        T a{(F_left - F_right) / (del_left - del_right)};
        diff = (del_left - F_left / a) / 2.;
        if (diff > 0)
        {
            left = mid;
            f_left = f_mid;
        }
        else if (diff < 0)
        {
            right = mid;
            f_right = f_mid;
        }
        else
        {
            left = mid + del_left / 2.;
            right = mid + del_right / 2.;
            f_left = func(left);
            f_right = func(right);
        }
        mid += diff;
        f_mid = func(mid);
        if (mid < left)
        {
            left = mid - 2. * (right - mid);
            f_left = func(left);
        }
        if (mid > right)
        {
            right = mid - 2. * (left - mid);
            f_left = func(right);
        }

        del_left = left - mid;
        del_right = right - mid;
    }
    if (std::abs(right - left) > tol_x * (tol_x + std::abs(left) + std::abs(right)) && !(f_left == f_mid && f_right == f_mid))
    {
        std::cerr << __FILE__ << ':' << __LINE__
                  << ": optimize1d_quad warning: might be an improper solution: left, right="
                  << left << ", " << right << std::endl;
    }
    return mid;
}

/* Kiefer method
 */
template <typename T>
T minimize1d_kiefer(std::function<T(T)> func, T a, T b,
                    const T tol = std::numeric_limits<T>::epsilon())
{
    constexpr T phi{2. / (std::sqrt(5.) + 1.)};
    T del{b - a};
    T x1{a + (1. - phi) * del};
    T func1{func(x1)};
    T x2{a + phi * del};
    T func2{func(x2)};
    T prev_func{func1};
    const T tol_x{std::sqrt(tol)};
    const T tol_y{tol};

    while (std::abs(del) > tol_x * std::max(std::abs(a), std::abs(b)))
    {
        if (func1 > func2)
        {
            a = x1;
            x1 = x2;
            prev_func = func1;
            func1 = func2;
            del = b - a;
            x2 = a + phi * del;
            func2 = func(x2);
        }
        else
        {
            b = x2;
            x2 = x1;
            prev_func = func2;
            func2 = func1;
            del = b - a;
            x1 = a + (1. - phi) * del;
            func1 = func(x1);
        }
        if (std::abs(func1 - func2) < tol_y * std::abs(func1) && std::abs(prev_func - std::max(func1, func2)) < tol_y * std::abs(func1))
        {
            break;
        }
    }
    return func1 > func2 ? x2 : x1;
}

// ===============================================================

// center of mass
template <typename T, size_t N>
std::array<T, N> _nelder_mead_cm(const std::multimap<T, std::array<T, N>> &fx_x);
// return (1. - coef) * ori + coef * to_add
template <typename T, size_t N>
std::array<T, N> _nelder_mead_vec_add(const std::array<T, N> &ori,
                                      const std::array<T, N> &to_add,
                                      T coef);

template <typename T, size_t N>
std::array<T, N> nelder_mead(std::function<T(const std::array<T, N> &)> func,
                             const std::array<T, N> &x0,
                             const std::array<T, N> &other_coef,
                             const size_t max_times = 10000,
                             const T tol = std::numeric_limits<T>::epsilon())
{
    constexpr T alpha{1.};
    constexpr T gamma{2.};
    constexpr T rho{.5};
    constexpr T sigma{.5};

    std::multimap<T, std::array<T, N>> fx_x{{func(x0), x0}}; // always sorted according to fx
    for (size_t i = 0; i < N; i++)
    {
        std::array<T, N> temp{x0};
        temp[i] += other_coef[i];
        T value{func(temp)};
        fx_x.emplace(std::make_pair(value, std::move(temp)));
    }

    T diff;
    T average_func;
    std::array<T, N> co_mass;
    for (size_t i = 0; i < max_times; i++)
    {
        co_mass = _nelder_mead_cm(fx_x);
        std::array<T, N> reflect{_nelder_mead_vec_add(co_mass, fx_x.crbegin()->second, -alpha)};
        T func_reflect{func(reflect)};
        if (fx_x.cbegin()->first <= func_reflect && func_reflect < std::next(fx_x.crbegin())->first)
        {
            fx_x.erase(std::prev(fx_x.end()));
            fx_x.emplace(std::make_pair(func_reflect, std::move(reflect)));
        }
        else if (func_reflect < fx_x.cbegin()->first)
        {
            std::array<T, N> expand{_nelder_mead_vec_add(co_mass, reflect, gamma)};
            T func_expand{func(expand)};
            fx_x.erase(std::prev(fx_x.end()));
            if (func_expand < func_reflect)
            {
                fx_x.emplace(std::make_pair(func_expand, std::move(expand)));
            }
            else
            {
                fx_x.emplace(std::make_pair(func_reflect, std::move(reflect)));
            }
        }
        else
        {
            std::array<T, N> contract{_nelder_mead_vec_add(co_mass, fx_x.crbegin()->second, rho)};
            T func_contract{func(contract)};
            if (func_contract < fx_x.crbegin()->first)
            {
                fx_x.erase(std::prev(fx_x.end()));
                fx_x.emplace(std::make_pair(func_contract, std::move(contract)));
            }
            else
            {
                std::multimap<T, std::array<T, N>> temp_fx_x{fx_x.cbegin(), std::next(fx_x.cbegin())};
                const std::array<T, N> &x0{temp_fx_x.cbegin()->second};
                for (auto it = std::next(fx_x.cbegin()); it != fx_x.cend(); it++)
                {
                    std::array<T, N> new_x{_nelder_mead_vec_add(x0, it->second, sigma)};
                    T value{func(new_x)};
                    temp_fx_x.emplace(std::make_pair(value, std::move(new_x)));
                }
                fx_x = std::move(temp_fx_x);
            }
        }

        average_func = 0.;
        for (const auto &pair : fx_x)
        {
            average_func += pair.first;
        }
        average_func /= fx_x.size();

        diff = 0.;
        for (const auto &pair : fx_x)
        {
            T del = pair.first - average_func;
            diff += del * del;
        }
        diff /= fx_x.size();
        diff = std::sqrt(diff);

        if (diff < tol * std::abs(average_func))
        {
            return co_mass;
        }
    }
    std::cerr << __FILE__ << ':' << __LINE__
              << ": NM simplex might not converge; diff, average_func: "
              << diff << ", " << average_func << std::endl;
    return co_mass;
}

template <size_t N, typename T, typename... X>
std::array<T, N> nelder_mead(std::function<T(X...)> func,
                             const std::array<T, N> &x0,
                             const std::array<T, N> &other_coef,
                             const size_t max_times = 10000,
                             const T tol = std::numeric_limits<T>::epsilon())
{
    std::function<T(const std::array<T, N> &)> new_func{
        [&func](const std::array<T, N> &x) {
            return std::apply(func, x);
        }};
    return nelder_mead(new_func, x0, other_coef, max_times, tol);
}

// center of mass
template <typename T, size_t N>
std::array<T, N> _nelder_mead_cm(const std::multimap<T, std::array<T, N>> &fx_x)
{
    // sum
    auto it{std::next(fx_x.crbegin())};
    std::array<T, N> res{it->second};
    for (it++; it != fx_x.crend(); it++)
    {
        const std::array<T, N> &to_add{it->second};
        for (size_t i = 0; i < N; i++)
        {
            res[i] += to_add[i];
        }
    }
    // get average
    for (size_t i = 0; i < N; i++)
    {
        res[i] /= fx_x.size() - 1;
    }
    return res;
}

template <typename T, size_t N>
std::array<T, N> _nelder_mead_vec_add(const std::array<T, N> &ori,
                                      const std::array<T, N> &to_add,
                                      T coef)
{
    std::array<T, N> res{ori};
    for (size_t i = 0; i < N; i++)
    {
        res[i] += coef * (to_add[i] - ori[i]);
    }
    return res;
}

// =============================================================================

template <typename T, size_t N>
std::array<T, N> grad(std::function<T(const std::array<T, N> &)> func,
                      std::array<T, N> x0,
                      const T step = std::cbrt(6. * std::numeric_limits<T>::epsilon()));

/*beg:cg*/

/* Multi-dimensional optimization
 * using Polak-Ribiere Algorithm
 * modify is used to modify x0's generated
 */
template <typename T, size_t N>
std::array<T, N> conj_grad(std::function<T(const std::array<T, N> &)> func,
                           std::function<std::array<T, N>(const std::array<T, N> &)> grad_f,
                           std::array<T, N> x0,
                           std::function<void(std::array<T, N> &)> modify,
                           const size_t max_times = 1000,
                           const T tol = std::numeric_limits<T>::epsilon())
{
    std::array<T, N> g{grad_f(x0)};
    T gg{g * g};
    constexpr T phi{(std::sqrt(5.) + 1.) / 2.}; // 1.618
    T alpha{1. / std::sqrt(gg)};
    T up_bound{alpha};
    std::array<T, N> d;
    for (size_t i = 0; i < N; i++)
    {
        d[i] = -g[i];
    }

    for (size_t i = 0; i < max_times; i++)
    {
        // find alpha
        std::function<T(T)> func_alpha{
            [&func, &x0, &d](double alpha) {
                std::array<T, N> x{x0};
                for (size_t i = 0; i < N; i++)
                {
                    x[i] += alpha * d[i];
                }
                return func(x);
            }};
        up_bound = alpha;
        do
        {
            up_bound *= phi;
            alpha = minimize1d_kiefer(func_alpha, 0., up_bound);
        } while (alpha >= up_bound / phi);
        // new x0
        for (size_t i = 0; i < N; i++)
        {
            x0[i] += alpha * d[i];
        }
        // stop or not
        modify(x0);

        if (alpha < tol * std::sqrt((x0 * x0) / (d * d)))
        {
            return x0;
        }
        // new g & beta
        std::array<T, N> new_g{grad_f(x0)};
        T new_gg{new_g * new_g};
        T beta{(new_gg - new_g * g) / gg};
        gg = new_gg;
        g = new_g;
        // new d
        for (size_t i = 0; i < N; i++)
        {
            d[i] *= beta;
            d[i] -= g[i];
        }
    }
    std::cerr << __FILE__ << ':' << __LINE__
              << ": CG might not converge; grad, func: "
              << std::sqrt(gg) << ", " << func(x0) << std::endl;
    return x0;
}
/*end:cg*/

template <typename T, size_t N>
std::array<T, N> conj_grad(std::function<T(const std::array<T, N> &)> func,
                           std::array<T, N> x0,
                           const size_t max_times = 1000,
                           const T tol = std::numeric_limits<T>::epsilon())
{
    std::function<void(std::array<T, N> &)> modify{
        [](std::array<T, N> &) {}};
    std::function<std::array<T, N>(const std::array<T, N> &)> grad_f{
        [&func](const std::array<T, N> &x0) {
            return grad(func, x0);
        }};
    return conj_grad(func, grad_f, std::move(x0), modify, max_times, tol);
}

template <typename T, size_t N>
std::array<T, N> conj_grad(std::function<T(const std::array<T, N> &)> func,
                           std::function<std::array<T, N>(const std::array<T, N> &)> grad_f,
                           std::array<T, N> x0,
                           const size_t max_times = 1000,
                           const T tol = std::numeric_limits<T>::epsilon())
{
    std::function<void(std::array<T, N> &)> modify{
        [](std::array<T, N> &) {}};
    return conj_grad(func, grad_f, std::move(x0), modify, max_times, tol);
}

template <typename T, size_t N>
std::array<T, N> conj_grad(std::function<T(const std::array<T, N> &)> func,
                           std::array<T, N> x0,
                           std::function<void(std::array<T, N> &)> modify,
                           const size_t max_times = 1000,
                           const T tol = std::numeric_limits<T>::epsilon())
{
    std::function<std::array<T, N>(const std::array<T, N> &)> grad_f{
        [&func](const std::array<T, N> &x0) {
            return grad(func, x0);
        }};
    return conj_grad(func, grad_f, std::move(x0), modify, max_times, tol);
}

template <typename T, size_t N>
std::array<T, N> grad(std::function<T(const std::array<T, N> &)> func,
                      std::array<T, N> x0,
                      const T step)
{
    std::array<T, N> res;

    for (size_t i = 0; i < N; i++)
    {
        x0[i] -= step;
        T f_left{func(x0)};
        x0[i] += 2. * step;
        T f_right{func(x0)};
        res[i] = (f_right - f_left) / 2. / step;
        x0[i] -= step;
    }
    return res;
}

} // namespace Misc

#endif // MISC_OPTIMIZE
