/* Basic integrations
 * ==================================
 *
 */

#ifndef MISC_INTEGRAL__INTEGRAL
#define MISC_INTEGRAL__INTEGRAL

#include <functional>
#include <map>
#include <iterator>
#include <cmath>
#include <random>
#include <iostream>
#include <tuple>

namespace Misc
{

// 1-D & multidimensional Monte-Carlo Integral
template <typename T>
inline T monte_carlo_integrate(std::function<T(T)> func, T a, T b, size_t rand_times = 1000);
// multi-dimensional Monte-Carlo Integral
template <typename T, typename... X>
inline T monte_carlo_integrate(std::function<T(X...)> func,
                               std::array<T, sizeof...(X)> out_bound_left,
                               std::array<T, sizeof...(X)> out_bound_right,
                               std::function<bool(X...)> in_bound, size_t rand_times = 1048576);

// generalized integration from 0 to b, with singularity at 0
// never directly call this function
template <typename T>
inline T _general_integrate(std::function<T(T)> func, T b);

// use Rumberg method to calculate integral
template <typename T>
inline T integrate(std::function<T(T)> func, T a, T b, size_t max_times = 50,
                   double abs_epsilon = 1e-10, double rel_epsilon = 1e-10)
{
    // optimization when a==b
    int reverse{1};
    if (a > b)
    {
        reverse = -1;
        std::swap(a, b);
    }
    else if (a == b)
    {
        return 0.;
    }
    // inf boundary
    T rem{0.};
    if (std::isinf(a) && std::isinf(b))
    {
        std::cerr << __FILE__ << ':' << __LINE__ << ": all boundary are inf; please use inf_integrate to "
                                                    "specify where to start general integration; using +-32."
                  << std::endl;
        rem += _general_integrate(
            std::function<T(T)>{
                [&func](T x) {
                    return func(-1. / x) / x / x;
                }},
            1. / 32);
        rem += _general_integrate(
            std::function<T(T)>{
                [&func](T x) {
                    return func(1. / x) / x / x;
                }},
            1. / 32);
        a = -32.;
        b = 32.;
    }
    else if (std::isinf(a))
    {
        a = -std::max(10 * std::abs(b), 32.);
        rem += _general_integrate(
            std::function<T(T)>{
                [&func](T x) {
                    return func(-1. / x) / x / x;
                }},
            -1. / a);
    }
    else if (std::isinf(b))
    {
        b = std::max(10 * std::abs(a), 32.);
        rem += _general_integrate(
            std::function<T(T)>{
                [&func](T x) {
                    return func(1. / x) / x / x;
                }},
            1. / b);
    }

    // initiate the x - f(x) map; 4*32 + 1 = 129 x's are selected
    std::map<T, T> xf{};
    T h{(b - a) / 128};
    for (size_t i = 0; i <= 128; i++)
    {
        T x{a + i * h};
        xf[x] = func(x);
    }
    // calc the first time
    T prev{0.};
    for (auto it = xf.begin(); std::next(it) != xf.end();)
    {
        prev += 7 * it++->second + 32 * it++->second + 12 * it++->second + 32 * it++->second + 7 * it->second;
    }
    prev *= 2 * h / 45;
    // iteration for max_times
    T current{0.};
    for (size_t i = 0; i < max_times; i++)
    {
        // double the grid
        h /= 2;
        for (auto it = std::next(xf.begin()); it != xf.end(); it++)
        {
            T prev_x{std::prev(it)->first};
            T curr_x{it->first};
            T x{prev_x + h};
            xf[x] = func(x);
        }
        // calculate the new integral
        for (auto it = xf.begin(); std::next(it) != xf.end();)
        {
            current += 7 * it++->second + 32 * it++->second + 12 * it++->second + 32 * it++->second + 7 * it->second;
        }
        current *= 2 * h / 45;
        // check the relative epsilon
        T delta{(current - prev) / 63};
        if (std::abs(delta) / std::abs(current) < rel_epsilon || std::abs(delta) < abs_epsilon)
        {
            return (current + delta + rem) * reverse;
        }
        // cleaning
        prev = current;
        current = 0.;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": integral not converge";
    return NAN;
}

// integrate from -inf to inf
// integrate ordinarily from a to b
template <typename T>
inline T inf_integrate(std::function<T(T)> func, T a, T b, size_t max_times = 50,
                       double abs_epsilon = 1e-10, double rel_epsilon = 1e-10)
{
    T res{0.};
    res += _general_integrate(
        std::function<T(T)>{
            [&func](T x) {
                return func(-1. / x) / x / x;
            }},
        -1. / a);
    res += _general_integrate(
        std::function<T(T)>{
            [&func](T x) {
                return func(1. / x) / x / x;
            }},
        1. / b);
    res += integrate(func, a, b, max_times, abs_epsilon, rel_epsilon);
    return res;
}

// use Rumberg method to calculate integral
//
template <typename T, typename It>
inline T integrate(std::function<T(T)> func, T a, T b, It b_sing, It e_sing, size_t max_times = 50,
                   double abs_epsilon = 1e-10, double rel_epsilon = 1e-10)
{
    int reverse{1};
    if (a > b)
    {
        reverse = -1;
        std::swap(a, b);
    }
    else if (a == b)
    {
        return 0.;
    }

    std::vector<T> sings(b_sing, e_sing);
    std::sort(sings.begin(), sings.end());

    auto first_sing{std::lower_bound(sings.begin(), sings.end(), a)};
    T res{0.};
    T buffer{0.};
    if (first_sing == sings.end() || *first_sing > b)
    {
        res = integrate(func, a, b);
    }
    else
    {
        buffer = std::min(1e-4, .01 * (*first_sing - a));
        res = *first_sing == a
                  ? 0.
                  : integrate(func, a, *first_sing - buffer) + _general_integrate(
                                                                   std::function<T(T)>{
                                                                       [&func, &first_sing](T x) {
                                                                           return func(*first_sing - x);
                                                                       }},
                                                                   buffer);
    }

    for (auto it = first_sing; it != sings.end() && !std::isnan(res) && *it < b; it++)
    {
        if (std::next(it) == sings.end() || *std::next(it) > b)
        {
            buffer = std::min(1e-4, .01 * (b - *it));
            res += _general_integrate(
                std::function<T(T)>{
                    [&func, &it](T x) {
                        return func(*it + x);
                    }},
                buffer);
            res += integrate(func, *it + buffer, b, max_times, abs_epsilon, rel_epsilon);
            break;
        }
        else
        {
            buffer = std::min(1e-4, .01 * (*std::next(it) - *it));
            res += _general_integrate(
                std::function<T(T)>{
                    [&func, &it](T x) {
                        return func(*it + x);
                    }},
                buffer);
            res += integrate(func, *it + buffer, *std::next(it) - buffer, max_times, abs_epsilon, rel_epsilon);
            res += _general_integrate(
                std::function<T(T)>{
                    [&func, &it](T x) {
                        return func(*std::next(it) - x);
                    }},
                buffer);
            break;
        }
    }
    return res * reverse;
}

// generalized integration from 0 to b, with singularity at 0
// never directly call this function
template <typename T>
inline T _general_integrate(std::function<T(T)> func, T b)
{
    // find the order of the singularity
    T p{b < 1e-4 ? std::log2(func(b / 2) / func(b)) : std::log2(func(5e-5) / func(1e-4))};
    if (p >= 1.)
    {
        return NAN;
    }
    T res;
    if (std::isnan(p))
    {
        res = p;
    }
    else
    {
        // the polynomial part: g(x) = f(x)*x^p = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4
        std::function<T(T)> g{[&func, &p](T x) { return func(x) * std::pow(x, p); }};
        std::array<T, 5> fs{g(1e-4 * b), g(.25 * b), g(.5 * b), g(.75 * b), g(b)};
        std::array<T, 5> as;
        as = Matrix<T, 5, 5>{
                 {1, 0, 0, 0, 0},
                 {-25 / 3, 16, -12, 16 / 3, -1},
                 {70 / 3, -208 / 3, 76, -112 / 3, 22 / 3},
                 {-80 / 3, 96, -128, 224 / 3, -16},
                 {32 / 3, -128 / 3, 64, -128 / 3, 32 / 3}} *
             fs;
        // std::cout << p << ' ' << as[0] << ' ' << as[1] << ' ' << as[2] << ' ' << as[3] << ' ' << as[4] << std::endl;

        res = std::pow(b, 1. - p) * (as[0] / (1. - p) + b * (as[1] / (2. - p) + b * (as[2] / (3. - p) + b * (as[3] / (4. - p) + b * as[4] / (5. - p)))));
        // the remain part
        std::function<T(T)> rem{[&func, &as, &p](T x) { return x < 1e-10 ? 0 : func(x) - (as[0] + x * (as[1] + x * (as[2] + x * (as[3] + x * as[4])))) / std::pow(x, p); }};
        // res += b / 36 * (5*rem((1-std::sqrt(.6))/4* b) + 8*rem(.25* b) + 5*rem((1+std::sqrt(.6))/4* b) + 5*rem((3-std::sqrt(.6))/4* b) + 8*rem(.75* b) + 5*rem((3+std::sqrt(.6))/4* b));
        res += integrate(rem, 0., b, 50, 1e-8, 1e-8);
        // if is nan, use Monte Carlo
    }
    if (std::isnan(res))
    {
        std::cerr << __FILE__ << ':' << __LINE__ << ": get nan in _general_integrate, trying Monte-Carlo" << std::endl;
        return monte_carlo_integrate(func, 0., b);
    }
    return res;
}

// 1-D & multidimensional Monte-Carlo Integral
template <typename T>
inline T monte_carlo_integrate(std::function<T(T)> func, T a, T b, size_t rand_times)
{
    static std::default_random_engine ran{};
    T res{0.};
    for (size_t i = 0; i < rand_times; i++)
    {
        res += func(std::uniform_real_distribution<>{a, b}(ran)) / rand_times;
    }
    return res * (b - a);
}

// multi-dimensional Monte-Carlo Integral
template <typename T, typename... X>
inline T monte_carlo_integrate(std::function<T(X...)> func,
                               std::array<T, sizeof...(X)> out_bound_left,
                               std::array<T, sizeof...(X)> out_bound_right,
                               std::function<bool(X...)> in_bound, size_t rand_times)
{
    static std::default_random_engine ran{};
    T res{0.};
    T vol{1.};
    for (size_t i = 0; i < sizeof...(X); i++)
    {
        vol *= out_bound_right[i] - out_bound_left[i];
    }
    if (vol == 0.)
    {
        return 0.;
    }
    size_t vol_count{0};
    std::array<T, sizeof...(X)> xs;
    for (size_t i = 0; i < rand_times; i++)
    {
        for (size_t i = 0; i < sizeof...(X); i++)
        {
            xs[i] = std::uniform_real_distribution<T>{out_bound_left[i], out_bound_right[i]}(ran);
        }
        if (std::apply(in_bound, xs))
        {
            res += std::apply(func, xs) / rand_times;
            vol_count++;
        }
    }
    return res * vol * vol_count / rand_times;
}

} // namespace Misc

#endif // MISC_INTEGRAL__INTEGRAL
