#ifndef MISC_INTEGRAL
#define MISC_INTEGRAL
/* Basic integrations
 * ==================================
 *
 */

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
inline T _improper_integrate(std::function<T(T)> func, T b);

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
    T delta;
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
        delta = (current - prev) / 63;
        if (std::abs(delta / current) < rel_epsilon || std::abs(delta) < abs_epsilon)
        {
            return (current + delta) * reverse;
        }
        // cleaning
        prev = current;
        current = 0.;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": integral not converge\n";
    std::cerr << '\t' << "rel_eps=" << std::abs(delta / current) << " abs_eps=" << std::abs(delta) << std::endl;
    return NAN;
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
        res = integrate_DE(func, a, b, max_times, abs_epsilon, rel_epsilon);
    }
    else
    {
        res = integrate_DE(func, a, *first_sing, max_times, abs_epsilon, rel_epsilon);
    }

    for (auto it = first_sing; it != sings.end() && !std::isnan(res) && *it < b; it++)
    {
        if (std::next(it) == sings.end() || *std::next(it) > b)
        {
            res += integrate_DE(func, *it, b, max_times, abs_epsilon, rel_epsilon);
            break;
        }
        else
        {
            res += integrate_DE(func, *it, *std::next(it), max_times, abs_epsilon, rel_epsilon);
        }
    }
    return res * reverse;
}

/* use double exponential method to calculate integral
 * I = \int f(x(t)) dx/dt dt, where
 *      x = (b+a)/2 + (b-a)/2 * tanh(sinh(t))
 *      dx/dt = (b-a)/2 * sech^2(sinh(t))*cosh(t)
 * using q = exp(-2*sinh(t)), then
 *      x = b - (b-a)*q/(1+q) := b - rel_x
 *      dx/dt = 2*(b-a) * q/(1+q)^2 * cosh(t) = 2*rel_x/(1+q) * cosh(t)
 * integrate from t=-h to h
 * also,
 *      x(t) + x(-t) = a + b
 *      dx/dt(at -t) = dx/dt(at t)
 */
template <typename T>
inline T integrate_DE(std::function<T(T)> func, T a, T b, T h = T{}, size_t max_times = 50,
                      double abs_epsilon = 1e-10, double rel_epsilon = 1e-10)
{
    // optimization when a==b
    if (a == b)
    {
        return 0.;
    }
    // if infty boundary
    if ((std::isinf(a) || std::isinf(b)) && h == T{})
    {
        h = 4.;
    }
    if (std::isinf(a) && std::isinf(b))
    {
        return integrate(
            std::function<T(T)>{[&func](T t) {
                return func(std::sinh(std::sinh(t))) * std::cosh(std::sinh(t)) * std::cosh(t);
            }},
            -h, h, max_times, abs_epsilon, rel_epsilon);
    }
    else if (std::isinf(a))
    {
        return integrate(
            std::function<T(T)>{[&func, &b](T t) {
                return func(b - std::exp(2 * std::sinh(t))) * 2 * std::exp(2 * std::sinh(t)) * std::cosh(t);
            }},
            -h, h, max_times, abs_epsilon, rel_epsilon);
    }
    else if (std::isinf(b))
    {
        return integrate(
            std::function<T(T)>{[&func, &a](T t) {
                return func(a + std::exp(2 * std::sinh(t))) * 2 * std::exp(2 * std::sinh(t)) * std::cosh(t);
            }},
            -h, h, max_times, abs_epsilon, rel_epsilon);
    }

    // find the minimum h which would not overflow
    if (h == T{})
    {
        h = std::log(std::numeric_limits<T>::epsilon());
        h -= std::log(std::abs(b - a));
        h += std::log(std::max(std::abs(a), std::abs(b)));
        h = std::asinh(-.5 * h);
    }

    // initiate the res
    T t{h};
    T q{std::exp(-2 * std::sinh(t))};
    T rel_x{(b - a) * q / (1 + q)};
    T prev{h * (func((a + b) / 2) * (b - a) / 2 + (func(b - rel_x) + func(a + rel_x)) * rel_x / (1 + q) * std::cosh(t))};

    // iteration for max_times
    T current{0.};
    T delta;
    T step{h};
    for (size_t i = 0; i < max_times; i++)
    {
        t = step / 2.;
        for (size_t i = 0; t < h; i++, t = step / 2. + i * step)
        {
            q = std::exp(-2 * std::sinh(t));
            rel_x = (b - a) * q / (1 + q);
            current += (func(b - rel_x) + func(a + rel_x)) * 2. * rel_x / (1 + q) * std::cosh(t);
        }
        current *= step;
        current += prev;
        current /= 2.;

        delta = current - prev;
        if (std::abs(delta / current) < rel_epsilon || std::abs(delta) < abs_epsilon)
        {
            return current;
        }
        prev = current;
        current = 0.;
        step /= 2.;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": integral not converge\n";
    std::cerr << '\t' << "rel_eps=" << std::abs(delta / current) << " abs_eps=" << std::abs(delta) << std::endl;
    return NAN;
}

// generalized integration from 0 to b, with singularity at 0
// never directly call this function
template <typename T>
inline T _improper_integrate(std::function<T(T)> func, T b)
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
        std::cerr << __FILE__ << ':' << __LINE__ << ": get nan in _improper_integrate, trying Monte-Carlo" << std::endl;
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

#endif // MISC_INTEGRAL
