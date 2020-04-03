/* 1D ROOT FINDING
 * =========================
 * graphical methods:
 *      bisection
 *      Dekker-Brent method
 * iterative methods:
 *      Steffenson acceleration
 * polynomial:
 *      see polynomial.h
 */

#ifndef MISC_ROOTFIND
#define MISC_ROOTFIND

#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace Misc
{

/* bisection
 */
template <typename T>
inline T bisection(std::function<T(T)> func, T a, T b, T val = T{}, T tol = std::numeric_limits<T>::epsilon())
{
    if (val != T{})
    {
        func = std::function<T(T)>{
            [func, &val](T x) {
                return func(x) - val;
            }};
    }

    T func_a(func(a));
    if (std::signbit(func_a) == std::signbit(func(b)))
    {
        std::cerr << __FILE__ << ':' << __LINE__ << ": bisection received: a=" << a
                  << ", func(a)=" << func_a << "; b=" << b << ", func(b)=" << func(b) << std::endl;
        return std::numeric_limits<T>::quiet_NaN();
    }
    T c;
    T func_c;
    while (std::abs(b - a) / 2. > tol * std::abs(b) && func_a)
    {
        // std::cerr << "a b: " << a << ' ' << b << ' ' << std::setprecision(16) << (b - a) / 2 << std::endl;
        c = a + (b - a) / 2.;
        func_c = func(c);
        if (func_c == 0.)
        {
            return c;
        }
        else if (std::signbit(func_a) != std::signbit(func_c))
        {
            b = c;
        }
        else
        {
            a = c;
            func_a = func_c;
        }
    }
    if (!(std::abs(func_a) < 1e-10))
    {
        std::cerr << __FILE__ << ':' << __LINE__
                  << ": bisection warning: might be a improper solution: func(a)=" << func_a << std::endl;
    }
    return a;
}

/* Dekker-Brent method
 */
template <typename T>
inline T dekker_brent(std::function<T(T)> func, T a, T b, T val = 0., T tol = std::numeric_limits<T>::epsilon())
{
    if (val != T{})
    {
        func = std::function<T(T)>{
            [func, &val](T x) {
                return func(x) - val;
            }};
    }

    T func_a{func(a)}; // the opposite point of b
    T func_b{func(b)}; // current best estimation
    if (std::signbit(func_a) == std::signbit(func_b))
    {
        std::cerr << __FILE__ << ':' << __LINE__ << ": bisection received: a=" << a
                  << ", func(a)=" << func_a << "; b=" << b << ", func(b)=" << func_b << std::endl;
        return std::numeric_limits<T>::quiet_NaN();
    }
    T c{b}; // previous solution
    T func_c{func_b};
    T s; // next solution
    T func_s;
    T inc{a - b}; // s - b
    T temp_inc;   // temp inc of interplotation

    T criteria{std::abs(b - a)}; // the difference of s and b should be larger than it
    // const T bisection_criteria{std::abs(b-a) * std::numeric_limits<T>::epsilon() * std::pow(2., 10.)}; // when

    while (std::abs(b - a) / 2 > tol * std::abs(b) && func_s)
    {
        // std::cerr << "s a b c |func_b|>|func_a|: " << s << ' ' << a << ' ' << b << ' ' << c << ' ' << (std::abs(func_b) > std::abs(func_a)) << ' ' << std::setprecision(16) << (b - a) / 2 << std::endl;
        if (std::abs(func_b) > std::abs(func_c))
        // condition 1: interval is small enough
        // condition 2: the result of quadric interpolation must fall out of [a, b]
        // so, use bisection
        {
            // make func_b the smaller
            if (std::abs(func_b) > std::abs(func_a))
            {
                std::swap(a, b);
                std::swap(func_a, func_b);
            }
            inc = (a - b) / 2.;
            criteria = std::abs(inc);
        }
        else
        {
            if (a == c || b == c)
            // only two available points for interpolation
            {
                // make func_b the smaller
                if (std::abs(func_b) > std::abs(func_a))
                {
                    std::swap(a, b);
                    std::swap(func_a, func_b);
                }
                temp_inc = -func_b * (b - a) / (func_b - func_a);
            }
            else
            {
                temp_inc = func_b / (func_a - func_c);
                temp_inc *= func_c * (a - b) / (func_a - func_b) + func_a * (c - b) / (func_b - func_c);
            }
            if (std::min(criteria / 2, .75 * std::abs(b - a)) < std::abs(temp_inc) || std::abs(temp_inc / b) < std::numeric_limits<T>::epsilon())
            // this means that sectant method does not converge rapidly
            // or the refinement is too small
            {
                inc = (a - b) / 2.;
                criteria = std::abs(inc);
            }
            else
            {
                criteria = std::abs(inc);
                inc = temp_inc;
            }
        }

        s = b + inc;
        func_s = func(s);
        c = b;
        func_c = func_b;
        if (std::signbit(func_a) != std::signbit(func_s))
        {
            b = s;
            func_b = func_s;
        }
        else
        {
            a = s;
            func_a = func_s;
        }
    }
    if (!(std::abs(func_s) < 1e-10))
    {
        std::cerr << __FILE__ << ':' << __LINE__
                  << ": dekker_brent warning: might be a improper solution: func(s)=" << func_s << std::endl;
    }
    return s;
}

/* Steffenson iterative
 * when original is true, solve f(x) = 0,
 * otherwise, x = f(x)
*/
template <typename T>
inline T steffenson(std::function<T(T)> func, T x0, bool original = false, size_t max_times = 1000, T tol = std::numeric_limits<T>::epsilon())
{
    if (original == false)
    {
        func = std::function<T(T)>{
            [func](T x) {
                return func(x) - x;
            }};
    }

    T diff;
    for (size_t i = 0; i < max_times; i++)
    {
        T func_x{func(x0)};
        T func_func_x{func(x0 + func(x0))};
        if (func_x){
            diff = func_x * func_x / (func_x - func_func_x);
        }
        else
        {
            return x0;
        }
        if (std::abs(diff) < tol * std::abs(x0))
        {
            return x0;
        }
        x0 += diff;
    }
    std::cerr << __FILE__ << ':' << __LINE__ << ": Steffenson method not converged: diff, x0=" << diff << ", " << x0 << std::endl;
    return std::numeric_limits<T>::quiet_NaN();
}

} // namespace Misc

#endif // MISC_ROOTFIND
