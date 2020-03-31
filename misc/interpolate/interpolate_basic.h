/* Lagrangian interpolations
 * ===================================
 * Neville_Interpolator<T>
 *      DO generate a interpolation POLYNOMIAL on its right;
 *      using Neville method, which might be numerically unstable
 * Newton_Interpolator<T>
 *      generate a Lagrangian interpolation using Newton's method
 *      cannot provide a polynomial format result, but is more stable
 */

#ifndef MISC_INTERPOLATE_BASIC
#define MISC_INTERPOLATE_BASIC

#include <vector>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <functional>

#include "../polynomial.h"

namespace Misc
{

template <typename T>
class Neville_Interpolator
{
    std::vector<Polynomial<T>> polys;
    std::vector<T> xs;

public:
    template <typename It1, typename It2>
    explicit Neville_Interpolator(It1 bx, It1 ex, It2 by, It2 ey)
        : polys{}, xs{}
    {
        for (; bx != ex && by != ey; bx++, by++)
        {
            add_point(*bx, *by);
        }
    }
    template <typename C1, typename C2>
    explicit Neville_Interpolator(const C1 &xs, const C2 &ys)
        : Neville_Interpolator(xs.begin(), xs.end(), ys.begin(), ys.end())
    {
    }
    explicit Neville_Interpolator(std::initializer_list<T> xs, std::initializer_list<T> ys)
        : Neville_Interpolator(xs.begin(), xs.end(), ys.begin(), ys.end())
    {
    }

    // ================= evaluation =======================
    T operator()(T x) const
    {
        return polys.back()(x);
    }
    // =========== utilities =================
    void add_point(T x, T y)
    {
        Polynomial<T> temp1{y};
        Polynomial<T> temp2;

        std::size_t i = size() - 1;
        for (auto it = polys.begin(); it != polys.end(); it++, i--)
        {
            temp2 = Polynomial<T>{1, -xs[i]} / (x - xs[i]) * temp1 - Polynomial<T>{1, -x} / (x - xs[i]) * (*it);
            std::swap(*it, temp1);
            std::swap(temp1, temp2);
        }

        xs.push_back(x);
        polys.push_back(temp1);
    }
    // ============= valid boundary ==================
    T left() const
    {
        return *std::min_element(xs.begin(), xs.end());
    }
    T right() const
    {
        return *std::max_element(xs.begin(), xs.end());
    }
    // ============= properties ======================
    const Polynomial<T> &poly() const
    {
        return polys.back();
    }
    const std::vector<Polynomial<T>> &polynomials() const
    {
        return polys;
    }
    const std::vector<T> &partitions() const
    {
        return xs;
    }
    size_t size() const
    {
        return xs.size();
    }
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Neville_Interpolator<T> &a)
{
    os << "Neville_" << a.poly() << " at [" << a.left() << ", " << a.right() << "]";
    return os;
}

// ====================================================================================

template <typename T>
class Newton_Interpolator
{
    std::vector<T> coefs;
    std::vector<T> xs;

public:
    template <typename It1, typename It2>
    explicit Newton_Interpolator(It1 bx, It1 ex, It2 by, It2 ey)
        : coefs{}, xs{}
    {
        for (; bx != ex && by != ey; bx++, by++)
        {
            add_point(*bx, *by);
        }
    }
    template <typename C1, typename C2>
    explicit Newton_Interpolator(const C1 &xs, const C2 &ys)
        : Newton_Interpolator(xs.begin(), xs.end(), ys.begin(), ys.end())
    {
    }
    explicit Newton_Interpolator(std::initializer_list<T> xs, std::initializer_list<T> ys)
        : Newton_Interpolator(xs.begin(), xs.end(), ys.begin(), ys.end())
    {
    }

    // ================= evaluation =======================
    T operator()(T x) const
    {
        T res{0.};
        for (auto it = coefs.crbegin(), itx = xs.cbegin(); it != coefs.crend(); it++, itx++)
        {
            res *= x - *itx;
            res += *it;
        }
        return res;
    }
    // =========== utilities =================
    void add_point(T x, T y)
    {
        T temp1{y};
        T temp2;

        std::size_t i = size() - 1;
        for (auto it = coefs.begin(); it != coefs.end(); it++, i--)
        {
            temp2 = (temp1 - *it) / (x - xs[i]);
            std::swap(*it, temp1);
            std::swap(temp1, temp2);
        }

        xs.push_back(x);
        coefs.push_back(temp1);
    }
    // ============= valid boundary ==================
    T left() const
    {
        return *std::min_element(xs.begin(), xs.end());
    }
    T right() const
    {
        return *std::max_element(xs.begin(), xs.end());
    }
    // ============= properties ======================
    const std::vector<T> &coefficients() const
    {
        return coefs;
    }
    const std::vector<T> &partitions() const
    {
        return xs;
    }
    size_t size() const
    {
        return xs.size();
    }
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Newton_Interpolator<T> &a)
{
    os << "Newton_Interpolator(";
    for (auto it = a.coefficients().crbegin(); it != a.coefficients().crend(); it++)
    {
        os << *it << ", ";
    }
    os << ") at [" << a.left() << ", " << a.right() << "]";
    return os;
}

} // namespace Misc

#endif // MISC_INTERPOLATE_BASIC
