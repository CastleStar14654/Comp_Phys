/* piecewise and spline interpolations
 * ===================================
 * Hermite_Spline<T>
 *      generate a Hermite cubic piecewise interpolation
 * Cubic_Spline<T>
 *      generate a cubic piecewise spline interpolation
 */

#ifndef MISC_INTERPOLATE_SPLINE
#define MISC_INTERPOLATE_SPLINE

#include <vector>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <functional>
#include <utility>

#include "../polynomial.h"
#include "../matrixlib.h"
#include "../linear_eq.h"

namespace Misc
{

template <typename T>
class Cubic_Spline;

template <typename T>
class Cubic_Hermite
{
protected:
    std::vector<T> xs;
    std::vector<T> ys;
    std::vector<T> dydxs;
    std::vector<Polynomial<T>> polys;

public:
    friend Cubic_Spline<T>;
    template <typename It1, typename It2, typename It3>
    explicit Cubic_Hermite(It1 bx, It1 ex, It2 by, It2 ey, It3 bdydx, It3 edydx)
        : xs{}, ys{}, dydxs{}, polys{}
    {
        xs.push_back(*bx++);
        ys.push_back(*by++);
        dydxs.push_back(*bdydx++);
        add_points(bx, ex, by, ey, bdydx, edydx);
    }
    template <typename C1, typename C2, typename C3>
    explicit Cubic_Hermite(const C1 &xs, const C2 &ys, const C3 &dydxs)
        : Cubic_Hermite(xs.begin(), xs.end(), ys.begin(), ys.end(), dydxs.begin(), dydxs.end())
    {
    }
    explicit Cubic_Hermite(std::initializer_list<T> xs, std::initializer_list<T> ys, std::initializer_list<T> dydxs)
        : Cubic_Hermite(xs.begin(), xs.end(), ys.begin(), ys.end(), dydxs.begin(), dydxs.end())
    {
    }

    // ================= evaluation =======================
    T operator()(T x) const
    {
        T res{0.};
        auto right{std::lower_bound(xs.begin(), xs.end(), x)};
        if (x == *right)
        {
            return *(ys.begin() + (right - xs.begin()));
        }
        else
        {
            T eta{x - *std::prev(right)};
            T del{*right - *std::prev(right)};
            return (*(polys.begin() + (std::prev(right) - xs.begin())))(eta / del);
        }
    }
    T deriv(T x) const
    {
        T res{0.};
        auto right{std::lower_bound(xs.begin(), xs.end(), x)};
        if (x == *right)
        {
            return *(dydxs.begin() + (right - xs.begin()));
        }
        else
        {
            T eta{x - *std::prev(right)};
            T del{*right - *std::prev(right)};
            return (polys.begin() + (std::prev(right) - xs.begin()))->deriv()(eta / del) / del;
        }
    }
    T deriv2(T x) const
    {
        T res{0.};
        auto right{std::upper_bound(xs.begin(), xs.end(), x)};
        if (right == xs.end())
        {
            right--;
        }
        T eta{x - *std::prev(right)};
        T del{*right - *std::prev(right)};
        return (polys.begin() + (std::prev(right) - xs.begin()))->deriv(2)(eta / del) / del / del;
    }
    // =========== utilities =================
    void add_point(T x, T y, T dydx)
    {
        T del{x - xs.back()};
        polys.emplace_back({2 * (ys.back() - y) + del * (dydxs.back() + dydx),
                            3 * (y - ys.back()) - del * (2 * dydxs.back() + dydx),
                            dydxs.back() * del,
                            ys.back()});
        xs.push_back(x);
        ys.push_back(y);
        dydxs.push_back(dydx);
    }

    template <typename It1, typename It2, typename It3>
    void add_points(It1 bx, It1 ex, It2 by, It2 ey, It3 bdydx, It3 edydx, bool back = true)
    {
        if (!std::is_sorted(bx, ex, std::less_equal<>()))
        {
            throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": xs are not strictly increasing");
        }
        if (back)
        {
            if (*bx <= xs.back())
            {
                throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": xs are not strictly increasing");
            }
            for (; bx != ex && by != ey && bdydx != edydx; bx++, by++, bdydx++)
            {
                add_point(*bx, *by, *bdydx);
            }
        }
        else
        {
            if (*std::prev(ex) >= xs.front())
            {
                throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": xs are not strictly increasing");
            }
            Cubic_Hermite temp{bx, ex, by, ey, bdydx, edydx};
            temp.add_point(xs.front(), ys.front(), dydxs.front());
            std::move(xs.begin() + 1, xs.end(), std::back_inserter(temp.xs));
            std::move(ys.begin() + 1, ys.end(), std::back_inserter(temp.ys));
            std::move(dydxs.begin() + 1, dydxs.end(), std::back_inserter(temp.dydxs));
            std::move(polys.begin(), polys.end(), std::back_inserter(temp.polys));
            std::swap(xs, temp.xs);
            std::swap(ys, temp.ys);
            std::swap(dydxs, temp.dydxs);
            std::swap(polys, temp.polys);
        }
    }
    template <typename C1, typename C2, typename C3>
    void add_points(const C1 &xs, const C2 &ys, const C3 &dydxs, bool back = true)
    {
        add_points(xs.begin(), xs.end(), ys.begin(), ys.end(), dydxs.begin(), dydxs.end(), back);
    }
    void add_points(std::initializer_list<T> xs, std::initializer_list<T> ys, std::initializer_list<T> dydxs, bool back = true)
    {
        add_points(xs.begin(), xs.end(), ys.begin(), ys.end(), dydxs.begin(), dydxs.end(), back);
    }
    // ============= valid boundary ==================
    T left() const
    {
        return xs.front();
    }
    T right() const
    {
        return xs.back();
    }
    // ============= properties ======================
    const std::vector<T> &values() const
    {
        return ys;
    }
    const std::vector<T> &derivatives() const
    {
        return dydxs;
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

// ===================================================================

template <typename T>
class Cubic_Spline
{
public:
    template <size_t N>
    struct Knots
    {
    };
    enum class BC
    {
        not_a_knot = -1,
        periodic = 0,
        clamped = 1,
        natural = 2,
        specified = 3
    };
    enum class BD
    {
        first = 1,
        second = 2
    };

private:
    Cubic_Hermite<T> base;

    using BD_Value = std::pair<std::pair<BD, T>, std::pair<BD, T>>;

    constexpr static BD_Value clamped{{BD::first, T{}}, {BD::first, T{}}};
    constexpr static BD_Value natural{{BD::second, T{}}, {BD::second, T{}}};

public:
    template <typename It1, typename It2, size_t N>
    explicit Cubic_Spline(It1 bx, It1 ex, It2 by, It2 ey,
                          Knots<N>, BC bc_type = BC::natural,
                          BD_Value bc_value = natural)
        : base{bx, bx + 2, by, by + 2, by, by + 2}
    {
        std::vector<T> dydxs{calc_dydxs<N>(bx, ex, by, ey, bc_type, bc_value)};
        base = Cubic_Hermite<T>{bx, ex, by, ey, dydxs.begin(), dydxs.end()};
    }
    template <typename C1, typename C2, size_t N>
    explicit Cubic_Spline(const C1 &xs, const C2 &ys, Knots<N>,
                          BC bc_type = BC::natural,
                          BD_Value bc_value = natural)
        : Cubic_Spline(xs.begin(), xs.end(), ys.begin(), ys.end(),
                       Knots<N>{}, bc_type, bc_value)
    {
    }
    template <size_t N>
    explicit Cubic_Spline(std::initializer_list<T> xs, std::initializer_list<T> ys,
                          Knots<N>, BC bc_type = BC::natural,
                          BD_Value bc_value = natural)
        : Cubic_Spline(xs.begin(), xs.end(), ys.begin(), ys.end(),
                       Knots<N>{}, bc_type, bc_value)
    {
    }

    // ================= evaluation =======================
    T operator()(T x) const
    {
        return base(x);
    }
    T deriv(T x) const
    {
        return base.deriv(x);
    }
    T deriv2(T x) const
    {
        return base.deriv2(x);
    }
    // =========== utilities =================
    // void add_point(T x, T y, T dydx)
    // {
    //     T del{x - base.xs.back()};
    //     polys.push_back(Polynomial<T>{
    //         2 * (base.ys.back() - y) + del * (base.dydxs.back() + dydx),
    //         3 * (y - base.ys.back()) - del * (2 * base.dydxs.back() + dydx),
    //         base.dydxs.back() * del,
    //         base.ys.back()});
    //     base.xs.push_back(x);
    //     base.ys.push_back(y);
    //     base.dydxs.push_back(dydx);
    // }

    // template <typename It1, typename It2, typename It3>
    // void add_points(It1 bx, It1 ex, It2 by, It2 ey, It3 bdydx, It3 edydx, bool back = true)
    // {
    //     if (!std::is_sorted(bx, ex, std::less_equal<>()))
    //     {
    //         throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": xs are not strictly increasing");
    //     }
    //     if (back)
    //     {
    //         if (*bx <= xs.back())
    //         {
    //             throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": xs are not strictly increasing");
    //         }
    //         for (; bx != ex && by != ey && bdydx != edydx; bx++, by++, bdydx++)
    //         {
    //             add_point(*bx, *by, *bdydx);
    //         }
    //     }
    //     else
    //     {
    //         if (*std::prev(ex) >= xs.front())
    //         {
    //             throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": xs are not strictly increasing");
    //         }
    //         Cubic_Spline temp{bx, ex, by, ey, bdydx, edydx};
    //         temp.add_point(xs.front(), ys.front(), dydxs.front());
    //         std::move(xs.begin() + 1, xs.end(), std::back_inserter(temp.xs));
    //         std::move(ys.begin() + 1, ys.end(), std::back_inserter(temp.ys));
    //         std::move(dydxs.begin() + 1, dydxs.end(), std::back_inserter(temp.dydxs));
    //         std::move(polys.begin(), polys.end(), std::back_inserter(temp.polys));
    //         std::swap(xs, temp.xs);
    //         std::swap(ys, temp.ys);
    //         std::swap(dydxs, temp.dydxs);
    //         std::swap(polys, temp.polys);
    //     }
    // }
    // template <typename C1, typename C2, typename C3>
    // void add_points(const C1 &xs, const C2 &ys, const C3 &dydxs, bool back = true)
    // {
    //     add_points(xs.begin(), xs.end(), ys.begin(), ys.end(), dydxs.begin(), dydxs.end(), back);
    // }
    // void add_points(std::initializer_list<T> xs, std::initializer_list<T> ys, std::initializer_list<T> dydxs, bool back = true)
    // {
    //     add_points(xs.begin(), xs.end(), ys.begin(), ys.end(), dydxs.begin(), dydxs.end(), back);
    // }
    // ============= valid boundary ==================
    T left() const
    {
        return base.left();
    }
    T right() const
    {
        return base.back();
    }
    // ============= properties ======================
    const std::vector<T> &values() const
    {
        return base.values();
    }
    const std::vector<T> &derivatives() const
    {
        return base.derivatives();
    }
    const std::vector<T> &partitions() const
    {
        return base.partitions();
    }
    size_t size() const
    {
        return base.size();
    }

private:
    template <size_t N, typename It1, typename It2>
    std::vector<T> calc_dydxs(It1 bx, It1 ex, It2 by, It2 ey,
                              BC bc_type = BC::natural,
                              BD_Value bc_value = natural) const
    {
        static_assert(N >= 5);

        if ((ex - bx != N) || (ey - by != N))
        {
            throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": wrong knots number");
        }
        // prepare the tridiagonal matrix
        std::array<T, N> mu_band_x;   // also work as the dydxs
        std::array<T, N> diag_band_b; // also work as the non-homogeneous vector
        std::array<T, N> lbd_band;
        diag_band_b.fill(2.);
        auto xl{bx};
        auto xm{std::next(bx)};
        auto xr{std::next(xm)};
        for (auto lbd = lbd_band.begin(), mu = std::next(mu_band_x.begin()); xr != ex; lbd++, mu++, xl++, xm++, xr++)
        {
            *lbd = (*xr - *xm) / (*xr - *xl);
            *mu = 1. - *lbd;
        }
        Band_Matrix<T, N, 1> mtx{mu_band_x, diag_band_b, lbd_band};
        // prepare the non-homogeneous value
        xm = std::next(bx);
        xr = std::next(xm);
        auto ity{by};
        T prev{(*std::next(ity) - *ity) / (*xm - *bx)};
        T current;
        ity++;
        for (auto itb = std::next(diag_band_b.begin()), mu = std::next(mu_band_x.begin()), lbd = lbd_band.begin();
             xr != ex; itb++, lbd++, mu++, xm++, xr++, ity++)
        {
            current = (*std::next(ity) - *ity) / (*xr - *xm);
            *itb = 3. * (*lbd * prev + *mu * current);
            prev = current;
        }
        // make more specification of the matrix and b
        switch (bc_type)
        {
        case BC::clamped:
            bc_value = clamped;
            goto specified;
            break;
        case BC::natural:
            bc_value = natural;
            goto specified;
            break;
        case BC::periodic:
        {
            if (*by != *std::prev(ey))
            {
                throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": not periodic");
            }
            mtx(N - 1, N - 2) = 0.;
            mtx(N - 1, N - 1) = 1.;
            mtx(N - 2, N - 1) = 0.;
            T h0{*std::next(bx) - *bx};
            T hn{*std::prev(ex) - *std::prev(std::prev(ex))};
            mtx(0, 1) = hn / (hn + h0);

            diag_band_b.back() = 0.;
            T dydx0{(*std::next(by) - *by) / h0};
            T dydxn{(*std::prev(ey) - *std::prev(std::prev(ey))) / hn};
            diag_band_b.front() = 3. * ((1. - mtx(0, 1)) * dydxn + mtx(0, 1) * dydx0);

            Sparse_Matrix<T, N, N> sp_mtx{mtx};
            sp_mtx(0, N - 2) = 1. - sp_mtx(0, 1);
            sp_mtx(N - 2, 0) = 1. - sp_mtx(N - 2, N - 3);
            if (jacobi(sp_mtx, diag_band_b, mu_band_x, true) != 0)
            {
                throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": not converged");
            }
            mu_band_x.back() = mu_band_x.front();
            return std::vector<T>(mu_band_x.begin(), mu_band_x.end());
            break;
        }
        case BC::not_a_knot:
        {
            mtx(0, 0) = mtx(1, 0);
            mtx(0, 1) = 1.;
            T fl{(*std::next(by) - *by) / (*std::next(bx) - *bx)};
            T fr{(*std::next(std::next(by)) - *std::next(by)) / (*std::next(std::next(bx)) - *std::next(bx))};
            diag_band_b[0] = (2 + mtx(1, 2)) * mtx(1, 0) * fl + mtx(1, 2) * mtx(1, 2) * fr;

            mtx(N - 1, N - 1) = mtx(N - 2, N - 1);
            mtx(N - 1, N - 2) = 1.;
            fl = (*std::prev(std::prev(ey)) - *std::prev(std::prev(std::prev(ey)))) / (*std::prev(std::prev(ex)) - *std::prev(std::prev(std::prev(ex))));
            fr = (*std::prev(ey) - *std::prev(std::prev(ey))) / (*std::prev(ex) - *std::prev(std::prev(ex)));
            diag_band_b[N - 1] = mtx(N - 2, N - 3) * mtx(N - 2, N - 3) * fl + mtx(N - 2, N - 1) * (2 + mtx(N - 2, N - 3)) * fr;
            break;
        }
        default:
        specified:
            if (bc_value.first.first == BD::first)
            {
                mtx(0, 0) = 1.;
                mtx(0, 1) = 0.;
                diag_band_b[0] = bc_value.first.second;
            }
            else
            {
                mtx(0, 1) = 1.;
                T h{*std::next(bx) - *bx};
                diag_band_b[0] = 3 * (*std::next(by) - *by) / h - h / 2 * bc_value.first.second;
            }
            if (bc_value.second.first == BD::first)
            {
                mtx(N - 1, N - 1) = 1.;
                mtx(N - 1, N - 2) = 0.;
                diag_band_b[N - 1] = bc_value.second.second;
            }
            else
            {
                mtx(N - 1, N - 2) = 1.;
                T h{*std::prev(ex) - *std::prev(std::prev(ex))};
                diag_band_b[N - 1] = 3 * (*std::prev(ey) - *std::prev(std::prev(ey))) / h + h / 2 * bc_value.second.second;
            }
            break;
        }
        // solve the matrix
        Up_Band_Matrix<T, N, 1> u;
        Low_Band_Matrix<T, N, 1> l;
        tri_factor(mtx, l, u);
        back_sub(l, diag_band_b, mu_band_x);
        back_sub(u, mu_band_x);
        // return
        return std::vector<T>(mu_band_x.begin(), mu_band_x.end());
    }
};

template <typename T>
constexpr typename Cubic_Spline<T>::BD_Value Cubic_Spline<T>::clamped;
template <typename T>
constexpr typename Cubic_Spline<T>::BD_Value Cubic_Spline<T>::natural;

} // namespace Misc

#endif // MISC_INTERPOLATE_SPLINE
