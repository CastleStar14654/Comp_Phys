#ifndef MISC_IVP
#define MISC_IVP

#include <vector>
#include <iostream>
#include <array>
#include <string>
#include <cmath>
#include <map>
#include <functional>
#include <complex>

#include "matrixlib.h"

namespace Misc
{

    /*beg:ivp_solver*/
    template <size_t N>
    class _IVP_solver
    {
    protected:
        std::function<std::array<double, N>(double, const std::array<double, N> &)> func;
        std::vector<double> ts;
        std::array<std::vector<double>, N> ys{};
        std::array<double, N> y;

    public:
        _IVP_solver(std::function<std::array<double, N>(double, const std::array<double, N> &)> function,
                    double t0, const std::array<double, N> &y0)
            : func{function}, ts{t0}, y{y0}
        {
            for (size_t i = 0; i < N; i++)
            {
                ys[i].emplace_back(y0[i]);
            }
        }
        virtual std::array<double, N> step_inc(double t, double h, const std::array<double, N> &y) const = 0;
        double last_t() const { return ts.back(); }
        const std::array<double, N> &last_y() const { return y; }
        void step(double h, const std::array<double, N> &y_inc)
        {
            ts.push_back(ts.back() + h);
            for (size_t i = 0; i < N; i++)
            {
                y[i] += y_inc[i];
                ys[i].emplace_back(y[i]);
            }
        }
        std::pair<std::vector<double>, std::array<std::vector<double>, N>> finish()
        {
            return std::make_pair(std::move(ts), std::move(ys));
        }
    };
    /*end:ivp_solver*/

    template <size_t N>
    class _Euler_forward : public _IVP_solver<N>
    {
    public:
        _Euler_forward(std::function<std::array<double, N>(double, const std::array<double, N> &)> function,
                       double t0, const std::array<double, N> &y0)
            : _IVP_solver<N>{function, t0, y0} {}
        std::array<double, N> step_inc(double t, double h, const std::array<double, N> &y) const override
        {
            std::array<double, N> res{this->func(t, y)};
            for (auto &y_p : res)
            {
                y_p *= h;
            }
            return res;
        }
    };

    /*beg:radau*/
    template <size_t N>
    class _Radau_5 : public _IVP_solver<N>
    {
    private:
        const Matrix<double, 3, 3> mat_a{
            {0.19681547722366044, -0.06553542585019838, 0.02377097434822015},
            {0.3944243147390873, 0.29207341166522843, -0.04154875212599792},
            {0.37640306270046725, 0.5124858261884216, 0.1111111111111111}};
        constexpr static std::array<double, 3> vec_b{
            0.37640306270046725, 0.5124858261884216, 0.1111111111111111};
        constexpr static std::array<double, 3> vec_c{
            0.15505102572168222, 0.6449489742783178, 1.};
        _Euler_forward<N> euler_solver;
        mutable Matrix<double, 3, N> mat_g_f{};

    public:
        _Radau_5(std::function<std::array<double, N>(double, const std::array<double, N> &)> function,
                 double t0, const std::array<double, N> &y0)
            : _IVP_solver<N>{function, t0, y0}, euler_solver{function, t0, y0} {}
        std::array<double, N> step_inc(double t, double h, const std::array<double, N> &y) const override
        {
            std::array<double, N> temp;
            const std::array<double, 3> t_s{t + vec_c[0] * h, t + vec_c[1] * h, t + vec_c[2] * h};

            temp = euler_solver.step_inc(t, h, y);
            for (size_t i = 0; i < 3; i++)
            {
                for (size_t j = 0; j < N; j++)
                {
                    mat_g_f(i, j) = y[j] + temp[j] * vec_c[i]; // g now
                }
            }
            std::array<double, N> euler_sol{temp};

            for (size_t i = 0; i < 2; i++)
            {
                for (size_t j = 0; j < 3; j++)
                {
                    std::copy(&mat_g_f(j, 0), &mat_g_f(j, N), temp.begin());
                    mat_g_f.row(j) = _IVP_solver<N>::func(t_s[j], temp); // f now
                }
                mat_g_f = mat_a * mat_g_f;
                for (size_t j = 0; j < 3; j++)
                    for (size_t k = 0; k < N; k++)
                    {
                        mat_g_f(j, k) *= h;
                        mat_g_f(j, k) += y[k]; // g now
                    }
            }

            for (size_t j = 0; j < 3; j++)
            {
                std::copy(&mat_g_f(j, 0), &mat_g_f(j, N), temp.begin());
                mat_g_f.row(j) = _IVP_solver<N>::func(t_s[j], temp); // f now
            }

            temp.fill(0);
            for (size_t j = 0; j < 3; j++)
                for (size_t k = 0; k < N; k++)
                {
                    temp[k] += vec_b[j] * mat_g_f(j, k);
                }
            for (auto &y_p : temp)
            {
                y_p *= h;
            }

            return temp;
        }
    };
    /*end:radau*/

    template <size_t N>
    constexpr std::array<double, 3> _Radau_5<N>::vec_b;
    template <size_t N>
    constexpr std::array<double, 3> _Radau_5<N>::vec_c;

    // ===========================================================

    template <size_t N>
    std::pair<std::vector<double>,
              std::array<std::vector<double>, N>>
    ivp_euler(std::function<std::array<double, N>(double, const std::array<double, N> &)> function,
              double t0, double t1, const std::array<double, N> &y0, double h0)
    {
        _Euler_forward<N> solver{function, t0, y0};
        double t{t0};
        while (t < t1)
        {
            auto temp{solver.step_inc(t, h0, solver.last_y())};
            solver.step(h0, temp);
            t = solver.last_t();
        }
        return solver.finish();
    }

    /*beg:radau_func*/
    template <size_t N>
    std::pair<std::vector<double>,
              std::array<std::vector<double>, N>>
    ivp_radau(std::function<std::array<double, N>(double, const std::array<double, N> &)> function,
              double t0, double t1, const std::array<double, N> &y0, double h0, bool auto_h = true,
              double h_min = 0., double h_max = std::numeric_limits<double>::infinity(),
              double rel_tol = 9.5367431640625e-07, double abs_tol = 9.313225746154785e-10)
    {
        _Radau_5<N> solver{function, t0, y0};
        double t{t0};

        std::array<double, N> inc_2h;
        while (t < t1)
        {
            if (auto_h)
            {
                inc_2h = solver.step_inc(t, 2 * h0, solver.last_y());
            }

            auto inc_h1{solver.step_inc(t, h0, solver.last_y())};
            solver.step(h0, inc_h1);
            t = solver.last_t();

            auto inc_h2{solver.step_inc(t, h0, solver.last_y())};
            solver.step(h0, inc_h2);
            t = solver.last_t();

            if (auto_h && h0 > 2. * h_min && h0 * 2 < h_max)
            {
                bool double_h{true};
                for (size_t i = 0; i < N; i++)
                {
                    double delta{std::abs(inc_2h[i] - inc_h1[i] - inc_h2[i])};
                    double criteria_7{7. * (abs_tol + rel_tol * std::abs(solver.last_y()[i]))}; // criteria * 7
                    if (delta > 2 * criteria_7)
                    {
                        h0 /= 2.;
                        break;
                    }
                    else if (delta * 8 > criteria_7)
                    {
                        double_h = false;
                    }
                }
                if (double_h)
                {
                    h0 *= 2;
                }
            }
        }
        return solver.finish();
    }
    /*end:radau_func*/

} // namespace Misc

#endif // MISC_IVP
