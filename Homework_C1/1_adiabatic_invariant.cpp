#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <functional>

#include "../misc/ivp.h"
#include "../misc/integral.h"
#include "../misc/interpolate.h"
#include "../misc/rootfind.h"

using namespace std;
using namespace Misc;

/*beg:hamiltonian*/
struct HamiltonianEq
{
    function<double(double)> g_func;
    function<double(double)> x_func;
    HamiltonianEq(function<double(double)> g,
                  function<double(double)> x)
        : g_func{g}, x_func{x} {}
    array<double, 2> operator()(double t, const array<double, 2> &y_p)
    {
        double x{x_func(t)};
        return {
            y_p[1],
            -g_func(t) - y_p[0] * (1. - 1. / sqrt(x * x + y_p[0] * y_p[0]))};
    };
};
/*end:hamiltonian*/

template <size_t N>
void solve_berry(double y0, double x_g_amplitude, int problem_number);

int main(int argc, const char **argv)
{
    /*beg:prob_2*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '2')
    {
        cout << "PROBLEM 2\n"
             << "===========================" << endl;
        auto func_g{
            [](double) -> double { return 0.; }};
        for (auto &&v : {1. / 4., 1. / 16., 1. / 64., 1. / 256.})
        {
            // solve ivp
            auto func_x{
                [v](double t) -> double { return 2. - v * t; }};
            HamiltonianEq func{func_g, func_x};
            auto res{ivp_radau(function(func), 0., 2. / v, {0.1, 0.}, 0.1, false)}; // t, (y, p)
            /*end:solve_ivp*/
            // to calc adiabatic invariant J
            vector<double> p_dots{};
            vector<double> p_y_dots{};
            vector<double> p_y_dot__dots{};
            size_t N{res.first.size()};
            for (size_t i = 0; i < N; i++)
            {
                double p_dot{func(res.first[i], {res.second[0][i], res.second[1][i]})[1]};
                p_dots.emplace_back(p_dot);
                p_y_dots.emplace_back(res.second[1][i] * res.second[1][i]); // p q_dot = p p
                p_y_dot__dots.emplace_back(2. * res.second[1][i] * p_dot);  // (p q_dot)_dot = 2 p p_dot
            }

            Cubic_Hermite<double> spline_p_y_dot{res.first, p_y_dots, p_y_dot__dots};
            /*end:spline_p_y_dot*/
            Cubic_Hermite<double> spline_p{res.first, res.second[1], p_dots};
            // find t of where p == 0
            vector<double> p_t_zeroes{0.};
            for (size_t i = 1; i < N - 1; i++)
            {
                if (signbit(res.second[1][i]) != signbit(res.second[1][i + 1]))
                {
                    p_t_zeroes.push_back(res.first[i]);
                }
            }
            for (auto &zero : p_t_zeroes)
            {
                zero = dekker_brent(function(spline_p), zero, zero + 0.2);
            }
            /*end:find_zero*/
            // get x - J
            vector<double> xs{};
            vector<double> Js{};
            if (p_t_zeroes.size() > 2)
                for (auto t_it = next(p_t_zeroes.begin()); t_it != prev(p_t_zeroes.end()); t_it++)
                {
                    xs.emplace_back(func_x(*t_it));
                    Js.emplace_back(integrate(function(spline_p_y_dot), *prev(t_it), *next(t_it)));
                }
            /*beg:2_output*/
            // output y, p, x, J
            ofstream ofs_y{"output/adiabatic_1.2_y" + to_string(int(1. / v)) + ".bin", ios_base::binary};
            ofstream ofs_p{"output/adiabatic_1.2_p" + to_string(int(1. / v)) + ".bin", ios_base::binary};
            ofstream ofs_x{"output/adiabatic_1.2_x" + to_string(int(1. / v)) + ".bin", ios_base::binary};
            ofstream ofs_J{"output/adiabatic_1.2_J" + to_string(int(1. / v)) + ".bin", ios_base::binary};

            ofs_y.write(static_cast<char *>((void *)(&res.second[0][0])), N * sizeof(double));
            ofs_p.write(static_cast<char *>((void *)(&res.second[1][0])), N * sizeof(double));
            ofs_x.write(static_cast<char *>((void *)(&xs[0])), xs.size() * sizeof(double));
            ofs_J.write(static_cast<char *>((void *)(&Js[0])), Js.size() * sizeof(double));
        }
    }
    /*beg:prob_3*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '3')
    {
        cout << "PROBLEM 3\n"
             << "===========================" << endl;
        auto func_x{
            [](double) -> double { return 0.2; }};
        for (auto &&nu : {1. / 4., 1. / 16., 1. / 64., 1. / 256.})
        {
            // solve ivp
            auto func_g{
                [nu](double t) -> double { return 2. * cos(6.283185307179586 * nu * t); }};
            HamiltonianEq func{func_g, func_x};
            auto res{ivp_radau(function(func), 0., .5 / nu, {-2., 0.}, 0.1, false)}; // t, (y, p)

            // to calc adiabatic invariant J
            vector<double> p_dots{};
            vector<double> p_y_dots{};
            vector<double> p_y_dot__dots{};
            size_t N{res.first.size()};
            for (size_t i = 0; i < N; i++)
            {
                double p_dot{func(res.first[i], {res.second[0][i], res.second[1][i]})[1]};
                p_dots.emplace_back(p_dot);
                p_y_dots.emplace_back(res.second[1][i] * res.second[1][i]); // p q_dot = p p
                p_y_dot__dots.emplace_back(2. * res.second[1][i] * p_dot);  // (p q_dot)_dot = 2 p p_dot
            }

            Cubic_Hermite<double> spline_p_y_dot{res.first, p_y_dots, p_y_dot__dots};
            Cubic_Hermite<double> spline_p{res.first, res.second[1], p_dots};
            // find t of where p == 0
            vector<double> p_t_zeroes{0.};
            for (size_t i = 1; i < N - 1; i++)
            {
                if (signbit(res.second[1][i]) != signbit(res.second[1][i + 1]))
                {
                    p_t_zeroes.push_back(res.first[i]);
                }
            }
            for (auto &zero : p_t_zeroes)
            {
                zero = dekker_brent(function(spline_p), zero, zero + 0.2);
            }
            // get x - J
            vector<double> gs{};
            vector<double> Js{};
            if (p_t_zeroes.size() > 2)
                for (auto t_it = next(p_t_zeroes.begin()); t_it != prev(p_t_zeroes.end()); t_it++)
                {
                    gs.emplace_back(func_g(*t_it));
                    Js.emplace_back(integrate(function(spline_p_y_dot), *prev(t_it), *next(t_it)));
                }
            // output y, p, x, J
            ofstream ofs_y{"output/adiabatic_1.3_y" + to_string(int(1. / nu)) + ".bin", ios_base::binary};
            ofstream ofs_p{"output/adiabatic_1.3_p" + to_string(int(1. / nu)) + ".bin", ios_base::binary};
            ofstream ofs_g{"output/adiabatic_1.3_g" + to_string(int(1. / nu)) + ".bin", ios_base::binary};
            ofstream ofs_J{"output/adiabatic_1.3_J" + to_string(int(1. / nu)) + ".bin", ios_base::binary};

            ofs_y.write(static_cast<char *>((void *)(&res.second[0][0])), N * sizeof(double));
            ofs_p.write(static_cast<char *>((void *)(&res.second[1][0])), N * sizeof(double));
            ofs_g.write(static_cast<char *>((void *)(&gs[0])), gs.size() * sizeof(double));
            ofs_J.write(static_cast<char *>((void *)(&Js[0])), Js.size() * sizeof(double));
        }
    }
    /*beg:prob_4*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '4')
    {
        cout << "PROBLEM 4\n"
             << "===========================" << endl;
        solve_berry<1>(-2.1, 2., 4);
        solve_berry<2>(-2.1, 2., 4);
        solve_berry<3>(-2.1, 2., 4);
        solve_berry<4>(-2.1, 2., 4);
        solve_berry<5>(-2.1, 2., 4);
        solve_berry<6>(-2.1, 2., 4);
    }
    /*beg:prob_5*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '5')
    {
        cout << "PROBLEM 5\n"
             << "===========================" << endl;
        solve_berry<1>(.32, .3, 5);
        solve_berry<2>(.32, .3, 5);
        solve_berry<3>(.32, .3, 5);
        solve_berry<4>(.32, .3, 5);
        solve_berry<5>(.32, .3, 5);
        solve_berry<6>(.32, .3, 5);
    }
    /*end:prob_5*/
}

/*beg:berry*/
template <size_t N>
void solve_berry(double y0, double x_g_amplitude, int problem_number)
{
    /*step1*/
    constexpr size_t unit{256};
    double nu{1. / (unit * N)};
    // solve ivp
    auto func_x{
        [nu, x_g_amplitude](double t) -> double {
            return x_g_amplitude * sin(6.283185307179586 * nu * t);
        }};
    auto func_g{
        [nu, x_g_amplitude](double t) -> double {
            return x_g_amplitude * cos(6.283185307179586 * nu * t);
        }};
    HamiltonianEq func{func_g, func_x};
    auto res{ivp_radau(function(func), 0., 1. / nu, {y0, 0.}, 0.1, false)}; // t, (y, p)
    // spline of p
    vector<double> p_dots{};
    for (size_t j = 0; j < res.first.size(); j++)
    {
        p_dots.emplace_back(func(res.first[j], {res.second[0][j], res.second[1][j]})[1]);
    }
    Cubic_Hermite<double> spline_y{res.first, res.second[0], res.second[1]};
    Cubic_Hermite<double> spline_p{res.first, res.second[1], p_dots};
    /*step2*/
    // get omega
    vector<double> ts{};
    vector<double> omegas{};
    double t0{};
    for (size_t i = 0; i < res.first.size(); i += 10)
    {
        // solve ivp
        ts.push_back(res.first[i]);
        double t{ts.back()};
        double x{func_x(t)};
        double g{func_g(t)};
        HamiltonianEq static_func{[g](double) { return g; }, [x](double) { return x; }};
        auto static_res{
            ivp_radau(function(static_func),
                      0., 20.,
                      {res.second[0][i], res.second[1][i]},
                      0.05, false)};
        size_t static_res_size{static_res.first.size()};
        // find zero candidates
        vector<double> p_t_zeroes{};
        for (size_t j = 1; j < static_res_size - 1; j++)
        {
            if (signbit(static_res.second[1][j]) != signbit(static_res.second[1][j + 1]))
            {
                p_t_zeroes.push_back(static_res.first[j]);
            }
        }
        // spline of p
        vector<double> static_p_dots{};
        for (size_t j = 0; j < static_res_size; j++)
        {
            static_p_dots.emplace_back(
                static_func(static_res.first[j],
                            {static_res.second[0][j], static_res.second[1][j]})[1]);
        }
        Cubic_Hermite<double> static_spline_p{static_res.first, static_res.second[1], static_p_dots};
        // find zero
        for (auto &&zero : p_t_zeroes)
        {
            zero = dekker_brent(function(static_spline_p), zero, zero + 0.05);
        }
        // got t for omega0
        if (i == 0)
        {
            Cubic_Hermite<double> static_spline_y{static_res.first, static_res.second[0], static_res.second[1]};
            t0 = dekker_brent(function(static_spline_y), 0.,
                              p_t_zeroes[0], spline_y(1. / nu));
            if (abs(static_spline_p(t0) - spline_p(1. / nu)) > 0.05)
            {
                t0 = dekker_brent(function(static_spline_y), p_t_zeroes[0],
                                  p_t_zeroes[1], spline_y(1. / nu));
            }
            double t1{dekker_brent(function(static_spline_p), t0 - 0.1,
                                   t0 + 0.1, spline_p(1. / nu))};
            t0 = (t1 + t0) / 2.;
        }
        // get omega
        omegas.push_back(6.283185307179586 / (p_t_zeroes[2] - p_t_zeroes[0]));
    }
    /*step3*/
    // interpolation
    Cubic_Spline<double> spline_omega{ts, omegas, Cubic_Spline<double>::Knots<unit * N + 1>{}};
    /*step4*/
    // calc phase
    double phase{spline_omega(0.) * t0};

    phase -= integrate(function(move(spline_omega)), 0., 1. / nu);
    phase = fmod(phase, 6.283185307179586);
    /*step_output*/
    cout << "t0 " << setprecision(16) << t0
         << "\t1/nu " << int(1. / nu)
         << "\tphase " << phase << endl;

    // output y, p
    ofstream ofs_t{"output/adiabatic_1." + to_string(problem_number) + "_t" + to_string(int(1. / nu)) + ".bin",
                   ios_base::binary};
    ofstream ofs_omega{"output/adiabatic_1." + to_string(problem_number) + "_omega" + to_string(int(1. / nu)) + ".bin",
                       ios_base::binary};
    ofstream ofs_y{"output/adiabatic_1." + to_string(problem_number) + "_y" + to_string(int(1. / nu)) + ".bin",
                   ios_base::binary};
    ofstream ofs_p{"output/adiabatic_1." + to_string(problem_number) + "_p" + to_string(int(1. / nu)) + ".bin",
                   ios_base::binary};

    ofs_t.write(static_cast<char *>((void *)(&ts[0])), ts.size() * sizeof(double));
    ofs_omega.write(static_cast<char *>((void *)(&omegas[0])), omegas.size() * sizeof(double));
    ofs_y.write(static_cast<char *>((void *)(&res.second[0][0])), res.first.size() * sizeof(double));
    ofs_p.write(static_cast<char *>((void *)(&res.second[1][0])), res.first.size() * sizeof(double));
}
/*end:berry*/
