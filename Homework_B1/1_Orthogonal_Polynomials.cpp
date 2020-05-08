#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <fstream>

#include "../misc/eigvals.h"

using namespace std;
using namespace Misc;

/*beg:Pn_dec*/

/* normalized Legendre polynomials.
 * recursive relation:
 *      (n+1) P_{n+1}/√(2(n+1)+1)
 *          = x*(2n+1) P_n/√(2n+1) - n P_{n-1}/√(2(n-1)-1)
 *      P_0 / √(2*0+1) = 1/√2
 *      P_1 / √(2*1+1) = x/√2
 * RETURN:
 *      P_n(x)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_legendre_norm(double x, size_t n);
/*end:Pn_dec*/
/*beg:Ln_dec*/

/* Laguerre polynomials.
 * recursive relation:
 *      (n+1) L_{n+1} = (2n+1-x) L_n - n L_{n-1}
 *      L_0 = 1
 *      L_1 = -x + 1
 * RETURN:
 *      L_n(x)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_laguerre(double x, size_t n);
/*end:Ln_dec*/
/*beg:hermite_dec*/

/* recurrence matrix of normalized Hermite polynomials
 * ORDINARY Hermite:
 *      2x H_n = H_{n+1} + 2n H_{n-1}
 *      H_0 = 1, H_1 = 2x
 * NORMALIZED:
 *      H'_n = H_n / √(2^n n! √π)
 *      √2 x H'_n = √(n+1) H'_{n+1} + √(n) H'_{n-1}
 * so, the matrix:
 *      0   √(1/2)
 *   √(1/2)  0  √(2/2)
 *         √(2/2)  0  √(3/2)
 *             √(3/2)  ....
 */
template <size_t N>
Symm_Band_Matrix<double, N, 1> hermite_mat();
/*end:hermite_dec*/
/*beg:chebyshev_dec*/

/* recurrence matrix of normalized Chebyshev polynomials
 * ORDINARY Chebyshev:
 *      x T_n = (T_{n+1} + T_{n-1}) / 2
 *      T_0 = 1, T_1 = x
 * NORMALIZED:
 *      T'_0 = T_0 / √π, T'_i = T_i / √(π/2)
 *      x T'_0 = T'_1 / √2
 *      x T'_1 = T'_0 / √2 + T'_2 / 2
 *      x T'_n = (T'_{n+1} + T'_{n-1}) / 2
 * so, the matrix:
 *      0   1/√2
 *    1/√2  0   .5
 *          .5  0   .5
 *              .5  ....
 */
template <size_t N>
Symm_Band_Matrix<double, N, 1> chebyshev_mat();
/*end:chebyshev_dec*/

/* 5 point derivative
 */
double poly_deriv(function<double(double, size_t)> func, double x, size_t n);

/*beg:hermite_dec*/

/* normalized Hermite polynomials.
 * recursive relation:
 *      √(n+1) H_{n+1} = √2 x H_n - √n H_{n-1}
 *      H_0 = 1 / √√π
 *      H_1 = x / √√π
 * RETURN:
 *      H_n(x)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_hermite(double x, size_t n);
/*end:hermite_dec*/

/*beg:hermite_weighted_dec*/

/* normalized Hermite polynomials with weight
 * RETURN:
 *      H_n(x) * exp(-x^2 / 2)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_hermite_weighted(double x, size_t n);
/*end:hermite_weighted_dec*/

/* normalized Hermite polynomials derivatives
 * using 5 points deri.
 * recursive relation:
 *      √(n+1) H_{n+1} = √2 x H_n - √n H_{n-1}
 *      H_0 = 1 / √√π
 *      H_1 = x / √√π
 * RETURN:
 *      H'_n(x)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_hermite_deri(double x, size_t n);
double my_hermite_deri_weighted(double x, size_t n);

/* normalized Chebyshev polynomials.
 * recursive relation:
 *      T_0 = T_0 / √π, T_i = T_i / √(π/2)
 *      x T_0 = T_1 / √2
 *      x T_1 = T_0 / √2 + T_2 / 2
 *      x T_n = (T_{n+1} + T_{n-1}) / 2
 * RETURN:
 *      T_n(x)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_chebyshev(double x, size_t n);

/* normalized Chebyshev polynomials derivatives
 * using 5 points deri.
 * recursive relation:
 *      T_0 = _T_0 / √π, T_i = _T_i / √(π/2)
 *      x T_0 = T_1 / √2
 *      x T_1 = T_0 / √2 + T_2 / 2
 *      x T_n = (T_{n+1} + T_{n-1}) / 2
 * RETURN:
 *      T'_n(x)
 * PARAMETERS:
 *      double x, size_t n
 */
double my_chebyshev_deri(double x, size_t n);

/*beg:legendre_mat_dec*/

/* recurrence matrix of normalized Legendre polynomials
 * recursive relation:
 *      (n+1) P_{n+1}/√(2(n+1)+1)
 *          = x*(2n+1) P_n/√(2n+1) - n P_{n-1}/√(2(n-1)-1)
 *      P_0 / √(2*0+1) = 1/√2
 *      P_1 / √(2*1+1) = x/√2
 * so α_n = (n+1) / √((2n+1)(2n+3))
 * so, the matrix:
 *      0  1/√3
 *    1/√3   0   2/√15
 *         2/√15   0    3/√35
 *               3/√35  ....
 */
template <size_t N>
Symm_Band_Matrix<double, N, 1> legendre_mat();
/*end:legendre_mat_dec*/

// ====================================================================================

int main(int argc, const char *argv[])
{
    cout << setprecision(16);
    /*beg:prob2.1_P*/
    // 2.1
    cout << "==========================" << '\n'
         << "PROBLEM 2.1" << '\n'
         << "--------------------------" << '\n'
         << "Legendre polynomial P_n(x) at x=0.5" << endl;
    for (size_t n : {2, 16, 128, 1024})
    {
        cout << "n=" << n << ",\tP_n(x)=" << my_legendre_norm(0.5, n) << endl;
    }
    /*end:prob2.1_P*/
    cout << "--------------------------" << '\n'
         << "Laguerre polynomial L_n(x) at x=0.5" << endl;
    for (size_t n : {2, 16, 128, 1024})
    {
        cout << "n=" << n << ",\tL_n(x)=" << my_laguerre(0.5, n) << endl;
    }
    /*end:prob2.1*/
    // 2.2
    cout << "==========================" << '\n'
         << "PROBLEM 2.2" << '\n'
         << "--------------------------" << '\n'
         << "j'th zero of Hermite polynomial of degree 1024 H_n(x)" << '\n'
         << "j=\t2\t\t\t16\t\t\t128\t\t\t1024" << endl;
    auto mat{hermite_mat<1024>()};
    auto res_hermite{eig_vals_symm(mat, Eig_Symm_Method::QR)};
    cout << "QR";
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << res_hermite[i - 1];
    }
    cout << endl;
    res_hermite = eig_vals_symm(mat, Eig_Symm_Method::bisection);
    cout << "bisec";
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << res_hermite[i - 1];
    }
    cout << endl;
    /*end:prob2.2_hermite*/
    cout << "--------------------------" << '\n'
         << "j'th zero of Chebyshev polynomial of degree 1024 T_n(x)" << '\n'
         << "j=\t2\t\t\t16\t\t\t128\t\t\t1024" << endl;
    mat = chebyshev_mat<1024>();
    auto res_cheby{eig_vals_symm(mat, Eig_Symm_Method::QR)};
    cout << "QR";
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << res_cheby[i - 1];
    }
    cout << endl;
    res_cheby = eig_vals_symm(mat, Eig_Symm_Method::bisection);
    cout << "bisec";
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << res_cheby[i - 1];
    }
    cout << endl;
    /*end:prob2.2*/
    if (argc > 1 && argv[1][0] == 't')
    {
        constexpr size_t qr_times{50};
        constexpr size_t bisec_times{5};
        cout << "--------------------------" << '\n'
             << "TIMING Chebyshev * " << qr_times
             << " for QR and * " << bisec_times << " for bisection" << endl;
        auto t1{chrono::steady_clock::now()};
        for (size_t i = 0; i < qr_times; i++)
        {
            eig_vals_symm(mat, Eig_Symm_Method::QR);
        }
        auto t2{chrono::steady_clock::now()};
        cout << "QR (mus per run)\t"
             << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() / double(qr_times)
             << endl;
        t1 = chrono::steady_clock::now();
        for (size_t i = 0; i < bisec_times; i++)
        {
            eig_vals_symm(mat, Eig_Symm_Method::bisection);
        }
        t2 = chrono::steady_clock::now();
        cout << "bisec (mus per run)\t"
             << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() / double(bisec_times)
             << endl;
    }
    /*end:timing*/
    cout << "==========================" << '\n'
         << "PROBLEM 2.3" << '\n'
         << "--------------------------" << '\n'
         << "j'th Gaussian integration weight of Hermite polynomial of degree 1024 H_n(x)" << '\n'
         << "j=\t2\t\t\t16\t\t\t128\t\t\t1024" << endl;
    auto func_hermite{
        [&](double xj) {
            constexpr size_t n{1024};
            double res{my_hermite(xj, n - 1) * my_hermite_deri(xj, n)};
            res *= std::sqrt(n / 2.);
            return 1. / res;
        }};
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << func_hermite(res_hermite[i - 1]);
    }
    /*end:prob2.3_hermite*/
    cout << "\n--------------------------" << '\n'
         << "j'th Gaussian integration weight times rho(x) of Hermite polynomial of degree 1024 H_n(x)" << '\n'
         << "j=\t2\t\t\t16\t\t\t128\t\t\t1024" << endl;
    auto func_hermite_weighted{
        [&](double xj) {
            constexpr size_t n{1024};
            double res{my_hermite_weighted(xj, n - 1) * my_hermite_deri_weighted(xj, n)};
            res *= std::sqrt(n / 2.);
            return 1. / res;
        }};
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << func_hermite_weighted(res_hermite[i - 1]);
    }
    /*end:prob2.3_hermite_weighted*/
    if (argc > 1 && argv[1][0] == 'o' && argv[1][1] == 'h')
    {
        cout << "\n--------------------------" << '\n'
             << "OUTPUT Hermite interpolation polys" << endl;
        for (size_t i : {2, 16, 128, 1024})
        {
            ofstream ofs{"figures/hermite_weighted_" + std::to_string(i)};
            ofs << setprecision(10);
            double coef{my_hermite_weighted(res_hermite[i - 1], 1023)};
            double x{res_hermite[i - 1] - 10};
            for (size_t j = 0; j < 2000; j++)
            {
                ofs << x << '\t'
                    << coef * my_hermite_weighted(x, 1024) / (x - res_hermite[i - 1]) << '\n';
                x += 0.01;
            }
            ofs.close();
        }
    }
    /*end:prob2.3_hermite_output*/
    cout << "\n--------------------------" << '\n'
         << "j'th Gaussian integration weight of Chebyshev polynomial of degree 1024 T_n(x)" << '\n'
         << "j=\t2\t\t\t16\t\t\t128\t\t\t1024" << endl;
    auto func_cheby{
        [&](double xj) {
            constexpr size_t n{1024};
            double res{my_chebyshev(xj, n - 1) * my_chebyshev_deri(xj, n)};
            res *= .5;
            return 1. / res;
        }};
    for (size_t i : {2, 16, 128, 1024})
    {
        cout << '\t' << func_cheby(res_cheby[i - 1]);
    }
    cout << endl;
    if (argc > 1 && argv[1][0] == 'o' && argv[1][1] == 'c')
    {
        cout << "\n--------------------------" << '\n'
             << "OUTPUT Chebyshev interpolation polys" << endl;
        for (size_t i : {2, 16, 128, 1024})
        {
            ofstream ofs{"figures/chebyshev_weighted_" + std::to_string(i)};
            ofs << setprecision(10);
            double coef{.5 * std::sqrt(1 - res_cheby[i - 1] * res_cheby[i - 1])};
            coef *= my_chebyshev(res_cheby[i - 1], 1023);
            double x{-.9996};
            for (size_t j = 0; j < 4999; j++)
            {
                ofs << x << '\t'
                    << coef * my_chebyshev(x, 1024) / ((x - res_cheby[i - 1]) * sqrt(1 - x * x)) << '\n';
                x += 0.0004;
            }
            ofs.close();
        }
    }
    cout << my_hermite(0.5, 2) << endl;
    cout << my_hermite_deri(0.5, 2) << endl;
    /*end:prob2.3*/
    cout << "==========================" << '\n'
         << "PROBLEM 3.1" << endl;
    if (argc > 1 && argv[1][0] == 'o' && argv[1][1] == 'l' && argv[1][2] == 'z')
    {
        cout << "--------------------------\n"
             << "OUTPUTing zeroes of P_n(x), n=16, 64, 256, 1024" << endl;
        auto legendre_mat_16{legendre_mat<16>()};
        auto legendre_mat_64{legendre_mat<64>()};
        auto legendre_mat_256{legendre_mat<256>()};
        auto legendre_mat_1024{legendre_mat<1024>()};
        auto legendre_zero_16{eig_vals_symm(legendre_mat_16)};
        auto legendre_zero_64{eig_vals_symm(legendre_mat_64)};
        auto legendre_zero_256{eig_vals_symm(legendre_mat_256)};
        auto legendre_zero_1024{eig_vals_symm(legendre_mat_1024)};
        ofstream ofs_16{"figures/legendre_zeroes_16"};
        ofstream ofs_64{"figures/legendre_zeroes_64"};
        ofstream ofs_256{"figures/legendre_zeroes_256"};
        ofstream ofs_1024{"figures/legendre_zeroes_1024"};
        ofs_16 << setprecision(10);
        ofs_64 << setprecision(10);
        ofs_256 << setprecision(10);
        ofs_1024 << setprecision(10);
        for (auto &&x : legendre_zero_16)
            ofs_16 << x << '\n';
        for (auto &&x : legendre_zero_64)
            ofs_64 << x << '\n';
        for (auto &&x : legendre_zero_256)
            ofs_256 << x << '\n';
        for (auto &&x : legendre_zero_1024)
            ofs_1024 << x << '\n';
        /*end:prob3.1_zeroes*/
        cout << "--------------------------\n"
             << "OUTPUTing some eigenfunction of x,\n"
             << "n=16, 64, 256, 1024, j=6, 28, 110, 400" << endl;
        size_t j{6};
        size_t n{16};
        ofstream ofs16{"figures/legendre_poly_16_" + to_string(j)};
        ofs16 << setprecision(10);
        double coef{n / sqrt((2. * n - 1.) * (2. * n + 1.))};
        coef *= my_legendre_norm(legendre_zero_16[j - 1], n - 1);
        double x{-.9996};
        for (size_t i = 0; i < 4999; i++)
        {
            ofs16 << x << '\t'
                  << coef * my_legendre_norm(x, n) / (x - legendre_zero_16[j - 1]) << '\n';
            x += 0.0004;
        }
        ofs16.close();
        // -------------------------------------
        j = 28;
        n = 64;
        ofstream ofs64{"figures/legendre_poly_64_" + to_string(j)};
        ofs64 << setprecision(10);
        coef = n / sqrt((2. * n - 1.) * (2. * n + 1.));
        coef *= my_legendre_norm(legendre_zero_64[j - 1], n - 1);
        x = -.9996;
        for (size_t i = 0; i < 4999; i++)
        {
            ofs64 << x << '\t'
                  << coef * my_legendre_norm(x, n) / (x - legendre_zero_64[j - 1]) << '\n';
            x += 0.0004;
        }
        ofs64.close();
        // -------------------------------------
        j = 110;
        n = 256;
        ofstream ofs256{"figures/legendre_poly_256_" + to_string(j)};
        ofs256 << setprecision(10);
        coef = n / sqrt((2. * n - 1.) * (2. * n + 1.));
        coef *= my_legendre_norm(legendre_zero_256[j - 1], n - 1);
        x = -.9996;
        for (size_t i = 0; i < 4999; i++)
        {
            ofs256 << x << '\t'
                   << coef * my_legendre_norm(x, n) / (x - legendre_zero_256[j - 1]) << '\n';
            x += 0.0004;
        }
        ofs256.close();
        // -------------------------------------
        j = 400;
        n = 1024;
        ofstream ofs1024{"figures/legendre_poly_1024_" + to_string(j)};
        ofs1024 << setprecision(10);
        coef = n / sqrt((2. * n - 1.) * (2. * n + 1.));
        coef *= my_legendre_norm(legendre_zero_1024[j - 1], n - 1);
        x = -.9996;
        for (size_t i = 0; i < 4999; i++)
        {
            ofs1024 << x << '\t'
                    << coef * my_legendre_norm(x, n) / (x - legendre_zero_1024[j - 1]) << '\n';
            x += 0.0004;
        }
        ofs1024.close();
    }
    /*end:prob3.1*/
}

// ====================================================================================

/*beg:Pn*/
double my_legendre_norm(double x, size_t n)
{
    if (n == 0)
    {
        return 0.7071067811865476; // 1/√2
    }
    else if (n == 1)
    {
        return 1.224744871391589 * x; // √(3/2) x
    }
    else
    {
        double prev_P{1.};   // P_0 √(2/(2*0+1))
        double current_P{x}; // P_1 √(2/(2*1+1))
        for (size_t count = 1; count < n; count++)
        {
            prev_P = x * (2 * count + 1) * current_P - count * prev_P;
            prev_P /= count + 1;
            std::swap(prev_P, current_P);
        }
        return current_P * std::sqrt(n + .5);
    }
}
/*end:Pn*/

/*beg:Ln*/
double my_laguerre(double x, size_t n)
{
    if (n == 0)
    {
        return 1.;
    }
    else if (n == 1)
    {
        return 1. - x;
    }
    else
    {
        double prev_P{1.};        // L_0
        double current_P{1. - x}; // L_1
        for (size_t count = 1; count < n; count++)
        {
            prev_P = (1. - x + 2. * count) * current_P - count * prev_P;
            prev_P /= count + 1;
            std::swap(prev_P, current_P);
        }
        return current_P;
    }
}
/*end:Ln*/

/*beg:hermite*/
template <size_t N>
Symm_Band_Matrix<double, N, 1> hermite_mat()
{
    Symm_Band_Matrix<double, N, 1> res{};
    for (size_t i = 1; i < N; i++)
    {
        res(i, i - 1) = std::sqrt(i / 2.);
    }
    return res;
}
/*end:hermite*/

/*beg:chebyshev*/
template <size_t N>
Symm_Band_Matrix<double, N, 1> chebyshev_mat()
{
    Symm_Band_Matrix<double, N, 1> res{};
    res(1, 0) = 0.7071067811865476;
    for (size_t i = 1; i < N - 1; i++)
    {
        res(i + 1, i) = .5;
    }
    return res;
}
/*end:chebyshev*/

/*beg:deriv*/
double poly_deriv(function<double(double, size_t)> func, double x, size_t n)
{
    constexpr double h{1. / (1 << 25)};        // 2^-25 ~ 3e-8
    constexpr double double_h{2. / (1 << 25)}; // 2^-24 ~ 6e-8
    constexpr double inv_12h{(1 << 25) / 12.}; // 1/12h
    double res{func(x + h, n) - func(x - h, n)};
    res *= 8.;
    res += func(x - double_h, n) - func(x + double_h, n);
    return res * inv_12h;
}
/*end:deriv*/

/*beg:hermite*/
double my_hermite(double x, size_t n)
{
    constexpr double inv_qqrt_pi{0.7511255444649425}; // 1 / √√π
    constexpr double sqrt_2{1.4142135623730951};
    if (n == 0)
    {
        return inv_qqrt_pi;
    }
    else if (n == 1)
    {
        return x * sqrt_2 * inv_qqrt_pi;
    }
    else
    {
        double prev_P{1.};            // H_0 * √√π
        double current_P{x * sqrt_2}; // H_1 * √√π
        double sqrt_count{1.};
        for (size_t count = 1; count < n; count++)
        {
            prev_P = x * sqrt_2 * current_P - sqrt_count * prev_P;
            sqrt_count = std::sqrt(count + 1);
            prev_P /= sqrt_count;
            std::swap(prev_P, current_P);
        }
        return current_P * inv_qqrt_pi;
    }
}
/*end:hermite*/

/*beg:hermite_weighted*/
double my_hermite_weighted(double x, size_t n)
{
    constexpr double inv_qqrt_pi{0.7511255444649425}; // 1 / √√π
    constexpr double sqrt_2{1.4142135623730951};
    double weight_every_time{n ? exp(-x * x / (2. * n)) : 1.};
    if (n == 0)
    {
        return inv_qqrt_pi * exp(-x * x / 2.);
    }
    else if (n == 1)
    {
        return x * inv_qqrt_pi * exp(-x * x / 2.);
    }
    else
    {
        double prev_P{weight_every_time};                 // H_0 * √√π
        double current_P{x * sqrt_2 * weight_every_time}; // H_1 * √√π
        double sqrt_count{1.};
        for (size_t count = 1; count < n; count++)
        {
            prev_P = x * sqrt_2 * current_P - sqrt_count * prev_P;
            sqrt_count = std::sqrt(count + 1);
            prev_P /= sqrt_count;
            prev_P *= weight_every_time;
            current_P *= weight_every_time;
            std::swap(prev_P, current_P);
        }
        return current_P * inv_qqrt_pi;
    }
}
/*end:hermite_weighted*/

double my_hermite_deri(double x, size_t n)
{
    return poly_deriv(my_hermite, x, n);
}

/*beg:hermite_deri_weighted*/
double my_hermite_deri_weighted(double x, size_t n)
{
    return x * my_hermite_weighted(x, n) + poly_deriv(my_hermite_weighted, x, n);
}
/*end:hermite_deri_weighted*/

double my_chebyshev(double x, size_t n)
{
    constexpr double inv_sqrt_pi{0.5641895835477563};      // 1 / √π
    constexpr double inv_sqrt_half_pi{0.7978845608028654}; // 1 / √(π/2)
    if (n == 0)
    {
        return inv_sqrt_pi;
    }
    else if (n == 1)
    {
        return x * inv_sqrt_half_pi;
    }
    else if (n == 2)
    {
        return (2. * x * x - 1.) * inv_sqrt_half_pi;
    }
    else
    {
        double prev_P{x};                  // T_1 * √(π/2)
        double current_P{2. * x * x - 1.}; // T_2 * √(π/2)
        for (size_t count = 2; count < n; count++)
        {
            prev_P = x * 2. * current_P - prev_P;
            std::swap(prev_P, current_P);
        }
        return current_P * inv_sqrt_half_pi;
    }
}

double my_chebyshev_deri(double x, size_t n)
{
    return poly_deriv(my_chebyshev, x, n);
}

/*beg:legendre_mat*/
template <size_t N>
Symm_Band_Matrix<double, N, 1> legendre_mat()
{
    Symm_Band_Matrix<double, N, 1> res{};
    for (size_t i = 1; i < N; i++)
    {
        res(i, i - 1) = i / std::sqrt(((i << 1) - 1.) * ((i << 1) + 1.));
    }
    return res;
}
/*end:legendre_mat*/
