#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>
#include <iomanip>

#include "../misc/linear_eq.h"

using namespace std;
using namespace Misc;

constexpr double pi{3.141592653589793};

// return the matrix for Laplacian operator
// template parameters:
//      M: partition in x       N: partition in y
// parameeters:
//      xrange                  yrange
template <size_t M, size_t N>
Sparse_Matrix<double, (M - 1) * (N - 1), (M - 1) * (N - 1)> laplacian(double xrange, double yrange);
// return the inhomogeneous item 2 * pi^2 * sin(pi*x) * sin(pi*y) for the Poisson eq
// template parameters:
//      M: partition in x       N: partition in y
// parameeters:
//      xa, xb, ya, yb
template <size_t M, size_t N>
array<double, (M - 1) * (N - 1)> sinxy(double xa, double xb, double ya, double yb);
// return the interpolation result of the Gaussian points
// parameeters:
//      grid_res: values at grid
template <size_t M, size_t N>
array<double, 4 * M * N> gauss_pts(const array<double, (M - 1) * (N - 1)> &grid_res);
/*beg:ana_sol*/
// analytical solution of this problem
double ana_sol(double x, double y)
{
    return sin(pi * x) * sin(pi * y);
}
/*end:ana_sol*/
/*beg:solve_decl*/
template <size_t N>
void solve(bool output = true, string prefix = "3_Poisson/");
/*end:solve_decl*/

int main()
{
    solve<4>();
    solve<5>();
    solve<6>();
    solve<7>();
}

/*beg:laplacian*/
template <size_t M, size_t N>
Sparse_Matrix<double, (M - 1) * (N - 1), (M - 1) * (N - 1)> laplacian(double xrange, double yrange)
{
    static_assert(M >= 1 && N >= 1);
    constexpr size_t mat_size{(M - 1) * (N - 1)};
    Sparse_Matrix<double, mat_size, mat_size> res{};
    double hx{xrange / M};
    double hy{yrange / N};
    double delx2{1 / hx / hx};
    double dely2{1 / hy / hy};

    double value{2 * (delx2 + dely2)};
    for (size_t i = 0; i < mat_size; i++)
    {
        res(i, i) = value;
    }

    value = -dely2;
    for (size_t i = 0; i < mat_size; i += N - 1)
        for (size_t j = i; j < i + N - 2; j++)
        {
            res(j, j + 1) = value;
            res(j + 1, j) = value;
        }

    value = -delx2;
    for (size_t i = 0; i < mat_size - (N - 1); i++)
    {
        res(i, i + N - 1) = value;
        res(i + N - 1, i) = value;
    }

    return res;
}
/*end:laplacian*/

/*beg:sinxy*/
template <size_t M, size_t N>
array<double, (M - 1) * (N - 1)> sinxy(double xa, double xb, double ya, double yb)
{
    static_assert(M >= 1 && N >= 1);
    constexpr size_t mat_size{(M - 1) * (N - 1)};
    array<double, mat_size> res;
    double hx{(xb - xa) / M};
    double hy{(yb - ya) / N};

    for (size_t i = 1, k = 0; i <= M - 1; i++)
    {
        double sinx{2 * pi * pi * sin(pi * (xa + i * hx))};
        for (size_t j = 1; j <= N - 1; j++, k++)
        {
            res[k] = sinx * sin(pi * (ya + j * hy));
        }
    }

    return res;
}
/*end:sinxy*/

/*beg:gauss_pts*/
template <size_t M, size_t N>
array<double, 4 * M * N> gauss_pts(const array<double, (M - 1) * (N - 1)> &grid_res)
{
    array<double, 4 * M * N> res;
    constexpr double _2_sqrt3{2. + sqrt(3.)};
    res.fill(0.);
    for (size_t i = 0; i < M - 1; i++)
    {
        size_t count{(i + 1) * (N - 1)};
        for (size_t value_idx = i * (N - 1), res_idx = 4 * i * N; value_idx < count; value_idx++, res_idx += 2)
        {
            res[res_idx] += grid_res[value_idx] / 6. / _2_sqrt3;
            res[res_idx + 1] += grid_res[value_idx] / 6.;
            res[res_idx + 2] += grid_res[value_idx] / 6.;
            res[res_idx + 3] += grid_res[value_idx] / 6. / _2_sqrt3;
        }
        for (size_t value_idx = i * (N - 1), res_idx = (4 * i + 2) * N; value_idx < count; value_idx++, res_idx += 2)
        {
            res[res_idx] += grid_res[value_idx] / 6.;
            res[res_idx + 1] += grid_res[value_idx] / 6. * _2_sqrt3;
            res[res_idx + 2] += grid_res[value_idx] / 6. * _2_sqrt3;
            res[res_idx + 3] += grid_res[value_idx] / 6.;
        }
        for (size_t value_idx = i * (N - 1), res_idx = (4 * i + 4) * N; value_idx < count; value_idx++, res_idx += 2)
        {
            res[res_idx] += grid_res[value_idx] / 6.;
            res[res_idx + 1] += grid_res[value_idx] / 6. * _2_sqrt3;
            res[res_idx + 2] += grid_res[value_idx] / 6. * _2_sqrt3;
            res[res_idx + 3] += grid_res[value_idx] / 6.;
        }
        for (size_t value_idx = i * (N - 1), res_idx = (4 * i + 6) * N; value_idx < count; value_idx++, res_idx += 2)
        {
            res[res_idx] += grid_res[value_idx] / 6. / _2_sqrt3;
            res[res_idx + 1] += grid_res[value_idx] / 6.;
            res[res_idx + 2] += grid_res[value_idx] / 6.;
            res[res_idx + 3] += grid_res[value_idx] / 6. / _2_sqrt3;
        }
    }
    return res;
}
/*end:gauss_pts*/

/*beg:solve_init*/
template <size_t N>
void solve(bool output, string prefix)
{
    ostringstream oss;
    oss << prefix << N << "_";
    prefix = oss.str();
    constexpr size_t partition{1 << N};
    constexpr size_t mat_size{(partition - 1) * (partition - 1)};

    constexpr double xa{0.};
    constexpr double xb{1.};
    constexpr double ya{0.};
    constexpr double yb{1.};
    /*end:solve_init*/
    // 1. solve the equation: op * solution = b
    auto op{laplacian<partition, partition>(xb - xa, yb - ya)};
    array<double, mat_size> b{sinxy<partition, partition>(xa, xb, ya, yb)};
    array<double, mat_size> solution;
    solution.fill(1.);
    suc_over_rel(op, b, solution, 1.95);
    if (output)
    {
        ofstream ofs{prefix + "1"};
        for (auto &&i : solution)
        {
            ofs << setprecision(16) << i << '\n';
        }
    }
    /*end:prob_1*/
    // 2. linear interpolation; get value at gauss points
    array<double, 4 * partition * partition> interp{gauss_pts<partition, partition>(solution)};
    if (output)
    {
        ofstream ofs{prefix + "2"};
        for (auto &&i : interp)
        {
            ofs << setprecision(16) << i << '\n';
        }
    }
    /*end:prob_2*/
    // 3. Gauss integration of each grid
    // 3.1 difference between exact value at gauss points; modify `interp` in situ
    constexpr double left_gauss{1. / (3. + sqrt(3.))};
    constexpr double hx{(xb - xa) / partition};
    constexpr double hy{(yb - ya) / partition};
    for (size_t i = 0; i < partition; i++)
    {
        double x{xa + i * hx + left_gauss * hx};
        for (size_t j = 0, idx = i * 4 * partition; j < partition; j++)
        {
            double y{ya + j * hy + left_gauss * hy};
            interp[idx] -= ana_sol(x, y);
            idx++;
            y += (1 - 2 * left_gauss) * hy;
            interp[idx] -= ana_sol(x, y);
            idx++;
        }
        x += (1 - 2 * left_gauss) * hx;
        for (size_t j = 0, idx = (i * 4 + 2) * partition; j < partition; j++)
        {
            double y{ya + j * hy + left_gauss * hy};
            interp[idx] -= ana_sol(x, y);
            idx++;
            y += (1 - 2 * left_gauss) * hy;
            interp[idx] -= ana_sol(x, y);
            idx++;
        }
    }
    /*end:prob_3.1*/
    // 3.2 integration in each grid
    array<double, partition * partition> integrals;
    for (size_t k = 0; k < partition * partition; k++)
    {
        size_t orig{k / partition * 4 * partition + k % partition * 2};
        integrals[k] = interp[orig] * interp[orig] + interp[orig + 1] * interp[orig + 1];
        orig += 2 * partition;
        integrals[k] += interp[orig] * interp[orig] + interp[orig + 1] * interp[orig + 1];
        integrals[k] *= hx * hy / 4.;
    }
    if (output)
    {
        ofstream ofs{prefix + "3"};
        for (auto &&i : integrals)
        {
            ofs << setprecision(16) << i << '\n';
        }
    }
    /*end:prob_3*/
    // 4. Error
    double error{};
    for (auto x : integrals)
    {
        error += x;
    }
    error = sqrt(error);
    if (output)
    {
        ofstream ofs{prefix + "4"};
        ofs << "partition number is 2^" << N << '\n';
        ofs << "rms error:\t" << setprecision(16) << error << endl;
    }
    cerr << "partition number is 2^" << N << '\n';
    cerr << "rms error:\t" << setprecision(16) << error << endl;
    /*end:prob_4*/
}
