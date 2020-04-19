#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <iomanip>
#include <random>
#include <chrono>
#include <map>
#include <utility>

#include "../misc/optimize.h"
#include "../misc/eigvals.h"

using namespace std;
using namespace Misc;

constexpr double pi{3.141592653589793};

/*beg:potential*/
// the first point is stick to theta=pi/2, phi=0
template <size_t N>
double potential(const array<double, N * 2 - 2> &angles)
{
    double res{};
    for (size_t j = 0; j < N * 2 - 2; j += 2)
    {
        res += 1. / sqrt(2. * (1. - sin(angles[j]) * cos(angles[j + 1])));
    }
    for (size_t i = 0; i < N * 2 - 2; i += 2)
    {
        double theta1{angles[i]};
        double phi1{angles[i + 1]};
        double sin_theta1{sin(theta1)};
        for (size_t j = i + 2; j < N * 2 - 2; j += 2)
        {
            double theta2{angles[j]};
            double temp{1. - cos(theta1 - theta2) + sin_theta1 * sin(theta2) * (1. - cos(phi1 - angles[j + 1]))};
            res += 1. / sqrt(2. * temp);
        }
    }
    return res;
}
/*end:potential*/

/*beg:distances*/
template <size_t N>
vector<double> distances(const array<double, N * 2 - 2> &angles, ostream &os = cout)
{
    vector<double> res{};
    for (size_t j = 0; j < N * 2 - 2; j += 2)
    {
        res.push_back(sqrt(2. * (1. - sin(angles[j]) * cos(angles[j + 1]))));
    }
    for (size_t i = 0; i < N * 2 - 2; i += 2)
    {
        double theta1{angles[i]};
        double phi1{angles[i + 1]};
        double sin_theta1{sin(theta1)};
        for (size_t j = i + 2; j < N * 2 - 2; j += 2)
        {
            double theta2{angles[j]};
            double temp{1. - cos(theta1 - theta2) + sin_theta1 * sin(theta2) * (1. - cos(phi1 - angles[j + 1]))};
            res.push_back(sqrt(2. * temp));
        }
    }
    sort(res.begin(), res.end());
    return res;
}
/*end:distances*/

/*beg:grad_potential*/
template <size_t N>
array<double, N * 2 - 2> grad_potential(const array<double, N * 2 - 2> &angles)
{
    array<double, N * 2 - 2> res;
    for (size_t j = 0; j < N * 2 - 2; j += 2)
    {
        double r3{2. * (1. - sin(angles[j]) * cos(angles[j + 1]))};
        r3 *= sqrt(r3);
        res[j] = cos(angles[j]) * cos(angles[j + 1]) / r3;
        res[j + 1] = -sin(angles[j]) * sin(angles[j + 1]) / r3;
    }
    for (size_t i = 0; i < N * 2 - 2; i += 2)
    {
        double theta1{angles[i]};
        double phi1{angles[i + 1]};
        double sin_theta1{sin(theta1)};
        double cos_theta1{cos(theta1)};
        for (size_t j = i + 2; j < N * 2 - 2; j += 2)
        {
            double theta2{angles[j]};
            double phi2{angles[j + 1]};
            double sin_theta2{sin(theta2)};
            double cos_phis{cos(phi1 - phi2)};                                                       // cos(phi1-phi2)
            double r3{2. * (1. - cos(theta1 - theta2) + sin_theta1 * sin_theta2 * (1. - cos_phis))}; // r^3
            r3 *= sqrt(r3);
            double sin_thetas{sin(theta1 - theta2)}; // sin(theta1-theta2)
            res[i] -= (sin_thetas + cos_theta1 * sin_theta2 * (1. - cos_phis)) / r3;
            double grad_phi{-sin_theta1 * sin_theta2 * sin(phi1 - phi2) / r3};
            res[i + 1] += grad_phi;
            res[j] -= (-sin_thetas + cos(theta2) * sin_theta1 * (1. - cos_phis)) / r3;
            res[j + 1] -= grad_phi;
        }
    }
    return res;
}
/*end:grad_potential*/

/*beg:hessian_potential*/
template <size_t N>
Matrix<double, N * 2, N * 2> hessian_potential(const array<double, N * 2> &angles)
{
    Matrix<double, N * 2, N * 2> res{};
    for (size_t i = 0; i < N * 2; i += 2)
    {
        double theta1{angles[i]};
        double phi1{angles[i + 1]};
        double sin_theta1{sin(theta1)};
        double cos_theta1{cos(theta1)};
        for (size_t j = i + 2; j < N * 2; j += 2)
        {
            double theta2{angles[j]};
            double phi2{angles[j + 1]};
            double sin_theta2{sin(theta2)};
            double cos_theta2{cos(theta2)};
            double cos_phis{cos(phi1 - phi2)};                                             // cos(phi1-phi2)
            double sin_phis{sin(phi1 - phi2)};                                             // cos(phi1-phi2)
            double sin_thetas{sin(theta1 - theta2)};                                       // cos(theta1-theta2)
            double cos_thetas{cos(theta1 - theta2)};                                       // cos(theta1-theta2)
            double r2{2. * (1. - cos_thetas + sin_theta1 * sin_theta2 * (1. - cos_phis))}; // r^2
            double r3{r2 * sqrt(r2)};                                                      // r^3
            double r5{r2 * r3};                                                            // r^5
            double Theta1{-sin_thetas - cos_theta1 * sin_theta2 * (1. - cos_phis)};
            double Phi1{-sin_theta1 * sin_theta2 * sin_phis};
            double Theta2{sin_thetas - sin_theta1 * cos_theta2 * (1. - cos_phis)};
            double temp{(-cos_thetas + sin_theta1 * sin_theta2 * (1. - cos_phis)) / r3};
            res(i, i) += temp + 3 * Theta1 * Theta1 / r5;
            res(j, j) += temp + 3 * Theta2 * Theta2 / r5;
            temp = -cos_theta1 * sin_theta2 * sin_phis / r3 + 3 * Phi1 * Theta1 / r5;
            res(i, i + 1) += temp;
            res(i, j + 1) -= temp;
            res(i, j) += (cos_thetas - cos_theta1 * cos_theta2 * (1. - cos_phis)) / r3 + 3 * Theta1 * Theta2 / r5;
            temp = -sin_theta1 * sin_theta2 * cos_phis / r3 + 3 * Phi1 * Phi1 / r5;
            res(i + 1, i + 1) += temp;
            res(i + 1, j + 1) -= temp;
            res(j + 1, j + 1) += temp;
            temp = -sin_theta1 * cos_theta2 * sin_phis / r3 + 3 * Phi1 * Theta2 / r5;
            res(i + 1, j) += temp;
            res(j, j + 1) -= temp;
        }
    }
    for (size_t i = 0; i < N * 2; i++)
        for (size_t j = 0; j < i; j++)
        {
            res(i, j) = res(j, i);
        }
    return res;
}
/*end:hessian_potential*/

/*beg:modify*/
template <size_t N>
void modify(array<double, N * 2 - 2> &angles)
{
    for (size_t i = 0; i < N * 2 - 2; i += 2)
    {
        angles[i] = fmod(angles[i], 2. * pi);
        angles[i + 1] = fmod(angles[i + 1], 2. * pi);
    }
}
/*end:modify*/

/*beg:solve*/
template <size_t N>
void solve(default_random_engine &ran, ostream &os = cout, size_t n = 5)
{
    array<double, 2 * N - 2> x0;
    vector<double> solutions;
    os << N << '\t';
    for (size_t i = 0; i < n; i++)
    {
        for (size_t i = 0; i < 2 * N - 2; i += 2)
        {
            x0[i] = acos(uniform_real_distribution<>{-1., 1.}(ran));
            x0[i + 1] = uniform_real_distribution<>{0., 2. * pi}(ran);
        }
        auto res{conj_grad<double, 2 * N - 2>(function(potential<N>), function(grad_potential<N>), x0, function(modify<N>))};
        solutions.push_back(potential<N>(res));
        os << solutions.back() << '\t';
    }
    double minimum(*min_element(solutions.begin(), solutions.end()));
    os << "min" << '\t' << minimum << endl;
    cerr << N << "\tmin" << '\t' << minimum << endl;
}
/*end:solve*/

int main()
{
    ofstream ofs{"4_Thomson_1"};
    if (!ofs)
    {
        throw runtime_error("cannot open the file");
    }
    std::default_random_engine ran{747579};
    ofs << setprecision(16);
    cout << setprecision(16);
    cerr << setprecision(16);
    /* ============== problem 1 ================ */
    solve<2>(ran, ofs);
    solve<3>(ran, ofs);
    solve<4>(ran, ofs);
    solve<5>(ran, ofs);
    solve<6>(ran, ofs);
    solve<7>(ran, ofs);
    solve<8>(ran, ofs);
    solve<9>(ran, ofs);
    solve<10>(ran, ofs);
    solve<11>(ran, ofs);
    solve<12>(ran, ofs);
    solve<13>(ran, ofs);
    solve<14>(ran, ofs);
    solve<15>(ran, ofs);
    solve<16>(ran, ofs);
    solve<17>(ran, ofs);
    solve<18>(ran, ofs);
    solve<19>(ran, ofs);
    solve<20>(ran, ofs);
    solve<21>(ran, ofs);
    solve<22>(ran, ofs);
    solve<23>(ran, ofs);
    solve<24>(ran, ofs);
    solve<25>(ran, ofs);
    solve<26>(ran, ofs);
    solve<27>(ran, ofs);
    solve<28>(ran, ofs);
    solve<29>(ran, ofs);
    solve<30>(ran, ofs);
    solve<31>(ran, ofs);
    solve<32>(ran, ofs);
    solve<33>(ran, ofs);
    solve<34>(ran, ofs);
    solve<35>(ran, ofs);
    solve<36>(ran, ofs);
    solve<37>(ran, ofs);
    solve<38>(ran, ofs);
    solve<39>(ran, ofs);
    solve<40>(ran, ofs);
    solve<41>(ran, ofs);
    solve<42>(ran, ofs);
    solve<43>(ran, ofs);
    solve<44>(ran, ofs);
    solve<45>(ran, ofs);
    solve<46>(ran, ofs);
    solve<47>(ran, ofs);
    solve<48>(ran, ofs);
    solve<49>(ran, ofs);
    solve<50>(ran, ofs);
    solve<51>(ran, ofs);
    solve<52>(ran, ofs);
    solve<53>(ran, ofs);
    solve<54>(ran, ofs);
    solve<55>(ran, ofs);
    solve<56>(ran, ofs);
    solve<57>(ran, ofs);
    solve<58>(ran, ofs);
    solve<59>(ran, ofs);
    solve<60>(ran, ofs);
    solve<61>(ran, ofs);
    solve<62>(ran, ofs);
    solve<63>(ran, ofs);
    solve<64>(ran, ofs);
    /* =============== append runs ================ */
    solve<47>(ran, ofs);
    solve<53>(ran, ofs);
    solve<56>(ran, ofs);
    solve<58>(ran, ofs);
    /*beg:12calc*/
    /* ============== problem 2 ================= */
    constexpr size_t N{12};
    array<double, 2 * N - 2> x0;
    for (size_t i = 0; i < 2 * N - 2; i += 2)
    {
        x0[i] = acos(uniform_real_distribution<>{-1., 1.}(ran));
        x0[i + 1] = uniform_real_distribution<>{0., 2. * pi}(ran);
    }

    auto res{conj_grad<double, 2 * N - 2>(function(potential<N>), function(grad_potential<N>), x0, function(modify<N>))};
    cout << "N=" << N << ", minimum potential: " << potential<N>(res) << endl;
    /*end:12calc*/
    // check whether they form a regular icosahedron
    auto dists{distances<N>(res)};
    cout << "minimum distance:\t" << dists[0] << endl;
    cout << "30th minimum distance:\t" << dists[30 - 1] << endl; // regular icosahedron have 30 edges
    /*end:icosahedron*/
    // add the fixed point
    array<double, 2 * N> complete_res{pi / 2., 0.};
    copy(res.begin(), res.end(), &complete_res[2]);
    cout << "\ntheta, phi pairs:\n";
    for (auto x: complete_res)
    {
        cout << x << ", ";
    }
    cout << endl;
    // compute Hessian (of V)
    auto hessian{hessian_potential<N>(complete_res)};
    /*end:calc_hessian*/
    // compute 1/sqrt(T)
    array<double, 2 * N> diag_line;
    for (size_t i = 0; i < 2 * N; i += 2)
    {
        diag_line[i] = 1.;
        diag_line[i + 1] = 1. / abs(sin(complete_res[i]));
    }
    Diag_Matrix<double, 2 * N> rev_sqrt_kinetic{diag_line.begin(), diag_line.end()};
    // make T and V diagonal simultaneously
    hessian = rev_sqrt_kinetic * hessian * rev_sqrt_kinetic;
    /*end:calc_V_prime*/
    // calc eigenvalues (i.e. omega^2)
    auto eig_pair{eig_vals(hessian)};
    if (eig_pair.second.size())
    {
        cerr << __FILE__ << ':' << __LINE__ << ": complex eigenvalues occurred, at:\n\t";
        for (auto i : eig_pair.first)
        {
            cerr << i << ", ";
        }
        cerr << endl;
    }
    else
    {
        sort(eig_pair.first.begin(), eig_pair.first.end());
    }
    cout << "\neigenvalues (i.e. omega^2):\n";
    for (auto x : eig_pair.first)
    {
        cout << x << ", ";
    }
    /*end:eigs*/
}
