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

using namespace std;
using namespace Misc;

constexpr double pi{3.141592653589793};

template <size_t N>
double potential(const array<double, N * 2 - 2> &angles)
{
    double res{};
    for (size_t j = 0; j < N * 2 - 2; j += 2)
    {
        res += .5 / abs(sin(angles[j] / 2.));
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

template <size_t N>
array<double, N * 2 - 2> grad_potential(const array<double, N * 2 - 2> &angles)
{
    array<double, N * 2 - 2> res;
    for (size_t j = 0; j < N * 2 - 2; j += 2)
    {
        double temp{abs(sin(angles[j] / 2.))};
        res[j] = -sin(angles[j]) / (8. * temp * temp * temp);
        res[j + 1] = 0.;
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

template <size_t N>
void modify(array<double, N * 2 - 2> &angles)
{
    for (size_t i = 0; i < N * 2 - 2; i += 2)
    {
        angles[i] = fmod(angles[i], 2. * pi);
        angles[i + 1] = fmod(angles[i + 1], 2. * pi);
    }
}

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
        os << setprecision(16) << solutions.back() << '\t';
    }
    double minimum(*min_element(solutions.begin(), solutions.end()));
    os << "min" << '\t' << setprecision(16) << minimum << endl;
    cerr << N << "\tmin" << '\t' << setprecision(16) << minimum << endl;
}

int main()
{
    ofstream ofs{"4_Thomson_1"};
    if (!ofs)
    {
        throw runtime_error("cannot open the file");
    }
    std::default_random_engine ran{};
    /* ============== first_run ================ */
    // solve<2>(ran, ofs);
    // solve<3>(ran, ofs);
    // solve<4>(ran, ofs);
    // solve<5>(ran, ofs);
    // solve<6>(ran, ofs);
    // solve<7>(ran, ofs);
    // solve<8>(ran, ofs);
    // solve<9>(ran, ofs);
    // solve<10>(ran, ofs);
    // solve<11>(ran, ofs);
    // solve<12>(ran, ofs);
    // solve<13>(ran, ofs);
    // solve<14>(ran, ofs);
    // solve<15>(ran, ofs);
    // solve<16>(ran, ofs);
    // solve<17>(ran, ofs);
    // solve<18>(ran, ofs);
    // solve<19>(ran, ofs);
    // solve<20>(ran, ofs);
    // solve<21>(ran, ofs);
    // solve<22>(ran, ofs);
    // solve<23>(ran, ofs);
    // solve<24>(ran, ofs);
    // solve<25>(ran, ofs);
    // solve<26>(ran, ofs);
    // solve<27>(ran, ofs);
    // solve<28>(ran, ofs);
    // solve<29>(ran, ofs);
    // solve<30>(ran, ofs);
    // solve<31>(ran, ofs);
    // solve<32>(ran, ofs);
    // solve<33>(ran, ofs);
    // solve<34>(ran, ofs);
    // solve<35>(ran, ofs);
    // solve<36>(ran, ofs);
    // solve<37>(ran, ofs);
    // solve<38>(ran, ofs);
    // solve<39>(ran, ofs);
    // solve<40>(ran, ofs);
    // solve<41>(ran, ofs);
    // solve<42>(ran, ofs);
    // solve<43>(ran, ofs);
    // solve<44>(ran, ofs);
    // solve<45>(ran, ofs);
    // solve<46>(ran, ofs);
    // solve<47>(ran, ofs);
    // solve<48>(ran, ofs);
    // solve<49>(ran, ofs);
    // solve<50>(ran, ofs);
    // solve<51>(ran, ofs);
    // solve<52>(ran, ofs);
    // solve<53>(ran, ofs);
    // solve<54>(ran, ofs);
    // solve<55>(ran, ofs);
    // solve<56>(ran, ofs);
    // solve<57>(ran, ofs);
    // solve<58>(ran, ofs);
    // solve<59>(ran, ofs);
    // solve<60>(ran, ofs);
    // solve<61>(ran, ofs);
    // solve<62>(ran, ofs);
    // solve<63>(ran, ofs);
    // solve<64>(ran, ofs);
    /* =============== second runs ================ */
    // solve<16>(ran, ofs, 10);
    // solve<36>(ran, ofs, 10);
    // solve<56>(ran, ofs, 10);
    // solve<60>(ran, ofs, 10);
    // solve<64>(ran, ofs, 10);
}
