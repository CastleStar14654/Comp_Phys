#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <map>

#include "../misc/Matrix_Catalogue.h"
#include "../misc/linear_eq_direct.h"
#include "../misc/linear_eq_iterative.h"

using namespace Misc;
using namespace std;

template <size_t N>
inline void calc_square();
template <size_t N>
inline Symm_Band_Matrix<double, N *(N + 2), N + 1> square_network();
template <size_t N>
inline size_t _square_index(size_t x, size_t y);
template <size_t N>
inline pair<size_t, size_t> _square_xy(size_t index);
template <size_t N>
inline vector<size_t> _square_connections(size_t index);

template <size_t N>
inline void calc_triangle();
template <size_t N>
inline Symm_Band_Matrix<double, N *(N + 3) / 2, N + 1> triangle_network();
template <size_t N>
inline size_t _triangle_index(size_t x, size_t y);
template <size_t N>
inline pair<size_t, size_t> _triangle_xy(size_t index);
template <size_t N>
inline vector<size_t> _triangle_connections(size_t index);

template <size_t N>
inline void calc_hex();
template <size_t N>
inline Symm_Band_Matrix<double, N *(N + 4), N + 1> hex_network();
template <size_t N>
inline size_t _hex_index(size_t x, size_t y);
template <size_t N>
inline pair<size_t, size_t> _hex_xy(size_t index);
template <size_t N>
inline vector<size_t> _hex_connections(size_t index);

enum class AC_Type
{
    resister,
    inductance,
    capacitance
};
template <size_t N>
inline void calc_triangle_ac(double omega);
template <size_t N>
inline Symm_Band_Matrix<complex<double>, N *(N + 3) / 2, N + 1> triangle_ac_network(double omega, bool hermite=false);
template <size_t N>
inline map<size_t, AC_Type> _triangle_ac_connections(size_t index);

int main()
{
    cout << "==== SQUARE ====" << endl;
    calc_square<1>();
    calc_square<4>();
    calc_square<16>();
    calc_square<64>();
    cout << "==== TRIANGLE ====" << endl;
    calc_triangle<1>();
    calc_triangle<4>();
    calc_triangle<16>();
    calc_triangle<64>();
    cout << "==== HEXAGON ====" << endl;
    calc_hex<4>();
    calc_hex<16>();
    calc_hex<64>();
    cout << "==== TRIANGLE AC ====" << endl;
    for (size_t i = 0; i < 50; i++)
    {
        static double omega;
        omega = 0.05 + 0.1*i;
        calc_triangle_ac<4>(omega);
    }
}

template <size_t N>
inline void calc_square()
{
    constexpr size_t mat_side{N * (N + 2)};
    constexpr size_t band_width{(N + 1)};
    cout << "===========================\n"
         << "N = " << N << endl;

    cout << "generating coefficient matrix... ";
    auto mat_g{square_network<N>()};
    cout << "finished" << endl;

    cout << "------- a, b ---------" << endl;

    cout << ">>> direct method --- ldl" << endl;

    Low_Band_Matrix<double, mat_side, band_width> l{};
    array<double, mat_side> d{};

    cout << "LDL factoring... ";
    ldl_factor(mat_g, l, d);
    cout << "finished ";

    cout << "generating b... ";
    array<double, mat_side> b{};
    size_t idx_b{_square_index<N>(N, 0)};
    b[idx_b] = 1;
    cout << "finished ";

    cout << "calculating x... ";
    array<double, mat_side> x1{};
    back_sub(l, b, x1);
    back_sub(d, x1);
    back_sub(l, x1, true);
    cout << "finished" << endl;

    cout << ">>> U_ba = " << setprecision(16) << x1[idx_b] << endl;

    cout << ">>> iterative method --- conjugate gradient" << endl;

    cout << "solving x... ";
    array<double, mat_side> x2{};
    if (conj_grad(mat_g, b, x2, true) == 0)
    {
        cout << "finished" << endl;

        cout << ">>> U_ba = " << setprecision(16) << x2[idx_b] << endl;
    }
    else
    {
        cerr << "not converged" << endl;
    }

    cout << "------- a, c ---------" << endl;

    cout << ">>> direct method --- ldl" << endl;

    cout << "changing b... ";
    b[idx_b] = 0;
    idx_b = _square_index<N>(N, N);
    b[idx_b] = 1;
    cout << "finished ";

    cout << "calculating x... ";
    back_sub(l, b, x1);
    back_sub(d, x1);
    back_sub(l, x1, true);
    cout << "finished" << endl;

    cout << ">>> U_ca = " << setprecision(16) << x1[idx_b] << endl;

    cout << ">>> iterative method --- conjugate gradient" << endl;

    cout << "solving x... ";
    if (conj_grad(mat_g, b, x2, true) == 0)
    {
        cout << "finished" << endl;

        cout << ">>> U_ca = " << setprecision(16) << x2[idx_b] << endl;
    }
    else
    {
        cerr << "not converged" << endl;
    }

    cout << "=========== END ===========" << endl;
}

// the N-sides*N-sides network is like:
//     N*(N+1)-1   N*(N+1)     ... ...     (N+1)^2 -2
//     ...         ...         ... ...     ...
//     N           N+1         ... ...     2*N
//     X           0           ... ...     N-1
// then, point X is set as 0V
// so, the return matrix is of (N^2+2*N) * (N^2+2*N)
template <size_t N>
inline Symm_Band_Matrix<double, N *(N + 2), N + 1> square_network()
{
    constexpr size_t mat_length{N * (N + 2)};
    Symm_Band_Matrix<double, mat_length, N + 1> res{};
    vector<size_t> connections;

    for (size_t i = 0; i < mat_length; i++)
    {
        connections = _square_connections<N>(i);
        res(i, i) += connections.size();
        for (const auto &j : connections)
        {
            res(j, j)++;
            res(j, i) = -1;
        }
    }

    // did not count point (0, 0)'s effect on point (1, 0) & (0, 1)
    // their diagonal item should ++
    for (const auto &i : {_square_index<N>(1, 0), _square_index<N>(0, 1)})
    {
        res(i, i)++;
    }

    return res;
}

// the N-sides*N-sides network is like:
// y ^  N*(N+1)-1   N*(N+1)     ... ...     (N+1)^2 -2
//   |  ...         ...         ... ...     ...
//   |  N           N+1         ... ...     2*N
//   |  X           0           ... ...     N-1
//   -----------------------------------------------> x
// also, point X is set as 0V, which is not indexed
// this function return the index of a point (x, y)
template <size_t N>
inline size_t _square_index(size_t x, size_t y)
{
    return y * (N + 1) + x - 1;
}

// the N-sides*N-sides network is like:
// y ^  N*(N+1)-1   N*(N+1)     ... ...     (N+1)^2 -2
//   |  ...         ...         ... ...     ...
//   |  N           N+1         ... ...     2*N
//   |  X           0           ... ...     N-1
//   -----------------------------------------------> x
// also, point X is set as 0V, which is not indexed
// this function return (x, y) of a point at index
template <size_t N>
inline pair<size_t, size_t> _square_xy(size_t index)
{
    return {(index + 1) % (N + 1), (index + 1) / (N + 1)};
}

// return connections with higher indices
template <size_t N>
inline vector<size_t> _square_connections(size_t index)
{
    auto xy{_square_xy<N>(index)};
    vector<size_t> res{};
    if (xy.first < N)
    {
        res.push_back(_square_index<N>(xy.first + 1, xy.second));
    }
    if (xy.second < N)
    {
        res.push_back(_square_index<N>(xy.first, xy.second + 1));
    }
    return res;
}

// ==============================================================

template <size_t N>
inline void calc_triangle()
{
    constexpr size_t mat_side{N * (N + 3) / 2};
    constexpr size_t band_width{(N + 1)};
    cout << "===========================\n"
         << "N = " << N << endl;

    cout << "generating coefficient matrix... ";
    auto mat_g{triangle_network<N>()};
    cout << "finished" << endl;

    cout << "------- a, b ---------" << endl;

    cout << ">>> direct method --- ldl" << endl;

    Low_Band_Matrix<double, mat_side, band_width> l{};
    array<double, mat_side> d{};

    cout << "LDL factoring... ";
    ldl_factor(mat_g, l, d);
    cout << "finished ";

    cout << "generating b... ";
    array<double, mat_side> b{};
    size_t idx_b{_triangle_index<N>(0, N)};
    b[idx_b] = 1;
    cout << "finished ";

    cout << "calculating x... ";
    array<double, mat_side> x1{};
    back_sub(l, b, x1);
    back_sub(d, x1);
    back_sub(l, x1, true);
    cout << "finished" << endl;

    cout << ">>> U_ba = " << setprecision(16) << x1[idx_b] << endl;

    cout << ">>> iterative method --- conjugate gradient" << endl;

    cout << "solving x... ";
    array<double, mat_side> x2{};
    if (conj_grad(mat_g, b, x2, true) == 0)
    {
        cout << "finished" << endl;

        cout << ">>> U_ba = " << setprecision(16) << x2[idx_b] << endl;
    }
    else
    {
        cerr << "not converged" << endl;
    }

    cout << "=========== END ===========" << endl;
}

// the N-sides triangular network is like:
//     N*(N+3)/2-1
//     N*(N+3)/2-3 N*(N+3)/2-2
//     ...         ...         ...
//     2*N         2*N+1       ...   3*N-2
//     N           N+1         ... ...     2*N-1
//     X           0           ... ...     N-2      N-1
// then, point X is set as 0V
// so, the return matrix is of (N*(N+3)/2) * (N*(N+3)/2)
template <size_t N>
inline Symm_Band_Matrix<double, N *(N + 3) / 2, N + 1> triangle_network()
{
    constexpr size_t mat_length{N * (N + 3) / 2};
    Symm_Band_Matrix<double, mat_length, N + 1> res{};
    vector<size_t> connections;

    for (size_t i = 0; i < mat_length; i++)
    {
        connections = _triangle_connections<N>(i);
        res(i, i) += connections.size();
        for (const auto &j : connections)
        {
            res(j, j)++;
            res(j, i) = -1;
        }
    }

    // did not count point (0, 0)'s effect on point (1, 0) & (0, 1)
    // their diagonal item should ++
    for (const auto &i : {_triangle_index<N>(1, 0), _triangle_index<N>(0, 1)})
    {
        res(i, i)++;
    }

    return res;
}

// the N-sides triangular network is like:
// y ^  N*(N+3)/2-1
//   |  N*(N+3)/2-3 N*(N+3)/2-2
//   |  ...         ...         ...
//   |  2*N         2*N+1       ...   3*N-2
//   |  N           N+1         ... ...     2*N-1
//   |  X           0           ... ...     N-2      N-1
//   -----------------------------------------------> x
// also, point X is set as 0V, which is not indexed
// this function return the index of a point (x, y)
template <size_t N>
inline size_t _triangle_index(size_t x, size_t y)
{
    return y * (2 * N + 3 - y) / 2 + x - 1;
}

// the N-sides triangular network is like:
// y ^  N*(N+3)/2-1
//   |  N*(N+3)/2-3 N*(N+3)/2-2
//   |  ...         ...         ...
//   |  2*N         2*N+1       ...   3*N-2
//   |  N           N+1         ... ...     2*N-1
//   |  X           0           ... ...     N-2      N-1
//   -----------------------------------------------> x
// also, point X is set as 0V, which is not indexed
// this function return (x, y) of a point at index
template <size_t N>
inline pair<size_t, size_t> _triangle_xy(size_t index)
{
    size_t x{index + 1};
    size_t y{0};
    for (; x > N - y; y++)
    {
        x -= N + 1 - y;
    }

    return {x, y};
}

// return connections with higher indices
template <size_t N>
inline vector<size_t> _triangle_connections(size_t index)
{
    auto xy{_triangle_xy<N>(index)};
    vector<size_t> res{};
    if (xy.first < N - xy.second)
    {
        res.push_back(_triangle_index<N>(xy.first + 1, xy.second));
    }
    if (xy.first > 0)
    {
        res.push_back(_triangle_index<N>(xy.first - 1, xy.second + 1));
    }
    if (xy.second < N && xy.first < N - xy.second)
    {
        res.push_back(_triangle_index<N>(xy.first, xy.second + 1));
    }
    return res;
}

// ==============================================================

template <size_t N>
inline void calc_hex()
{
    constexpr size_t mat_side{N * (N + 4)};
    constexpr size_t band_width{N + 1};
    cout << "===========================\n"
         << "N = " << N << endl;

    cout << "generating coefficient matrix... ";
    auto mat_g{hex_network<N>()};
    cout << "finished" << endl;

    cout << "------- a, b ---------" << endl;

    cout << ">>> direct method --- ldl" << endl;

    Low_Band_Matrix<double, mat_side, band_width> l{};
    array<double, mat_side> d{};

    cout << "LDL factoring... ";
    ldl_factor(mat_g, l, d);
    cout << "finished ";

    cout << "generating b... ";
    array<double, mat_side> b{};
    size_t idx_b{_hex_index<N>(0, 2*N-3)};
    b[idx_b] = 1;
    cout << "finished ";

    cout << "calculating x... ";
    array<double, mat_side> x1{};
    back_sub(l, b, x1);
    back_sub(d, x1);
    back_sub(l, x1, true);
    cout << "finished" << endl;

    cout << ">>> U_ba = " << setprecision(16) << x1[idx_b] << endl;

    cout << ">>> iterative method --- conjugate gradient" << endl;

    cout << "solving x... ";
    array<double, mat_side> x2{};
    if (conj_grad(mat_g, b, x2, true) == 0)
    {
        cout << "finished" << endl;

        cout << ">>> U_ba = " << setprecision(16) << x2[idx_b] << endl;
    }
    else
    {
        cerr << "not converged" << endl;
    }

    cout << "=========== END ===========" << endl;
}

// the N-sides hexagonal network is like:
//                        X
//                      0   1
//                      2   3
//                    4   5   6
//                    7   8   9
//                  10  11  12  13
//                ... ... ... ... ...
//        ... ... ... ... ... ... ...
// N(N+1)-2  N(N+1)-1  ... ... ...  N(N+2)-3  N(N+2)-2
// N(N+2)-1  N(N+2) ... ... ... ... N(N+3)-2  N(N+3)-1
//      N(N+3)  ... ... ... ... ... ... N(N+4)-1
// then, point X is set as 0V
// so, the return matrix is of (N*(N+4)) * (N*(N+4))
template <size_t N>
inline Symm_Band_Matrix<double, N *(N + 4), N + 1> hex_network()
{
    constexpr size_t mat_length{N *(N + 4)};
    Symm_Band_Matrix<double, mat_length, N + 1> res{};
    vector<size_t> connections;

    for (size_t i = 0; i < mat_length; i++)
    {
        connections = _hex_connections<N>(i);
        res(i, i) += connections.size();
        for (const auto &j : connections)
        {
            res(j, j)++;
            res(j, i) = -1;
        }
    }

    // did not count point (-1, 0)'s effect on point (0, 0) & (1, 0)
    // their diagonal item should ++
    for (const auto &i : {_hex_index<N>(0, 0), _hex_index<N>(1, 0)})
    {
        res(i, i)++;
    }

    return res;
}

// the N-sides hexagonal network is like:
//                ------------------------------------------> x
//               /                  X
//            0 /                 0   1
//           1 /                  2   3
//          2 /                 4   5   6
//         3 /                  7   8   9
//        4 /                 10  11  12  13
//       . /                ... ... ... ... ...
//      . /         ... ... ... ... ... ... ...
//  2N-2 /   N(N+1)-2  N(N+1)-1  ... ... ...  N(N+2)-3  N(N+2)-2
// 2N-1 /    N(N+2)-1  N(N+2) ... ... ... ... N(N+3)-2  N(N+3)-1
//  2N /          N(N+3)  ... ... ... ... ... ... N(N+4)-1
//    V y
// also, point X is set as 0V, which is not indexed
// this function return the index of a point (x, y)
template <size_t N>
inline size_t _hex_index(size_t x, size_t y)
{
    size_t l{y / 2};
    return (l + y % 2) * (l + 2) + l + x;
}

// the N-sides hexagonal network is like:
//                ------------------------------------------> x
//               /                  X
//            0 /                 0   1
//           1 /                  2   3
//          2 /                 4   5   6
//         3 /                  7   8   9
//        4 /                 10  11  12  13
//       . /                ... ... ... ... ...
//      . /         ... ... ... ... ... ... ...
//  2N-2 /   N(N+1)-2  N(N+1)-1  ... ... ...  N(N+2)-3  N(N+2)-2
// 2N-1 /    N(N+2)-1  N(N+2) ... ... ... ... N(N+3)-2  N(N+3)-1
//  2N /          N(N+3)  ... ... ... ... ... ... N(N+4)-1
//    V y
// also, point X is set as 0V, which is not indexed
// this function return (x, y) of a point at index
template <size_t N>
inline pair<size_t, size_t> _hex_xy(size_t index)
{
    size_t x{index};
    size_t y{0};
    for (; x > 3+y; y+=2)
    {
        x -= 4+y;
    }
    if (x > y/2+1)
    {
        x -= y/2+2;
        y++;
    }
    if (y == 2*N)
    {
        x--;
    }
    return {x, y};
}

// return connections with higher indices
template <size_t N>
inline vector<size_t> _hex_connections(size_t index)
{
    auto xy{_hex_xy<N>(index)};
    vector<size_t> res{};
    if (index < N*(N+3)-1)
    {
        res.push_back(_hex_index<N>(xy.first, xy.second + 1));
    }

    if (xy.second == 2*N-1)
    {
        if (xy.first > 0)
        {
            res.push_back(_hex_index<N>(xy.first-1, xy.second + 1));
        }
    }
    else if (xy.second%2)
    {
        res.push_back(_hex_index<N>(xy.first+1, xy.second + 1));
    }
    return res;
}

// =============================================================

template <size_t N>
inline void calc_triangle_ac(double omega)
{
    constexpr size_t mat_side{N * (N + 3) / 2};
    constexpr size_t band_width{(N + 1)};
    cout << "===========================\n"
         << "N = " << N
         << "\nomega  = " << omega << endl;

    cout << "generating coefficient matrix... ";
    auto mat_g{triangle_ac_network<N>(omega)};
    auto mat_g_dagger{triangle_ac_network<N>(omega, true)};
    Hermite_Matrix<complex<double>, mat_side> mat_a;
    {
        auto temp {mat_g_dagger*mat_g};
        for (size_t i = 0; i < mat_side; i++)
            for (size_t j = 0; j <= i; j++)
            {
                if (temp(i, j) != std::conj(temp(j, i)))
                {
                    cerr << "not hermite: " << temp(i, j) << ' ' << temp(j, i);
                }
                mat_a(i, j) = temp(i, j);
            }
    }
    cout << "finished" << endl;

    cout << "------- a, b ---------" << endl;

    cout << ">>> direct method --- inv" << endl;

    cout << "inversing... ";
    auto inv_g {inv(mat_g)};
    cout << "finished ";

    cout << "generating b... ";
    array<complex<double>, mat_side> b{};
    size_t idx_b{_triangle_index<N>(0, N)};
    b[idx_b] = 1.;
    // b = mat_g_dagger*b;
    cout << "finished ";

    cout << "calculating x... ";
    array<complex<double>, mat_side> x1{
        inv_g*b};
    cout << "finished" << endl;

    cout << ">>> U_ba = " << setprecision(16) << x1[idx_b] << endl;

    cout << ">>> iterative method --- conjugate gradient" << endl;

    cout << "generating b... ";
    b = mat_g_dagger*b;
    cout << "finished ";

    cout << "solving x... ";
    array<complex<double>, mat_side> x2{};
    if (conj_grad(mat_a, b, x2, true) == 0)
    {
        cout << "finished" << endl;

        cout << ">>> U_ba = " << setprecision(16) << x2[idx_b] << endl;
    }
    else
    {
        cerr << "not converged" << endl;
    }

    cout << "=========== END ===========" << endl;
}

// the N-sides triangular network is like:
//     N*(N+3)/2-1
//     N*(N+3)/2-3 N*(N+3)/2-2
//     ...         ...         ...
//     2*N         2*N+1       ...   3*N-2
//     N           N+1         ... ...     2*N-1
//     X           0           ... ...     N-2      N-1
// then, point X is set as 0V
// so, the return matrix is of (N*(N+3)/2) * (N*(N+3)/2)
template <size_t N>
inline Symm_Band_Matrix<complex<double>, N *(N + 3) / 2, N + 1> triangle_ac_network(double omega, bool hermite)
{
    constexpr size_t mat_length{N * (N + 3) / 2};
    Symm_Band_Matrix<complex<double>, mat_length, N + 1> res{};
    map<size_t, AC_Type> connections;

    for (size_t i = 0; i < mat_length; i++)
    {
        connections = _triangle_ac_connections<N>(i);
        complex<double> g;
        for (const auto &j : connections)
        {
            switch (j.second)
            {
            case AC_Type::resister:
                g = 1;
                break;
            case AC_Type::inductance:
                g = hermite?1i/omega:-1i/omega;
                break;
            case AC_Type::capacitance:
                g = hermite?-1i*omega:1i*omega;
                break;
            default:
                break;
            }
            res(i, i) += g;
            res(j.first, j.first) += g;
            res(j.first, i) -= g;
        }
    }

    // did not count point (0, 0)'s effect on point (1, 0) & (0, 1)
    // their diagonal item should ++
    {
        auto i {_triangle_index<N>(1, 0)};
        res(i, i) += hermite?-1i*omega:1i*omega;
    }
    {
        auto i {_triangle_index<N>(0, 1)};
        res(i, i) += 1;
    }
    return res;
}

// return connections with higher indices
template <size_t N>
inline map<size_t, AC_Type> _triangle_ac_connections(size_t index)
{
    auto xy{_triangle_xy<N>(index)};
    map<size_t, AC_Type> res{};
    if (xy.first < N - xy.second)
    {
        res[_triangle_index<N>(xy.first + 1, xy.second)] = AC_Type::capacitance;
    }
    if (xy.first > 0)
    {
        res[_triangle_index<N>(xy.first - 1, xy.second + 1)] = AC_Type::inductance;
    }
    if (xy.second < N && xy.first < N - xy.second)
    {
        res[_triangle_index<N>(xy.first, xy.second + 1)] = AC_Type::resister;
    }
    return res;
}
