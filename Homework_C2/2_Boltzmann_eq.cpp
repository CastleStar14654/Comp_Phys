#include <iostream>
#include <fstream>
#include <array>
#include <cmath>

using namespace std;

/*beg:grid*/
constexpr double MAX_P{20.}; // - MAX_P <= p <= MAX_P
constexpr double MAX_Z{1.};  // - MAX_Z <= z <= MAX_Z
constexpr size_t N_P{200};
constexpr size_t N_Z{200};
constexpr double H_P{2. * MAX_P / N_P};
constexpr double INV_H_P{1. / H_P};
constexpr double H_Z{2. * MAX_Z / N_Z};
constexpr double INV_H_Z{1. / H_Z};
constexpr double TAU{0.0004};
constexpr double COLL_NU{100.};
constexpr double G{10.};

array<double, N_P> p_s;
array<double, N_Z> z_s;
/*end:grid*/

/*beg:step_dec*/

/* get data from `from`, write next frame into `to`
 * f0: Boltzmann distribution
 * g:  gravity
 */
void step(const array<array<double, N_P>, N_Z> &from,
          array<array<double, N_P>, N_Z> &to,
          const array<array<double, N_P>, N_Z> &f0, double g = 0.);
/*end:step_dec*/

void output(const array<array<double, N_P>, N_Z> &from,
            ofstream &ofs);

/*beg:equi_dec*/

/* generate the equilibrium data of T
 * from, to: range of Z index
 */
void equilibrium(array<array<double, N_P>, N_Z> &data, double T,
                 double g = 0., size_t from = 0, size_t to = N_Z);
/*end:equi_dec*/

double H_status(const array<array<double, N_P>, N_Z> &data);

/*beg:init*/
int main()
{
    for (size_t j = 0; j < N_P; j++)
        p_s[j] = -MAX_P + H_P * (j + .5);
    for (size_t i = 0; i < N_Z; i++)
        z_s[i] = -MAX_Z + H_Z * (i + .5);
    /*end:init*/
    auto p_equi_array{new array<array<double, N_P>, N_Z>};
    equilibrium(*p_equi_array, 30.);

    /*problem_1*/
    auto current{new array<array<double, N_P>, N_Z>};
    equilibrium(*current, 10., 0., 0, N_Z / 2);
    equilibrium(*current, 20., 0., N_Z / 2, N_Z);
    {
        cout << "PROBLEM 1" << endl;
        ofstream ofs{"output/1.bin", ios_base::binary};
        output(*current, ofs);
    }

    /*problem_2*/
    auto next_frame{new array<array<double, N_P>, N_Z>};
    {
        cout << "PROBLEM 2" << endl;
        ofstream ofs_phase{"output/2_phase.bin", ios_base::binary};
        ofstream ofs_H{"output/2_H.bin", ios_base::binary};
        for (size_t i = 0; i < 201; i++)
        {
            step(*current, *next_frame, *p_equi_array);
            if (!(i % 10) && i < 150)
            {
                cout << i << ' ';
                output(*current, ofs_phase);
            }
            double H{H_status(*current)};
            ofs_H.write(static_cast<char *>((void *)(&H)), sizeof(double));
            swap(current, next_frame);
        }
    }
    /*end:problem_2*/

    /*problem_3*/
    equilibrium(*p_equi_array, 30., G);
    equilibrium(*current, 10., G, 0, N_Z / 2);
    equilibrium(*current, 20., G, N_Z / 2, N_Z);
    {
        cout << "\nPROBLEM 3" << endl;
        ofstream ofs_phase{"output/3_phase.bin", ios_base::binary};
        ofstream ofs_H{"output/3_H.bin", ios_base::binary};
        for (size_t i = 0; i < 201; i++)
        {
            step(*current, *next_frame, *p_equi_array, G);
            if (!(i % 10) && i < 150)
            {
                cout << i << ' ';
                output(*current, ofs_phase);
            }
            double H{H_status(*current)};
            ofs_H.write(static_cast<char *>((void *)(&H)), sizeof(double));
            swap(current, next_frame);
        }
    }

    delete p_equi_array;
    delete current;
    delete next_frame;
}

/*beg:step*/
void step(const array<array<double, N_P>, N_Z> &from,
          array<array<double, N_P>, N_Z> &to,
          const array<array<double, N_P>, N_Z> &f0, double g)
{
    g /= H_P;
    for (size_t i = 0; i < N_Z; i++)
        for (size_t j = 0; j < N_P; j++)
        {
            to[i][j] = TAU * (COLL_NU * f0[i][j] - (COLL_NU + INV_H_Z * abs(p_s[j])) * from[i][j]);
        }

    for (size_t i = 0; i < N_Z; i++)
        for (size_t j = 0; j < N_P - 1; j++)
        {
            to[i][j] += TAU * g * (from[i][j + 1] - from[i][j]);
        }

    for (size_t j = 0; j < N_P / 2; j++)
        to[0][j] -= TAU * INV_H_Z * p_s[j] * from[1][j];
    for (size_t j = N_P / 2; j < N_P; j++)
        to[0][j] += TAU * INV_H_Z * p_s[j] * from[0][N_P - 1 - j];
    for (size_t i = 1; i < N_Z - 1; i++)
    {
        for (size_t j = 0; j < N_P / 2; j++)
            to[i][j] -= TAU * INV_H_Z * p_s[j] * from[i + 1][j];
        for (size_t j = N_P / 2; j < N_P; j++)
            to[i][j] += TAU * INV_H_Z * p_s[j] * from[i - 1][j];
    }
    for (size_t j = 0; j < N_P / 2; j++)
        to[N_Z - 1][j] -= TAU * INV_H_Z * p_s[j] * from[N_Z - 1][N_P - 1 - j];
    for (size_t j = N_P / 2; j < N_P; j++)
        to[N_Z - 1][j] += TAU * +INV_H_Z * p_s[j] * from[N_Z - 2][j];

    for (size_t i = 0; i < N_Z; i++)
        for (size_t j = 0; j < N_P; j++)
        {
            to[i][j] += from[i][j];
        }
}
/*end:step*/

void output(const array<array<double, N_P>, N_Z> &from,
            ofstream &ofs)
{
    ofs.write(static_cast<char *>((void *)(&from[0][0])), N_P * N_Z * sizeof(double));
}

/*beg:equi*/
void equilibrium(array<array<double, N_P>, N_Z> &data, double T,
                 double g, size_t from, size_t to)
{
    double half_beta{.5 / T};
    double coef{1. / sqrt(6.283185307179586 * T)};
    if (g)
    {
        g /= T;
        if (from != 0 || to != N_Z)
        {
            double half_g{.5 * g};
            coef *= half_g / sinh(half_g);
            coef *= from != 0 ? exp(half_g) : exp(-half_g);
        }
        else
        {
            coef *= g / sinh(g);
        }
    }
    for (size_t j = 0; j < N_P; j++)
        for (size_t i = from; i < to; i++)
        {
            data[i][j] = coef * exp(-p_s[j] * p_s[j] * half_beta - g * z_s[i]);
        }
}
/*end:equi*/

/*beg:H*/
double H_status(const array<array<double, N_P>, N_Z> &data)
{
    double res{};
    for (size_t i = 0; i < N_Z; i++)
        for (size_t j = 0; j < N_P; j++)
        {
            res += data[i][j] * log(data[i][j]);
        }
    res *= H_P * H_Z;
    return res;
}
/*end:H*/
