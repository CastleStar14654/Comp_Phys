#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <string>
#include <iomanip>

#include "../misc/svd.h"

using namespace std;
using namespace Misc;

/* Decentralize `data` in situ.  Return the data centre
 */
template <size_t R, size_t C>
array<double, C> decentralize(Matrix<double, R, C> &data);

/* Load matrix from file.
 * Shape of the matrix should be given explicitly.
 */
template <size_t R, size_t C>
Matrix<double, R, C> load(string file_name);

/* Save matrix to file as plain text.
 * Data will be output in `fixed` form with precision set to 4.
 */
template <size_t R, size_t C>
void save(string file_name, const Matrix<double, R, C> &data);

/* Save array to file as plain text.
 * Data will be output in `fixed` form with precision set to 4.
 */
template <size_t N>
void save(string file_name, const array<double, N> &data);

int main()
{
    auto data{load<10000, 6>("data.txt")};
    auto centre{decentralize(data)};
    save("output/centre.txt", centre);

    Matrix<double, 10000, 6> u_mat{};
    Matrix<double, 6, 6> v_mat{};
    auto singular_values{svd(data, &u_mat, &v_mat)};

    save("output/singular_values.txt", singular_values);
    save("output/v_mat.txt", v_mat);
    save("output/u_mat.txt", u_mat);

    u_mat = data * v_mat;
    save("output/answer_svd.txt", u_mat);
}

template <size_t R, size_t C>
array<double, C> decentralize(Matrix<double, R, C> &data)
{
    array<double, C> res;
    for (size_t i = 0; i < R; i++)
        for (size_t j = 0; j < C; j++)
        {
            res[j] += data(i, j);
        }
    for (size_t j = 0; j < C; j++)
    {
        res[j] *= 1. / R;
    }
    for (size_t i = 0; i < R; i++)
        for (size_t j = 0; j < C; j++)
        {
            data(i, j) -= res[j];
        }
    return res;
}

template <size_t R, size_t C>
Matrix<double, R, C> load(string file_name)
{
    ifstream ifs{file_name};
    Matrix<double, R, C> res{};
    for (size_t i = 0; i < R; i++)
        for (size_t j = 0; j < C; j++)
        {
            ifs >> (res(i, j));
        }
    return res;
}

template <size_t R, size_t C>
void save(string file_name, const Matrix<double, R, C> &data)
{
    ofstream ofs{file_name};
    ofs << fixed << setprecision(4);
    for (size_t i = 0; i < R; i++)
    {
        if (C)
        {
            ofs << data(i, 0);
        }
        for (size_t j = 1; j < C; j++)
        {
            ofs << " " << data(i, j);
        }
        ofs << '\n';
    }
}

template <size_t N>
void save(string file_name, const array<double, N> &data)
{
    ofstream ofs{file_name};
    ofs << fixed << setprecision(4);
    if (N)
    {
        ofs << data[0];
    }
    for (size_t j = 1; j < N; j++)
    {
        ofs << " " << data[j];
    }
    ofs << '\n';
}
