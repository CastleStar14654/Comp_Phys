// Basic definitions of Matrices for Computational Physics of 2020 Spring
// Author: Lin Xuchen, 29 Feb 2020
// Benefit a lot from Matrix.h by Bjarne Stroustrup

#ifndef MISC_MATRIX
#define MISC_MATRIX

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>

#include "Base_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Matrix with every item specified
template <typename T, size_t R, size_t C>
class Matrix : public Base_Matrix<T, R, C>
{
private:
    using Base_Matrix<T, R, C>::elem;
    using Base_Matrix<T, R, C>::data_ln;

public:
    using typename Base_Matrix<T, R, C>::size_type;

    explicit Matrix(T deft = T{}) : Base_Matrix<T, R, C>{R, deft} {}
    Matrix(const Matrix &mat) = default;
    Matrix(Matrix &&mat) = default;
    Matrix(const Base_Matrix<T, R, C> &mat)
        : Base_Matrix<T, R, C>{R, new T[R][C]}
    {
        for (size_type i = 0; i < R; i++)
            for (size_type j = 0; j < C; j++)
            {
                (*this)(i, j) = mat(i, j);
            }
    }
    Matrix(Base_Matrix<T, R, C> &&mat)
        : Base_Matrix<T, R, C>{R, new T[R][C]}
    {
        for (size_type i = 0; i < R; i++)
            for (size_type j = 0; j < C; j++)
            {
                (*this)(i, j) = std::move(mat(i, j));
            }
    }

    // Matrix<double> mat {
    //     {1, 2, 4},
    //     {6, 3, 2},
    //     {9, 5, -2},
    //     {3, 3, 4}
    // };
    // would be
    // [
    //     [1 2 4]
    //     [6 3 2]
    //     [9 5 -2]
    //     [3 3 4]
    // ]
    explicit Matrix(std::initializer_list<std::initializer_list<T>> ini)
        : Base_Matrix<T, R, C>{R, nullptr}
    {
        if (ini.size() != R)
        {
            throw std::invalid_argument((__FILE__ ":") + std::to_string(__LINE__) + ": wrong row number");
        }
        for (auto i = ini.begin(); i != ini.end(); i++)
            if (i->size() != C)
            {
                throw std::invalid_argument((__FILE__ ":") + std::to_string(__LINE__) + ": non-uniform column number");
            }
        elem = new T[R][C];
        auto it{ini.begin()};
        for (std::size_t i = 0; i < R; i++, it++)
        {
            std::move(it->begin(), it->end(), elem[i]);
        }
    }

    Matrix &operator=(const Matrix &mat) = default;
    Matrix &operator=(Matrix &&mat) = default;
    Matrix &operator=(const Base_Matrix<T, R, C> &mat)
    {
        Matrix<T, R, C> temp{mat};
        *this = std::move(temp);
        return *this;
    }
    Matrix &operator=(Base_Matrix<T, R, C> &&mat)
    {
        Matrix<T, R, C> temp{mat};
        *this = std::move(temp);
        return *this;
    }

    T &operator()(size_type row, size_type col) override { return elem[row][col]; }
    const T &operator()(size_type row, size_type col) const override { return elem[row][col]; }
};

// -------------------------------------------------------------------------

template <typename T, size_t R, size_t N, size_t C>
Matrix<T, R, C> operator*(const Base_Matrix<T, R, N> &a, const Base_Matrix<T, N, C> &b)
{
    Matrix<T, R, C> res{};
    // This should be ten times faster.
    for (std::size_t i = 0; i < R; i++)
        for (std::size_t k = 0; k < N; k++)
        {
            T r = a(i, k);
            for (std::size_t j = 0; j < C; j++)
            {
                res(i, j) += r * b(k, j);
            }
        }
    return res;
}

template <typename T, size_t R, size_t C, size_t A, size_t B>
Matrix<T, R, C> operator*(const Column<T, R, A> &a, const Row<T, B, C> &b)
{
    Matrix<T, R, C> res{};
    for (std::size_t i = 0; i < R; i++)
        for (std::size_t j = 0; j < C; j++)
        {
            res(i, j) = a[i] * b[j];
        }
    return res;
}

template <typename T, size_t R, size_t N, size_t C>
Matrix<T, 1, C> operator*(const Row<T, R, N> &a, const Base_Matrix<T, N, C> &b)
{
    Matrix<T, 1, C> res{};
    for (std::size_t i = 0; i < C; i++)
    {
        res(0, i) = a * b.column(i);
    }
    return res;
}

template <typename T, size_t R, size_t N, size_t C>
Matrix<T, R, 1> operator*(const Base_Matrix<T, R, N> &a, const Column<T, N, C> &b)
{
    Matrix<T, R, 1> res{};
    for (std::size_t i = 0; i < R; i++)
    {
        res(i, 0) = a.row(i) * b;
    }
    return res;
}

} // namespace Misc

#endif // MISC_MATRIX
