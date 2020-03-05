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
template <typename T>
class Matrix : public Base_Matrix<T>
{
private:
    using Base_Matrix<T>::rs;
    using Base_Matrix<T>::cs;
    using Base_Matrix<T>::elem;
    using Base_Matrix<T>::data_sz;

public:
    using typename Base_Matrix<T>::size_type;

    explicit Matrix(size_type r, size_type c, T deft = T{});
    Matrix(const Matrix &mat) = default;
    Matrix(Matrix &&mat) = default;
    Matrix(const Base_Matrix<T> &mat);
    Matrix(Base_Matrix<T> &&mat);
    explicit Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Matrix &operator=(const Matrix &mat) = default;
    Matrix &operator=(Matrix &&mat) = default;
    Matrix &operator=(const Base_Matrix<T> &mat);
    Matrix &operator=(Base_Matrix<T> &&mat);

    T &operator()(size_type row, size_type col) override { return elem[row * cs + col]; }
    const T &operator()(size_type row, size_type col) const override { return elem[row * cs + col]; }
};

// -------------------------------------------------------------------------

template <typename T>
Matrix<T> operator*(const Base_Matrix<T> &a, const Base_Matrix<T> &b)
{
    if (a.cols() != b.rows())
    {
        throw std::domain_error("Matrix::operator*(): invalid shapes.");
    }

    Matrix<T> res{a.rows(), b.cols()};
    // This should be ten times faster.
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t k = 0; k < a.cols(); k++)
        {
            T r = a(i, k);
            for (std::size_t j = 0; j < res.cols(); j++)
            {
                res(i, j) += r * b(k, j);
            }
        }
    return res;
}

template <typename T>
Matrix<T> operator*(const Column<T>& a, const Row<T>& b)
{
    Matrix<T> res {a.size(), b.size()};
    for (std::size_t i = 0; i < a.size(); i++)
        for (std::size_t j = 0; j < b.size(); j++)
        {
            res(i, j) = a[i]*b[j];
        }
    return res;
}

// =====================Matrix===============================

template <typename T>
Matrix<T>::Matrix(size_type r, size_type c, T deft)
    : Base_Matrix<T>{r, c, r * c, deft}
{
}

template <typename T>
Matrix<T>::Matrix(const Base_Matrix<T> &mat)
    : Base_Matrix<T>{mat.rows(), mat.cols(), mat.size(), new T[mat.size()]}
{
    for (size_type i = 0; i < rs; i++)
        for (size_type j = 0; j < cs; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
}

template <typename T>
Matrix<T>::Matrix(Base_Matrix<T> &&mat)
    : Base_Matrix<T>{mat.rows(), mat.cols(), mat.size(), new T[mat.size()]}
{
    for (size_type i = 0; i < rs; i++)
        for (size_type j = 0; j < cs; j++)
        {
            this->operator()(i, j) = std::move(mat(i, j));
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
template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Matrix<T>{ini.size(), ini.begin()->size(),
                     ini.size() * ini.begin()->size(), nullptr}
{
    for (auto i = ini.begin(); i != ini.end(); i++)
        if (i->size() != cs)
        {
            throw std::invalid_argument("Matrix::Matrix: non-uniform column number");
        }
    elem = new T[data_sz];
    for (std::size_t i = 0; i < rs; i++)
    {
        std::move((ini.begin() + i)->begin(), (ini.begin() + i)->end(), &(*this)(i, 0));
    }
}

// ------------------- Matrix operator= -----------------------------

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Base_Matrix<T> &mat)
{
    if (this->shape() != mat.shape())
    {
        throw std::runtime_error("Matrix assignment: non-uniform shape.");
    }

    Matrix<T> temp{mat};
    *this = std::move(temp);
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Base_Matrix<T> &&mat)
{
    if (this->shape() != mat.shape())
    {
        throw std::runtime_error("Matrix assignment: different shapes.");
    }

    Matrix<T> temp{mat};
    *this = std::move(temp);
    return *this;
}

} // namespace Misc

#endif // MISC_MATRIX
