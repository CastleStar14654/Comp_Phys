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
    Matrix(const Base_Matrix<T> &mat);
    Matrix(const Matrix &mat);
    Matrix(Base_Matrix<T> &&mat);
    Matrix(Matrix &&mat);
    explicit Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Matrix &operator=(const Base_Matrix<T> &mat);
    Matrix &operator=(const Matrix &mat);
    Matrix &operator=(Base_Matrix<T> &&mat);
    Matrix &operator=(Matrix &&mat);

    Row<T> row(size_type pos) const override;
    Column<T> column(size_type pos) const override;
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
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = 0; j < res.cols(); j++)
            for (std::size_t k = 0; k < a.cols(); k++)
            {
                res(i, j) += a(i, k) * b(k, j);
            }
    return res;
}

// =====================Matrix===============================

template <typename T>
Matrix<T>::Matrix(size_type r, size_type c, T deft)
    : Base_Matrix<T>{r, c, r * c, new T[r * c]{}}
{
    if (deft != T{})
        for (std::size_t i = 0; i < rs * cs; i++)
        {
            elem[i] = deft;
        }
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
Matrix<T>::Matrix(const Matrix &mat)
    : Base_Matrix<T>{mat.rs, mat.cs, mat.data_sz, new T[mat.data_sz]}
{
    std::copy(mat.elem, mat.elem + data_sz, elem);
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

template <typename T>
Matrix<T>::Matrix(Matrix &&mat)
    : Base_Matrix<T>{mat.rs, mat.cs, mat.data_sz, new T[mat.data_sz]}
{
    mat.elem = nullptr;
}

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
    Matrix<T> temp{mat};
    *this = std::move(temp);
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix &mat)
{
    T *temp = new T[mat.size()];
    std::copy(mat.elem, mat.elem + rs * cs, temp);
    delete[] elem;
    elem = temp;
    rs = mat.rs;
    cs = mat.cs;
    data_sz = mat.data_sz;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Base_Matrix<T> &&mat)
{
    Matrix<T> temp{mat};
    *this = std::move(temp);
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix &&mat)
{
    delete[] elem;
    elem = mat.elem;
    mat.elem = nullptr;
    rs = mat.rs;
    cs = mat.cs;
    data_sz = mat.data_sz;
    return *this;
}

// ------------------------- Matrix row() & column() -----------------------------

template <typename T>
Row<T> Matrix<T>::row(size_type pos) const
{
    Row<T> res(&elem[pos*cs], &elem[(pos+1)*cs]);
    // for (std::size_t c = 0; c < cs; c++)
    // {
    //     res[c] = this->operator()(pos, c);
    // }
    return res;
}

template <typename T>
Column<T> Matrix<T>::column(size_type pos) const
{
    Column<T> res(rs);
    for (std::size_t r = 0; r < rs; r++)
    {
        res[r] = this->operator()(r, pos);
    }
    return res;
}

// ===============================================================================

} // namespace Misc

#endif // MISC_MATRIX
