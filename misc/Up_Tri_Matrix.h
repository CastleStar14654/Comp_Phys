#ifndef MISC_UP_TRI_MATRIX
#define MISC_UP_TRI_MATRIX

#include <iterator>
#include <initializer_list>
#include <iostream>

#include "Base_Matrix.h"
// #include "Matrix.h"
#include "Base_Tri_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Upper Triangular Matrix; only half of the items are stored
template <typename T, size_t N>
class Up_Tri_Matrix : public Base_Tri_Matrix<T, N>
{
private:
    // using Base_Tri_Matrix<T, N>::rs;
    // using Base_Tri_Matrix<T, N>::cs;
    using Base_Tri_Matrix<T, N>::elem;
    using Base_Tri_Matrix<T, N>::data_ln;

public:
    using typename Base_Tri_Matrix<T, N>::size_type;
    using Base_Tri_Matrix<T, N>::Base_Tri_Matrix;

    Up_Tri_Matrix(const Up_Tri_Matrix &mat) = default;
    Up_Tri_Matrix(Up_Tri_Matrix &&mat) = default;
    Up_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Up_Tri_Matrix &operator=(const Up_Tri_Matrix &mat) = default;
    Up_Tri_Matrix &operator=(Up_Tri_Matrix &&mat) = default;

    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;
};

// -------------------------------------------------------------------------

template <typename T, size_t N>
Up_Tri_Matrix<T, N> operator*(const Up_Tri_Matrix<T, N> &a, const Up_Tri_Matrix<T, N> &b)
{
    Up_Tri_Matrix<T, N> res{};
    for (std::size_t j = 0; j < N; j++)
        for (std::size_t k = 0; k <= j; k++)
        {
            T temp {b(k, j)};
            for (std::size_t i = k; i <= j; i++)
            {
                res(i, j) += a(i, k) * temp;
            }
        }
    return res;
}

template <typename T, size_t N>
Up_Tri_Matrix<T, N> operator*(const Up_Tri_Matrix<T, N> &a, const Diag_Matrix<T, N> &b)
{
    Up_Tri_Matrix<T, N> res{};
    for (std::size_t j = 0; j < N; j++)
        for (std::size_t i = 0; i <= j; i++)
        {
            res(i, j) += a(i, j) * b(j);
        }
    return res;
}

template <typename T, size_t N>
Up_Tri_Matrix<T, N> operator*(const Diag_Matrix<T, N> &a, const Up_Tri_Matrix<T, N> &b)
{
    Up_Tri_Matrix<T, N> res{};
    for (std::size_t j = 0; j < N; j++)
        for (std::size_t i = 0; i <= j; i++)
        {
            res(i, j) += a(i) * b(i, j);
        }
    return res;
}

// ========================== Up_Tri_Matrix =================================

template <typename T, size_t N>
Up_Tri_Matrix<T, N>::Up_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Tri_Matrix<T, N>{}
{
    int count{1};
    for (auto i = std::rbegin(ini); i != std::rend(ini); i++)
    {
        if (i->size() != count)
        {
            throw std::invalid_argument("Matrix::Matrix: wrong column number");
        }
        ++count;
    }
    // elem = new T[data_sz];

    auto ini_r{ini.begin()};
    for (std::size_t r = 0; r < N; r++)
    {
        auto ini_c{ini_r->begin()};
        for (std::size_t c = r; c < N; c++)
        {
            (*this)(r, c) = *ini_c;
            ++ini_c;
        }
        ++ini_r;
    }
}

// -------------------- Up_Tri_Matrix: row & column ----------------------------

template <typename T, size_t N>
T &Up_Tri_Matrix<T, N>::operator()(size_type row, size_type col)
{
    if (col < row)
    {
        throw std::out_of_range("Up_Tri_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return Base_Tri_Matrix<T, N>::operator()(col, row);
    }
}

template <typename T, size_t N>
const T &Up_Tri_Matrix<T, N>::operator()(size_type row, size_type col) const
{
    if (col < row)
    {
        return Base_Matrix<T, N, N>::zero;
    }
    else
    {
        return Base_Tri_Matrix<T, N>::operator()(col, row);
    }
}

} // namespace Misc

#endif // MISC_UP_TRI_MATRIX
