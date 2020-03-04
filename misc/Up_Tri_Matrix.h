#ifndef MISC_LOW_TRI_MATRIX
#define MISC_LOW_TRI_MATRIX

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
template <typename T>
class Up_Tri_Matrix : public Base_Tri_Matrix<T>
{
private:
    using Base_Tri_Matrix<T>::rs;
    using Base_Tri_Matrix<T>::cs;
    using Base_Tri_Matrix<T>::elem;
    using Base_Tri_Matrix<T>::data_sz;

public:
    using typename Base_Tri_Matrix<T>::size_type;
    using Base_Tri_Matrix<T>::Base_Tri_Matrix;

    Up_Tri_Matrix(const Up_Tri_Matrix &mat)=default;
    Up_Tri_Matrix(Up_Tri_Matrix &&mat)=default;
    Up_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Up_Tri_Matrix &operator=(const Up_Tri_Matrix &mat);
    Up_Tri_Matrix &operator=(Up_Tri_Matrix &&mat);

    Row<T> row(size_type pos) const override;
    Column<T> column(size_type pos) const override;
    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;
};

// -------------------------------------------------------------------------

template <typename T>
Up_Tri_Matrix<T> operator*(const Up_Tri_Matrix<T> &a, const Up_Tri_Matrix<T> &b)
{
    if (a.shape() != b.shape())
    {
        throw std::domain_error("Up_Tri_Matrix::operator*(): invalid shapes.");
    }

    Up_Tri_Matrix<T> res{a.rows()};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = i; j < res.cols(); j++)
            for (std::size_t k = i; k <= j; k++)
            {
                res(i, j) += a(i, k) * b(k, j);
            }
    return res;
}

template <typename T>
Up_Tri_Matrix<T> operator*(const Up_Tri_Matrix<T> &a, const Diag_Matrix<T> &b)
{
    return a * (*reinterpret_cast<const Up_Tri_Matrix<T> *>(&b));
}

template <typename T>
Up_Tri_Matrix<T> operator*(const Diag_Matrix<T> &a, const Up_Tri_Matrix<T> &b)
{
    return (*reinterpret_cast<const Up_Tri_Matrix<T> *>(&a)) * b;
}


// ========================== Up_Tri_Matrix =================================

template <typename T>
Up_Tri_Matrix<T>::Up_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Tri_Matrix<T>{ini.size()}
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

    auto ini_r {ini.begin()};
    for (std::size_t r = 0; r < rs; r++)
    {
        auto ini_c {ini_r->begin()};
        for (std::size_t c = r; c < cs; c++)
        {
            (*this)(r, c) = *ini_c;
            ++ini_c;
        }
        ++ini_r;
    }
}

// -------------------- Up_Tri_Matrix: row & column ----------------------------

template <typename T>
Row<T> Up_Tri_Matrix<T>::row(size_type pos) const
{
    Row<T> res(pos, Base_Tri_Matrix<T>::zero);
    Base_Tri_Matrix<T>::insert_column(pos, std::back_inserter(res));
    return res;
}

template <typename T>
Column<T> Up_Tri_Matrix<T>::column(size_type pos) const
{
    Column<T> res{Base_Tri_Matrix<T>::row(pos)};
    res.resize(rs, Base_Tri_Matrix<T>::zero);
    return res;
}

template <typename T>
T &Up_Tri_Matrix<T>::operator()(size_type row, size_type col)
{
    if (col < row)
    {
        throw std::out_of_range("Up_Tri_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return Base_Tri_Matrix<T>::operator()(col, row);
    }
}

template <typename T>
const T &Up_Tri_Matrix<T>::operator()(size_type row, size_type col) const
{
    if (col < row)
    {
        return Base_Matrix<T>::zero;
    }
    else
    {
        return Base_Tri_Matrix<T>::operator()(col, row);
    }
}

} // namespace Misc

#endif // MISC_LOW_TRI_MATRIX
