#ifndef MISC_BASE_TRI_MATRIX
#define MISC_BASE_TRI_MATRIX

#include <stdexcept>
#include <initializer_list>
#include <utility>

#include "Base_Matrix.h"
#include "Matrix.h"
#include "Diag_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Base Triangular Matrix; only half of the items are stored
// being inherited by Symm_Matrix, Up_Tri_Matrix, Low_Tri_Matrix
template <typename T>
class Base_Tri_Matrix : public Base_Matrix<T>
{
public:
    using typename Base_Matrix<T>::size_type;

    explicit Base_Tri_Matrix(size_type n, T deft = T{});
    Base_Tri_Matrix(const Base_Tri_Matrix &mat) = default;
    Base_Tri_Matrix(Base_Tri_Matrix &&mat) = default;
    Base_Tri_Matrix(const Diag_Matrix<T> &mat);
    Base_Tri_Matrix(Diag_Matrix<T> &&mat);
    Base_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Base_Tri_Matrix &operator=(const Base_Tri_Matrix &mat) = default;
    Base_Tri_Matrix &operator=(Base_Tri_Matrix &&mat) = default;
    Base_Tri_Matrix &operator=(const Diag_Matrix<T> &mat);
    Base_Tri_Matrix &operator=(Diag_Matrix<T> &&mat);

protected:
    using Base_Matrix<T>::rs;
    using Base_Matrix<T>::cs;
    using Base_Matrix<T>::elem;
    using Base_Matrix<T>::data_sz;

    // protected, only return PART of the column/row
    virtual Row<T> row(size_type pos) const = 0;
    // virtual Column<T> column(size_type pos) const override;
    virtual T &operator()(size_type row, size_type col) override;
    virtual const T &operator()(size_type row, size_type col) const override;

    template <typename OutputIt>
    OutputIt insert_column(size_type pos, OutputIt out) const;
};

// ========================== Base_Tri_Matrix =================================

template <typename T>
Base_Tri_Matrix<T>::Base_Tri_Matrix(size_type n, T deft)
    : Base_Matrix<T>{n, n, n * (n + 1) / 2, deft}
{
}

template <typename T>
Base_Tri_Matrix<T>::Base_Tri_Matrix(const Diag_Matrix<T> &mat)
    : Base_Matrix<T>{mat.rows(), mat.cols(),
                     mat.rows() * (mat.rows() + 1) / 2,
                     T{}}
{
    for (std::size_t i = 0; i < rs; i++)
    {
        (*this)(i, i) = mat(i, i);
    }
}

template <typename T>
Base_Tri_Matrix<T>::Base_Tri_Matrix(Diag_Matrix<T> &&mat)
    : Base_Matrix<T>{mat.rows(), mat.cols(),
                     mat.rows() * (mat.rows() + 1) / 2,
                     T{}}
{
    for (std::size_t i = 0; i < rs; i++)
    {
        (*this)(i, i) = std::move(mat(i, i));
    }
}

template <typename T>
Base_Tri_Matrix<T>::Base_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Matrix<T>{ini.size(), ini.size(),
                     ini.size() * (ini.size() + 1) / 2, nullptr}
{
    int count{1};
    for (auto i = ini.begin(); i != ini.end(); i++)
    {
        if (i->size() != count)
        {
            throw std::invalid_argument("Matrix::Matrix: wrong column number");
        }
        ++count;
    }
    elem = new T[data_sz];
    for (std::size_t i = 0; i < rs; i++)
    {
        std::move((ini.begin() + i)->begin(), (ini.begin() + i)->end(), &(*this)(i, 0));
    }
}

// -------------------- Base_Tri_Matrix: operator= ----------------------------

template <typename T>
Base_Tri_Matrix<T> &Base_Tri_Matrix<T>::operator=(const Diag_Matrix<T> &mat)
{
    if (this->shape() != mat.shape())
    {
        throw std::runtime_error("Matrix assignment: non-uniform shape.");
    }

    T* temp = new T[data_sz]{};
    delete[] elem;
    elem = temp;
    for (std::size_t i = 0; i < rs; i++)
    {
        (*this)(i, i) = mat(i, i);
    }
    return *this;
}

template <typename T>
Base_Tri_Matrix<T> &Base_Tri_Matrix<T>::operator=(Diag_Matrix<T> &&mat)
{
    if (this->shape() != mat.shape())
    {
        throw std::runtime_error("Matrix assignment: non-uniform shape.");
    }

    T* temp = new T[data_sz]{};
    delete[] elem;
    elem = temp;
    for (std::size_t i = 0; i < rs; i++)
    {
        (*this)(i, i) = std::move(mat(i, i));
    }
    return *this;
}

// -------------------- Base_Tri_Matrix: row & column ----------------------------
// protected, only return PART of the column/row

template <typename T>
Row<T> Base_Tri_Matrix<T>::row(size_type pos) const
{
    Row<T> res(&elem[pos * (pos + 1) / 2], &elem[(pos + 1) * (pos + 2) / 2]);
    return res;
}

template <typename T>
T &Base_Tri_Matrix<T>::operator()(size_type row, size_type col)
{
    if (row < col)
    {
        throw std::out_of_range("Base_Tri_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return elem[row * (row + 1) / 2 + col];
    }
}

template <typename T>
const T &Base_Tri_Matrix<T>::operator()(size_type row, size_type col) const
{
    if (row < col)
    {
        throw std::out_of_range("Base_Tri_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return elem[row * (row + 1) / 2 + col];
    }
}

template <typename T>
template <typename OutputIt>
OutputIt Base_Tri_Matrix<T>::insert_column(size_type pos, OutputIt out) const
{
    for (std::size_t i = pos; i < rs; i++)
    {
        *out = (*this)(i, pos);
        ++out;
    }
    return out;
}

} // namespace Misc

#endif // MISC_BASE_TRI_MATRIX
