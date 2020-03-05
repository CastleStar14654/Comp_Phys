// Basic definitions of Base_Matrix for Computational Physics of 2020 Spring
// Author: Lin Xuchen, 1 Mar 2020
// Benefit a lot from Matrix.h by Bjarne Stroustrup

#ifndef MISC_BASE_MATRIX
#define MISC_BASE_MATRIX

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>

#include "Matrix.h"

// the namespace miscellany
namespace Misc
{
template <typename T>
class Base_Matrix;

template <typename T>
class Matrix;


// Row and Column for convenience
// A deputy for columns & rows in a matrix
// NO data is copied
// template <typename T>
// class Base_Vector
// {
// public:
//     using size_type = std::size_t;
//     using value_type = T;

//     size_type size() const {return m.cols();}

// protected:
//     Base_Vector(Base_Matrix<T>& mat) : m{mat} {}

//     Base_Matrix<T>& m;
// };

template <typename T>
class Row
{
public:
    using size_type = std::size_t;
    using value_type = T;

    size_type size() const {return m.cols();}

    Row(Base_Matrix<T>& mat, size_type row) : m{mat}, r{row} {}

    T& operator[](size_type col) {return m(r, col);}
    const T& operator[](size_type col) const {
        return const_cast<const Base_Matrix<T>&>(m)(r, col);
    }
    Row& operator=(const Base_Matrix<T>& mat_row);
    Row& operator=(Base_Matrix<T>&& mat_row);
private:
    Base_Matrix<T>& m;
    size_type r;
};

// int a[3]{1, 2, 3};
// int(&b)[3] = a;
// void f(const int(&)[]);
// void f(const int(&a)[]) {
//     a;
// }

template <typename T>
class Column
{
public:
    using size_type = std::size_t;
    using value_type = T;

    size_type size() const {return m.rows();}

    Column(Base_Matrix<T>& mat, size_type col) : m{mat}, c{col} {}

    T& operator[](size_type row) {return m(row, c);}
    const T& operator[](size_type row) const {
        return const_cast<const Base_Matrix<T>&>(m)(row, c);
    }
    Column& operator=(const Base_Matrix<T>& mat_col);
    Column& operator=(Base_Matrix<T>&& mat_col);
private:
    Base_Matrix<T>& m;
    size_type c;
};

template <typename T>
T operator*(const Row<T>& a, const Column<T>& b)
{
    if (a.size() != b.size())
    {
        throw std::domain_error("Row*Column: invalid shapes.");
    }

    T res {};
    for (std::size_t i = 0; i < a.size(); i++)
    {
        res += a[i]*b[i];
    }
    return res;
}

template <typename T>
Matrix<T> operator*(const Row<T>& a, const Base_Matrix<T>& b)
{
    if (a.size() != b.rows())
    {
        throw std::domain_error("Row*Base_Matrix: invalid shapes.");
    }

    Matrix<T> res {1, a.size()};
    for (std::size_t i = 0; i < a.size(); i++)
    {
        res(0, i) = a*b.column(i);
    }
    return res;
}

template <typename T>
Matrix<T> operator*(const Base_Matrix<T>& a, const Column<T>& b)
{
    if (a.cols() != b.size())
    {
        throw std::domain_error("Base_Matrix*Column: invalid shapes.");
    }

    Matrix<T> res {b.size(), 1};
    for (std::size_t i = 0; i < a.size(); i++)
    {
        res(i, 0) = a.row(i)*b;
    }
    return res;
}



// ================================================================================

// Never initiate a Base_Matrix
template <typename T>
class Base_Matrix
{
public:
    using size_type = std::size_t;
    using value_type = T;

protected:
    const size_type rs;
    const size_type cs;
    const size_type data_sz;
    T *elem;
    constexpr static T zero{};

public:
    Base_Matrix(size_type r, size_type c, size_type dt_sz, T *elm)
        : rs{r}, cs{c}, data_sz{dt_sz}, elem{elm} {}
    Base_Matrix(size_type r, size_type c, size_type dt_sz, T deft=T{});
    Base_Matrix(const Base_Matrix &mat);
    Base_Matrix(Base_Matrix &&mat);

    virtual ~Base_Matrix() { delete[] elem; }

    Base_Matrix &operator=(const Base_Matrix &mat);
    Base_Matrix &operator=(Base_Matrix &&mat);

    Row<T> row(size_type pos) { return Row<T>(*this, pos);}
    const Row<T> row(size_type pos) const { return Row<T>(*this, pos);}
    Row<T> operator[](size_type pos) { return Row<T>(*this, pos);}
    const Row<T> operator[](size_type pos) const { return Row<T>(*this, pos);}
    Column<T> column(size_type pos) { return Column<T>(*this, pos);};
    const Column<T> column(size_type pos) const { return Column<T>(*this, pos);};

    virtual T &operator()(size_type row, size_type col) = 0;
    virtual const T &operator()(size_type row, size_type col) const = 0;

    T *data() { return elem; }
    const T *data() const { return elem; }
    size_type size() const { return rs * cs; }
    size_type data_size() const { return data_sz; }
    size_type rows() const { return rs; }
    size_type cols() const { return cs; }
    std::pair<size_type, size_type> shape() const { return {rs, cs}; }
};

template <typename T>
constexpr T Base_Matrix<T>::zero;

template <typename T>
std::ostream &operator<<(std::ostream &os, const Base_Matrix<T> &mat)
{
    os << "[\n";
    for (std::size_t r = 0; r < mat.rows(); r++)
    {
        os << "[";
        if (mat.cols())
        {
            os << mat(r, 0);
        }
        for (std::size_t c = 1; c < mat.cols(); c++)
        {
            os << "\t" << mat(r, c);
        }
        os << "]\n";
    }
    os << ']';
    return os;
}

// =========================== Base_Matrix =============================

template <typename T>
Base_Matrix<T>::Base_Matrix(size_type r, size_type c, size_type dt_sz, T deft)
    : rs{r}, cs{c}, data_sz{dt_sz}, elem{new T[dt_sz]{}}
{
    if (deft != T{})
        for (std::size_t i = 0; i < rs * cs; i++)
        {
            elem[i] = deft;
        }
}

template <typename T>
Base_Matrix<T>::Base_Matrix(const Base_Matrix<T> &mat)
    : rs{mat.rs}, cs{mat.cs}, data_sz{mat.data_sz}, elem{new T[mat.data_sz]}
{
    std::copy(mat.elem, mat.elem+data_sz, elem);
}

template <typename T>
Base_Matrix<T>::Base_Matrix(Base_Matrix<T> &&mat)
    : rs{mat.rs}, cs{mat.cs}, data_sz{mat.data_sz}, elem{mat.elem}
{
    mat.elem = nullptr;
}

// ========================= operator= ===================================

template <typename T>
Base_Matrix<T>& Base_Matrix<T>::operator=(const Base_Matrix<T> &mat)
{
    if (this != &mat)
    {
        if (this->shape() != mat.shape())
        {
            throw std::runtime_error("Base_Matrix assignment: different shapes");
        }
        if (this->data_sz != mat.data_sz)
        {
            throw std::runtime_error("Base_Matrix assignment: different data_size");
        }
        std::copy(mat.elem, mat.elem+data_sz, elem);
    }
    return *this;
}

template <typename T>
Base_Matrix<T>& Base_Matrix<T>::operator=(Base_Matrix<T> &&mat)
{
    if (this != &mat)
    {
        if (this->shape() != mat.shape())
        {
            throw std::runtime_error("Base_Matrix assignment: different shapes");
        }
        if (this->data_sz != mat.data_sz)
        {
            throw std::runtime_error("Base_Matrix assignment: different data_size");
        }
        delete[] elem;
        elem = std::exchange(mat.elem, nullptr);
    }
    return *this;
}


} // namespace Misc

#endif // MISC_BASE_MATRIX
