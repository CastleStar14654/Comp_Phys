// Basic definitions of Base_Matrix for Computational Physics of 2020 Spring
// Author: Lin Xuchen, 1 Mar 2020
// Benefit a lot from Matrix.h by Bjarne Stroustrup

#ifndef MISC_BASE_MATRIX
#define MISC_BASE_MATRIX

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>

// the namespace miscellany
namespace Misc
{

// Alias for convenience
template <typename T>
using Row = std::vector<T>;
template <typename T>
using Column = std::vector<T>;

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

    virtual Row<T> row(size_type pos) const = 0;
    virtual Column<T> column(size_type pos) const = 0;

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
