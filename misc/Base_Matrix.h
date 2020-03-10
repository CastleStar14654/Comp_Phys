// Basic definitions of Base_Matrix for Computational Physics of 2020 Spring
// Author: Lin Xuchen, 1 Mar 2020
// Benefit a lot from Matrix.h by Bjarne Stroustrup

#ifndef MISC_BASE_MATRIX
#define MISC_BASE_MATRIX

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>

// the namespace miscellany
namespace Misc
{
template <typename T, size_t R, size_t C>
class Base_Matrix;

// ==================== Row =============================

template <typename T, size_t R, size_t C>
class Row
{
public:
    using size_type = std::size_t;
    using value_type = T;

    constexpr size_type size() const { return C; }

    Row(Base_Matrix<T, R, C> &mat, size_type row) : m{mat}, r{row} {}

    T &operator[](size_type col) { return m(r, col); }
    const T &operator[](size_type col) const
    {
        return const_cast<const Base_Matrix<T, R, C> &>(m)(r, col);
    }
    Row &operator=(const Base_Matrix<T, 1, C> &mat_row)
    {
        for (size_t i = 0; i < C; i++)
        {
            (*this)[i] = mat_row(0, i);
        }
        return *this;
    }
    Row &operator=(Base_Matrix<T, 1, C> &&mat_row)
    {
        for (size_t i = 0; i < C; i++)
        {
            (*this)[i] = std::move(mat_row(0, i));
        }
        return *this;
    }
    template <size_t _R>
    Row &operator=(const Row<T, _R, C> &other)
    {
        for (size_t i = 0; i < C; i++)
        {
            (*this)[i] = other[i];
        }
        return *this;
    }
    template <size_t _R>
    Row &operator=(Row<T, _R, C> &&other)
    {
        for (size_t i = 0; i < C; i++)
        {
            (*this)[i] = other[i];
        }
        return *this;
    }

private:
    Base_Matrix<T, R, C> &m;
    size_type r;
};

// ====================== Column ======================

template <typename T, size_t R, size_t C>
class Column
{
public:
    using size_type = std::size_t;
    using value_type = T;

    constexpr size_type size() const { return R; }

    Column(Base_Matrix<T, R, C> &mat, size_type col) : m{mat}, c{col} {}

    T &operator[](size_type row) { return m(row, c); }
    const T &operator[](size_type row) const
    {
        return const_cast<const Base_Matrix<T, R, C> &>(m)(row, c);
    }
    Column &operator=(const Base_Matrix<T, R, 1> &mat_col)
    {
        for (size_t i = 0; i < R; i++)
        {
            (*this)[i] = mat_row(i, 0);
        }
        return *this;
    }
    Column &operator=(Base_Matrix<T, R, 1> &&mat_col)
    {
        for (size_t i = 0; i < R; i++)
        {
            (*this)[i] = std::move(mat_row(i, 0));
        }
        return *this;
    }
    template <size_t _C>
    Column &operator=(const Column<T, R, _C> &other)
    {
        for (size_t i = 0; i < R; i++)
        {
            (*this)[i] = other[i];
        }
        return *this;
    }
    template <size_t _C>
    Column &operator=(Column<T, R, _C> &&other)
    {
        for (size_t i = 0; i < R; i++)
        {
            (*this)[i] = other[i];
        }
        return *this;
    }

private:
    Base_Matrix<T, R, C> &m;
    size_type c;
};

// ------------------------- operator * --------------------------

template <typename T, size_t R, size_t N, size_t C>
T operator*(const Row<T, R, N> &a, const Column<T, N, C> &b)
{
    T res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

// ================================================================================

// Never initiate a Base_Matrix
template <typename T, size_t R, size_t C>
class Base_Matrix
{
public:
    using size_type = std::size_t;
    using value_type = T;

protected:
    const size_type data_ln;
    using p_ta = T (*)[C];
    p_ta elem;
    constexpr static T zero{};

public:
    Base_Matrix(size_type dt_ln, T (*elm)[C])
        : data_ln{dt_ln}, elem{elm} {}
    Base_Matrix(size_type dt_ln, T deft = T{})
        : data_ln{dt_ln}, elem{new T[dt_ln][C]{}}
    {
        if (deft != T{})
            for (std::size_t i = 0; i < R; i++)
                for (std::size_t j = 0; j < C; j++)
                {
                    elem[i][j] = deft;
                }
    }
    Base_Matrix(const Base_Matrix &mat)
        : data_ln{mat.data_ln}, elem{new T[mat.data_ln][C]}
    {
        std::copy(mat.elem[0], mat.elem[data_ln], elem[0]);
    }
    Base_Matrix(Base_Matrix &&mat)
        : data_ln{mat.data_ln}, elem{mat.elem}
    {
        mat.elem = nullptr;
    }

    virtual ~Base_Matrix() { delete[] elem; }

    Base_Matrix &operator=(const Base_Matrix &mat)
    {
        if (this != &mat)
        {
            if (this->data_ln != mat.data_ln)
            {
                throw std::runtime_error("Base_Matrix assignment: different data_size");
            }
            std::copy(mat.elem, mat.elem + data_ln, elem);
        }
        return *this;
    }
    Base_Matrix &operator=(Base_Matrix &&mat)
    {
        if (this != &mat)
        {
            if (this->data_ln != mat.data_ln)
            {
                throw std::runtime_error("Base_Matrix assignment: different data_size");
            }
            delete[] elem;
            elem = std::exchange(mat.elem, nullptr);
        }
        return *this;
    }

    Row<T, R, C> row(size_type pos) { return Row<T, R, C>(*this, pos); }
    const Row<T, R, C> row(size_type pos) const { return Row<T, R, C>(*this, pos); }
    Row<T, R, C> operator[](size_type pos) { return Row<T, R, C>(*this, pos); }
    const Row<T, R, C> operator[](size_type pos) const { return Row<T, R, C>(*this, pos); }
    Column<T, R, C> column(size_type pos) { return Column<T, R, C>(*this, pos); };
    const Column<T, R, C> column(size_type pos) const { return Column<T, R, C>(*this, pos); };

    virtual T &operator()(size_type row, size_type col) = 0;
    virtual const T &operator()(size_type row, size_type col) const = 0;

    p_ta data() { return elem; }
    const p_ta data() const { return elem; }
    constexpr size_type size() const { return R * C; }
    size_type data_size() const { return data_ln * C; }
    size_type data_lines() const { return data_ln; }
    constexpr size_type rows() const { return R; }
    constexpr size_type cols() const { return C; }
    std::pair<size_type, size_type> shape() const { return {R, C}; }
};

template <typename T, size_t R, size_t C>
constexpr T Base_Matrix<T, R, C>::zero;

template <typename T, size_t R, size_t C>
std::ostream &operator<<(std::ostream &os, const Base_Matrix<T, R, C> &mat)
{
    os << "[\n";
    for (std::size_t r = 0; r < R; r++)
    {
        os << "[";
        if (C)
        {
            os << mat(r, 0);
        }
        for (std::size_t c = 1; c < C; c++)
        {
            os << "\t" << mat(r, c);
        }
        os << "]\n";
    }
    os << ']';
    return os;
}

} // namespace Misc

#endif // MISC_BASE_MATRIX
