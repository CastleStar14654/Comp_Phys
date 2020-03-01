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
    Base_Matrix(const Base_Matrix &mat) = delete;
    Base_Matrix(Base_Matrix &&mat) = delete;

    virtual ~Base_Matrix() { delete[] elem; }

    Base_Matrix &operator=(const Base_Matrix &mat) = delete;
    Base_Matrix &operator=(Base_Matrix &&mat) = delete;

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

// ======================================================================

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
    Matrix(std::initializer_list<std::initializer_list<T>> ini);

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

// ======================================================================

// Diagonal Matrix; only diagonal items are stored
template <typename T>
class Diag_Matrix : public Base_Matrix<T>
{
private:
    using Base_Matrix<T>::rs;
    using Base_Matrix<T>::cs;
    using Base_Matrix<T>::elem;
    using Base_Matrix<T>::data_sz;

public:
    using typename Base_Matrix<T>::size_type;

    explicit Diag_Matrix(size_type n, T deft = T{});
    Diag_Matrix(const Diag_Matrix &mat);
    Diag_Matrix(Diag_Matrix &&mat);
    template <typename It>
    explicit Diag_Matrix(It b, It e);

    Diag_Matrix &operator=(const Diag_Matrix &mat);
    Diag_Matrix &operator=(Diag_Matrix &&mat);

    Row<T> row(size_type pos) const override;
    Column<T> column(size_type pos) const override;
    T &operator()(size_type row, size_type col) override;
    T &operator()(size_type n) { return elem[n]; }
    const T &operator()(size_type row, size_type col) const override;
    const T &operator()(size_type n) const { return elem[n]; }
};

// -------------------------------------------------------------------------

template <typename T>
Diag_Matrix<T> operator*(const Diag_Matrix<T> &a, const Diag_Matrix<T> &b)
{
    if (a.cols() != b.rows())
    {
        throw std::domain_error("Matrix::operator*(): invalid shapes.");
    }

    Matrix<T> res{a.rows(), b.cols()};
    for (std::size_t i = 0; i < res.rows(); i++)
    {
        res(i, i) = a(i, i) * b(i, i);
    }
    return res;
}

template <typename T>
Matrix<T> operator*(const Diag_Matrix<T> &a, const Base_Matrix<T> &b)
{
    if (a.cols() != b.rows())
    {
        throw std::domain_error("Matrix::operator*(): invalid shapes.");
    }

    Matrix<T> res{a.rows(), b.cols()};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = 0; j < res.cols(); j++)
        {
            res(i, j) = a(i, i) * b(i, j);
        }
    return res;
}

template <typename T>
Matrix<T> operator*(const Base_Matrix<T> &a, const Diag_Matrix<T> &b)
{
    if (a.cols() != b.rows())
    {
        throw std::domain_error("Matrix::operator*(): invalid shapes.");
    }

    Matrix<T> res{a.rows(), b.cols()};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = 0; j < res.cols(); j++)
        {
            res(i, j) = a(i, j) * b(j, j);
        }
    return res;
}

// =====================Matrix===============================

template <typename T>
Matrix<T>::Matrix(size_type r, size_type c, T deft)
    : Base_Matrix<T>{r, c, r * c, new T[r * c]{}}
{
    if (deft != T{})
        for (size_t i = 0; i < rs * cs; i++)
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
    Row<T> res(cs);
    for (size_t c = 0; c < cs; c++)
    {
        res[c] = this->operator()(pos, c);
    }
    return res;
}

template <typename T>
Column<T> Matrix<T>::column(size_type pos) const
{
    Column<T> res(rs);
    for (size_t r = 0; r < rs; r++)
    {
        res[r] = this->operator()(r, pos);
    }
    return res;
}

// ===========================Diag_Matrix====================================

template <typename T>
Diag_Matrix<T>::Diag_Matrix(size_type n, T deft)
    : Base_Matrix<T>{n, n, n, new T[n]{}}
{
    if (deft != T{})
        for (size_t i = 0; i < data_sz; i++)
        {
            elem[i] = deft;
        }
}

template <typename T>
Diag_Matrix<T>::Diag_Matrix(const Diag_Matrix &mat)
    : Base_Matrix<T>{mat.rs, mat.cs, mat.data_sz, new T[mat.data_sz]}
{
    std::copy(mat.elem, mat.elem + data_sz, elem);
}

template <typename T>
Diag_Matrix<T>::Diag_Matrix(Diag_Matrix &&mat)
    : Base_Matrix<T>{mat.rs, mat.cs, mat.data_sz, mat.elem}
{
    mat.elem = nullptr;
}

template <typename T>
template <typename It>
Diag_Matrix<T>::Diag_Matrix(It b, It e)
    : Base_Matrix<T>(e - b, e - b, e - b, new T[e - b])
{
    std::copy(b, e, elem);
}

// ----------------------- Diag_Matrix operator= --------------------------

template <typename T>
Diag_Matrix<T> &Diag_Matrix<T>::operator=(const Diag_Matrix &mat)
{
    T *temp = new T[mat.data_size()];
    std::copy(mat.elem, mat.elem + data_sz, temp);
    delete[] elem;
    elem = temp;
    rs = mat.rs;
    cs = mat.cs;
    data_sz = mat.data_sz;
    return *this;
}

template <typename T>
Diag_Matrix<T> &Diag_Matrix<T>::operator=(Diag_Matrix &&mat)
{
    delete[] elem;
    elem = mat.elem;
    mat.elem = nullptr;
    rs = mat.rs;
    cs = mat.cs;
    data_sz = mat.data_sz;
    return *this;
}

// ------------------------Diag_Matrix row() & column() ------------------------------

template <typename T>
Row<T> Diag_Matrix<T>::row(size_type pos) const
{
    Row<T> res(cs, Base_Matrix<T>::zero);
    res[pos] = elem[pos];
    return res;
}

template <typename T>
Column<T> Diag_Matrix<T>::column(size_type pos) const
{
    Column<T> res(rs, Base_Matrix<T>::zero);
    res[pos] = elem[pos];
    return res;
}

template <typename T>
T &Diag_Matrix<T>::operator()(size_type row, size_type col)
{
    if (row == col)
    {
        return (*this)(row);
    }
    else
    {
        throw std::out_of_range("Diag_Matrix::operator(): trying to access 0.");
    }
}

template <typename T>
const T &Diag_Matrix<T>::operator()(size_type row, size_type col) const
{
    if (row == col)
    {
        return (*this)(row);
    }
    else
    {
        return Base_Matrix<T>::zero;
    }
}

} // namespace Misc

#endif // MISC
