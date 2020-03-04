#ifndef MISC_DIAG_MATRIX
#define MISC_DIAG_MATRIX

#include "Base_Matrix.h"
#include "Matrix.h"

// the namespace miscellany
namespace Misc
{

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
    Diag_Matrix(const Diag_Matrix &mat)=default;
    Diag_Matrix(Diag_Matrix &&mat)=default;
    template <typename It>
    explicit Diag_Matrix(It b, It e);
    explicit Diag_Matrix(std::initializer_list<T> ini);

    Diag_Matrix &operator=(const Diag_Matrix &mat)=default;
    Diag_Matrix &operator=(Diag_Matrix &&mat)=default;

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
        throw std::domain_error("Diag_Matrix::operator*(): invalid shapes.");
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
        throw std::domain_error("Diag_Matrix::operator*(): invalid shapes.");
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
        throw std::domain_error("Diag_Matrix::operator*(): invalid shapes.");
    }

    Matrix<T> res{a.rows(), b.cols()};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = 0; j < res.cols(); j++)
        {
            res(i, j) = a(i, j) * b(j, j);
        }
    return res;
}

// ===========================Diag_Matrix====================================

template <typename T>
Diag_Matrix<T>::Diag_Matrix(size_type n, T deft)
    : Base_Matrix<T>{n, n, n, deft}
{
}

template <typename T>
template <typename It>
Diag_Matrix<T>::Diag_Matrix(It b, It e)
    : Base_Matrix<T>(e - b, e - b, e - b, new T[e - b])
{
    std::copy(b, e, elem);
}

template <typename T>
Diag_Matrix<T>::Diag_Matrix(std::initializer_list<T> ini)
    : Base_Matrix<T>(ini.size(), ini.size(), ini.size(), new T[ini.size()])
{
    std::move(ini.begin(), ini.end(), elem);
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

#endif // MISC_DIAG_MATRIX
