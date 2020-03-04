#ifndef MISC_LOW_TRI_MATRIX
#define MISC_LOW_TRI_MATRIX

#include <iterator>

#include "Base_Matrix.h"
#include "Diag_Matrix.h"
#include "Base_Tri_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Lower Triangular Matrix; only half of the items are stored
template <typename T>
class Low_Tri_Matrix : public Base_Tri_Matrix<T>
{
private:
    using Base_Tri_Matrix<T>::rs;
    using Base_Tri_Matrix<T>::cs;
    using Base_Tri_Matrix<T>::elem;
    using Base_Tri_Matrix<T>::data_sz;

public:
    using typename Base_Tri_Matrix<T>::size_type;
    using Base_Tri_Matrix<T>::Base_Tri_Matrix;
    using Base_Tri_Matrix<T>::operator=;
    using Base_Tri_Matrix<T>::operator();

    Low_Tri_Matrix(const Low_Tri_Matrix &mat)=default;
    Low_Tri_Matrix(Low_Tri_Matrix &&mat)=default;

    Low_Tri_Matrix &operator=(const Low_Tri_Matrix &mat) = default;
    Low_Tri_Matrix &operator=(Low_Tri_Matrix &&mat) = default;

    Row<T> row(size_type pos) const override;
    Column<T> column(size_type pos) const override;
    const T &operator()(size_type row, size_type col) const override;
};

// -------------------------------------------------------------------------

template <typename T>
Low_Tri_Matrix<T> operator*(const Low_Tri_Matrix<T> &a, const Low_Tri_Matrix<T> &b)
{
    if (a.shape() != b.shape())
    {
        throw std::domain_error("Low_Tri_Matrix::operator*(): invalid shapes.");
    }

    Low_Tri_Matrix<T> res{a.rows()};
    for (std::size_t i = 0; i < res.rows(); i++)
        for (std::size_t j = 0; j <= i; j++)
            for (std::size_t k = j; k <= i; k++)
            {
                res(i, j) += a(i, k) * b(k, j);
            }
    return res;
}

// ---------------------------------------------------------------------------

template <typename T>
Low_Tri_Matrix<T> operator*(const Low_Tri_Matrix<T> &a, const Diag_Matrix<T> &b)
{
    return a * (*reinterpret_cast<const Low_Tri_Matrix<T> *>(&b));
}

template <typename T>
Low_Tri_Matrix<T> operator*(const Diag_Matrix<T> &a, const Low_Tri_Matrix<T> &b)
{
    return (*reinterpret_cast<const Low_Tri_Matrix<T> *>(&a)) * b;
}

// ========================== Low_Tri_Matrix =================================

// -------------------- Low_Tri_Matrix: row & column ----------------------------

template <typename T>
Row<T> Low_Tri_Matrix<T>::row(size_type pos) const
{
    Row<T> res{Base_Tri_Matrix<T>::row(pos)};
    res.resize(cs, Base_Tri_Matrix<T>::zero);
    return res;
}

template <typename T>
Column<T> Low_Tri_Matrix<T>::column(size_type pos) const
{
    Column<T> res(pos, Base_Tri_Matrix<T>::zero);
    Base_Tri_Matrix<T>::insert_column(pos, std::back_inserter(res));
    return res;
}

template <typename T>
const T &Low_Tri_Matrix<T>::operator()(size_type row, size_type col) const
{
    if (row < col)
    {
        return Base_Matrix<T>::zero;
    }
    else
    {
        return Base_Tri_Matrix<T>::operator()(row, col);
    }
}

} // namespace Misc

#endif // MISC_LOW_TRI_MATRIX
