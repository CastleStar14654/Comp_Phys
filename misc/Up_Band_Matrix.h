#ifndef MISC_SYMM_BAND_MATRIX
#define MISC_SYMM_BAND_MATRIX

#include <stdexcept>
#include <initializer_list>
#include <utility>

#include "Base_Matrix.h"
#include "Matrix.h"
#include "Diag_Matrix.h"
#include "Base_Half_Band_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Band Matrix
template <typename T, size_t N, size_t M>
class Up_Band_Matrix : public Base_Half_Band_Matrix<T, N, M>
{
public:
    using typename Base_Half_Band_Matrix<T, N, M>::size_type;

    explicit Up_Band_Matrix(T deft = T{});
    Up_Band_Matrix(const Up_Band_Matrix &mat) = default;
    Up_Band_Matrix(Up_Band_Matrix &&mat) = default;

    /* For example, the Matrix
        [[2 6 2 0 0 0 0]
         [0 5 3 7 0 0 0]
         [0 0 5 4 9 0 0]
         [0 0 0 5 4 2 0]
         [0 0 0 0 2 3 5]
         [0 0 0 0 0 1 6]
         [0 0 0 0 0 0 7]
        ]
      should be input as
        {{2, 5, 5, 5, 2, 1, 7},
         {6, 3, 4, 4, 3, 6},
         {2, 7, 9, 2, 5}
        }
     */
    Up_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Up_Band_Matrix &operator=(const Up_Band_Matrix &mat) = default;
    Up_Band_Matrix &operator=(Up_Band_Matrix &&mat) = default;
    template <size_t OLD_M>
    Up_Band_Matrix &operator=(const Up_Band_Matrix<T, N, OLD_M> &mat);
    template <size_t OLD_M>
    Up_Band_Matrix &operator=(Up_Band_Matrix<T, N, OLD_M> &&mat);

    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;

private:
    using Base_Half_Band_Matrix<T, N, M>::elem;
    using Base_Half_Band_Matrix<T, N, M>::data_ln;

    // size_type elem_line_len() const { return 2 * rs - h_bd_w - 1; }
    // std::pair<size_type, size_type> row_col(size_type idx) const;
};

// ========================== Up_Band_Matrix =================================

template <typename T, size_t N, size_t M>
Up_Band_Matrix<T, N, M>::Up_Band_Matrix(T deft)
    : Base_Half_Band_Matrix<T, N, M>{deft}
{
}

template <typename T, size_t N, size_t M>
Up_Band_Matrix<T, N, M>::Up_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Half_Band_Matrix<T, N, M>{ini}
{
}

// -------------------- Up_Band_Matrix: operator= ----------------------------

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Up_Band_Matrix<T, N, M> &Up_Band_Matrix<T, N, M>::operator=(const Up_Band_Matrix<T, N, OLD_M> &mat)
{
    Base_Half_Band_Matrix<T, N, M>::operator=(mat);
    return *this;
}

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Up_Band_Matrix<T, N, M> &Up_Band_Matrix<T, N, M>::operator=(Up_Band_Matrix<T, N, OLD_M> &&mat)
{
    Base_Half_Band_Matrix<T, N, M>::operator=(mat);
    return *this;
}

// -------------------- Up_Band_Matrix: row & column ----------------------------

template <typename T, size_t N, size_t M>
T &Up_Band_Matrix<T, N, M>::operator()(size_type row, size_type col)
{
    if (col < row)
    {
        throw std::out_of_range("Up_Band_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return Base_Half_Band_Matrix<T, N, M>::operator()(col, row);
    }
}

template <typename T, size_t N, size_t M>
const T &Up_Band_Matrix<T, N, M>::operator()(size_type row, size_type col) const
{
    if (col < row)
    {
        return Base_Matrix<T, N, N>::zero;
    }
    else
    {
        return Base_Half_Band_Matrix<T, N, M>::operator()(col, row);
    }
}

} // namespace Misc

#endif // MISC_SYMM_BAND_MATRIX
