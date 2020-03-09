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
class Low_Band_Matrix : public Base_Half_Band_Matrix<T, N, M>
{
public:
    using typename Base_Half_Band_Matrix<T, N, M>::size_type;
    using Base_Half_Band_Matrix<T, N, M>::operator();

    explicit Low_Band_Matrix(T deft = T{});
    Low_Band_Matrix(const Low_Band_Matrix &mat) = default;
    Low_Band_Matrix(Low_Band_Matrix &&mat) = default;

    /* For example, the Matrix
        [[2 0 0 0 0 0 0],
         [6 5 0 0 0 0 0],
         [2 3 5 0 0 0 0],
         [0 7 4 5 0 0 0],
         [0 0 9 4 2 0 0],
         [0 0 0 2 3 1 0],
         [0 0 0 0 5 6 7]]
      should be input as
        {{2, 5, 5, 5, 2, 1, 7},
         {6, 3, 4, 4, 3, 6},
         {2, 7, 9, 2, 5}
        }
     */
    Low_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Low_Band_Matrix &operator=(const Low_Band_Matrix &mat) = default;
    Low_Band_Matrix &operator=(Low_Band_Matrix &&mat) = default;
    template <size_t OLD_M>
    Low_Band_Matrix &operator=(const Low_Band_Matrix<T, N, OLD_M> &mat);
    template <size_t OLD_M>
    Low_Band_Matrix &operator=(Low_Band_Matrix<T, N, OLD_M> &&mat);

    // T &operator()(size_type row, size_type col) override;
    // const T &operator()(size_type row, size_type col) const override;

private:
    // using Base_Half_Band_Matrix<T, N, M>::rs;
    // using Base_Half_Band_Matrix<T, N, M>::cs;
    using Base_Half_Band_Matrix<T, N, M>::elem;
    using Base_Half_Band_Matrix<T, N, M>::data_ln;
    // const size_type h_bd_w; // half_band_width

    // size_type elem_line_len() const { return 2 * rs - h_bd_w - 1; }
    // std::pair<size_type, size_type> row_col(size_type idx) const;
};

// ========================== Low_Band_Matrix =================================

template <typename T, size_t N, size_t M>
Low_Band_Matrix<T, N, M>::Low_Band_Matrix(T deft)
    : Base_Half_Band_Matrix<T, N, M>{deft}
{
}

template <typename T, size_t N, size_t M>
Low_Band_Matrix<T, N, M>::Low_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Half_Band_Matrix<T, N, M>{ini}
{
}

// -------------------- Low_Band_Matrix: operator= ----------------------------

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Low_Band_Matrix<T, N, M> &Low_Band_Matrix<T, N, M>::operator=(const Low_Band_Matrix<T, N, OLD_M> &mat)
{
    Base_Half_Band_Matrix<T, N, M>::operator=(mat);
    return *this;
}

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Low_Band_Matrix<T, N, M> &Low_Band_Matrix<T, N, M>::operator=(Low_Band_Matrix<T, N, OLD_M> &&mat)
{
    Base_Half_Band_Matrix<T, N, M>::operator=(mat);
    return *this;
}

} // namespace Misc

#endif // MISC_SYMM_BAND_MATRIX
