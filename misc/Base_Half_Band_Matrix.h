#ifndef MISC_BASE_HALF_BAND_MATRIX
#define MISC_BASE_HALF_BAND_MATRIX

#include <stdexcept>
#include <initializer_list>
#include <utility>

#include "Base_Matrix.h"
#include "Matrix.h"
#include "Diag_Matrix.h"
#include "Band_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Band Matrix
template <typename T, size_t N, size_t M>
class Base_Half_Band_Matrix : public Band_Matrix<T, N, M>
{
public:
    using typename Band_Matrix<T, N, M>::size_type;

    Base_Half_Band_Matrix &operator=(const Base_Half_Band_Matrix &mat) = default;
    Base_Half_Band_Matrix &operator=(Base_Half_Band_Matrix &&mat) = default;
    template <size_t OLD_M>
    Base_Half_Band_Matrix &operator=(const Base_Half_Band_Matrix<T, N, OLD_M> &mat);
    template <size_t OLD_M>
    Base_Half_Band_Matrix &operator=(Base_Half_Band_Matrix<T, N, OLD_M> &&mat);

    virtual T &operator()(size_type row, size_type col) override;
    virtual const T &operator()(size_type row, size_type col) const override;

protected:
    // using Band_Matrix<T, N, M>::rs;
    // using Band_Matrix<T, N, M>::cs;
    using Band_Matrix<T, N, M>::elem;
    using Band_Matrix<T, N, M>::data_ln;
    // const size_type h_bd_w; // half_band_width

    explicit Base_Half_Band_Matrix(T deft = T{});
    Base_Half_Band_Matrix(const Base_Half_Band_Matrix &mat) = default;
    Base_Half_Band_Matrix(Base_Half_Band_Matrix &&mat) = default;

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
    Base_Half_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini);

private:
    // size_type elem_line_len() const { return 2 * rs - h_bd_w - 1; }
    // get an item's place in `elem' from row & col
    std::pair<size_type, size_type> index(size_type row, size_type col) const override
    {
        return (col > row || row > M + col)
            ? std::pair<size_type, size_type>{-1, -1}
            : std::pair<size_type, size_type>{row - col, col};
    }
    // std::pair<size_type, size_type> row_col(size_type idx) const;
};

// ========================== Base_Half_Band_Matrix =================================

template <typename T, size_t N, size_t M>
Base_Half_Band_Matrix<T, N, M>::Base_Half_Band_Matrix(T deft)
    : Band_Matrix<T, N, M>{true, deft}
{
}

template <typename T, size_t N, size_t M>
Base_Half_Band_Matrix<T, N, M>::Base_Half_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Band_Matrix<T, N, M>{true, ini}
{
}

// -------------------- Base_Half_Band_Matrix: operator= ----------------------------

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Base_Half_Band_Matrix<T, N, M> &Base_Half_Band_Matrix<T, N, M>::operator=(const Base_Half_Band_Matrix<T, N, OLD_M> &mat)
{
    static_assert(M > OLD_M, "Base_Half_Band_Matrix: too wide.");
    std::copy(mat.data()[0], mat.data()[mat.data_lines()], elem[0]);
    for (std::size_t i = OLD_M + 1; i < data_ln; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            elem[i][j] = T{};
        }
    return *this;
}

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Base_Half_Band_Matrix<T, N, M> &Base_Half_Band_Matrix<T, N, M>::operator=(Base_Half_Band_Matrix<T, N, OLD_M> &&mat)
{
    static_assert(M > OLD_M, "Base_Half_Band_Matrix: too wide.");
    std::move(mat.data()[0], mat.data()[mat.data_lines()], elem[0]);
    for (std::size_t i = OLD_M + 1; i < data_ln; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            elem[i][j] = T{};
        }
    return *this;
}

// -------------------- Base_Half_Band_Matrix: row & column ----------------------------
template <typename T, size_t N, size_t M>
T &Base_Half_Band_Matrix<T, N, M>::operator()(size_type row, size_type col)
{
    auto idx{index(row, col)};
    if (idx.first == -1)
    {
        throw std::out_of_range("Base_Half_Band_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return elem[idx.first][idx.second];
    }
}

template <typename T, size_t N, size_t M>
const T &Base_Half_Band_Matrix<T, N, M>::operator()(size_type row, size_type col) const
{
    auto idx{index(row, col)};
    if (idx.first == -1)
    {
        return Base_Matrix<T, N, N>::zero;
    }
    else
    {
        return elem[idx.first][idx.second];
    }
}

// // get an item's row & col from index in `elem'
// template <typename T, size_t N, size_t M>
// std::pair<typename Base_Half_Band_Matrix<T, N, M>::size_type, typename Base_Half_Band_Matrix<T, N, M>::size_type>
// Base_Half_Band_Matrix<T, N, M>::row_col(size_type idx) const
// {
//     size_type line_no{idx / rs};
//     size_type rel_loc{idx % rs};
//     if (line_no > h_bd_w)
//     {
//         return {line_no - h_bd_w + rel_loc, rel_loc};
//     }
//     else
//     {
//         return {rel_loc, h_bd_w - line_no + rel_loc};
//     }
// }

} // namespace Misc

#endif // MISC_BASE_HALF_BAND_MATRIX
