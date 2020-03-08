#ifndef MISC_BAND_MATRIX
#define MISC_BAND_MATRIX

#include <stdexcept>
#include <initializer_list>
#include <utility>

#include "Base_Matrix.h"
#include "Matrix.h"
#include "Diag_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Band Matrix
template <typename T, size_t N, size_t M>
class Band_Matrix : public Base_Matrix<T, N, N>
{
public:
    using typename Base_Matrix<T, N, N>::size_type;

    explicit Band_Matrix(T deft = T{});
    Band_Matrix(const Band_Matrix &mat) = default;
    Band_Matrix(Band_Matrix &&mat) = default;

    /* For example, the Matrix
        [[2 3 4 0 0 0 0],
         [6 5 4 6 0 0 0],
         [2 3 5 7 8 0 0],
         [0 7 4 5 7 3 0],
         [0 0 9 4 2 4 7],
         [0 0 0 2 3 1 4],
         [0 0 0 0 5 6 7]]
      should be input as
        {{4, 6, 8, 3, 7},
         {3, 4, 7, 7, 4, 4},
         {2, 5, 5, 5, 2, 1, 7},
         {6, 3, 4, 4, 3, 6},
         {2, 7, 9, 2, 5}
        }
     */
    Band_Matrix(std::initializer_list<std::initializer_list<T>> ini);

    Band_Matrix &operator=(const Band_Matrix &mat) = default;
    Band_Matrix &operator=(Band_Matrix &&mat) = default;
    template <size_t OLD_M>
    Band_Matrix &operator=(const Band_Matrix<T, N, OLD_M> &mat);
    template <size_t OLD_M>
    Band_Matrix &operator=(Band_Matrix<T, N, OLD_M> &&mat);

    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;

private:
    // using Base_Matrix<T, N, N>::rs;
    // using Base_Matrix<T, N, N>::cs;
    using Base_Matrix<T, N, N>::elem;
    using Base_Matrix<T, N, N>::data_ln;
    // const size_type h_bd_w; // half_band_width

    // size_type elem_line_len() const { return 2 * rs - h_bd_w - 1; }
    // get an item's place in `elem' from row & col
    std::pair<size_type, size_type> index(size_type row, size_type col) const
    {
        return (col > M + row || row > M + col)
            ? std::pair<size_type, size_type>{-1, -1}
            : std::pair<size_type, size_type>{M + row - col, (row > col) ? col : row};
    }
    // std::pair<size_type, size_type> row_col(size_type idx) const;
};

// ========================== Band_Matrix =================================

template <typename T, size_t N, size_t M>
Band_Matrix<T, N, M>::Band_Matrix(T deft)
    : Base_Matrix<T, N, N>{2 * M + 1, deft}
{
}

template <typename T, size_t N, size_t M>
Band_Matrix<T, N, M>::Band_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Matrix<T, N, N>{2 * M + 1, nullptr}
{
    if (!(ini.size() % 2))
    {
        throw std::invalid_argument("Band_Matrix::Band_Matrix: wrong band width");
    }
    auto it = ini.begin();
    for (size_type count = N - M; count < N; count++)
    {
        if (it->size() != count)
        {
            throw std::invalid_argument("Band_Matrix::Band_Matrix: wrong column number");
        }
        ++it;
    }
    for (size_type count = N; count >= N - M; count--)
    {
        if (it->size() != count)
        {
            throw std::invalid_argument("Band_Matrix::Band_Matrix: wrong column number");
        }
        ++it;
    }

    elem = new T[data_ln][N]{};
    std::size_t count{0};
    for (auto it = ini.begin(); it != ini.end(); it++)
    {
        std::move(it->begin(), it->end(), elem[count]);
        ++count;
    }
}

// -------------------- Band_Matrix: operator= ----------------------------

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Band_Matrix<T, N, M> &Band_Matrix<T, N, M>::operator=(const Band_Matrix<T, N, OLD_M> &mat)
{
    static_assert(M > OLD_M, "Band_Matrix: too wide.");
    for (std::size_t i = 0; i < M - OLD_M; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            elem[i][j] = T{};
        }
    std::copy(mat.data()[0], mat.data()[mat.data_lines()], elem[M - OLD_M]);
    for (std::size_t i = M + OLD_M + 1; i < data_ln; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            elem[i][j] = T{};
        }
    return *this;
}

template <typename T, size_t N, size_t M>
template <size_t OLD_M>
Band_Matrix<T, N, M> &Band_Matrix<T, N, M>::operator=(Band_Matrix<T, N, OLD_M> &&mat)
{
    static_assert(M > OLD_M, "Band_Matrix: too wide.");
    for (std::size_t i = 0; i < M - OLD_M; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            elem[i][j] = T{};
        }
    std::move(mat.data()[0], mat.data()[mat.data_lines()], elem[M - OLD_M]);
    for (std::size_t i = M + OLD_M + 1; i < data_ln; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            elem[i][j] = T{};
        }
    return *this;
}

// -------------------- Band_Matrix: row & column ----------------------------
template <typename T, size_t N, size_t M>
T &Band_Matrix<T, N, M>::operator()(size_type row, size_type col)
{
    auto idx{index(row, col)};
    if (idx.first == -1)
    {
        throw std::out_of_range("Band_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return elem[idx.first][idx.second];
    }
}

template <typename T, size_t N, size_t M>
const T &Band_Matrix<T, N, M>::operator()(size_type row, size_type col) const
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
// std::pair<typename Band_Matrix<T, N, M>::size_type, typename Band_Matrix<T, N, M>::size_type>
// Band_Matrix<T, N, M>::row_col(size_type idx) const
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

#endif // MISC_BAND_MATRIX
