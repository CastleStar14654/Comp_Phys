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
template <typename T>
class Band_Matrix : public Base_Matrix<T>
{
public:
    using typename Base_Matrix<T>::size_type;

    explicit Band_Matrix(size_type n, size_type m, T deft = T{});
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

    Band_Matrix &operator=(const Band_Matrix &mat);
    Band_Matrix &operator=(Band_Matrix &&mat);

    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;

private:
    using Base_Matrix<T>::rs;
    using Base_Matrix<T>::cs;
    using Base_Matrix<T>::elem;
    using Base_Matrix<T>::data_sz;
    const size_type h_bd_w; // half_band_width

    // size_type elem_line_len() const { return 2 * rs - h_bd_w - 1; }
    size_type index(size_type row, size_type col) const;
    std::pair<size_type, size_type> row_col(size_type idx) const;
};

// ========================== Band_Matrix =================================

template <typename T>
Band_Matrix<T>::Band_Matrix(size_type n, size_type m, T deft)
    : Base_Matrix<T>{n, n, n * (2 * m + 1), deft}, h_bd_w{m}
{
}

template <typename T>
Band_Matrix<T>::Band_Matrix(std::initializer_list<std::initializer_list<T>> ini)
    : Base_Matrix<T>{ini.begin()->size() + ini.size() / 2,
                     ini.begin()->size() + ini.size() / 2,
                     ini.size() * (ini.begin()->size() + ini.size() / 2),
                     nullptr},
      h_bd_w{ini.size() / 2}
{
    if (!(ini.size() % 2))
    {
        throw std::invalid_argument("Band_Matrix::Band_Matrix: wrong band width");
    }
    auto it = ini.begin();
    for (size_type count = cs - h_bd_w; count < cs; count++)
    {
        if (it->size() != count)
        {
            throw std::invalid_argument("Band_Matrix::Band_Matrix: wrong column number");
        }
        ++it;
    }
    for (size_type count = cs; count >= cs - h_bd_w; count--)
    {
        if (it->size() != count)
        {
            throw std::invalid_argument("Band_Matrix::Band_Matrix: wrong column number");
        }
        ++it;
    }

    elem = new T[data_sz]{};
    std::size_t count{0};
    for (auto it = ini.begin(); it != ini.end(); it++)
    {
        std::move(it->begin(), it->end(), &elem[count * rs]);
        ++count;
    }
}

// -------------------- Band_Matrix: operator= ----------------------------

template <typename T>
Band_Matrix<T> &Band_Matrix<T>::operator=(const Band_Matrix<T> &mat)
{
    if (this != &mat)
    {
        if (h_bd_w < mat.h_bd_w)
        {
            throw std::runtime_error("Base_Matrix assignment: different half_band_width");
        }
        else if (h_bd_w == mat.h_bd_w)
        {
            Base_Matrix<T>::operator=(mat);
        }
        else
        {
            for (std::size_t i = 0; i < (h_bd_w - mat.h_bd_w) * cs; i++)
            {
                elem[i] = T{};
            }
            std::copy(mat.elem, &mat.elem[mat.data_sz], &elem[(h_bd_w - mat.h_bd_w) * cs]);
            for (std::size_t i = (h_bd_w + 1 + mat.h_bd_w) * cs; i < data_sz; i++)
            {
                elem[i] = T{};
            }
        }
    }
    return *this;
}

template <typename T>
Band_Matrix<T> &Band_Matrix<T>::operator=(Band_Matrix<T> &&mat)
{
    if (this != &mat)
    {
        if (this->h_bd_w < mat.h_bd_w)
        {
            throw std::runtime_error("Base_Matrix assignment: different half_band_width");
        }
        else if (this->h_bd_w == mat.h_bd_w)
        {
            Base_Matrix<T>::operator=(mat);
        }
        else
        {
            for (std::size_t i = 0; i < (h_bd_w - mat.h_bd_w) * cs; i++)
            {
                elem[i] = T{};
            }
            std::move(mat.elem, &mat.elem[mat.data_sz], &elem[(h_bd_w - mat.h_bd_w) * cs]);
            for (std::size_t i = (h_bd_w + 1 + mat.h_bd_w) * cs; i < data_sz; i++)
            {
                elem[i] = T{};
            }
        }
    }
    return *this;
}

// -------------------- Band_Matrix: row & column ----------------------------

template <typename T>
T &Band_Matrix<T>::operator()(size_type row, size_type col)
{
    size_type idx{index(row, col)};
    if (idx == -1)
    {
        throw std::out_of_range("Band_Matrix::operator(): trying to access empty area.");
    }
    else
    {
        return elem[idx];
    }
}

template <typename T>
const T &Band_Matrix<T>::operator()(size_type row, size_type col) const
{
    size_type idx{index(row, col)};
    if (idx == -1)
    {
        return Base_Matrix<T>::zero;
    }
    else
    {
        return elem[idx];
    }
}

// get an item's index in `elem' from row & col
template <typename T>
typename Band_Matrix<T>::size_type Band_Matrix<T>::index(size_type row, size_type col) const
{
    if (col > h_bd_w + row || row > h_bd_w + col)
    {
        return -1;
    }
    else
    {
        size_type line_no{h_bd_w + row - col};
        if (row > col)
        {
            return line_no * rs + col;
        }
        else
        {
            return line_no * rs + row;
        }
    }
}

// get an item's row & col from index in `elem'
template <typename T>
std::pair<typename Band_Matrix<T>::size_type, typename Band_Matrix<T>::size_type>
Band_Matrix<T>::row_col(size_type idx) const
{
    size_type line_no{idx / rs};
    size_type rel_loc{idx % rs};
    if (line_no > h_bd_w)
    {
        return {line_no - h_bd_w + rel_loc, rel_loc};
    }
    else
    {
        return {rel_loc, h_bd_w - line_no + rel_loc};
    }
}

} // namespace Misc

#endif // MISC_BAND_MATRIX
