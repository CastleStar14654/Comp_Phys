#ifndef MISC_SPARSE_MATRIX
#define MISC_SPARSE_MATRIX

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>
#include <map>
#include <array>

#include "Base_Matrix.h"

// the namespace miscellany
namespace Misc
{

template <typename T, size_t R, size_t C>
class Sparse_Matrix : public Base_Matrix<T, R, C>, public std::array<std::map<size_t, T>, R>
{
public:
    using typename Base_Matrix<T, R, C>::size_type;
    using std::array<std::map<size_t, T>, R>::operator[];

    Sparse_Matrix()
        : Base_Matrix<T, R, C>{0, nullptr}, std::array<std::map<size_t, T>, R>{} {}
    Sparse_Matrix(const Base_Matrix<T, R, C> &mat)
        : Base_Matrix<T, R, C>{0, nullptr}, std::array<std::map<size_t, T>, R>{}
    {
        for (size_t i = 0; i < R; i++)
            for (size_t j = 0; j < C; j++)
                if (mat(i, j) != Base_Matrix<T, R, C>::zero)
                {
                    std::array<std::map<size_t, T>, R>::operator[](i)[j] = mat(i, j);
                }
    }
    Sparse_Matrix(Base_Matrix<T, R, C> &&mat)
        : Base_Matrix<T, R, C>{0, nullptr}, std::array<std::map<size_t, T>, R>{}
    {
        for (size_t i = 0; i < R; i++)
            for (size_t j = 0; j < C; j++)
                if (mat(i, j) != Base_Matrix<T, R, C>::zero)
                {
                    std::array<std::map<size_t, T>, R>::operator[](i)[j] = std::move(mat(i, j));
                }
    }
    template <size_t M>
    Sparse_Matrix(const Band_Matrix<T, R, M> &mat)
        : Base_Matrix<T, R, C>{0, nullptr}, std::array<std::map<size_t, T>, R>{}
    {
        static_assert(R == C);
        for (size_t i = 0; i < R; i++)
            for (size_t j = i > M ? i - M : 0; j < std::min(i + M + 1, R); j++)
                if (mat(i, j) != Base_Matrix<T, R, C>::zero)
                {
                    std::array<std::map<size_t, T>, R>::operator[](i)[j] = mat(i, j);
                }
    }
    template <size_t M>
    Sparse_Matrix(Band_Matrix<T, R, M> &&mat)
        : Base_Matrix<T, R, C>{0, nullptr}, std::array<std::map<size_t, T>, R>{}
    {
        static_assert(R == C);
        for (size_t i = 0; i < R; i++)
            for (size_t j = i > M ? i - M : 0; j < std::min(i + M + 1, R); j++)
                if (mat(i, j) != Base_Matrix<T, R, C>::zero)
                {
                    std::array<std::map<size_t, T>, R>::operator[](i)[j] = std::move(mat(i, j));
                }
    }

    Sparse_Matrix &operator=(const Base_Matrix<T, R, C> &mat)
    {
        for (size_t i = 0; i < R; i++)
            for (size_t j = 0; j < C; j++)
                if (mat(i, j) != Base_Matrix<T, R, C>::zero)
                {
                    std::array<std::map<size_t, T>, R>::operator[](i)[j] = mat(i, j);
                }
        return *this;
    }
    Sparse_Matrix &operator=(Base_Matrix<T, R, C> &&mat)
    {
        for (size_t i = 0; i < R; i++)
            for (size_t j = 0; j < C; j++)
                if (mat(i, j) != Base_Matrix<T, R, C>::zero)
                {
                    std::array<std::map<size_t, T>, R>::operator[](i)[j] = std::move(mat(i, j));
                }
        return *this;
    }

    T &operator()(size_type row, size_type col) override
    {
        return std::array<std::map<size_t, T>, R>::operator[](row)[col];
    }
    const T &operator()(size_type row, size_type col) const override
    {
        auto it {std::array<std::map<size_t, T>, R>::operator[](row).find(col)};
        if (it == std::array<std::map<size_t, T>, R>::operator[](row).cend())
        {
            return Base_Matrix<T, R, C>::zero;
        }
        else
        {
            return it->second;
        }
    }
};

} // namespace Misc

#endif // MISC_SPARSE_MATRIX
