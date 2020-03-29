#ifndef MISC_BASE_TRI_MATRIX
#define MISC_BASE_TRI_MATRIX

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

// Base Triangular Matrix; only half of the items are stored
// being inherited by Symm_Matrix, Up_Tri_Matrix, Low_Tri_Matrix
template <typename T, size_t N>
class Base_Tri_Matrix : public Base_Matrix<T, N, N>
{
public:
    using typename Base_Matrix<T, N, N>::size_type;

    explicit Base_Tri_Matrix(T deft = T{})
        : Base_Matrix<T, N, N>{N / 2 + 1, deft} {}
    Base_Tri_Matrix(const Base_Tri_Matrix &mat) = default;
    Base_Tri_Matrix(Base_Tri_Matrix &&mat) = default;
    Base_Tri_Matrix(const Diag_Matrix<T, N> &mat)
        : Base_Matrix<T, N, N>{N / 2 + 1, T{}}
    {
        for (std::size_t i = 0; i < N; i++)
        {
            (*this)(i, i) = mat(i);
        }
    }
    Base_Tri_Matrix(Diag_Matrix<T, N> &&mat)
        : Base_Matrix<T, N, N>{N / 2 + 1, T{}}
    {
        for (std::size_t i = 0; i < N; i++)
        {
            (*this)(i, i) = std::move(mat(i));
        }
    }
    template <size_t M>
    Base_Tri_Matrix(const Base_Half_Band_Matrix<T, N, M> &mat)
        : Base_Matrix<T, N, N>{N / 2 + 1, T{}}
    {
        for (size_t del = 0; del <= M; del++)
            for (size_t j = 0; j < N - del; j++)
            {
                (*this)(j + del, j) = mat.Base_Half_Band_Matrix<T, N, M>::operator()(j + del, j);
            }
    }
    template <size_t M>
    Base_Tri_Matrix(Base_Half_Band_Matrix<T, N, M> &&mat)
        : Base_Matrix<T, N, N>{N / 2 + 1, T{}}
    {
        for (size_t del = 0; del <= M; del++)
            for (size_t j = 0; j < N - del; j++)
            {
                (*this)(j + del, j) = std::move(mat.Base_Half_Band_Matrix<T, N, M>::operator()(j + del, j));
            }
    }
    Base_Tri_Matrix(std::initializer_list<std::initializer_list<T>> ini)
        : Base_Matrix<T, N, N>{N / 2 + 1, nullptr}
    {
        int count{1};
        for (auto i = ini.begin(); i != ini.end(); i++)
        {
            if (i->size() != count)
            {
                throw std::invalid_argument(__FILE__ + ":" + __LINE__ + ": wrong column number");
            }
            ++count;
        }

        elem = new T[data_ln][N];
        auto it = ini.begin();
        for (std::size_t i = 0; i < N; i++)
        {
            std::move(it->begin(), it->end(), &(*this)(i, 0));
            it++;
        }
    }

    Base_Tri_Matrix &operator=(const Base_Tri_Matrix &mat) = default;
    Base_Tri_Matrix &operator=(Base_Tri_Matrix &&mat) = default;
    Base_Tri_Matrix &operator=(const Diag_Matrix<T, N> &mat)
    {
        T(*temp)
        [N] = new T[data_ln][N]{};
        delete[] elem;
        elem = temp;
        for (std::size_t i = 0; i < N; i++)
        {
            (*this)(i, i) = mat(i, i);
        }
        return *this;
    }
    Base_Tri_Matrix &operator=(Diag_Matrix<T, N> &&mat)
    {
        T(*temp)
        [N] = new T[data_ln][N]{};
        delete[] elem;
        elem = temp;
        for (std::size_t i = 0; i < N; i++)
        {
            (*this)(i, i) = std::move(mat(i, i));
        }
        return *this;
    }
    template <size_t M>
    Base_Tri_Matrix &operator=(const Base_Half_Band_Matrix<T, N, M> &mat)
    {
        T(*temp)
        [N] = new T[data_ln][N]{};
        delete[] elem;
        elem = temp;
        for (size_t del = 0; del <= M; del++)
            for (size_t j = 0; j < N - del; j++)
            {
                (*this)(j + del, j) = mat.Base_Half_Band_Matrix<T, N, M>::operator()(j + del, j);
            }
        return *this;
    }
    template <size_t M>
    Base_Tri_Matrix &operator=(Base_Half_Band_Matrix<T, N, M> &&mat)
    {
        T(*temp)
        [N] = new T[data_ln][N]{};
        delete[] elem;
        elem = temp;
        for (size_t del = 0; del <= M; del++)
            for (size_t j = 0; j < N - del; j++)
            {
                (*this)(j + del, j) = std::move(mat.Base_Half_Band_Matrix<T, N, M>::operator()(j + del, j));
            }
        return *this;
    }

protected:
    using Base_Matrix<T, N, N>::elem;
    using Base_Matrix<T, N, N>::data_ln;

    // return the indices of the first item of a row in the elem
    constexpr size_type head_i(size_type row) const { return (row >= (N - 1) / 2) ? (row - (N - 1) / 2) : (N / 2 - 1 - row); }
    constexpr size_type head_j(size_type row) const { return (row >= (N - 1) / 2) ? (0) : (N - 1 - row); }

    virtual T &operator()(size_type row, size_type col) override
    {
        if (row < col)
        {
            throw std::out_of_range(__FILE__ + ":" + __LINE__ + ": trying to access empty area.");
        }
        else
        {
            return elem[head_i(row)][head_j(row) + col];
        }
    }
    virtual const T &operator()(size_type row, size_type col) const override
    {
        if (row < col)
        {
            throw std::out_of_range(__FILE__ + ":" + __LINE__ + ": trying to access empty area.");
        }
        else
        {
            return elem[head_i(row)][head_j(row) + col];
        }
    }
};

} // namespace Misc

#endif // MISC_BASE_TRI_MATRIX
