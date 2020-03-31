#ifndef MISC_SYMM_MATRIX
#define MISC_SYMM_MATRIX

#include <iterator>

#include "Base_Matrix.h"
#include "Base_Tri_Matrix.h"
#include "Symm_Band_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Symmetrical Matrix; only half of the items are stored
template <typename T, size_t N>
class Symm_Matrix : public Base_Tri_Matrix<T, N>
{
private:
    using Base_Tri_Matrix<T, N>::elem;
    using Base_Tri_Matrix<T, N>::data_ln;

public:
    using typename Base_Tri_Matrix<T, N>::size_type;
    using Base_Tri_Matrix<T, N>::Base_Tri_Matrix;

    Symm_Matrix(const Symm_Matrix &mat) = default;
    Symm_Matrix(Symm_Matrix &&mat) = default;
    template <size_t M>
    Symm_Matrix(const Symm_Band_Matrix<T, N, M> &mat)
        : Base_Tri_Matrix<T, N>{mat} {}
    template <size_t M>
    Symm_Matrix(Symm_Band_Matrix<T, N, M> &&mat)
        : Base_Tri_Matrix<T, N>{mat} {}

    Symm_Matrix &operator=(const Symm_Matrix &mat) = default;
    Symm_Matrix &operator=(Symm_Matrix &&mat) = default;
    Symm_Matrix &operator=(const Diag_Matrix<T, N> &mat)
    {
        Base_Tri_Matrix<T, N>::operator=(mat);
        return *this;
    }
    Symm_Matrix &operator=(Diag_Matrix<T, N> &&mat)
    {
        Base_Tri_Matrix<T, N>::operator=(mat);
        return *this;
    }
    template <size_t M>
    Symm_Matrix &operator=(const Symm_Band_Matrix<T, N, M> &mat)
    {
        Base_Tri_Matrix<T, N>::operator=(mat);
        return *this;
    }
    template <size_t M>
    Symm_Matrix &operator=(Symm_Band_Matrix<T, N, M> &&mat)
    {
        Base_Tri_Matrix<T, N>::operator=(mat);
        return *this;
    }

    T &operator()(size_type row, size_type col) override
    {
        if (row < col)
        {
            return Base_Tri_Matrix<T, N>::operator()(col, row);
        }
        else
        {
            return Base_Tri_Matrix<T, N>::operator()(row, col);
        }
    }
    const T &operator()(size_type row, size_type col) const override
    {
        if (row < col)
        {
            return Base_Tri_Matrix<T, N>::operator()(col, row);
        }
        else
        {
            return Base_Tri_Matrix<T, N>::operator()(row, col);
        }
    }
};

} // namespace Misc

#endif // MISC_SYMM_MATRIX
