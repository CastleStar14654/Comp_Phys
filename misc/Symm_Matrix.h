#ifndef MISC_SYMM_MATRIX
#define MISC_SYMM_MATRIX

#include <iterator>

#include "Base_Matrix.h"
// #include "Matrix.h"
#include "Base_Tri_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Symmetrical Matrix; only half of the items are stored
template <typename T, size_t N>
class Symm_Matrix : public Base_Tri_Matrix<T, N>
{
private:
    // using Base_Tri_Matrix<T, N>::rs;
    // using Base_Tri_Matrix<T, N>::cs;
    using Base_Tri_Matrix<T, N>::elem;
    using Base_Tri_Matrix<T, N>::data_ln;

public:
    using typename Base_Tri_Matrix<T, N>::size_type;
    using Base_Tri_Matrix<T, N>::Base_Tri_Matrix;

    Symm_Matrix(const Symm_Matrix &mat) = default;
    Symm_Matrix(Symm_Matrix &&mat) = default;

    Symm_Matrix &operator=(const Symm_Matrix &mat) = default;
    Symm_Matrix &operator=(Symm_Matrix &&mat) = default;

    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;
};

// ========================== Symm_Matrix =================================

// -------------------- Symm_Matrix: row & column ----------------------------

template <typename T, size_t N>
T &Symm_Matrix<T, N>::operator()(size_type row, size_type col)
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

template <typename T, size_t N>
const T &Symm_Matrix<T, N>::operator()(size_type row, size_type col) const
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

} // namespace Misc

#endif // MISC_SYMM_MATRIX
