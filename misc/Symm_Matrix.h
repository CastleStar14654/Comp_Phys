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
template <typename T>
class Symm_Matrix : public Base_Tri_Matrix<T>
{
private:
    using Base_Tri_Matrix<T>::rs;
    using Base_Tri_Matrix<T>::cs;
    using Base_Tri_Matrix<T>::elem;
    using Base_Tri_Matrix<T>::data_sz;

public:
    using typename Base_Tri_Matrix<T>::size_type;
    using Base_Tri_Matrix<T>::Base_Tri_Matrix;

    Symm_Matrix(const Symm_Matrix &mat) = default;
    Symm_Matrix(Symm_Matrix &&mat) = default;

    Symm_Matrix &operator=(const Symm_Matrix &mat) = default;
    Symm_Matrix &operator=(Symm_Matrix &&mat) = default;

    T &operator()(size_type row, size_type col) override;
    const T &operator()(size_type row, size_type col) const override;
};

// ========================== Symm_Matrix =================================

// -------------------- Symm_Matrix: row & column ----------------------------

template <typename T>
T &Symm_Matrix<T>::operator()(size_type row, size_type col)
{
    if (row < col)
    {
        return Base_Tri_Matrix<T>::operator()(col, row);
    }
    else
    {
        return Base_Tri_Matrix<T>::operator()(row, col);
    }
}

template <typename T>
const T &Symm_Matrix<T>::operator()(size_type row, size_type col) const
{
    if (row < col)
    {
        return Base_Tri_Matrix<T>::operator()(col, row);
    }
    else
    {
        return Base_Tri_Matrix<T>::operator()(row, col);
    }
}

} // namespace Misc

#endif // MISC_SYMM_MATRIX
