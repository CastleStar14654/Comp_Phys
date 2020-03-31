#ifndef MISC_HERMITE_MATRIX
#define MISC_HERMITE_MATRIX

#include <iterator>
#include <complex>

#include "Base_Matrix.h"
#include "Base_Tri_Matrix.h"

// the namespace miscellany
namespace Misc
{

// Symmetrical Matrix; only half of the items are stored
template <typename T, size_t N>
class Hermite_Matrix : public Base_Tri_Matrix<T, N>
{
private:
    using Base_Tri_Matrix<T, N>::elem;
    using Base_Tri_Matrix<T, N>::data_ln;
    mutable T temp_return{};

public:
    using typename Base_Tri_Matrix<T, N>::size_type;
    using Base_Tri_Matrix<T, N>::Base_Tri_Matrix;

    Hermite_Matrix(const Hermite_Matrix &mat) = default;
    Hermite_Matrix(Hermite_Matrix &&mat) = default;
    // template <size_t M>
    // Hermite_Matrix(const Symm_Band_Matrix<T, N, M> &mat)
    //     : Base_Tri_Matrix<T, N>{mat} {}
    // template <size_t M>
    // Hermite_Matrix(Symm_Band_Matrix<T, N, M> &&mat)
    //     : Base_Tri_Matrix<T, N>{mat} {}

    Hermite_Matrix &operator=(const Hermite_Matrix &mat) = default;
    Hermite_Matrix &operator=(Hermite_Matrix &&mat) = default;
    // Hermite_Matrix &operator=(const Diag_Matrix<T, N> &mat)
    // {
    //     Base_Tri_Matrix<T, N>::operator=(mat);
    //     return *this;
    // }
    // Hermite_Matrix &operator=(Diag_Matrix<T, N> &&mat)
    // {
    //     Base_Tri_Matrix<T, N>::operator=(mat);
    //     return *this;
    // }
    // template <size_t M>
    // Hermite_Matrix &operator=(const Symm_Band_Matrix<T, N, M> &mat)
    // {
    //     Base_Tri_Matrix<T, N>::operator=(mat);
    //     return *this;
    // }
    // template <size_t M>
    // Hermite_Matrix &operator=(Symm_Band_Matrix<T, N, M> &&mat)
    // {
    //     Base_Tri_Matrix<T, N>::operator=(mat);
    //     return *this;
    // }

    T &operator()(size_type row, size_type col) override
    {
        if (row < col)
        {
            throw std::out_of_range((__FILE__ ":") + std::to_string(__LINE__) + ": trying to access conjugate area.");
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
            temp_return = std::conj(Base_Tri_Matrix<T, N>::operator()(col, row));
            return temp_return;
        }
        else
        {
            return Base_Tri_Matrix<T, N>::operator()(row, col);
        }
    }
};

} // namespace Misc

#endif // MISC_HERMITE_MATRIX
