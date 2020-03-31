#ifndef MISC_LOW_BAND_MATRIX
#define MISC_LOW_BAND_MATRIX

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

    explicit Low_Band_Matrix(T deft = T{}) : Base_Half_Band_Matrix<T, N, M>{deft} {}
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
    Low_Band_Matrix(std::initializer_list<std::initializer_list<T>> ini)
        : Base_Half_Band_Matrix<T, N, M>{ini} {}

    Low_Band_Matrix &operator=(const Low_Band_Matrix &mat) = default;
    Low_Band_Matrix &operator=(Low_Band_Matrix &&mat) = default;
    template <size_t OLD_M>
    Low_Band_Matrix &operator=(const Low_Band_Matrix<T, N, OLD_M> &mat)
    {
        Base_Half_Band_Matrix<T, N, M>::operator=(mat);
        return *this;
    }
    template <size_t OLD_M>
    Low_Band_Matrix &operator=(Low_Band_Matrix<T, N, OLD_M> &&mat)
    {
        Base_Half_Band_Matrix<T, N, M>::operator=(mat);
        return *this;
    }

private:
    using Base_Half_Band_Matrix<T, N, M>::elem;
    using Base_Half_Band_Matrix<T, N, M>::data_ln;
    // std::pair<size_type, size_type> row_col(size_type idx) const;
};

} // namespace Misc

#endif // MISC_LOW_BAND_MATRIX
