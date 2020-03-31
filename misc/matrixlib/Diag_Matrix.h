#ifndef MISC_DIAG_MATRIX
#define MISC_DIAG_MATRIX

#include "Base_Matrix.h"
#include "Matrix.h"

// the namespace miscellany
namespace Misc
{

// Diagonal Matrix; only diagonal items are stored
template <typename T, size_t N>
class Diag_Matrix : public Base_Matrix<T, N, N>
{
private:
    using Base_Matrix<T, N, N>::elem;
    using Base_Matrix<T, N, N>::data_ln;

public:
    using typename Base_Matrix<T, N, N>::size_type;

    explicit Diag_Matrix(T deft = T{}) : Base_Matrix<T, N, N>{1, deft} {}
    Diag_Matrix(const Diag_Matrix &mat) = default;
    Diag_Matrix(Diag_Matrix &&mat) = default;
    template <typename It>
    explicit Diag_Matrix(It b, It e)
        : Base_Matrix<T, N, N>(1, nullptr)
    {
        if (e - b != N)
        {
            throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": wrong iterators.");
        }
        elem = new T[1][N];
        std::copy(b, e, elem[0]);
    }
    explicit Diag_Matrix(std::initializer_list<T> ini)
        : Base_Matrix<T, N, N>(1, nullptr)
    {
        if (ini.size() != N)
        {
            throw std::runtime_error((__FILE__ ":") + std::to_string(__LINE__) + ": wrong iterators.");
        }
        elem = new T[1][N];
        std::move(ini.begin(), ini.end(), elem[0]);
    }

    Diag_Matrix &operator=(const Diag_Matrix &mat) = default;
    Diag_Matrix &operator=(Diag_Matrix &&mat) = default;

    T &operator()(size_type row, size_type col) override
    {
        if (row == col)
        {
            return (*this)(row);
        }
        else
        {
            throw std::out_of_range((__FILE__ ":") + std::to_string(__LINE__) + ": trying to access 0.");
        }
    }
    T &operator()(size_type n) { return elem[0][n]; }
    const T &operator()(size_type row, size_type col) const override
    {
        if (row == col)
        {
            return (*this)(row);
        }
        else
        {
            return Base_Matrix<T, N, N>::zero;
        }
    }
    const T &operator()(size_type n) const { return elem[0][n]; }
};

// -------------------------------------------------------------------------

template <typename T, size_t N>
inline Diag_Matrix<T, N> operator*(const Diag_Matrix<T, N> &a, const Diag_Matrix<T, N> &b)
{
    Diag_Matrix<T, N> res{};
    for (std::size_t i = 0; i < N; i++)
    {
        res(i) = a(i) * b(i);
    }
    return res;
}

template <typename T, size_t N>
inline Matrix<T, N, N> operator*(const Diag_Matrix<T, N> &a, const Base_Matrix<T, N, N> &b)
{
    Matrix<T, N, N> res{};
    for (std::size_t i = 0; i < N; i++)
    {
        T temp{a(i)};
        for (std::size_t j = 0; j < N; j++)
        {
            res(i, j) = temp * b(i, j);
        }
    }
    return res;
}

template <typename T, size_t N>
inline Matrix<T, N, N> operator*(const Base_Matrix<T, N, N> &a, const Diag_Matrix<T, N> &b)
{
    Matrix<T, N, N> res{};
    for (std::size_t i = 0; i < N; i++)
        for (std::size_t j = 0; j < N; j++)
        {
            res(i, j) = a(i, j) * b(j);
        }
    return res;
}

} // namespace Misc

#endif // MISC_DIAG_MATRIX
