#ifndef MISC_POLYNOMIAL
#define MISC_POLYNOMIAL

#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <iterator>
#include <initializer_list>

#include "linear_eq_direct.h"

namespace Misc
{

template <typename T>
class Polynomial : public std::vector<T>
{
public:
    // using std::vector<T>::size;
    // using typename std::vector<T>::value_type;
    // using std::vector<T>::operator[];

    explicit Polynomial(std::initializer_list<T> ini)
        : std::vector<T>(std::rbegin(ini), std::rend(ini))
    {
        cut_zero();
    }
    template <typename... Ts>
    explicit Polynomial(Ts... args, ...)
        : std::vector<T>(args...)
    {
        cut_zero();
    }
    // Polynomial(const Polynomial &other) = default;
    // Polynomial(Polynomial &&other) = default;

    T operator()(T x) const
    {
        T res{0.};
        for (auto it = this->crbegin(); it != this->crend(); it++)
        {
            res *= x;
            res += *it;
        }
        return res;
    }

    Polynomial deriv() const
    {
        Polynomial res{};
        for (size_t i = 1; i < this->size(); i++)
        {
            res.push_back((*this)[i] * i);
        }
        return res;
    }
    Polynomial deriv(size_t m) const
    {
        Polynomial res{};
        size_t coefficient{1};
        for (size_t i = 1; i < m; i++)
        {
            coefficient *= i;
        }
        for (size_t i = m; i < this->size(); i++)
        {
            coefficient *= i;
            res.push_back((*this)[i] * coefficient);
            coefficient /= i - m + 1;
        }
        return res;
    }

    Polynomial integ() const
    {
        std::vector<T> res{0.};
        for (size_t i = 0; i < this->size(); i++)
        {
            res.push_back((*this)[i] / (i + 1));
        }
        return Polynomial(std::move(res));
    }
    Polynomial integ(size_t m) const
    {
        std::vector<T> res(m, 0.);
        size_t coefficient{1};
        for (size_t i = 1; i < m; i++)
        {
            coefficient *= i;
        }
        for (size_t i = 0; i < this->size(); i++)
        {
            coefficient *= i + m;
            res.push_back((*this)[i] / coefficient);
            coefficient /= i + 1;
        }
        return Polynomial(std::move(res));
    }

    Polynomial &operator+=(const Polynomial<T> &other)
    {
        if (other.size() > this->size())
        {
            this->resize(other.size(), 0.);
        }
        for (size_t i = 0; i < other.size(); i++)
        {
            (*this)[i] += other[i];
        }
        cut_zero();
        return *this;
    }

    Polynomial &operator-=(const Polynomial<T> &other)
    {
        if (other.size() > this->size())
        {
            this->resize(other.size(), 0.);
        }
        for (size_t i = 0; i < other.size(); i++)
        {
            (*this)[i] -= other[i];
        }
        cut_zero();
        return *this;
    }

    std::pair<Polynomial, Polynomial> div(const Polynomial<T> &other) const
    {
        Polynomial remain{*this};
        if (this->size() < other.size())
        {
            return std::make_pair(Polynomial{}, remain);
        }
        else
        {
            size_t res_size{remain.size() - other.size()};
            std::vector<T> res(this->size() - other.size() + 1);
            T divider{other.back()};
            for (size_t i = 0; i < res.size(); i++)
            {
                static T alpha;
                static size_t remain_start;
                res[res.size() - i - 1] = remain[remain.size() - i - 1] / divider;
                alpha = res[res.size() - i - 1];
                remain_start = res_size - i;
                for (size_t j = 0; j < other.size() - 1; j++)
                {
                    remain[remain_start + j] -= alpha * other[j];
                }
            }
            return std::make_pair(res, Polynomial(remain.begin(), remain.begin() + (other.size() - 1)));
        }
    }

    std::vector<T> coef() const
    {
        return static_cast<std::vector<T>>(*this);
    }

private:
    void cut_zero()
    {
        for (size_t i = this->size() - 1; i != -1 && (*this)[i] == 0.; i--)
        {
            this->pop_back();
        }
    }
};

// ================= binary operator ============================

template <typename T>
inline Polynomial<T> operator+(const Polynomial<T> &a, const Polynomial<T> &b)
{
    Polynomial<T> res{a};
    res += b;
    return res;
}
template <typename T>
inline Polynomial<T> operator+(const Polynomial<T> &a, const T &b)
{
    return a + Polynomial<T>{b};
}
template <typename T>
inline Polynomial<T> operator+(const T &a, const Polynomial<T> &b)
{
    return Polynomial<T>{a} + b;
}

template <typename T>
inline Polynomial<T> operator-(const Polynomial<T> &a, const Polynomial<T> &b)
{
    Polynomial<T> res{a};
    res -= b;
    return res;
}
template <typename T>
inline Polynomial<T> operator-(const Polynomial<T> &a, const T &b)
{
    return a - Polynomial<T>{b};
}
template <typename T>
inline Polynomial<T> operator-(const T &a, const Polynomial<T> &b)
{
    return Polynomial<T>{a} - b;
}

template <typename T>
inline Polynomial<T> operator-(const Polynomial<T> &a)
{
    Polynomial<T> res{};
    res -= a;
    return res;
}

template <typename T>
inline Polynomial<T> operator*(const Polynomial<T> &a, const Polynomial<T> &b)
{
    if (a.size() and b.size())
    {
        std::vector<T> res(a.size() + b.size() - 1);
        for (size_t i = 0; i < a.size(); i++)
        {
            for (size_t j = 0; j < b.size(); j++)
            {
                res[i + j] += a[i] * b[j];
            }
        }
        return Polynomial<T>(std::move(res));
    }
    else
    {
        return Polynomial<T>{};
    }
}
template <typename T>
inline Polynomial<T> operator*(const Polynomial<T> &a, const T &b)
{
    return a * Polynomial<T>{b};
}
template <typename T>
inline Polynomial<T> operator*(const T &a, const Polynomial<T> &b)
{
    return Polynomial<T>{a} * b;
}

template <typename T>
inline Polynomial<T> operator/(const Polynomial<T> &a, const Polynomial<T> &b)
{
    return a.div(b).first;
}
template <typename T>
inline Polynomial<T> operator/(const Polynomial<T> &a, const T &b)
{
    return a / Polynomial<T>{b};
}
template <typename T>
inline Polynomial<T> operator/(const T &a, const Polynomial<T> &b)
{
    return Polynomial<T>{a} / b;
}

template <typename T>
inline Polynomial<T> operator%(const Polynomial<T> &a, const Polynomial<T> &b)
{
    return a.div(b).second;
}
template <typename T>
inline Polynomial<T> operator%(const Polynomial<T> &a, const T &b)
{
    return a % Polynomial<T>{b};
}
template <typename T>
inline Polynomial<T> operator%(const T &a, const Polynomial<T> &b)
{
    return Polynomial<T>{a} % b;
}

// ========================= IO ============================

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Polynomial<T> &a)
{
    os << "Polynomial(";
    if (a.size())
        for (auto it = a.crbegin(); it != a.crend(); it++)
        {
            os << *it << ", ";
        }
    else
    {
        os << 0. << ", ";
    }
    os << ")";
    return os;
}

} // namespace Misc

#endif // MISC_POLYNOMIAL
