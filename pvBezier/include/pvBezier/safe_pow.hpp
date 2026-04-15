#pragma once

#include <cmath>
#include <iostream>

inline static constexpr double powi(double base, uint32_t exp) {
    double result = 1.0;
    while (exp > 0) {
        if (exp & 1) result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

/**
 * @brief Computes the i-th power of t, with 0^0=1 and t^i=0 for i<0
 * @param t Value to raise.
 * @param i Exponent to raise to.
 * @return Returns t^i with some safety checks to avoid undefined values.
 */
inline static constexpr double safePow(double t, int i) noexcept
{
    if (i < 0)
        return 0;
    if (t == 0 && i == 0)
        return 1;
    return powi(t, i);
}