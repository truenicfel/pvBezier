#pragma once

#include "binomial_coefficient.hpp"

#include <array>

/**
 * @brief Class that stores binomial coefficients up to a certain value.
 *
 * @tparam TMax Maximum value for which the binomial coefficients are stored.
 */
template <int TMax>
class BinomialTable
{
public:
    using TableType = std::array<std::array<unsigned long long, TMax + 1>, TMax + 1>;

    /**
     * @brief Delete constructor.
     */
    BinomialTable() = delete;

    static constexpr std::array<std::array<unsigned long long, TMax + 1>, TMax + 1> table = []
    {
        std::array<std::array<unsigned long long, TMax + 1>, TMax + 1> t{};

        for (int n = 0; n <= TMax; ++n)
        {
            t[n][0] = t[n][n] = 1;
            for (int k = 1; k < n; ++k)
            {
                t[n][k] = t[n - 1][k - 1] + t[n - 1][k];
            }
        }
        return t;
    }();

    /**
     * @brief Gets the binomial coefficient from the table above or falls back to runtime method. Prefer to use the table directly.
     *
     * @param n Number of objects to select from.
     * @param i Number of objects to select.
     *
     * @return Binomial coefficient.
     */
    static constexpr unsigned long long get(int n, int i)
    {
        if (n <= TMax && i >= 0 && i <= n)
        {
            return table[n][i];
        }
        else
        {
            // fallback to runtime method
            return BinomialCoefficient::binomial(n, i);
        }
    }
};
