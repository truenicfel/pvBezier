#pragma once

/**
 * @brief Class that evaluates binomial coefficients.
 */
class BinomialCoefficient
{
public:
    /**
     * @brief Delete constructor.
     */
    BinomialCoefficient() = delete;

    /**
     * @brief Computes the binomial coefficient recursively at compile time. Returns the number of ways to select a sequence of i distinct objects, retaining the order of selection, from a set of n objects. Note that 0<=i<=n, otherwise the function returns 0.
     * @param n Number of objects to select from.
     * @param i Number of objects to select.
     * @return Binomial coefficient.
     */
    [[nodiscard]] static constexpr unsigned long long binomial(int n, int i) noexcept
    {

        if (i < 0 || i > n)
            return 0;
        if (i == 0 || i == n)
            return 1;

        // exploit symmetry: C(n, i) = C(n, n-i)
        if (i > n - i)
            i = n - i;

        unsigned long long result = 1;
        for (int j = 1; j <= i; ++j)
        {
            result = result * (n - j + 1) / j;
        }
        return result;
    }
};
