#include "pvBezier/bernstein_basis.hpp"

#include "pvBezier/binomial_coefficient.hpp"
#include "pvBezier/binomial_table.hpp"
#include "pvBezier/safe_pow.hpp"

#include <algorithm>
#include <cmath>


    static constexpr int factorial(int n) noexcept
    {
        if (n <= 1)
            return 1;
        return n * factorial(n - 1);
    }

    double BernsteinBasis::sample(double t, int i, int n) noexcept
    {
        //return BinomialCoefficient::binomial(n, i) * safePow(1 - t, n - i) * safePow(t, i);
        return BinomialTable<10>::get(n,i) * safePow(1 - t, n - i) * safePow(t, i);
        //return BinomialTable<10>::table[n][i] * safePow(1 - t, n - i) * safePow(t, i);
    }

    double BernsteinBasis::samplePartial(double t, int i, int n, uint8_t k) noexcept
    {
        // hard-coded versions for common cases
        switch (k)
        {
        case 1:
            return BinomialCoefficient::binomial(n, i) * (i * safePow(1 - t, n - i) * safePow(t, i - 1) - (n - i) * safePow(1 - t, n - i - 1) * safePow(t, i));
        case 2:
            return BinomialCoefficient::binomial(n, i) * ((n - i - 1) * (n - i) * safePow(1 - t, n - i - 2) * safePow(t, i) - 2 * i * (n - i) * safePow(1 - t, n - i - 1) * safePow(t, i - 1) + (i - 1) * i * safePow(1 - t, n - i) * safePow(t, i - 2));
        }

        // general formula for all derivatives
        double result = 0;
        for (int j = std::max(0, i + k - n); j <= std::min(i, (int)k); ++j)
        {
            result += safePow(-1, j + k) * BinomialCoefficient::binomial(k, j) * sample(t, i - j, n - k);
        }
        return factorial(n) / factorial(n - k) * result;
    }

