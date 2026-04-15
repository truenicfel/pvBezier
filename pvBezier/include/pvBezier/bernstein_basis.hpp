#pragma once

#include <cstdint>

/**
 * @brief Class that evaluates Bernstein basis polynomials.
 */
class BernsteinBasis
{
public:
    /**
     * @brief Delete constructor.
     */
    BernsteinBasis() = delete;

    /**
     * @brief Evaluates the Bernstein basis polynomial B_i^n(t) = (n i) * (1-t)^(n-i) * t^i
     * @param t Parameter value in [0,1]
     * @param i Index of the basis.
     * @param n Number of basis functions.
     * @return Value of the basis function.
     */
    [[nodiscard]] static double sample(double t, int i, int n) noexcept;

    /**
     * @brief Evaluates the kth-order derivative of the Bernstein polynomial.
     * @param t Parameter value in [0,1]
     * @param i Index of the basis.
     * @param n Number of basis functions.
     * @param k Degree of the derivative to compute.
     * @return kth-order derivative of the basis function.
     */
    [[nodiscard]] static double samplePartial(double t, int i, int n, uint8_t k) noexcept;
};
