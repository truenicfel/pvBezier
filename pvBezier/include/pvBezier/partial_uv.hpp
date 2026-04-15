#pragma once

#include <array>
#include <cstdint>

/**
 * @brief Partial derivative specification for a parametric surface.
 */
union PartialUV
{
    /**
     * @brief Number of dimensions.
     */
    static const int Dimensions = 2;

    /**
     * @brief Constructor from hash.
     * @param _hash Hash that specifies the derivative degree for all dimensions.
     */
    PartialUV(uint64_t _hash)
        : hash(_hash)
    {
    }

    /**
     * @brief Constructor from an array.
     * @param _degree Array that specifies the derivative degree for all dimensions.
     */
    PartialUV(std::array<uint8_t, Dimensions> _degree)
        : hash(0)
    {
        for (int i = 0; i < Dimensions; ++i)
            degree[i] = _degree[i];
    }

    /**
     * @brief Constructor with individual partial derivative degrees.
     * @param _du Degree of the u-partial derivative.
     * @param _dv Degree of the v-partial derivative.
     */
    PartialUV(uint8_t _du, uint8_t _dv)
        : hash(0)
    {
        degreeDu = _du;
        degreeDv = _dv;
    }

    uint64_t hash; /* Hash that can be used for comparison and switch cases. */
    struct
    {
        uint8_t degreeDu; /* Specifies the degree of the u-partial derivative. */
        uint8_t degreeDv; /* Specifies the degree of the v-partial derivative. */
    };

    static constexpr uint64_t c    = 0;               /* Hash of the function (no partial). */
    static constexpr uint64_t du   = 1 << 0;          /* Hash of the du partial derivative. */
    static constexpr uint64_t dv   = 1 << 8;          /* Hash of the dv partial derivative. */
    static constexpr uint64_t duu  = 2 << 0;          /* Hash of the duu partial derivative. */
    static constexpr uint64_t duv  = 1 << 0 | 1 << 8; /* Hash of the duv partial derivative. */
    static constexpr uint64_t dvv  = 2 << 8;          /* Hash of the dvv partial derivative. */
    static constexpr uint64_t duuu = 3 << 0;          /* Hash of the duuu partial derivative. */
    static constexpr uint64_t duuv = 2 << 0 | 1 << 8; /* Hash of the duuv partial derivative. */
    static constexpr uint64_t duvv = 1 << 0 | 2 << 8; /* Hash of the duvv partial derivative. */
    static constexpr uint64_t dvvv = 3 << 8;          /* Hash of the dvvv partial derivative. */

    /**
     * @brief Gets the derivative degree for a given dimension.
     * @param dimension Dimension to get the derivative degree for.
     * @return Derivative degree of given dimension.
     */
    [[nodiscard]] uint8_t operator[](uint8_t dimension) const
    {
        return degree[dimension];
    }

    /**
     * @brief Gets the derivative degree for a given dimension.
     * @param dimension Dimension to get the derivative degree for.
     * @return Derivative degree of given dimension.
     */
    [[nodiscard]] uint8_t& operator[](uint8_t dimension)
    {
        return degree[dimension];
    }

private:
    /**
     * @brief Array of derivative degrees for up to 8 dimensions.
     */
    uint8_t degree[8];
};
