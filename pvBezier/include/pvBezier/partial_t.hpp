#pragma once

/**
 * @brief Partial derivative specification for a parametric curve.
 */
union PartialT
{
    /**
     * @brief Number of dimensions.
     */
    static const int Dimensions = 1;

    /**
     * @brief Constructor from hash.
     * @param _hash Hash that specifies the derivative degree for all dimensions.
     */
    PartialT(uint64_t _hash)
        : hash(_hash)
    {
    }

    /**
     * @brief Constructor from an array.
     * @param _degree Array that specifies the derivative degree for all dimensions.
     */
    PartialT(std::array<uint8_t, Dimensions> _degree)
        : hash(0)
    {
        for (int i = 0; i < Dimensions; ++i)
            degree[i] = _degree[i];
    }

    /**
     * @brief Constructor with individual partial derivative degrees.
     * @param _dt Degree of the t-partial derivative.
     */
    PartialT(uint8_t _dt)
        : hash(0)
    {
        degreeDt = _dt;
    }

    uint64_t hash; /* Hash that can be used for comparison and switch cases. */
    struct
    {
        uint8_t degreeDt; /* Specifies the degree of the t-partial derivative. */
    };

    static constexpr uint64_t c    = 0;      /* Hash of the function (no partial). */
    static constexpr uint64_t dt   = 1 << 0; /* Hash of the dt partial derivative. */
    static constexpr uint64_t dtt  = 2 << 0; /* Hash of the dtt partial derivative. */
    static constexpr uint64_t dttt = 3 << 0; /* Hash of the dttt partial derivative. */

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
