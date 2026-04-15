#pragma once

#include "Eigen/Eigen"

/**
 * @brief Represents a plane in 3D.
 */
class Plane3d
{
public:
    using VectorType = Eigen::Vector3d;

    /**
     * @brief Constructor from two vectors. The origin of the plane is always in (0, 0, 0).
     *
     * @param vector1 First vector defining the plane.
     * @param vector2 Second vector defining the plane.
     */
    Plane3d(const VectorType& vector1, const VectorType& vector2) noexcept
        : mNormal(vector1.cross(vector2).normalized())
    {
    }

    /**
     * @brief Constructor from two vectors. The origin of the plane is always in (0, 0, 0).
     *
     * @param normal The normal vector of the plane.
     */
    Plane3d(const VectorType& normal) noexcept
        : mNormal(normal.normalized())
    {
    }

    /**
     * @brief Computes the signed distance from the given point to the plane.
     *
     * @param point the point to compute the signed distance to.
     *
     * @return The signed distance.
     */
    [[nodiscard]] double distance(const VectorType& point) const noexcept
    {

        return mNormal.dot(point);
    }

    [[nodiscard]] const VectorType& getNormal() const noexcept
    {
        return mNormal;
    }

private:

    /**
     * @brief The normal vector of the plane.
     */
    VectorType mNormal;
};
