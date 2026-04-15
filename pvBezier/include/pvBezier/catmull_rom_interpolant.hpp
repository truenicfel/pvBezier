#pragma once

#include "bezier_curve.hpp"

/**
 * @brief Computes the tri-cubic C1-continuous interpolant of the central cell in a 4x4x4 grid.
 */
class CatmullRomInterpolant
{
public:
    template <typename TValue>
    using Cells = std::array<std::array<std::array<TValue, 4>, 4>, 4>;

    template <typename TValue>
    static void convert(const Cells<TValue>& grid, Cells<TValue>& bezier)
    {
        // Matrix M (Eq. 16)
        Eigen::Matrix4d transform;
        transform << 0, 1, 0, 0,
            -1. / 6., 1, 1. / 6., 0,
            0, 1. / 6., 1, -1. / 6.,
            0, 0, 1, 0;

        Cells<TValue> temp1, temp2;

        // Tensor Contraction (Eq. 17 and Additional Material Algorithm 2)
        // transform along x-axis
        for (int y = 0; y <= 3; ++y)
            for (int z = 0; z <= 3; ++z)
            {
                Eigen::Matrix<TValue, 4, 1> col;
                for (int x = 0; x <= 3; ++x)
                    col(x) = grid[x][y][z];
                Eigen::Matrix<TValue, 4, 1> result;
                for (int i = 0; i <= 3; ++i)
                {
                    result(i) = TValue::Zero();
                    for (int j = 0; j <= 3; ++j)
                        result(i) += transform(i, j) * col(j);
                }
                for (int x = 0; x <= 3; ++x)
                    temp1[x][y][z] = result(x);
            }

        // transform along y-axis
        for (int x = 0; x <= 3; ++x)
            for (int z = 0; z <= 3; ++z)
            {
                Eigen::Matrix<TValue, 4, 1> col;
                for (int y = 0; y <= 3; ++y)
                    col(y) = temp1[x][y][z];
                Eigen::Matrix<TValue, 4, 1> result;
                for (int i = 0; i <= 3; ++i)
                {
                    result(i) = TValue::Zero();
                    for (int j = 0; j <= 3; ++j)
                        result(i) += transform(i, j) * col(j);
                }
                for (int y = 0; y <= 3; ++y)
                    temp2[x][y][z] = result(y);
            }

        // transform along z-axis
        for (int x = 0; x <= 3; ++x)
            for (int y = 0; y <= 3; ++y)
            {
                Eigen::Matrix<TValue, 4, 1> col;
                for (int z = 0; z <= 3; ++z)
                    col(z) = temp2[x][y][z];
                Eigen::Matrix<TValue, 4, 1> result;
                for (int i = 0; i <= 3; ++i)
                {
                    result(i) = TValue::Zero();
                    for (int j = 0; j <= 3; ++j)
                        result(i) += transform(i, j) * col(j);
                }
                for (int z = 0; z <= 3; ++z)
                    bezier[x][y][z] = result(z);
            }
    }
};
