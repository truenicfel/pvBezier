#pragma once

#include <Eigen/Eigen>

/**
 * @brief Class for tensor product surfaces with Bernstein basis functions.
 * @tparam TValue Type of values stored at the control points.
 * @tparam TN Degree of the surface in u direction.
 * @tparam TM Degree of the surface in v direction.
 */
template <int TN, int TM, typename TValue>
class BezierSurface;

/**
 * @brief Biquadratic tensor product surface with Bernstein basis functions.
 * @tparam TValue Coordinate type of the curve.
 */
template <typename TValue>
using BilinearBezierSurface = BezierSurface<1, 1, TValue>;

/**
 * @brief Biilinear uni-variate tensor product surface with Bernstein basis functions.
 */
using BilinearBezierSurface1d = BilinearBezierSurface<Eigen::Vector<double, 1>>;

/**
 * @brief Bilinear bi-variate tensor product surface with Bernstein basis functions.
 */
using BilinearBezierSurface2d = BilinearBezierSurface<Eigen::Vector2d>;

/**
 * @brief Bilinear tri-variate tensor product surface with Bernstein basis functions.
 */
using BilinearBezierSurface3d = BilinearBezierSurface<Eigen::Vector3d>;

/**
 * @brief Bilinear tetra-variate tensor product surface with Bernstein basis functions.
 */
using BilinearBezierSurface4d = BilinearBezierSurface<Eigen::Vector4d>;

/**
 * @brief Biquadratic tensor product surface with Bernstein basis functions.
 * @tparam TValue Coordinate type of the curve.
 */
template <typename TValue>
using BiquadraticBezierSurface = BezierSurface<2, 2, TValue>;

/**
 * @brief Biquadratic uni-variate tensor product surface with Bernstein basis functions.
 */
using BiquadraticBezierSurface1d = BiquadraticBezierSurface<Eigen::Vector<double, 1>>;

/**
 * @brief Biquadratic bi-variate tensor product surface with Bernstein basis functions.
 */
using BiquadraticBezierSurface2d = BiquadraticBezierSurface<Eigen::Vector2d>;

/**
 * @brief Biquadratic tri-variate tensor product surface with Bernstein basis functions.
 */
using BiquadraticBezierSurface3d = BiquadraticBezierSurface<Eigen::Vector3d>;

/**
 * @brief Biquadratic tetra-variate tensor product surface with Bernstein basis functions.
 */
using BiquadraticBezierSurface4d = BiquadraticBezierSurface<Eigen::Vector4d>;

/**
 * @brief Bicubic tensor product surface with Bernstein basis functions.
 * @tparam TValue Coordinate type of the curve.
 */
template <typename TValue>
using BicubicBezierSurface = BezierSurface<3, 3, TValue>;

/**
 * @brief Bicubic uni-variate tensor product surface with Bernstein basis functions.
 */
using BicubicBezierSurface1d = BicubicBezierSurface<Eigen::Vector<double, 1>>;

/**
 * @brief Bicubic bi-variate tensor product surface with Bernstein basis functions.
 */
using BicubicBezierSurface2d = BicubicBezierSurface<Eigen::Vector2d>;

/**
 * @brief Bicubic tri-variate tensor product surface with Bernstein basis functions.
 */
using BicubicBezierSurface3d = BicubicBezierSurface<Eigen::Vector3d>;

/**
 * @brief Bicubic tetra-variate tensor product surface with Bernstein basis functions.
 */
using BicubicBezierSurface4d = BicubicBezierSurface<Eigen::Vector4d>;
