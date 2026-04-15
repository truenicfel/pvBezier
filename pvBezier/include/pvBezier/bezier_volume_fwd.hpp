#pragma once

#include <Eigen/Eigen>

/**
 * @brief Class for tensor product volumes with Bernstein basis functions.
 * @tparam TValue Type of values stored at the control points.
 * @tparam TN Degree of the volume in u direction.
 * @tparam TM Degree of the volume in v direction.
 * @tparam TO Degree of the volume in w direction.
 */
template <int TN, int TM, int TO, typename TValue>
class BezierVolume;

/**
 * @brief Trilinear tensor product volume with Bernstein basis functions.
 * @tparam TValue Coordinate type of the volume.
 */
template <typename TValue>
using TrilinearBezierVolume = BezierVolume<1, 1, 1, TValue>;

/**
 * @brief Trilinear uni-variate tensor product volume with Bernstein basis functions.
 */
using TrilinearBezierVolume1d = TrilinearBezierVolume<Eigen::Vector<double, 1>>;

/**
 * @brief Trilinear bi-variate tensor product volume with Bernstein basis functions.
 */
using TrilinearBezierVolume2d = TrilinearBezierVolume<Eigen::Vector2d>;

/**
 * @brief Trilinear tri-variate tensor product volume with Bernstein basis functions.
 */
using TrilinearBezierVolume3d = TrilinearBezierVolume<Eigen::Vector3d>;

/**
 * @brief Trilinear tetra-variate tensor product volume with Bernstein basis functions.
 */
using TrilinearBezierVolume4d = TrilinearBezierVolume<Eigen::Vector4d>;

/**
 * @brief Triquadratic tensor product volume with Bernstein basis functions.
 * @tparam TValue Coordinate type of the volume.
 */
template <typename TValue>
using TriquadraticBezierVolume = BezierVolume<2, 2, 2, TValue>;

/**
 * @brief Triquadratic uni-variate tensor product volume with Bernstein basis functions.
 */
using TriquadraticBezierVolume1d = TriquadraticBezierVolume<Eigen::Vector<double, 1>>;

/**
 * @brief Triquadratic bi-variate tensor product volume with Bernstein basis functions.
 */
using TriquadraticBezierVolume2d = TriquadraticBezierVolume<Eigen::Vector2d>;

/**
 * @brief Triquadratic tri-variate tensor product volume with Bernstein basis functions.
 */
using TriquadraticBezierVolume3d = TriquadraticBezierVolume<Eigen::Vector3d>;

/**
 * @brief Triquadratic tetra-variate tensor product volume with Bernstein basis functions.
 */
using TriquadraticBezierVolume4d = TriquadraticBezierVolume<Eigen::Vector4d>;

/**
 * @brief Tricubic tensor product volume with Bernstein basis functions.
 * @tparam TValue Coordinate type of the volume.
 */
template <typename TValue>
using TricubicBezierVolume = BezierVolume<3, 3, 3, TValue>;

/**
 * @brief Tricubic uni-variate tensor product volume with Bernstein basis functions.
 */
using TricubicBezierVolume1d = TricubicBezierVolume<Eigen::Vector<double, 1>>;

/**
 * @brief Tricubic bi-variate tensor product volume with Bernstein basis functions.
 */
using TricubicBezierVolume2d = TricubicBezierVolume<Eigen::Vector2d>;

/**
 * @brief Tricubic tri-variate tensor product volume with Bernstein basis functions.
 */
using TricubicBezierVolume3d = TricubicBezierVolume<Eigen::Vector3d>;

/**
 * @brief Tricubic tetra-variate tensor product volume with Bernstein basis functions.
 */
using TricubicBezierVolume4d = TricubicBezierVolume<Eigen::Vector4d>;
