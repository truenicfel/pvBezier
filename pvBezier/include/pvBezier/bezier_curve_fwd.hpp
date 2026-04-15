#pragma once

#include <Eigen/Eigen>

/**
 * @brief Class that stores and manipulates bezier curves of general degree.
 * @tparam TDegree Degree of the curve.
 * @tparam TValue Type of values stored at the control points.
 */
template <int TDegree, typename TValue>
class BezierCurve;

/**
 * @brief Linear bezier curve.
 * @tparam TValue Coordinate type of the curve.
 */
template <typename TValue>
using LinearBezierCurve = BezierCurve<1, TValue>;

/**
 * @brief Linear uni-variate polynomial with Bernstein basis functions.
 */
using LinearBezierCurve1d = LinearBezierCurve<Eigen::Vector<double, 1>>;

/**
 * @brief Linear bi-variate polynomial with Bernstein basis functions.
 */
using LinearBezierCurve2d = LinearBezierCurve<Eigen::Vector2d>;

/**
 * @brief Linear tri-variate polynomial with Bernstein basis functions.
 */
using LinearBezierCurve3d = LinearBezierCurve<Eigen::Vector3d>;

/**
 * @brief Linear tetra-variate polynomial with Bernstein basis functions.
 */
using LinearBezierCurve4d = LinearBezierCurve<Eigen::Vector4d>;

/**
 * @brief Quadratic bezier curve.
 * @tparam TValue Coordinate type of the curve.
 */
template <typename TValue>
using QuadraticBezierCurve = BezierCurve<2, TValue>;

/**
 * @brief Quadratic uni-variate polynomial with Bernstein basis functions.
 */
using QuadraticBezierCurve1d = QuadraticBezierCurve<Eigen::Vector<double, 1>>;

/**
 * @brief Quadratic bi-variate polynomial with Bernstein basis functions.
 */
using QuadraticBezierCurve2d = QuadraticBezierCurve<Eigen::Vector2d>;

/**
 * @brief Quadratic tri-variate polynomial with Bernstein basis functions.
 */
using QuadraticBezierCurve3d = QuadraticBezierCurve<Eigen::Vector3d>;

/**
 * @brief Quadratic tetra-variate polynomial with Bernstein basis functions.
 */
using QuadraticBezierCurve4d = QuadraticBezierCurve<Eigen::Vector4d>;

/**
 * @brief Cubic bezier curve.
 * @tparam TValue Coordinate type of the curve.
 */
template <typename TValue>
using CubicBezierCurve = BezierCurve<3, TValue>;

/**
 * @brief Cubic uni-variate polynomial with Bernstein basis functions.
 */
using CubicBezierCurve1d = CubicBezierCurve<Eigen::Vector<double, 1>>;

/**
 * @brief Cubic bi-variate polynomial with Bernstein basis functions.
 */
using CubicBezierCurve2d = CubicBezierCurve<Eigen::Vector2d>;

/**
 * @brief Cubic tri-variate polynomial with Bernstein basis functions.
 */
using CubicBezierCurve3d = CubicBezierCurve<Eigen::Vector3d>;

/**
 * @brief Cubic tetra-variate polynomial with Bernstein basis functions.
 */
using CubicBezierCurve4d = CubicBezierCurve<Eigen::Vector4d>;
