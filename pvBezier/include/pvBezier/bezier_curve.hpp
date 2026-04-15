#pragma once

#include "bernstein_basis.hpp"
#include "bezier_curve_fwd.hpp"
#include "partial_t.hpp"

/**
 * @brief Class that stores and manipulates bezier curves of general degree.
 * @tparam TDegree Degree of the curve.
 * @tparam TValue Type of values stored at the control points.
 */
template <int TDegree, typename TValue>
class BezierCurve
{
public:
    /**
     * @brief Degree of the curve.
     */
    static constexpr int Degree = TDegree;

    /**
     * @brief Type of coordinate that specifies the location in the domain.
     */
    using DomainCoord = double;

    /**
     * @brief Type of domain bounding box.
     */
    using BoundingBox = Eigen::AlignedBox<double, 1>;

    /**
     * @brief Type of the values returned by the function.
     */
    using Value = TValue;

    /**
     * @brief Type of the partial derivative specification.
     */
    using Partial = PartialT;

    /**
     * @brief Type of the control points.
     */
    using ControlPointsType = std::array<TValue, Degree + 1>;

    using Vector1d = Eigen::Vector<double, 1>;

    /**
     * @brief Constructor with zero initialization.
     */
    BezierCurve() noexcept
        : mDomain(Vector1d(0), Vector1d(1))
        , mBoundingBox()
        , mControlPoints()
    {
        for (int i = 0; i <= Degree; ++i)
            mControlPoints[i] = TValue::Zero();
        recomputeBoundingBox();
    }

    /**
     * @brief Constructor with initial weights.
     * @param initial Initial set of control points.
     */
    explicit BezierCurve(const std::array<TValue, Degree + 1>& initial) noexcept
        : mDomain(Vector1d(0), Vector1d(1))
        , mBoundingBox()
        , mControlPoints(initial)
    {
        recomputeBoundingBox();
    }

    /**
     * @brief Evaluates the parametric curve algebraically at a given t coordinate.
     * @param t Domain location where to sample the curve.
     * @return Position of the curve.
     */
    [[nodiscard]] TValue sample(const DomainCoord& t) const noexcept
    {
        double tmin = mDomain.min().x(), tmax = mDomain.max().x();
        double trel   = (t - tmin) / (tmax - tmin);
        TValue result = TValue::Zero();
#if 0
        // Direct evaluation of Bernstein Bezier functions
        for (int i = 0; i <= Degree; ++i)
            result += mControlPoints[i] * BernsteinBasis::sample(trel, i, Degree);
#else
        // Volk-Schumaker Evaluation
        // Handle edge case: t == 1.0
        double Bu[Degree + 1] = { 0 };
        if (std::abs(trel - 1.0) < std::numeric_limits<double>::epsilon())
            Bu[Degree] = 1.0;
        else
        {
            Bu[0] = std::pow(1.0 - trel, Degree);
            for (int i = 1; i <= Degree; ++i)
                Bu[i] = Bu[i - 1] * trel * (Degree - i + 1) / ((1.0 - trel) * i);
        }

        // Evaluate curve
        for (int i = 0; i <= Degree; ++i)
            result += mControlPoints[i] * Bu[i];
#endif
        return result;
    }

    /**
     * @brief Evaluates a partial derivative at a given domain coordinate.
     * @param t Domain location to sample the field derivative at.
     * @param partial Specifies the desired partial derivative of each dimension.
     * @return Partial derivative at domain location.
     */
    [[nodiscard]] Value samplePartial(const DomainCoord& t, const Partial& partial) const noexcept
    {
        double tmin = mDomain.min().x(), tmax = mDomain.max().x();
        double trel   = (t - tmin) / (tmax - tmin);
        TValue result = TValue::Zero();
        for (int i = 0; i <= Degree; ++i)
            result += mControlPoints[i] * BernsteinBasis::samplePartial(trel, i, Degree, partial.degreeDt);
        return result / std::pow(tmax - tmin, partial.degreeDt);
    }

    /**
     * @brief Subdivides a Bezier curve into two curves at a specified parameter value t.
     * @param output1 First part of the split curve.
     * @param output2 Second part of the split curve.
     * @param t Parameter location where to split.
     */
    void subdivide(
        BezierCurve& output1,
        BezierCurve& output2,
        DomainCoord t = 0.5) const noexcept
    {

        double tmin = mDomain.min().x(), tmax = mDomain.max().x();

        BezierCurve::subdivide<BezierCurve>(getControlPoints(), output1.mControlPoints, output2.mControlPoints, tmin, tmax, t);

        output1.recomputeBoundingBox();
        output2.recomputeBoundingBox();
        output1.setDomain(Eigen::AlignedBox1d(Vector1d(tmin), Vector1d(t)));
        output2.setDomain(Eigen::AlignedBox1d(Vector1d(t), Vector1d(tmax)));
    }

    /**
     * @brief Subdivides a Bezier curve into two curves (just the control points) at a specified parameter value t.
     *
     * @tparam TBezierCurve Type of the bezier curve (also determines the type of the control points)
     *
     * @param controlPoints The control poitns of the bezier curve to split.
     * @param controlPointsOut1 The control points of the curve "left of t".
     * @param controlPointsOut2 The control points of the curve "right of t".
     * @param t where to split the curve
     */
    template <typename TBezierCurve>
    static void subdivide(
        const TBezierCurve::ControlPointsType& controlPoints,
        TBezierCurve::ControlPointsType& controlPointsOut1,
        TBezierCurve::ControlPointsType& controlPointsOut2,
        const DomainCoord& tmin,
        const DomainCoord& tmax,
        const DomainCoord& t) noexcept
    {
        const double trel = (t - tmin) / (tmax - tmin);

        // start point of output1 is the same as the start point of this curve
        controlPointsOut1[0] = controlPoints[0];
        // end point of output2 is the same as the end point of this curve
        controlPointsOut2[TDegree] = controlPoints[TDegree];

        std::array<std::array<TValue, TDegree + 1>, 2> controlPointsBuffer{};
        controlPointsBuffer[0] = controlPoints;
        int write              = 1;

        // Example Degree 5:
        // 1 2 3 4 5
        // the outer loop iterates over the control points of the 2 new curves
        // since the first control point of output1 and the last control point of output2 are already set,
        // this loop will start at 1
        for (int n = 1; n <= TDegree; ++n)
        {
            TValue& left = controlPointsBuffer[1 - write][0];
            // Example Degree 5:
            // 0 1 2 3 4 (5-1)
            // 0 1 2 3 (5-2)
            // 0 1 2 (5-3)
            // 0 1 (5-4)
            // 0 (5-5)
            // this inner loop iterates the intermediate control points which are always on less than the last set of control points
            for (int i = 0; i <= TDegree - n; ++i)
            {
                const TValue& right           = controlPointsBuffer[1 - write][i + 1];
                controlPointsBuffer[write][i] = (1 - trel) * left + trel * right;
                left                          = right;
            }
            // the next batch of intermediate control points is written so we can add them to the final curves
            // for output1 it is always the leftmost one in the current write buffer
            // Example Degree 5:
            // 1
            // 2
            // 3
            // 4
            // 5
            controlPointsOut1[n] = controlPointsBuffer[write][0];
            // for output2 it is always the rightmost one in the current write buffer
            // Example Degree 5:
            // 5-1 = 4
            // 5-2 = 3
            // 5-3 = 2
            // 5-4 = 1
            // 5-5 = 0
            controlPointsOut2[TDegree - n] = controlPointsBuffer[write][TDegree - n];
            // switch write buffer
            write = 1 - write;
        }
    }

    /**
     * @brief Clips away the part of the curve that is left of a specified parameter value t.
     * @param output Right part of the split curve.
     * @param t Parameter location where to split.
     */
    void clipLeft(
        BezierCurve& output,
        DomainCoord t = 0.5) const noexcept
    {
        double tmin = mDomain.min().x(), tmax = mDomain.max().x();
        double trel = (t - tmin) / (tmax - tmin);
        for (int n = 0; n <= TDegree; ++n)
        {
            output.mControlPoints[n] = intermediateControlPoint(n, TDegree - n, trel);
        }
        output.recomputeBoundingBox();
        output.setDomain(Eigen::AlignedBox1d(Vector1d(t), Vector1d(tmax)));
    }

    /**
     * @brief Clips away the part of the curve that is right of a specified parameter value t.
     * @param output Left part of the split curve.
     * @param t Parameter location where to split.
     */
    void clipRight(
        BezierCurve& output,
        DomainCoord t = 0.5) const noexcept
    {
        double tmin = mDomain.min().x(), tmax = mDomain.max().x();
        double trel = (t - tmin) / (tmax - tmin);
        for (int n = 0; n <= TDegree; ++n)
        {
            output.mControlPoints[n] = intermediateControlPoint(0, n, trel);
        }
        output.recomputeBoundingBox();
        output.setDomain(Eigen::AlignedBox1d(Vector1d(tmin), Vector1d(t)));
    }

    /**
     * @brief Recomputes the bounding box of the Bezier curve.
     */
    void recomputeBoundingBox() noexcept
    {
        mBoundingBox.setEmpty();
        for (int i = 0; i <= TDegree; ++i)
            mBoundingBox.extend(mControlPoints[i]);
    }

    /**
     * @brief Gets the bounding box of the vertices. Note that recomputeBoundingBox has to be called first!
     * @return Bounding box of this geometry.
     */
    [[nodiscard]] const Eigen::AlignedBox<double, TValue::scalar_count>& getBoundingBox() const noexcept { return mBoundingBox; }

    /**
     * @brief Gets the array of Degree+1 control points.
     * @return Control points.
     */
    [[nodiscard]] const std::array<TValue, Degree + 1>& getControlPoints() const noexcept { return mControlPoints; }

    /**
     * @brief Gets the array of Degree+1 control points.
     * @return Control points.
     */
    [[nodiscard]] std::array<TValue, Degree + 1>& getControlPoints() noexcept { return mControlPoints; }

    /**
     * @brief Sets the array of Degree+1 control points.
     * @param controlPoints Control points.
     * @param computeBoundingBox Optional flag that recomputes the bounding box according to the new control points.
     */
    void setControlPoints(const std::array<TValue, Degree + 1>& controlPoints, bool computeBoundingBox = true) noexcept
    {
        mControlPoints = controlPoints;
        if (computeBoundingBox)
            recomputeBoundingBox();
    }

    /**
     * @brief Sets the parameter domain over which the function is defined.
     * @param domain Parameter domain.
     */
    virtual void setDomain(const BoundingBox& domain)
    {
        mDomain = domain;
    }

private:
    /**
     * @brief Calculates the intermediate control point of the de Casteljau algorithm.
     * @param i Index of the control point.
     * @param r Degree of the intermediate point.
     * @param t Parameter location where to evaluate.
     * @return Position of the intermediate control point.
     */
    [[nodiscard]] TValue intermediateControlPoint(int i, int r, DomainCoord t) const noexcept
    {
        TValue value = TValue::Zero();
        for (int j = 0; j <= r; ++j)
            value += BernsteinBasis::sample(t, j, r) * mControlPoints[i + j];
        return value;
    }

    /**
     * @brief Parameter domain over which the function is defined.
     */
    BoundingBox mDomain;

    /**
     * @brief Bounding box of the curve.
     */
    Eigen::AlignedBox<double, TValue::scalar_count> mBoundingBox;

    /**
     * @brief Array of Degree+1 control points.
     */
    ControlPointsType mControlPoints;
};
