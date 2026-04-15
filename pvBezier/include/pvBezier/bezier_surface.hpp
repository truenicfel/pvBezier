#pragma once

#include "bernstein_basis.hpp"
#include "bezier_curve.hpp"
#include "bezier_surface_fwd.hpp"
#include "partial_uv.hpp"

/**
 * @brief Class for tensor product surfaces with Bernstein basis functions.
 * @tparam TValue Type of values stored at the control points.
 * @tparam TN Degree of the surface in u direction.
 * @tparam TM Degree of the surface in v direction.
 */
template <int TN, int TM, typename TValue>
class BezierSurface
{
public:
    /**
     * @brief Degree of the curve in u direction.
     */
    static constexpr int N = TN;

    /**
     * @brief Degree of the curve in v direction.
     */
    static constexpr int M = TM;

    /**
     * @brief Type of coordinate that specifies the location in the domain.
     */
    using DomainCoord = Eigen::Vector2d;

    /**
     * @brief Type of domain bounding box.
     */
    using BoundingBox = Eigen::AlignedBox<double, 2>;

    /**
     * @brief Type of the values returned by the function.
     */
    using Value = TValue;

    /**
     * @brief Type of the partial derivative specification.
     */
    using Partial = PartialUV;

    /**
     * @brief Type of the control points.
     */
    using ControlPointsType = std::array<std::array<TValue, M + 1>, N + 1>;

    /**
     * @brief Constructor with zero initialization.
     */
    BezierSurface() noexcept
        : mDomain(Eigen::Vector2d::Zero(), Eigen::Vector2d::Ones())
        , mBoundingBox()
        , mControlPoints()
    {
        for (int j = 0; j <= M; ++j)
            for (int i = 0; i <= N; ++i)
                mControlPoints[i][j] = TValue::Zero();
        recomputeBoundingBox();
    }

    /**
     * @brief Constructor with initial weights.
     * @param initial Initial set of control points.
     */
    explicit BezierSurface(const std::array<std::array<TValue, M + 1>, N + 1>& initial) noexcept
        : mDomain(Eigen::Vector2d::Zero(), Eigen::Vector2d::Ones())
        , mBoundingBox()
        , mControlPoints(initial)
    {
        recomputeBoundingBox();
    }

    /**
     * @brief Evaluates the tensor-product surface algebraically at a given uv coordinate.
     * @param uv Domain location where to evaluate the surface.
     * @return Position on the surface.
     */
    [[nodiscard]] TValue sample(const DomainCoord& uv) const noexcept
    {
        DomainCoord trel = (uv - mDomain.min()).cwiseQuotient(mDomain.max() - mDomain.min());
        const double u   = trel.x();
        const double v   = trel.y();
        TValue result    = TValue::Zero();
#if 0
            // Direct evaluation of Bernstein Bezier functions
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += mControlPoints[i][j] * BernsteinBasis::sample(u, i, N) * BernsteinBasis::sample(v, j, M);
#else
        // Volk-Schumaker Tensor Product Evaluation
        // Precompute basis for u and v using recurrence
        double Bu[N + 1] = { 0 };
        double Bv[M + 1] = { 0 };
        if (std::abs(u - 1.0) < std::numeric_limits<double>::epsilon())
            Bu[N] = 1.0;
        else
        {
            Bu[0] = std::pow(1 - u, N);
            for (int i = 1; i <= N; ++i)
                Bu[i] = Bu[i - 1] * u * (N - i + 1) / ((1 - u) * i);
        }

        if (std::abs(v - 1.0) < std::numeric_limits<double>::epsilon())
            Bv[M] = 1.0;
        else
        {
            Bv[0] = std::pow(1 - v, M);
            for (int j = 1; j <= M; ++j)
                Bv[j] = Bv[j - 1] * v * (M - j + 1) / ((1 - v) * j);
        }

        // Evaluate surface
        for (int j = 0; j <= M; ++j)
            for (int i = 0; i <= N; ++i)
                result += mControlPoints[i][j] * Bu[i] * Bv[j];
#endif
        return result;
    }

    /**
     * @brief Evaluates a partial derivative at a given domain coordinate.
     * @param uv Domain location to sample the derivative at.
     * @param partial Specifies the desired partial derivative of each dimension.
     * @return Partial derivative at domain location.
     */
    [[nodiscard]] Value samplePartial(const DomainCoord& uv, const Partial& partial) const noexcept
    {
        DomainCoord trel = (uv - mDomain.min()).cwiseQuotient(mDomain.max() - mDomain.min());
        TValue result    = TValue::Zero();
        for (int j = 0; j <= M; ++j)
            for (int i = 0; i <= N; ++i)
                result += mControlPoints[i][j] * BernsteinBasis::samplePartial(trel.x(), i, N, partial.degreeDu) * BernsteinBasis::samplePartial(trel.y(), j, M, partial.degreeDv);
        return result;
    }

    /**
     * @brief Generates a tensor product with degree [N, N] for the specified partial. Only enabled if both N and M are
     * equal and only works for non-mixed order on partials.
     *
     * @details Usually the tensor product will have its degree reduced by one for the dimension the partial is computed for
     * but after computing the partial that dimension is lifted to the degree N.
     *
     * @tparam DegreeU unused except for compile-time checks.
     * @tparam DegreeV unused except for compile-time checks.
     *
     * @param partial The partial to generate the tensor product for.
     *
     * @return A surface containing the partial tensor product.
     */
    template <int DegreeU = N, int DegreeV = M, typename = std::enable_if_t<(DegreeU == DegreeV)>>
    BezierSurface generatePartialTensorProduct(const PartialUV& partial) const
    {
        if (partial.hash != PartialUV::du && partial.hash != PartialUV::dv)
        {
            throw std::runtime_error("Partials with order larger than one are not supported yet!");
        }

        static constexpr int Degree = N;

        // compute the degrees for the individual dimensions of the partial tensor product.
        // when there is a partial for the given dimension, the degree reduces by one.
        // note that this only works as long as we only use first order and non-mixed partials.
        const int n1 = Degree - partial.degreeDu;
        const int n2 = Degree - partial.degreeDv;

        // this will use the control points type even though one of the dimensions is reduced by one
        ControlPointsType rawResult;

        // iterate all dimensions
        for (int i = 0; i <= n1; ++i)
        {
            for (int j = 0; j <= n2; ++j)
            {
                // increment the running index by one for the dimension we are currently computing the partial for
                // again note that this only works as long as we only use first order and non-mixed partials
                rawResult[i][j] = Degree * (mControlPoints[i + partial.degreeDu][j + partial.degreeDv] - mControlPoints[i][j]);
            }
        }

        // perform degree elevation
        // this result will be filled
        ControlPointsType result;

        // for u partial
        if (partial.hash == PartialUV::du)
        {
            for (int i = 0; i <= Degree; ++i)
            {
                for (int j = 0; j <= Degree; ++j)
                {
                    if (i == 0)
                    {
                        result[i][j] = rawResult[0][j] * (Degree - i) / Degree;
                    }
                    else if (i == Degree)
                    {
                        result[i][j] = rawResult[Degree - 1][j] * i / Degree;
                    }
                    else
                    {
                        result[i][j] = rawResult[i][j] * (Degree - i) / Degree + rawResult[i - 1][j] * i / Degree;
                    }
                }
            }
        }
        else
        {
            for (int i = 0; i <= Degree; ++i)
            {
                for (int j = 0; j <= Degree; ++j)
                {
                    if (j == 0)
                    {
                        result[i][j] = rawResult[i][0] * (Degree - j) / Degree;
                    }
                    else if (j == Degree)
                    {
                        result[i][j] = rawResult[i][Degree - 1] * j / Degree;
                    }
                    else
                    {
                        result[i][j] = rawResult[i][j] * (Degree - j) / Degree + rawResult[i][j - 1] * j / Degree;
                    }
                }
            }
        }

        return BezierSurface(result);
    }

    /**
     * @brief Subdivides a Bezier tensor product surface at a specified domain coordinate into four patches.
     * @param output00 First part of the split surface.
     * @param output01 Second part of the split surface.
     * @param output10 Third part of the split surface.
     * @param output11 Fourth part of the split surface.
     * @param uv Parameter location where to split.
     */
    void subdivide(
        BezierSurface& output00,
        BezierSurface& output01,
        BezierSurface& output10,
        BezierSurface& output11,
        DomainCoord uv = DomainCoord(0.5, 0.5)) const noexcept
    {
        double umin = mDomain.min().x(), umax = mDomain.max().x();
        double vmin = mDomain.min().y(), vmax = mDomain.max().y();
        Eigen::Vector2d uvrel(
            (uv.x() - umin) / (umax - umin),
            (uv.y() - vmin) / (vmax - vmin));
        for (int i = 0; i <= TN; ++i)
            for (int j = 0; j <= TM; ++j)
            {
                output00.mControlPoints[i][j] = intermediateControlPoint(0, 0, i, j, uvrel);
                output01.mControlPoints[i][j] = intermediateControlPoint(0, j, i, TM - j, uvrel);
                output10.mControlPoints[i][j] = intermediateControlPoint(i, 0, TN - i, j, uvrel);
                output11.mControlPoints[i][j] = intermediateControlPoint(i, j, TN - i, TM - j, uvrel);
            }

        output00.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(mDomain.min().x(), mDomain.min().y()),
            Eigen::Vector2d(uv.x(), uv.y())));
        output01.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(mDomain.min().x(), uv.y()),
            Eigen::Vector2d(uv.x(), mDomain.max().y())));
        output10.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(uv.x(), mDomain.min().y()),
            Eigen::Vector2d(mDomain.max().x(), uv.y())));
        output11.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(uv.x(), uv.y()),
            Eigen::Vector2d(mDomain.max().x(), mDomain.max().y())));

        output00.recomputeBoundingBox();
        output01.recomputeBoundingBox();
        output10.recomputeBoundingBox();
        output11.recomputeBoundingBox();
    }

    /**
     * @brief Subdivides a Bezier tensor product surface at a specified u domain coordinate into two patches.
     * @param output0 First part of the split surface.
     * @param output1 Second part of the split surface.
     * @param u Parameter location where to split.
     */
    void subdivideU(
        BezierSurface& output0,
        BezierSurface& output1,
        BezierCurve<TM, TValue>::DomainCoord u = typename BezierCurve<TM, TValue>::DomainCoord(0.5)) const noexcept
    {
        const double umin = mDomain.min().x();
        const double umax = mDomain.max().x();

        ControlPointsType& controlPointsOutput0 = output0.getControlPoints();
        ControlPointsType& controlPointsOutput1 = output1.getControlPoints();

        // when subdividing along u parameter, we subdivide each bezier curve for each v parameter separately.
        // from those curves we assemble the control nets for the two resulting surfaces.

        // iterate over all v indices
        for (int j = 0; j <= M; ++j)
        {
            // use the subdivision routine from BezierCurve to subdivide the bezier curve in each v row.
            // this writes the control points directly in the correct array of the output surfaces.
            // unfortunately, we must copy the u's for each v row because of the layout of the control points.
            // there is the option to instead use a view, but this would require changing the subdivision
            // routine in BezierCurve (and implementing the view of course).
            std::array<Value, N + 1> row;
            for (int i = 0; i <= N; ++i)
            {
                row[i] = mControlPoints[i][j];
            }
            std::array<Value, N + 1> rowOutput0;
            std::array<Value, N + 1> rowOutput1;
            BezierCurve<M, Value>::template subdivide<BezierCurve<M, Value>>(
                row,
                rowOutput0,
                rowOutput1,
                umin, umax, u);
            for (int i = 0; i <= N; ++i)
            {
                controlPointsOutput0[i][j] = rowOutput0[i];
                controlPointsOutput1[i][j] = rowOutput1[i];
            }
        }

        // set the domains
        output0.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(mDomain.min().x(), mDomain.min().y()),
            Eigen::Vector2d(u, mDomain.max().y())));
        output1.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(u, mDomain.min().y()),
            Eigen::Vector2d(mDomain.max().x(), mDomain.max().y())));
        // recompute the bounding boxes
        output0.recomputeBoundingBox();
        output1.recomputeBoundingBox();
    }

    /**
     * @brief Clips (cuts off) the part that is left/right of the specified location u.
     * @tparam Left The part to the left of u is clipped away if true, otherwise the right part is clipped away.
     * @param u Parameter location where to clip.
     */
    template <bool Left>
    BezierSurface<TN, TM, TValue> clipU(typename BezierCurve<TM, TValue>::DomainCoord u = typename BezierCurve<TM, TValue>::DomainCoord(0.5)) const noexcept
    {
        // declare
        BezierSurface<TN, TM, TValue> left;
        BezierSurface<TN, TM, TValue> right;

        // compute
        subdivideU(left, right, u);

        // return
        if constexpr (Left)
        {
            return right;
        }
        else
        {
            return left;
        }
    }

    /**
     * @brief Subdivides a Bezier tensor product surface at a specified v domain coordinate into two patches.
     * @param output0 First part of the split surface.
     * @param output1 Second part of the split surface.
     * @param v Parameter location where to split.
     */
    void subdivideV(
        BezierSurface& output0,
        BezierSurface& output1,
        BezierCurve<TN, TValue>::DomainCoord v = typename BezierCurve<TN, TValue>::DomainCoord(0.5)) const noexcept
    {
        const double vmin = mDomain.min().y();
        const double vmax = mDomain.max().y();

        // when subdividing along v parameter, we subdivide each bezier curve for each u parameter separately.
        // from those curves we assemble the control nets for the two resulting surfaces.

        // iterate over all u indices
        for (int i = 0; i <= N; ++i)
        {
            // use the subdivision routine from BezierCurve to subdivide the bezier curve in each u column
            // this writes the control points directly in the correct array of the output surfaces
            BezierCurve<M, Value>::template subdivide<BezierCurve<M, Value>>(
                mControlPoints[i],
                output0.getControlPoints()[i],
                output1.getControlPoints()[i],
                vmin, vmax, v);
        }

        // set the domains
        output0.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(mDomain.min().x(), mDomain.min().y()),
            Eigen::Vector2d(mDomain.max().x(), v)));
        output1.setDomain(Eigen::AlignedBox2d(
            Eigen::Vector2d(mDomain.min().x(), v),
            Eigen::Vector2d(mDomain.max().x(), mDomain.max().y())));
        // recompute the bounding boxes
        output0.recomputeBoundingBox();
        output1.recomputeBoundingBox();
    }

    /**
     * @brief Clips (cuts off) the part that is left/right of the specified location v.
     * @tparam Left The part to the left of v is clipped away if true, otherwise the right part is clipped away.
     * @param v Parameter location where to clip.
     */
    template <bool Left>
    BezierSurface clipV(BezierCurve<TN, TValue>::DomainCoord v = typename BezierCurve<TN, TValue>::DomainCoord(0.5)) const noexcept
    {
        // declare
        BezierSurface left;
        BezierSurface right;

        // compute
        subdivideV(left, right, v);

        // return
        if constexpr (Left)
        {
            return right;
        }
        else
        {
            return left;
        }
    }

    /**
     * @brief Recomputes the bounding box of the Bezier surface.
     */
    void recomputeBoundingBox() noexcept
    {
        mBoundingBox.setEmpty();
        for (int i = 0; i <= TN; ++i)
            for (int j = 0; j <= TM; ++j)
                mBoundingBox.extend(mControlPoints[i][j]);
    }

    /**
     * @brief Gets the bounding box of the Bezier surface.
     * @return Bounding box of this geometry.
     */
    [[nodiscard]] const Eigen::AlignedBox<double, TValue::scalar_count>& getBoundingBox() const noexcept { return mBoundingBox; }

    /**
     * @brief Gets the two-dimensional array of control points.
     * @return Two-dimensional array of control points.
     */
    [[nodiscard]] const std::array<std::array<TValue, M + 1>, N + 1>& getControlPoints() const noexcept { return mControlPoints; }

    /**
     * @brief Gets the two-dimensional array of control points.
     * @return Two-dimensional array of control points.
     */
    [[nodiscard]] std::array<std::array<TValue, M + 1>, N + 1>& getControlPoints() noexcept { return mControlPoints; }

    /**
     * @brief Sets the two-dimensional array of control points.
     * @param controlPoints Two-dimensional array of control points.
     * @param computeBoundingBox Optional flag that recomputes the bounding box according to the new control points.
     */
    void setControlPoints(const std::array<std::array<TValue, M + 1>, N + 1>& controlPoints, bool computeBoundingBox = true) noexcept
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

    /**
     * @brief Gets the parameter domain over which the function is defined.
     * @return Parameter domain.
     */
    [[nodiscard]] const BoundingBox& getDomain() const noexcept
    {
        return mDomain;
    }

private:
    /**
     * @brief Calculates the intermediate control point of the de Casteljau algorithm.
     * @param i Index of the control point in u-direction.
     * @param j Index of the control point in v-direction.
     * @param r Degree of the intermediate point in u-direction.
     * @param s Degree of the intermediate point in v-direction.
     * @param uv Parameter location where to evaluate.
     * @return Position of the intermediate control point.
     */
    [[nodiscard]] TValue intermediateControlPoint(int i, int j, int r, int s, const DomainCoord& uv) const noexcept
    {
        std::vector<double> basisTwoValues(s + 1);
        for (int jj = 0; jj <= s; ++jj)
        {
            basisTwoValues[jj] = BernsteinBasis::sample(uv.y(), jj, s);
        }

        TValue value = TValue::Zero();
        for (int ii = 0; ii <= r; ++ii)
        {
            const double basisOne = BernsteinBasis::sample(uv.x(), ii, r);
            for (int jj = 0; jj <= s; ++jj)
            {
                value += basisOne * basisTwoValues[jj] * mControlPoints[i + ii][j + jj];
            }
        }
        return value;
    }

    /**
     * @brief Parameter domain over which the function is defined.
     */
    BoundingBox mDomain;

    /**
     * @brief Bounding box of the surface.
     */
    Eigen::AlignedBox<double, TValue::scalar_count> mBoundingBox;

    /**
     * @brief Two-dimensional array of control points.
     */
    std::array<std::array<TValue, M + 1>, N + 1> mControlPoints;
};
