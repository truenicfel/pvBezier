#pragma once

#include "bernstein_basis.hpp"
#include "bezier_surface.hpp"
#include "bezier_volume_fwd.hpp"
#include "partial_uvw.hpp"

namespace BezierVolumeUtil
{
    /**
     * @brief The faces of the volume.
     *
     * @note To illustrate the conventions, we assume that u is point "right", v is pointing "up" and
     * w is pointing "back". When we look at the volume from a certain side to determine the face we
     * are working with, we always aim to have the corner in the lower left to be the new (0, 0)
     * coordinate.
     */
    enum EFaces : int32_t
    {
        /**
         * @brief The face where the u coordinate is the minimum (zero) and & w vary. The convention
         * is that v will stay the same and w will be flipped and become the new u direction of the
         * resulting bezier surface. We pretend to look at the volume from the left.
         */
        UMin = 0,

        /**
         * @brief The face where the u coordinate is maximum (N) and v & w vary. The convention is that
         * v will stay the same and w will be used as the new u direction of the resulting bezier surface.
         * We pretend to look at the volume from the right.
         */
        UMax = 1,

        /**
         * @brief The face where the v coordinate is minimum (0) and u & w vary. The convention is that
         * u will stay the same and w will be flipped and used as the new v direction of the resulting
         * bezier surface. We pretend to look at the volume from the bottom.
         */
        VMin = 2,

        /**
         * @brief The face where the v coordinate is maximum (M) and u & w vary. The convention is that
         * u will stay the same and w will be used as the new v direction of the resulting bezier surface.
         * We pretend to look at the volume from the top.
         */
        VMax = 3,

        /**
         * @brief The face where the w coordinate is minimum (0) and u & v vary. The convention is that
         * u will stay the same and v will stay the same. We pretend to look at the volume from the front.
         */
        WMin = 4,

        /**
         * @brief The face where the w coordinate is maximum (O) and u & v vary. The convention is that
         * v will stay the same and u will be flipped. We pretend to look at the volume from the back.
         */
        WMax = 5
    };

    static constexpr int FACE_COUNT = 6;

}

/**
 * @brief Class for tensor product volumes with Bernstein basis functions.
 * @tparam TValue Type of values stored at the control points.
 * @tparam TN Degree of the volume in u direction.
 * @tparam TM Degree of the volume in v direction.
 * @tparam TO Degree of the volume in w direction.
 *
 */
template <int TN, int TM, int TO, typename TValue>
class BezierVolume
{
public:
    /**
     * @brief Degree of the volume in u direction.
     */
    static constexpr int N = TN;

    /**
     * @brief Degree of the volume in v direction.
     */
    static constexpr int M = TM;

    /**
     * @brief Degree of the volume in w direction.
     */
    static constexpr int O = TO;

    /**
     * @brief Type of coordinate that specifies the location in the domain.
     */
    using DomainCoord = Eigen::Vector3d;

    /**
     * @brief Type of domain bounding box.
     */
    using BoundingBox = Eigen::AlignedBox<double, 3>;

    /**
     * @brief Type of the values returned by the function.
     */
    using Value = TValue;

    /**
     * @brief Type of the partial derivative specification.
     */
    using Partial = PartialUVW;

    /**
     * @brief Type for the control points.
     */
    using ControlPointsType = std::array<std::array<std::array<Value, O + 1>, M + 1>, N + 1>;

    template <BezierVolumeUtil::EFaces Face>
    using BoundarySurfaceType =
        std::conditional_t<
            (Face == BezierVolumeUtil::EFaces::UMin || Face == BezierVolumeUtil::EFaces::UMax),
            BezierSurface<TO, TM, TValue>,
            std::conditional_t<
                (Face == BezierVolumeUtil::EFaces::VMin || Face == BezierVolumeUtil::EFaces::VMax),
                BezierSurface<TN, TO, TValue>,
                BezierSurface<TN, TM, TValue>>>;

    /**
     * @brief Constructor with zero initialization.
     */
    BezierVolume() noexcept
        : mDomain(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones())
        , mBoundingBox()
        , mControlPoints()
    {
        for (int k = 0; k <= O; ++k)
        {
            for (int j = 0; j <= M; ++j)
            {
                for (int i = 0; i <= N; ++i)
                {
                    mControlPoints[i][j][k] = Value::Zero();
                }
            }
        }
        recomputeBoundingBox();
    }

    /**
     * @brief Constructor with initial weights.
     * @param initial Initial set of control points.
     */
    explicit BezierVolume(const ControlPointsType& initial) noexcept
        : mDomain(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones())
        , mBoundingBox()
        , mControlPoints(initial)
    {
        recomputeBoundingBox();
    }

    /**
     * @brief Evaluates the tensor-product volume algebraically at a given uvw coordinate.
     * @param uvw Domain location where to evaluate the volume.
     * @return Position on the volume.
     */
    [[nodiscard]] Value sample(const DomainCoord& uvw) const noexcept
    {
        DomainCoord trel = (uvw - mDomain.min()).cwiseQuotient(mDomain.max() - mDomain.min());
        const double u   = trel.x();
        const double v   = trel.y();
        const double w   = trel.z();
        Value result     = Value::Zero();
#if 0
            // Direct evaluation of Bernstein Bezier functions
            for (int k = 0; k <= O; ++k)
                for (int j = 0; j <= M; ++j)
                    for (int i = 0; i <= N; ++i)
                        result += mControlPoints[i][j][k] * BernsteinBasis::sample(u, i, N) * BernsteinBasis::sample(v, j, M) * BernsteinBasis::sample(w, k, O);
#else
        // Volk-Schumaker Tensor Product Evaluation
        // Precompute basis for u, v, and w using recurrence
        double Bu[N + 1] = { 0 };
        double Bv[M + 1] = { 0 };
        double Bw[O + 1] = { 0 };
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

        if (std::abs(w - 1.0) < std::numeric_limits<double>::epsilon())
            Bw[O] = 1.0;
        else
        {
            Bw[0] = std::pow(1 - w, O);
            for (int k = 1; k <= O; ++k)
                Bw[k] = Bw[k - 1] * w * (O - k + 1) / ((1 - w) * k);
        }

        // Evaluate volume
        for (int k = 0; k <= O; ++k)
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += mControlPoints[i][j][k] * Bu[i] * Bv[j] * Bw[k];
#endif
        return result;
    }

    /**
     * @brief Evaluates a partial derivative at a given domain coordinate.
     * @param uvw Domain location to sample the derivative at.
     * @param partial Specifies the desired partial derivative of each dimension.
     * @return Partial derivative at domain location.
     */
    [[nodiscard]] Value samplePartial(const DomainCoord& uvw, const Partial& partial) const noexcept
    {
        throw std::runtime_error("Not implemented!");
    }

    /**
     * @brief Generates a tensor product with degree [N, N, N] for the specified partial. Only enabled if all N, M and O are
     * equal and only works for non-mixed order on partials.
     *
     * @details Usually the tensor product will have its degree reduced by one for the dimension the partial is computed for
     * but after computing the partial that dimension is lifted to the degree N.
     *
     * @tparam DegreeU unused except for compile-time checks.
     * @tparam DegreeV unused except for compile-time checks.
     * @tparam DegreeW unused except for compile-time checks.
     *
     * @param partial The partial to generate the tensor product for.
     *
     * @return A volume containing the partial tensor product.
     */
    template <int DegreeU = N, int DegreeV = M, int DegreeW = O, typename = std::enable_if_t<(DegreeU == DegreeV && DegreeV == DegreeW)>>
    BezierVolume generatePartialTensorProduct(const PartialUVW& partial)
    {
        if (partial.hash != PartialUVW::du && partial.hash != PartialUVW::dv && partial.hash != PartialUVW::dw)
        {
            throw std::runtime_error("Partials with order larger than one are not supported yet!");
        }

        static constexpr int Degree = N;

        // compute the degrees for the individual dimensions of the partial tensor product.
        // when there is a partial for the given dimension, the degree reduces by one.
        // note that this only works as long as we only use first order and non-mixed partials.
        const int n1 = Degree - partial.degreeDu;
        const int n2 = Degree - partial.degreeDv;
        const int n3 = Degree - partial.degreeDw;

        // this will use the control points type even though one of the dimensions is reduced by one
        ControlPointsType rawResult;

        // iterate all dimensions
        for (int i = 0; i <= n1; ++i)
        {
            for (int j = 0; j <= n2; ++j)
            {
                for (int k = 0; k <= n3; ++k)
                {
                    // increment the running index by one for the dimension we are currently computing the partial for
                    // again note that this only works as long as we only use first order and non-mixed partials
                    // Eq. 22/23/24 depending on the supplied Degree
                    rawResult[i][j][k] = Degree * (mControlPoints[i + partial.degreeDu][j + partial.degreeDv][k + partial.degreeDw] - mControlPoints[i][j][k]);
                }
            }
        }

        // perform degree elevation
        // Eq. 25
        // this result will be filled
        ControlPointsType result;

        // for u partial
        if (partial.hash == PartialUVW::du)
        {
            for (int i = 0; i <= Degree; ++i)
            {
                for (int j = 0; j <= Degree; ++j)
                {
                    for (int k = 0; k <= Degree; ++k)
                    {
                        if (i == 0)
                        {
                            result[i][j][k] = rawResult[0][j][k] * (Degree - i) / Degree;
                        }
                        else if (i == Degree)
                        {
                            result[i][j][k] = rawResult[Degree - 1][j][k] * i / Degree;
                        }
                        else
                        {
                            result[i][j][k] = rawResult[i][j][k] * (Degree - i) / Degree + rawResult[i - 1][j][k] * i / Degree;
                        }
                    }
                }
            }
        }
        else if (partial.hash == PartialUVW::dv)
        {
            for (int i = 0; i <= Degree; ++i)
            {
                for (int j = 0; j <= Degree; ++j)
                {
                    for (int k = 0; k <= Degree; ++k)
                    {
                        if (j == 0)
                        {
                            result[i][j][k] = rawResult[i][0][k] * (Degree - j) / Degree;
                        }
                        else if (j == Degree)
                        {
                            result[i][j][k] = rawResult[i][Degree - 1][k] * j / Degree;
                        }
                        else
                        {
                            result[i][j][k] = rawResult[i][j][k] * (Degree - j) / Degree + rawResult[i][j - 1][k] * j / Degree;
                        }
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
                    for (int k = 0; k <= Degree; ++k)
                    {
                        if (k == 0)
                        {
                            result[i][j][k] = rawResult[i][j][0] * (Degree - k) / Degree;
                        }
                        else if (k == Degree)
                        {
                            result[i][j][k] = rawResult[i][j][Degree - 1] * k / Degree;
                        }
                        else
                        {
                            result[i][j][k] = rawResult[i][j][k] * (Degree - k) / Degree + rawResult[i][j][k - 1] * k / Degree;
                        }
                    }
                }
            }
        }

        return BezierVolume(result);
    }

    /**
     * @brief Subdivides a Bezier tensor product volume at a specified domain coordinate into eight patches.
     * @param output000 left, lower, front part of the split volume
     * @param output100 right, lower, front part of the split volume
     * @param output010 left, upper, front of the split volume
     * @param output110 right, upper, front part of the split volume
     * @param output001 left, lower, back part of the split volume
     * @param output101 right, lower, back part of the split volume
     * @param output011 left, upper, back part of the split volume
     * @param output111 right, upper, back part of the split volume
     * @param uvw Parameter location where to split.
     */
    void subdivide(
        BezierVolume& output000, // left, lower, front
        BezierVolume& output100, // right, lower, front
        BezierVolume& output010, // left, upper, front
        BezierVolume& output110, // right, upper, front
        BezierVolume& output001, // left, lower, back
        BezierVolume& output101, // right, lower, back
        BezierVolume& output011, // left, upper, back
        BezierVolume& output111, // right, upper, back
        DomainCoord uvw = DomainCoord(0.5, 0.5, 0.5)) const noexcept
    {
        double umin = mDomain.min().x(), umax = mDomain.max().x();
        double vmin = mDomain.min().y(), vmax = mDomain.max().y();
        double wmin = mDomain.min().z(), wmax = mDomain.max().z();
        DomainCoord uvwrel(
            (uvw.x() - umin) / (umax - umin),
            (uvw.y() - vmin) / (vmax - vmin),
            (uvw.z() - wmin) / (wmax - wmin));

        for (int i = 0; i <= TN; ++i)
        {
            for (int j = 0; j <= TM; ++j)
            {
                for (int k = 0; k <= TO; ++k)
                {
                    output000.mControlPoints[i][j][k] = intermediateControlPoint(
                        0, 0, 0,
                        i, j, k,
                        uvwrel);
                    output100.mControlPoints[i][j][k] = intermediateControlPoint(
                        i, 0, 0,
                        TN - i, j, k,
                        uvwrel);
                    output010.mControlPoints[i][j][k] = intermediateControlPoint(0, j, 0,
                                                                                 i, TM - j, k,
                                                                                 uvwrel);
                    output110.mControlPoints[i][j][k] = intermediateControlPoint(i, j, 0,
                                                                                 TN - i, TM - j, k,
                                                                                 uvwrel);
                    output001.mControlPoints[i][j][k] = intermediateControlPoint(0, 0, k,
                                                                                 i, j, TO - k,
                                                                                 uvwrel);
                    output101.mControlPoints[i][j][k] = intermediateControlPoint(i, 0, k,
                                                                                 TN - i, j, TO - k,
                                                                                 uvwrel);
                    output011.mControlPoints[i][j][k] = intermediateControlPoint(0, j, k,
                                                                                 i, TM - j, TO - k,
                                                                                 uvwrel);
                    output111.mControlPoints[i][j][k] = intermediateControlPoint(i, j, k,
                                                                                 TN - i, TM - j, TO - k,
                                                                                 uvwrel);
                }
            }
        }
        DomainCoord min = mDomain.min();
        DomainCoord max = mDomain.max();

        output000.setDomain(BoundingBox(
            DomainCoord(min),
            DomainCoord(uvw)));
        output100.setDomain(BoundingBox(
            DomainCoord(uvw.x(), min.y(), min.z()),
            DomainCoord(max.x(), uvw.y(), uvw.z())));
        output010.setDomain(BoundingBox(
            DomainCoord(min.x(), uvw.y(), min.z()),
            DomainCoord(uvw.x(), max.y(), uvw.z())));
        output110.setDomain(BoundingBox(
            DomainCoord(uvw.x(), uvw.y(), min.z()),
            DomainCoord(max.x(), max.y(), uvw.z())));
        output001.setDomain(BoundingBox(
            DomainCoord(min.x(), min.y(), uvw.z()),
            DomainCoord(uvw.x(), uvw.y(), max.z())));
        output101.setDomain(BoundingBox(
            DomainCoord(uvw.x(), min.y(), uvw.z()),
            DomainCoord(max.x(), uvw.y(), max.z())));
        output011.setDomain(BoundingBox(
            DomainCoord(min.x(), uvw.y(), uvw.z()),
            DomainCoord(uvw.x(), max.y(), max.z())));
        output111.setDomain(BoundingBox(
            DomainCoord(uvw.x(), uvw.y(), uvw.z()),
            DomainCoord(max.x(), max.y(), max.z())));

        output000.recomputeBoundingBox();
        output100.recomputeBoundingBox();
        output010.recomputeBoundingBox();
        output110.recomputeBoundingBox();
        output001.recomputeBoundingBox();
        output101.recomputeBoundingBox();
        output011.recomputeBoundingBox();
        output111.recomputeBoundingBox();
    }

    /**
     * @brief Recomputes the bounding box of the Bezier volume.
     */
    void recomputeBoundingBox() noexcept
    {
        mBoundingBox.setEmpty();
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                for (int k = 0; k <= O; ++k)
                {
                    mBoundingBox.extend(mControlPoints[i][j][k]);
                }
            }
        }
    }

    /**
     * @brief Gets the bounding box of the Bezier surface.
     * @return Bounding box of this geometry.
     */
    [[nodiscard]] const Eigen::AlignedBox<double, Value::scalar_count>& getBoundingBox() const noexcept { return mBoundingBox; }

    /**
     * @brief Gets the three-dimensional array of control points.
     * @return three-dimensional array of control points.
     */
    [[nodiscard]] const ControlPointsType& getControlPoints() const noexcept { return mControlPoints; }

    /**
     * @brief Gets the three-dimensional array of control points.
     * @return three-dimensional array of control points.
     */
    [[nodiscard]] ControlPointsType& getControlPoints() noexcept { return mControlPoints; }

    /**
     * @brief Sets the three-dimensional array of control points.
     * @param controlPoints three-dimensional array of control points.
     * @param computeBoundingBox Optional flag that recomputes the bounding box according to the new control points.
     */
    void setControlPoints(const ControlPointsType& controlPoints, bool computeBoundingBox = true) noexcept
    {
        mControlPoints = controlPoints;
        if (computeBoundingBox)
            recomputeBoundingBox();
    }

    template <BezierVolumeUtil::EFaces Face>
    void getBoundarySurface(BoundarySurfaceType<Face>& surface) const
    {

        if constexpr (Face == BezierVolumeUtil::EFaces::UMin)
        {
            getBoundarySurfaceUMin(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::VMin)
        {
            getBoundarySurfaceVMin(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::WMin)
        {
            getBoundarySurfaceWMin(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::UMax)
        {
            getBoundarySurfaceUMax(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::VMax)
        {
            getBoundarySurfaceVMax(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::WMax)
        {
            getBoundarySurfaceWMax(surface);
        }
        else
        {
            static_assert(false, "Invalid face type");
        }
    }

    template <BezierVolumeUtil::EFaces Face>
    BoundarySurfaceType<Face> getBoundarySurface() const
    {

        BoundarySurfaceType<Face> surface;

        if constexpr (Face == BezierVolumeUtil::EFaces::UMin)
        {
            getBoundarySurfaceUMin(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::VMin)
        {
            getBoundarySurfaceVMin(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::WMin)
        {
            getBoundarySurfaceWMin(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::UMax)
        {
            getBoundarySurfaceUMax(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::VMax)
        {
            getBoundarySurfaceVMax(surface);
        }
        else if constexpr (Face == BezierVolumeUtil::EFaces::WMax)
        {
            getBoundarySurfaceWMax(surface);
        }
        else
        {
            static_assert(false, "Invalid face type");
        }

        return surface;
    }

    /**
     * @brief Gets the boundary surface along u, v (from 0->1) for w=0.
     *
     * @see EFaces::WMin
     *
     * @param surface Resulting Bezier surface.
     */
    void getBoundarySurfaceWMin(BezierSurface<N, M, Eigen::Vector3d>& surface) const
    {

        typename BezierSurface<N, M, Eigen::Vector3d>::ControlPointsType& surfaceControlPoints = surface.getControlPoints();
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                surfaceControlPoints[i][j] = mControlPoints[i][j][0];
            }
        }
        surface.setDomain(Eigen::AlignedBox2d(mDomain.min().xy(), mDomain.max().xy()));
        surface.recomputeBoundingBox();
    }

    /**
     * @brief Gets the boundary surface along u, w (from 0->1) for v=0.
     *
     * @see EFaces::VMin
     *
     * @param surface Resulting Bezier surface.
     */
    void getBoundarySurfaceVMin(BezierSurface<N, O, Eigen::Vector3d>& surface) const
    {

        typename BezierSurface<N, O, Eigen::Vector3d>::ControlPointsType& surfaceControlPoints = surface.getControlPoints();
        for (int i = 0; i <= N; ++i)
        {
            for (int k = 0; k <= O; ++k)
            {
                surfaceControlPoints[i][k] = mControlPoints[i][0][O - k];
            }
        }
        const Eigen::Vector3d& min = mDomain.min();
        const Eigen::Vector3d& max = mDomain.max();
        surface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(min.x(), min.z()), Eigen::Vector2d(max.x(), max.z())));
        surface.recomputeBoundingBox();
    }

    /**
     * @brief Gets the boundary surface along v, w (from 0->1) for u=0.
     *
     * @see EFaces::UMin
     *
     * @param surface Resulting Bezier surface.
     */
    void getBoundarySurfaceUMin(BezierSurface<O, M, Eigen::Vector3d>& surface) const
    {
        typename BezierSurface<O, M, Eigen::Vector3d>::ControlPointsType& surfaceControlPoints = surface.getControlPoints();
        for (int j = 0; j <= M; ++j)
        {
            for (int k = 0; k <= O; ++k)
            {
                surfaceControlPoints[k][j] = mControlPoints[0][j][O - k];
            }
        }
        const Eigen::Vector3d& min = mDomain.min();
        const Eigen::Vector3d& max = mDomain.max();
        surface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(min.z(), min.y()), Eigen::Vector2d(max.z(), max.y())));
        surface.recomputeBoundingBox();
    }

    /**
     * @brief Gets the boundary surface along u, v (from 0->1) for w=1.
     *
     * @see EFaces::WMax
     *
     * @param surface Resulting Bezier surface.
     */
    void getBoundarySurfaceWMax(BezierSurface<N, M, Eigen::Vector3d>& surface) const
    {
        typename BezierSurface<N, M, Eigen::Vector3d>::ControlPointsType& surfaceControlPoints = surface.getControlPoints();
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                surfaceControlPoints[i][j] = mControlPoints[N - i][j][O];
            }
        }
        const Eigen::Vector3d& min = mDomain.min();
        const Eigen::Vector3d& max = mDomain.max();
        surface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(min.x(), min.y()), Eigen::Vector2d(max.x(), max.y())));
        surface.recomputeBoundingBox();
    }

    /**
     * @brief Gets the boundary surface along u, w (from 0->1) for v=1.
     *
     * @see EFaces::VMax
     *
     * @param surface Resulting Bezier surface.
     */
    void getBoundarySurfaceVMax(BezierSurface<N, O, Eigen::Vector3d>& surface) const
    {
        typename BezierSurface<N, O, Eigen::Vector3d>::ControlPointsType& surfaceControlPoints = surface.getControlPoints();
        for (int i = 0; i <= N; ++i)
        {
            for (int k = 0; k <= O; ++k)
            {
                surfaceControlPoints[i][k] = mControlPoints[i][M][k];
            }
        }
        const Eigen::Vector3d& min = mDomain.min();
        const Eigen::Vector3d& max = mDomain.max();
        surface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(min.x(), min.z()),
                                              Eigen::Vector2d(max.x(), max.z())));
        surface.recomputeBoundingBox();
    }

    /**
     * @brief Gets the boundary surface along v, w (from 0->1) for u=1.
     *
     * @see EFaces::UMax
     *
     * @param surface Resulting Bezier surface.
     */
    void getBoundarySurfaceUMax(BezierSurface<O, M, Eigen::Vector3d>& surface) const
    {
        typename BezierSurface<O, M, Eigen::Vector3d>::ControlPointsType& surfaceControlPoints = surface.getControlPoints();
        for (int j = 0; j <= M; ++j)
        {
            for (int k = 0; k <= O; ++k)
            {
                surfaceControlPoints[k][j] = mControlPoints[N][j][k];
            }
        }
        const Eigen::Vector3d& min = mDomain.min();
        const Eigen::Vector3d& max = mDomain.max();
        surface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(min.z(), min.y()),
                                              Eigen::Vector2d(max.z(), max.y())));
        surface.recomputeBoundingBox();
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
    ControlPointsType mControlPoints;

    /**
     * @brief Calculates the intermediate control point of the de Casteljau algorithm.
     * @param i Index of the control point in u-direction.
     * @param j Index of the control point in v-direction.
     * @param k Index of the control point in w-direction.
     * @param r Degree of the intermediate point in u-direction.
     * @param s Degree of the intermediate point in v-direction.
     * @param t Degree of the intermediate point in w-direction.
     * @param uvw Parameter location where to evaluate.
     * @return Position of the intermediate control point.
     */
    [[nodiscard]] TValue intermediateControlPoint(int i, int j, int k, int r, int s, int t, const DomainCoord& uvw) const noexcept
    {
        std::vector<double> basisThreeValues(t + 1);
        for (int kk = 0; kk <= t; ++kk)
        {
            basisThreeValues[kk] = BernsteinBasis::sample(uvw.z(), kk, t);
        }

        std::vector<double> basisTwoValues(s + 1);
        for (int jj = 0; jj <= s; ++jj)
        {
            basisTwoValues[jj] = BernsteinBasis::sample(uvw.y(), jj, s);
        }

        TValue value = TValue::Zero();
        for (int ii = 0; ii <= r; ++ii)
        {
            const double basisOne = BernsteinBasis::sample(uvw.x(), ii, r);
            for (int jj = 0; jj <= s; ++jj)
            {
                const double basisTwo = basisOne * basisTwoValues[jj];
                for (int kk = 0; kk <= t; ++kk)
                {
                    value += basisTwo * basisThreeValues[kk] * mControlPoints[i + ii][j + jj][k + kk];
                }
            }
        }
        return value;
    }
};
