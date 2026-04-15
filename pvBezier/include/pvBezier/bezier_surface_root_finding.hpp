#pragma once

#include "eclipping_strategy.hpp"
#include "implicit_line.hpp"
#include "partial_uv.hpp"
#include "plane3d.hpp"

#include <cstdint>
#include <vector>

#include <Eigen/Eigen>

/**
 * @brief Implements root finding strategies on bezier surfaces.
 *
 * @details
 */
template <typename TBezierSurface>
class BezierSurfaceRootFinding
{
public:
    /**
     * @brief The domain coordinate used by the bezier surface (usually Eigen::Vector2d for surfaces).
     */
    using DomainCoord = TBezierSurface::DomainCoord;

    /**
     * @brief The degree of the bezier surface in the first direction (index increments the slowest).
     */
    static constexpr int N = TBezierSurface::N;

    /**
     * @brief The degree of the bezier surface in the second direction (index increments the fastest).
     */
    static constexpr int M = TBezierSurface::M;

    /**
     * @brief The values returned by the bezier surface when sampling (also the type of the control points).
     */
    using Value = TBezierSurface::Value;

    /**
     * @brief The number of components the values that are returned have.
     */
    static constexpr int Components = Value::SizeAtCompileTime;
    static_assert(Components == 3, "Components must always be 3!");

    /**
     * @brief The type of the control points used by the bezier surface.
     */
    using ControlPointsType = TBezierSurface::ControlPointsType;

    /**
     * @brief The type that specifies the possible clipping range.
     */
    using ClipRange = std::pair<double, double>;

    /**
     * @brief Constructor.
     *
     * @param zeroEpsilon When the function reaches this value in its center, the value is considered to be zero. Defaults to 1e-8.
     * @param requiredDomainSize The function domain must reach a size smaller than this for the root location to be
     * precise enough. Defaults to 0.001% (1e-5) of a normalized domain [0, 1] in both u and v.
     * @param duplicateEpsilon If two roots are closer than this epsilon, the new root will be discarded. Defaults to 1e-3.
     * @param maximumDepth If the recursion reaches a depth larger than this, it will be stopped. Defaults to 256.
     * @param newtonStartSize If a hybrid approach is chosen, the newton iterations will take over when the domain is smaller
     * than this in both dimensions. Defaults to 1e-3.
     * @param maxNewtonIterations If a hybrid approach is chosen, this gives an upper bound for the number of newton steps that
     * can be done. Defaults to 2048.
     */
    explicit BezierSurfaceRootFinding(const double& zeroEpsilon         = 1e-8,
                                      const double& requiredDomainSize  = 1e-5,
                                      const double& duplicateEpsilon    = 1e-3,
                                      const uint64_t& maximumDepth      = 256,
                                      const double& newtonStartSize     = 1e-3,
                                      const double& maxNewtonIterations = 2048)

        : mZeroEpsilon(zeroEpsilon)
        , mRequiredDomainSize(requiredDomainSize)
        , mDuplicateEpsilon(duplicateEpsilon)
        , mMaximumDepth(maximumDepth)
        , mNewtonStartSize(newtonStartSize)
        , mMaxNewtonIterations(maxNewtonIterations)
    {
    }

    /**
     * @brief Solve (find roots) the given bezier surface using bezier clipping.
     *
     * @tparam ClippingStrategy The strategy to compute the clipping ranges.
     *
     * @param surface The surface where the roots are searched in.
     * @param uvs The uv coordinates of the roots.
     */
    template <EClippingStrategy ClippingStrategy = ProjectionBased>
    void clippingSolver(const TBezierSurface& surface, std::vector<Eigen::Vector2d>& uvs)
    {
        // Necessary for Eq. 41
        Eigen::Vector2d rangeU(surface.getDomain().min()[0], surface.getDomain().max()[0]);
        Eigen::Vector2d rangeV(surface.getDomain().min()[1], surface.getDomain().max()[1]);
        solveRecursiveClipping<ClippingStrategy, false>(surface, surface, rangeU, rangeV, mMaximumDepth, uvs);
    }

    /**
     * @brief Solve (find roots) the given bezier surface using bezier clipping with subsequent newton iterations.
     *
     * @tparam ClippingStrategy The strategy to compute the clipping ranges.
     *
     * @param surface The surface where the roots are searched in.
     * @param uvs The uv coordinates of the roots.
     */
    template <EClippingStrategy ClippingStrategy = ProjectionBased>
    void hybridClippingSolver(const TBezierSurface& surface, std::vector<Eigen::Vector2d>& uvs)
    {
        Eigen::Vector2d rangeU(surface.getDomain().min()[0], surface.getDomain().max()[0]);
        Eigen::Vector2d rangeV(surface.getDomain().min()[1], surface.getDomain().max()[1]);
        solveRecursiveClipping<ClippingStrategy, true>(surface, surface, rangeU, rangeV, mMaximumDepth, uvs);
    }

    /**
     * @brief Solve (find roots) the given bezier surface using bezier subdivision.
     *
     * @param surface The surface where the roots are searched in.
     * @param uvs The uv coordinates of the roots.
     */
    void bisectionSolver(const TBezierSurface& surface, std::vector<Eigen::Vector2d>& uvs)
    {
        solveRecursiveBisection<false>(surface, mMaximumDepth, uvs);
    }

    /**
     * @brief Solve (find roots) the given bezier surface using bezier subdivision subsequent newton iterations.
     *
     * @param surface The surface where the roots are searched in.
     * @param uvs The uv coordinates of the roots.
     */
    void hybridBisectionSolver(const TBezierSurface& surface, std::vector<Eigen::Vector2d>& uvs)
    {
        solveRecursiveBisection<true>(surface, mMaximumDepth, uvs);
    }

    /**
     * @brief Solve (find roots) the given bezier surface using newton iterations.
     *
     * @param surface The surface where the roots are searched in.
     * @param uvs The uv coordinates of the roots.
     */
    void newtonSolver(const TBezierSurface& surface, std::vector<Eigen::Vector2d>& uvs, const DomainCoord& startPoint = DomainCoord(0.5, 0.5))
    {
        DomainCoord candidate;

        // compute du and dv derivatives
        const TBezierSurface dU = surface.generatePartialTensorProduct(PartialUV::du);
        const TBezierSurface dV = surface.generatePartialTensorProduct(PartialUV::dv);
        candidate               = newtonRefinement(surface, dU, dV, startPoint);

        const bool accepted = candidate.x() >= 0 && candidate.x() <= 1 && candidate.y() >= 0 && candidate.y() <= 1;

        if (accepted)
        {
            const Eigen::AlignedBox2d& domain = surface.getDomain();
            // then we translate that to a point to the actual domain
            candidate.x() = (domain.max().x() - domain.min().x()) * candidate.x() + domain.min().x();
            candidate.y() = (domain.max().y() - domain.min().y()) * candidate.y() + domain.min().y();
            if (!domainCoordinateContained(candidate, uvs))
            {
                // add to uvs
                uvs.push_back(candidate);
            }
        }
    }

private:
    /**
     * @brief When the function reaches this value in its center the value is considered to be zero. Defaults to 1e-3.
     */
    const double mZeroEpsilon;

    /**
     * @brief The function domain must reach a size smaller than this for the root location to be
     * precise enough. Defaults to 0.001% (1e-5) of a normalized domain [0, 1] in both u and v.
     */
    const double mRequiredDomainSize;

    /**
     * @brief Epsilon to determine if a point is a duplicate or not (defaults to 1e-3).
     */
    const double mDuplicateEpsilon;

    /**
     * @brief The maximum number of recursions the algorithm is allowed to perform. Defaults to 256.
     */
    const uint64_t mMaximumDepth;

    /**
     * @brief If a hybrid approach is chosen, the newton iterations will take over when the domain is smaller
     * than this in both dimensions. Defaults to 1e-3.
     */
    const double mNewtonStartSize;

    /**
     * @brief If a hybrid approach is chosen, this gives an upper bound for the number of newton steps that
     * can be done. Defaults to 2048.
     */
    const double mMaxNewtonIterations;

    /**
     * @brief Checks whether a domain coordinate is contained in a list of domain coordinates.
     *
     * @param domainCoordinate the point to check against the list of domainCoordinates.
     * @param domainCoordinates the list of domainCoordinates.
     * @return true if the domainCoordinate is already in the list.
     */
    bool domainCoordinateContained(const DomainCoord& domainCoordinate, const std::vector<DomainCoord>& domainCoordinates)
    {
        // go through the list of point indices to check
        for (const DomainCoord& candidate : domainCoordinates)
        {
            // if an existing point is closer than the duplicate distance, return yes
            if ((candidate - domainCoordinate).squaredNorm() < mDuplicateEpsilon * mDuplicateEpsilon)
                return true;
        }
        // no other point is in duplicate distance. seems to be a new point
        return false;
    }

    /**
     * @brief Checks if the given control points have a zero crossing.
     *
     * @param controlPoints the control points to check.
     *
     * @return true if there is a zero crossing in ANY of the components.
     */
    bool hasZeroCrossing(const ControlPointsType& controlPoints)
    {
        std::array<bool, Components> anyAbove{}; // init to false
        std::array<bool, Components> anyBelow{}; // init to false
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                const Value& value = controlPoints[i][j];
                for (int componentIndex = 0; componentIndex < Components; ++componentIndex)
                {
                    if ((!anyAbove[componentIndex] || !anyBelow[componentIndex]) && value[componentIndex] == 0.0)
                    {
                        anyAbove[componentIndex] = true;
                        anyBelow[componentIndex] = true;
                    }
                    if (!anyAbove[componentIndex] && value[componentIndex] > 0)
                    {
                        anyAbove[componentIndex] = true;
                    }
                    if (!anyBelow[componentIndex] && value[componentIndex] < 0)
                    {
                        anyBelow[componentIndex] = true;
                    }
                }
            }
        }
        for (int componentIndex = 0; componentIndex < Components; ++componentIndex)
        {
            if (!(anyAbove[componentIndex] && anyBelow[componentIndex]))
            {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Returns the clipping range for the given control points.
     *
     * @tparam ClipU set this to true if clipping should happen in u direction and false for clipping in v direction.
     * @tparam Component The component the clipping is done for (must be smaller than Components of TBezierSurface).
     *
     * @param controlPoints the control points to do the clipping for.
     *
     * @return The clipping range [umin,umax]/[vmin, vmax] to discard parts of the surface that do not have a zero intersection.
     */
    template <bool ClipU, int Component>
    [[nodiscard]] ClipRange clip(const ControlPointsType& controlPoints)
    {

        static_assert(Component < Components, "Component must be smaller than Components");
        static_assert(Component >= 0, "Component must be larger than 0");

        // reset clip range
        double smin = 1;
        double smax = 0;

        // traverse the rows of the control polygon and collect the smallest and largest distance for each column.
        constexpr int Columns = ClipU ? N + 1 : M + 1;
        std::array<double, Columns> min;
        std::fill(min.begin(), min.end(), std::numeric_limits<double>::max());
        std::array<double, Columns> max;
        std::fill(max.begin(), max.end(), std::numeric_limits<double>::lowest());
        // the columns change depending on the clipping direction (U or V). either N or M determines the number of columns.
        if constexpr (ClipU)
        {
            for (int i = 0; i <= N; ++i)
            {
                for (int j = 0; j <= M; ++j)
                {
                    min[i] = std::min(min[i], controlPoints[i][j](Component));
                    max[i] = std::max(max[i], controlPoints[i][j](Component));
                }
            }
        }
        else
        {
            for (int i = 0; i <= N; ++i)
            {
                for (int j = 0; j <= M; ++j)
                {
                    min[j] = std::min(min[j], controlPoints[i][j](Component));
                    max[j] = std::max(max[j], controlPoints[i][j](Component));
                }
            }
        }

        // if the first column has a positive and negative distance in it, then this left boundary of the convex hull intersects the x-axis at 0.
        if (min.front() <= 0.0 && 0.0 <= max.front())
        {
            smin = 0;
        }

        // same for the last column
        if (min.back() <= 0.0 && 0.0 <= max.back())
        {
            smax = 1;
        }

        // Let's call the entries of controlPoints = {e_ij}
        // The columns/rows (depending on ClipU) of controlPoints are mapped to t-values[t0, t1, t2, t3] = [0, 1/3, 2/3, 1]
        // This gives a coordinate for each point in controlPoints.
        // for ClipU: (t0,e_00), (t0,e_01), ..., (tN,e_NM)
        // for ClipV: (t0,e_00), (t1,e_10), ..., (tM,e_NM)
        // We need to intersect the convex hull of these (N+1)*(M+1) points with the x-axis at y=0.
        // Instead of finding the convex hull, we get the lowest and highest values per column, and test each connecting edge among the lowest values, and the connecting edge among the highest values per column.
        // This will include the edges of the convex hull, but it is also contains a few unnecessary checks.
        for (int i = 0; i < Columns - 1; ++i)
            for (int j = i + 1; j < Columns; ++j)
            {
                // the two points of the segment are(t1, e1) and (t2, e2)
                double t1 = i / static_cast<double>(Columns - 1);
                double t2 = j / static_cast<double>(Columns - 1);
                double e1 = min[i];
                double e2 = min[j];
                // horizontal segment cannot intersect with the line
                if (e1 != e2)
                {
                    // does the segment intersect with 0?
                    const double u = (0 - e1) / (e2 - e1);
                    if (0 <= u && u <= 1)
                    {
                        const double t = u * (t2 - t1) + t1;
                        smin           = std::min(t, smin);
                        smax           = std::max(t, smax);
                    }
                }
                e1 = max[i];
                e2 = max[j];
                // horizontal segment cannot intersect with the line dmin or the line dmax
                if (e1 != e2)
                {
                    // does the segment intersect with 0?
                    const double u = (0 - e1) / (e2 - e1);
                    if (0 <= u && u <= 1)
                    {
                        const double t = u * (t2 - t1) + t1;
                        smin           = std::min(t, smin);
                        smax           = std::max(t, smax);
                    }
                }
            }
        return std::make_pair(smin, smax);
    }

    /**
     * @brief Returns the clipping range for the given distances to control points for a bezier surface of degree NxM.
     *
     * @details This function accepts the bezier representation of the signed distances (the control points are the signed distances
     * of the bezier surface we are clipping on). For simplicity it just accepts the control points of that signed distance bezier
     * surface. The major difference here is that the bezier surface we are clipping on is now a uni-variate function (scalar valued).
     *
     * @tparam ClipU Set this to true if clipping should happen in u direction and false for clipping in v direction.
     *
     * @param controlPoints the control points to do the clipping for.
     *
     * @return The clipping range [umin,umax]/[vmin, vmax] to discard parts of the surface that do not have a zero intersection.
     */
    template <bool ClipU>
    [[nodiscard]] ClipRange clip(const std::array<std::array<double, M + 1>, N + 1>& controlPoints)
    {
        // initialize clip range
        double smin = 1;
        double smax = 0;

        // traverse the rows of the control polygon and collect the smallest and largest distance for each column.
        constexpr int Columns = ClipU ? N + 1 : M + 1;
        std::array<double, Columns> min;
        std::fill(min.begin(), min.end(), std::numeric_limits<double>::max());
        std::array<double, Columns> max;
        std::fill(max.begin(), max.end(), std::numeric_limits<double>::lowest());
        // the columns change depending on the clipping direction (U or V). either N or M determines the number of columns.
        if constexpr (ClipU)
        {
            for (int i = 0; i <= N; ++i)
            {
                for (int j = 0; j <= M; ++j)
                {
                    min[i] = std::min(min[i], controlPoints[i][j]);
                    max[i] = std::max(max[i], controlPoints[i][j]);
                }
            }
        }
        else
        {
            for (int i = 0; i <= N; ++i)
            {
                for (int j = 0; j <= M; ++j)
                {
                    min[j] = std::min(min[j], controlPoints[i][j]);
                    max[j] = std::max(max[j], controlPoints[i][j]);
                }
            }
        }

        // if the first column has a positive and negative distance in it,
        // then this left boundary of the convex hull intersects the x-axis at 0.
        if (min.front() < 0 && 0 < max.front())
        {
            smin = 0;
        }

        // same for the last column
        if (min.back() < 0 && 0 < max.back())
        {
            smax = 1;
        }

        // Let's call the entries of controlPoints = {e_ij}
        // The columns/rows (depending on ClipU) of controlPoints are mapped to t-values[t0, t1, t2, t3] = [0, 1/3, 2/3, 1]
        // This gives a coordinate for each point in controlPoints.
        // for ClipU: (t0,e_00), (t0,e_01), ..., (tN,e_NM)
        // for ClipV: (t0,e_00), (t1,e_10), ..., (tM,e_NM)
        // We need to intersect the convex hull of these (N+1)*(M+1) points with the x-axis at y=0.
        // Instead of finding the convex hull, we get the lowest and highest values per column, and test each connecting edge among the lowest values, and the connecting edge among the highest values per column.
        // This will include the edges of the convex hull, but it is also contains a few unnecessary checks.
        for (int i = 0; i < Columns - 1; ++i)
            for (int j = i + 1; j < Columns; ++j)
            {
                // the two points of the segment are(t1, e1) and (t2, e2)
                double t1 = i / static_cast<double>(Columns - 1);
                double t2 = j / static_cast<double>(Columns - 1);
                double e1 = min[i];
                double e2 = min[j];
                // horizontal segment cannot intersect with the line
                if (e1 != e2)
                {
                    // does the segment intersect with 0?
                    const double u = (0 - e1) / (e2 - e1);
                    if (0 <= u && u <= 1)
                    {
                        const double t = u * (t2 - t1) + t1;
                        smin           = std::min(t, smin);
                        smax           = std::max(t, smax);
                    }
                }
                e1 = max[i];
                e2 = max[j];
                // horizontal segment cannot intersect with the line dmin or the line dmax
                if (e1 != e2)
                {
                    // does the segment intersect with 0?
                    const double u = (0 - e1) / (e2 - e1);
                    if (0 <= u && u <= 1)
                    {
                        const double t = u * (t2 - t1) + t1;
                        smin           = std::min(t, smin);
                        smax           = std::max(t, smax);
                    }
                }
            }
        return std::make_pair(smin, smax);
    }

    /**
     * @brief Returns the clipping range for the given control points for all components (minimums are maxed and maximums are "min'd")
     * @tparam ClipU set this to true if clipping should happen in u direction and false for clipping in v direction.
     * @tparam Is The indices for the components.
     *
     * @param controlPoints the control points to do the clipping for.k
     *
     * @return The aggregated clipping ranges (all mins are maxed and all maxes are "min'd").
     */
    template <bool ClipU, int... Is>
    ClipRange clipAll(const ControlPointsType& controlPoints, std::integer_sequence<int, Is...> sequence)
    {
        std::array<ClipRange, Components> ranges = { clip<ClipU, Is>(controlPoints)... };

        double maxOfMin = std::numeric_limits<double>::lowest();
        double minOfMax = std::numeric_limits<double>::max();

        // from all min values, find the largest one
        for (auto& [min, max] : ranges)
        {
            maxOfMin = std::max(maxOfMin, min);
            minOfMax = std::min(minOfMax, max);
        }

        return std::make_pair(maxOfMin, minOfMax);
    }

    template <bool ClipU>
    ClipRange clipAll(const ControlPointsType& controlPoints)
    {
        return clipAll<ClipU>(controlPoints, std::make_integer_sequence<int, Components>{});
    }

    /**
     * @brief Compute the clipping ranges for a 3 component bezier surface in both u and v direction using a plane.
     *
     * @param controlPoints the control points of the current surface.
     * @param initialSurfaceControlPoints the control points of the original surface (as a fallback).
     *
     * @return a pair of optional pairs containing the clipping ranges in u and v [[umin, umax], [vmin, vmax]].
     */
    std::pair<std::optional<ClipRange>, std::optional<ClipRange>> compute3ComponentDistanceClipUV(const std::array<std::array<Eigen::Vector3d, M + 1>, N + 1>& controlPoints, const std::array<std::array<Eigen::Vector3d, M + 1>, N + 1>& initialSurfaceControlPoints)
    {
        // Eq. 34 & 35
        // See Eq.(12) in Alexander Efremov, Vlastimil Havran, and Hans-Peter Seidel. 2005. Robust and numerically stable Bezier clipping method for ray tracing NURBS surfaces. In Proceedings of the 21st Spring Conference on Computer Graphics (SCCG '05). Association for Computing Machinery, New York, NY, USA, 127�135. https://doi.org/10.1145/1090122.1090144
        // vertical line: subdivide u
        Value edgeUOne = controlPoints[N][0] - controlPoints[0][0];
        Value edgeUTwo = controlPoints[N][M] - controlPoints[0][M];
        Value edgeU    = edgeUOne + edgeUTwo;
        if (edgeU.norm() < 1E-10)
        {
            // fall back to initial directions in case the surface already collapsed too far upon convergence
            edgeUOne = initialSurfaceControlPoints[N][0] - initialSurfaceControlPoints[0][0];
            edgeUTwo = initialSurfaceControlPoints[N][M] - initialSurfaceControlPoints[0][M];
            edgeU    = edgeUOne + edgeUTwo;
        }

        Value edgeVOne = controlPoints[0][M] - controlPoints[0][0];
        Value edgeVTwo = controlPoints[N][M] - controlPoints[N][0];
        Value edgeV    = edgeVOne + edgeVTwo;
        if (edgeV.norm() < 1E-10)
        {
            // fall back to initial directions in case the surface already collapsed too far upon convergence
            edgeVOne = initialSurfaceControlPoints[0][M] - initialSurfaceControlPoints[0][0];
            edgeVTwo = initialSurfaceControlPoints[N][M] - initialSurfaceControlPoints[N][0];
            edgeV    = edgeVOne + edgeVTwo;
        }

        std::array<std::array<double, M + 1>, N + 1> signedDistancesClipU;
        std::array<std::array<double, M + 1>, N + 1> signedDistancesClipV;

        // Eq. 36
        edgeU.normalize();
        edgeV.normalize();
        Value c = edgeU.cross(edgeV);
        c.normalize();
        Value normalClipU;
        Value normalClipV;
        normalClipU = edgeV.cross(c);
        normalClipV = c.cross(edgeU);

        // Eg. 37
        Plane3d planeClipU(normalClipU);
        Plane3d planeClipV(normalClipV);

        bool anyAboveU = false;
        bool anyBelowU = false;
        bool anyAboveV = false;
        bool anyBelowV = false;
        // fill the signed distances based on the control points
        // and compute the above/below criteria
        // Eq. 38 & 39
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                signedDistancesClipU[i][j] = planeClipU.distance(controlPoints[i][j]);
                signedDistancesClipV[i][j] = planeClipV.distance(controlPoints[i][j]);

                if (!anyAboveU && signedDistancesClipU[i][j] >= 0)
                {
                    anyAboveU = true;
                }
                if (!anyBelowU && signedDistancesClipU[i][j] <= 0)
                {
                    anyBelowU = true;
                }
                if (!anyAboveV && signedDistancesClipV[i][j] >= 0)
                {
                    anyAboveV = true;
                }
                if (!anyBelowV && signedDistancesClipV[i][j] <= 0)
                {
                    anyBelowV = true;
                }
            }
        }

        // combine the conditions
        const bool clipUPossible = anyAboveU && anyBelowU;
        const bool clipVPossible = anyAboveV && anyBelowV;

        return std::make_pair(
            clipUPossible ? std::optional<ClipRange>(clip<true>(signedDistancesClipU)) : std::nullopt,
            clipVPossible ? std::optional<ClipRange>(clip<false>(signedDistancesClipV)) : std::nullopt);
    }

    /**
     * @brief Compute the clipping ranges for a 3 component bezier surface in u direction a plane.
     *
     * @param controlPoints the control points of the current surface.
     * @param initialSurfaceControlPoints the control points of the original surface (as a fallback).
     *
     * @return a pair of optional pairs containing the clipping ranges in u [umin, umax].
     */
    std::optional<ClipRange> compute3ComponentDistanceClipU(const std::array<std::array<Eigen::Vector3d, M + 1>, N + 1>& controlPoints, const std::array<std::array<Eigen::Vector3d, M + 1>, N + 1>& initialSurfaceControlPoints)
    {
        // Eq. 34
        // See Eq.(12) in Alexander Efremov, Vlastimil Havran, and Hans-Peter Seidel. 2005. Robust and numerically stable Bezier clipping method for ray tracing NURBS surfaces. In Proceedings of the 21st Spring Conference on Computer Graphics (SCCG '05). Association for Computing Machinery, New York, NY, USA, 127�135. https://doi.org/10.1145/1090122.1090144
        // vertical line: subdivide u
        Value edgeUOne = controlPoints[N][0] - controlPoints[0][0];
        Value edgeUTwo = controlPoints[N][M] - controlPoints[0][M];
        Value edgeU    = edgeUOne + edgeUTwo;
        if (edgeU.norm() < 1E-10)
        {
            // fall back to initial directions in case the surface already collapsed too far upon convergence
            edgeUOne = initialSurfaceControlPoints[N][0] - initialSurfaceControlPoints[0][0];
            edgeUTwo = initialSurfaceControlPoints[N][M] - initialSurfaceControlPoints[0][M];
            edgeU    = edgeUOne + edgeUTwo;
        }

        Value edgeVOne = controlPoints[0][M] - controlPoints[0][0];
        Value edgeVTwo = controlPoints[N][M] - controlPoints[N][0];
        Value edgeV    = edgeVOne + edgeVTwo;
        if (edgeV.norm() < 1E-10)
        {
            // fall back to initial directions in case the surface already collapsed too far upon convergence
            edgeVOne = initialSurfaceControlPoints[0][M] - initialSurfaceControlPoints[0][0];
            edgeVTwo = initialSurfaceControlPoints[N][M] - initialSurfaceControlPoints[N][0];
            edgeV    = edgeVOne + edgeVTwo;
        }

        std::array<std::array<double, M + 1>, N + 1> signedDistancesClipU;

        // Eq. 36
        edgeU.normalize();
        edgeV.normalize();
        Value c = edgeU.cross(edgeV);
        c.normalize();
        // Eq. 37
        const Plane3d planeClipU(edgeV.cross(c));

        bool anyAboveU = false;
        bool anyBelowU = false;
        // fill the signed distances based on the control points
        // and compute the above/below criteria
        // Eq. 38
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                signedDistancesClipU[i][j] = planeClipU.distance(controlPoints[i][j]);

                if (!anyAboveU && signedDistancesClipU[i][j] >= 0)
                {
                    anyAboveU = true;
                }
                if (!anyBelowU && signedDistancesClipU[i][j] <= 0)
                {
                    anyBelowU = true;
                }
            }
        }

        // combine the conditions
        const bool clipUPossible = anyAboveU && anyBelowU;

        return clipUPossible ? std::optional<ClipRange>(clip<true>(signedDistancesClipU)) : std::nullopt;
    }

    /**
     * @brief Compute the clipping ranges for a 3 component bezier surface in v direction a plane.
     *
     * @param controlPoints the control points of the current surface.
     * @param initialSurfaceControlPoints the control points of the original surface (as a fallback).
     *
     * @return a pair of optional pairs containing the clipping ranges in v [vmin, vmax].
     */
    std::optional<ClipRange> compute3ComponentDistanceClipV(const std::array<std::array<Eigen::Vector3d, M + 1>, N + 1>& controlPoints, const std::array<std::array<Eigen::Vector3d, M + 1>, N + 1>& initialSurfaceControlPoints)
    {
        // Eq. 35
        // See Eq.(12) in Alexander Efremov, Vlastimil Havran, and Hans-Peter Seidel. 2005. Robust and numerically stable Bezier clipping method for ray tracing NURBS surfaces. In Proceedings of the 21st Spring Conference on Computer Graphics (SCCG '05). Association for Computing Machinery, New York, NY, USA, 127�135. https://doi.org/10.1145/1090122.1090144
        // vertical line: subdivide u
        Value edgeUOne = controlPoints[N][0] - controlPoints[0][0];
        Value edgeUTwo = controlPoints[N][M] - controlPoints[0][M];
        Value edgeU    = edgeUOne + edgeUTwo;
        if (edgeU.norm() < 1E-10)
        {
            // fall back to initial directions in case the surface already collapsed too far upon convergence
            edgeUOne = initialSurfaceControlPoints[N][0] - initialSurfaceControlPoints[0][0];
            edgeUTwo = initialSurfaceControlPoints[N][M] - initialSurfaceControlPoints[0][M];
            edgeU    = edgeUOne + edgeUTwo;
        }

        Value edgeVOne = controlPoints[0][M] - controlPoints[0][0];
        Value edgeVTwo = controlPoints[N][M] - controlPoints[N][0];
        Value edgeV    = edgeVOne + edgeVTwo;
        if (edgeV.norm() < 1E-10)
        {
            // fall back to initial directions in case the surface already collapsed too far upon convergence
            edgeVOne = initialSurfaceControlPoints[0][M] - initialSurfaceControlPoints[0][0];
            edgeVTwo = initialSurfaceControlPoints[N][M] - initialSurfaceControlPoints[N][0];
            edgeV    = edgeVOne + edgeVTwo;
        }

        std::array<std::array<double, M + 1>, N + 1> signedDistancesClipV;

        // Eq. 36
        edgeU.normalize();
        edgeV.normalize();
        Value c = edgeU.cross(edgeV);
        c.normalize();
        // Eq. 37
        Plane3d planeClipV(c.cross(edgeU));

        bool anyAboveV = false;
        bool anyBelowV = false;
        // fill the signed distances based on the control points
        // and compute the above/below criteria
        // Eq. 39
        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= M; ++j)
            {
                signedDistancesClipV[i][j] = planeClipV.distance(controlPoints[i][j]);

                if (!anyAboveV && signedDistancesClipV[i][j] >= 0)
                {
                    anyAboveV = true;
                }
                if (!anyBelowV && signedDistancesClipV[i][j] <= 0)
                {
                    anyBelowV = true;
                }
            }
        }

        // combine the conditions
        const bool clipVPossible = anyAboveV && anyBelowV;

        return clipVPossible ? std::optional<ClipRange>(clip<false>(signedDistancesClipV)) : std::nullopt;
    }

    /**
     * @brief Compute the clipping ranges for a bezier surface in both u and v direction.
     *
     * @param controlPoints the control points of the current surface.
     *
     * @return a pair of optional pairs containing the clipping ranges in u and v [[umin, umax], [vmin, vmax]].
     */
    std::pair<std::optional<ClipRange>, std::optional<ClipRange>> computeComponentClipUV(const ControlPointsType& controlPoints)
    {
        if (hasZeroCrossing(controlPoints))
        {
            const ClipRange clipRangeU = clipAll<true>(controlPoints);
            const ClipRange clipRangeV = clipAll<false>(controlPoints);

            // if first (min) is larger than second (max) then no root is possible
            const std::optional<ClipRange> resultU = clipRangeU.first > clipRangeU.second ? std::nullopt : std::optional<ClipRange>(clipRangeU);
            const std::optional<ClipRange> resultV = clipRangeV.first > clipRangeV.second ? std::nullopt : std::optional<ClipRange>(clipRangeV);

            return std::make_pair(resultU, resultV);
        }

        return std::make_pair(std::nullopt, std::nullopt);
    }

    /**
     * @brief Compute the clipping ranges for a bezier surface in u direction.
     *
     * @param controlPoints the control points of the current surface.
     *
     * @return a pair of optional pairs containing the clipping ranges in u [umin, umax].
     */
    std::optional<ClipRange> computeComponentClipU(const ControlPointsType& controlPoints)
    {
        if (hasZeroCrossing(controlPoints))
        {
            const ClipRange clipRangeV = clipAll<true>(controlPoints);
            // if first (min) is larger than second (max) then no root is possible
            return clipRangeV.first > clipRangeV.second ? std::nullopt : std::optional<ClipRange>(clipRangeV);
        }
        return std::nullopt;
    }

    /**
     * @brief Compute the clipping ranges for a bezier surface in v direction.
     *
     * @param controlPoints the control points of the current surface.
     *
     * @return a pair of optional pairs containing the clipping ranges in v [vmin, vmax].
     */
    std::optional<ClipRange> computeComponentClipV(const ControlPointsType& controlPoints)
    {
        if (hasZeroCrossing(controlPoints))
        {
            const ClipRange clipRangeV = clipAll<false>(controlPoints);
            // if first (min) is larger than second (max) then no root is possible
            return clipRangeV.first > clipRangeV.second ? std::nullopt : std::optional<ClipRange>(clipRangeV);
        }
        return std::nullopt;
    }

    /**
     * @brief Searches the given bezier surface recursively for roots using bezier clipping.
     *
     * @tparam ClippingStrategy The strategy to compute the clipping ranges.
     * @tparam Hybrid Set to true if after clipping down to mNewtonStartSize a newton solver should take over. If set to false clipping will happen until mRefineEpsilon.
     *
     * @param surface The bezier surface to solve for roots.
     * @param initialSurface The initial bezier surface, which is used to find the direction of the implicit line in case of numerical precisions problems.
     * @param rangeU Currently tested range of u-coordinates.
     * @param rangeV Currently tested range of v-coordinates.
     * @param depth The current depth of the recursion (when zero the recursion is left).
     * @param uvs List of uv coordinates found (roots).
     */
    template <EClippingStrategy ClippingStrategy = ProjectionBased,
              bool Hybrid                        = false>
    void solveRecursiveClipping(TBezierSurface surface, const TBezierSurface& initialSurface, Eigen::Vector2d rangeU, Eigen::Vector2d rangeV, uint64_t depth, std::vector<Eigen::Vector2d>& uvs)
    {

        static_assert(ClippingStrategy == ProjectionBased || ClippingStrategy == ComponentWise, "Clipping strategy not supported.");

        const ControlPointsType& controlPoints = surface.getControlPoints();

        // have we finished yet in the u or v directions? (Eq. 41)
        const bool narrowedDownEnoughU = std::abs(rangeU.y() - rangeU.x()) < (Hybrid ? mNewtonStartSize : mRequiredDomainSize);
        const bool narrowedDownEnoughV = std::abs(rangeV.y() - rangeV.x()) < (Hybrid ? mNewtonStartSize : mRequiredDomainSize);
        // Eq. 42
        const bool zeroEpsilonReached = surface.sample(Eigen::Vector2d(rangeU.mean(), rangeV.mean())).norm() < mZeroEpsilon;

        // narrowed down far enough? if we are using a non-hybrid approach, also the function value in the center must be small enough
        if (narrowedDownEnoughU && narrowedDownEnoughV && (Hybrid || zeroEpsilonReached))
        {
            Eigen::Vector2d candidate;
            bool accepted;
            if constexpr (Hybrid)
            {
                // for hybrid solvers we first do a newton refinement from the mean of the current normalized domain.
                const TBezierSurface dU = surface.generatePartialTensorProduct(PartialUV::du);
                const TBezierSurface dV = surface.generatePartialTensorProduct(PartialUV::dv);
                candidate               = newtonRefinement(surface, dU, dV);
                // the result must be in the normalized domain
                accepted = candidate.x() >= 0 && candidate.x() <= 1 && candidate.y() >= 0 && candidate.y() <= 1;
                if (accepted)
                {
                    // then we translate that to a point to the actual domain and check if it is too close to existing roots.
                    candidate.x() = (rangeU.y() - rangeU.x()) * candidate.x() + rangeU.x();
                    candidate.y() = (rangeV.y() - rangeV.x()) * candidate.y() + rangeV.x();
                    accepted      = !domainCoordinateContained(candidate, uvs);
                }
            }
            else
            {
                // for non-hybrid we use the center of the current range
                candidate.x() = rangeU.mean();
                candidate.y() = rangeV.mean();
                // and acceptance is determined by not being too close to existing roots
                accepted = !domainCoordinateContained(candidate, uvs);
            }

            if (accepted)
            {
                uvs.push_back(candidate);
            }
            return;
        }

        // not successful; recursion went too deep
        if (depth == 0)
        {
            return;
        }

        // next we must compute the new clipping ranges
        double smin;
        double smax;
        // is the clipping done in u?
        bool clipU;

        // compute U/V if they are not narrowed down enough
        // in rare cases both are narrowed down enough but we did not reach the required zero epsilon in the center. We then only clip on the longer edge.
        bool computeU = !narrowedDownEnoughU || (narrowedDownEnoughU && narrowedDownEnoughV && rangeU.y() - rangeU.x() > rangeV.y() - rangeV.x());
        bool computeV = !narrowedDownEnoughV || (narrowedDownEnoughU && narrowedDownEnoughV && rangeV.y() - rangeV.x() > rangeU.y() - rangeU.x());

        const ControlPointsType& initialSurfaceControlPoints = initialSurface.getControlPoints();

        // to setup the clipping process we must define the bezier representation of the signed distances
        // to an implicit line or a plane depending on a 3d or 2d output of our bezier surface
        // we simply use a nested array (2d) of control points containing the signed distances
        double sminU = -1.0;
        double smaxU = 2.0;
        double sminV = -1.0;
        double smaxV = 2.0;
        // check which ones we have to compute (both or just one of them?)
        // all of them essentially just compute the clipping range by calling the appropriate method
        if (computeU && computeV)
        {
            // U and V
            std::pair<std::optional<ClipRange>, std::optional<ClipRange>> computedClippingRanges;

            // Sec. 3.4 "Component-wise Clipping"
            if constexpr (ClippingStrategy == ComponentWise)
            {
                computedClippingRanges = computeComponentClipUV(controlPoints);
            }
            // Sec. 3.4 "Projection-based Clipping"
            else
            {
                computedClippingRanges = compute3ComponentDistanceClipUV(controlPoints, initialSurfaceControlPoints);
            }

            if (computedClippingRanges.first.has_value())
            {
                sminU = computedClippingRanges.first->first;
                smaxU = computedClippingRanges.first->second;
            }
            if (computedClippingRanges.second.has_value())
            {
                sminV = computedClippingRanges.second->first;
                smaxV = computedClippingRanges.second->second;
            }
        }
        else if (computeU)
        {
            // U
            std::optional<ClipRange> clippingRange;
            // Sec. 3.4 "Component-wise Clipping"
            if constexpr (ClippingStrategy == ComponentWise)
            {
                clippingRange = computeComponentClipU(controlPoints);
            }
            // Sec. 3.4 "Projection-based Clipping"
            else
            {
                clippingRange = compute3ComponentDistanceClipU(controlPoints, initialSurfaceControlPoints);
            }
            if (clippingRange.has_value())
            {
                sminU = clippingRange->first;
                smaxU = clippingRange->second;
            }
        }
        else if (computeV)
        {
            // V
            std::optional<ClipRange> clippingRange;
            // Sec. 3.4 "Component-wise Clipping"
            if constexpr (ClippingStrategy == ComponentWise)
            {
                clippingRange = computeComponentClipV(controlPoints);
            }
            // Sec. 3.4 "Projection-based Clipping"
            else
            {
                clippingRange = compute3ComponentDistanceClipV(controlPoints, initialSurfaceControlPoints);
            }
            if (clippingRange.has_value())
            {
                sminV = clippingRange->first;
                smaxV = clippingRange->second;
            }
        }
        else
        {
            return;
        }

        // if I wanted to have one of them but the clipping was not successful, then the clip failed, and there cannot be a zero crossing.
        // it is important to check first if I wanted that clip, or not.
        // < 0 and > 1 explicitly refer to the default values set for sminU, smaxU, sminV and smaxV. the clipping itself will never return
        // a value < 0 or > 1.
        const bool clipInUFailed = computeU && (sminU < 0.0 || smaxU > 1.0);
        const bool clipInVFailed = computeV && (sminV < 0.0 || smaxV > 1.0);

        // there cannot be a zero crossing anymore
        if (clipInUFailed || clipInVFailed)
        {
            return;
        }

        // clipU is true if it has the larger clip range
        // Eq. 40
        clipU = smaxU - sminU < smaxV - sminV;

        if (clipU)
        {
            smin = sminU;
            smax = smaxU;
        }
        else
        {
            smin = sminV;
            smax = smaxV;
        }

        if (smax - smin > 0.9) // reduction is less than 10% --> subdivide the surface
        {
            if (clipU)
            {
                // subdivide the surface and try recursively on the four pieces
                TBezierSurface clip0, clip1;
                surface.subdivideU(clip0, clip1, 0.5);
                clip0.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)));
                clip1.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)));
                solveRecursiveClipping<ClippingStrategy, Hybrid>(clip0, initialSurface, Eigen::Vector2d(rangeU.x(), rangeU.mean()), Eigen::Vector2d(rangeV.x(), rangeV.y()), depth - 1, uvs);
                solveRecursiveClipping<ClippingStrategy, Hybrid>(clip1, initialSurface, Eigen::Vector2d(rangeU.mean(), rangeU.y()), Eigen::Vector2d(rangeV.x(), rangeV.y()), depth - 1, uvs);
            }
            else
            {
                // subdivide the surface and try recursively on the four pieces
                TBezierSurface clip0, clip1;
                surface.subdivideV(clip0, clip1, 0.5);
                clip0.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)));
                clip1.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)));
                solveRecursiveClipping<ClippingStrategy, Hybrid>(clip0, initialSurface, Eigen::Vector2d(rangeU.x(), rangeU.y()), Eigen::Vector2d(rangeV.x(), rangeV.mean()), depth - 1, uvs);
                solveRecursiveClipping<ClippingStrategy, Hybrid>(clip1, initialSurface, Eigen::Vector2d(rangeU.x(), rangeU.y()), Eigen::Vector2d(rangeV.mean(), rangeV.y()), depth - 1, uvs);
            }
        }
        else
        {
            // adjust global range, since [smin,smax] are in [0,1] range
            if (clipU)
            {
                rangeU = Eigen::Vector2d(
                    rangeU[0] + (rangeU[1] - rangeU[0]) * smin,
                    rangeU[0] + (rangeU[1] - rangeU[0]) * smax);
            }
            else
            {
                rangeV = Eigen::Vector2d(
                    rangeV[0] + (rangeV[1] - rangeV[0]) * smin,
                    rangeV[0] + (rangeV[1] - rangeV[0]) * smax);
            }

            // subdivide if necessary
            if (smin > 0)
            {
                if (clipU)
                    surface = surface.template clipU<true>(smin);
                else
                    surface = surface.template clipV<true>(smin);
            }
            if (smax < 1)
            {
                if (clipU)
                    surface = surface.template clipU<false>(smax);
                else
                    surface = surface.template clipV<false>(smax);
            }
            surface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)));
            solveRecursiveClipping<ClippingStrategy, Hybrid>(surface, initialSurface, rangeU, rangeV, depth - 1, uvs);
        }
    }

    /**
     * @brief Searches the given bezier surface recursively for roots using bezier subdivision.
     *
     * @tparam Hybrid Set to true if after subdivision down to mNewtonStartSize a newton solver should take over. If set
     * to false subdivision will happen until mRefineEpsilon.
     *
     * @param surface The surface the roots are searched in.
     * @param depth The current recursion depth (if zero is reached the recursion stops).
     * @param uvs A list of the the uv coordinates where roots are found.
     */
    template <bool Hybrid = false>
    void solveRecursiveBisection(const TBezierSurface& surface, int depth, std::vector<Eigen::Vector2d>& uvs)
    {
        const ControlPointsType& controlPoints = surface.getControlPoints();

        if (!hasZeroCrossing(controlPoints))
        {
            return;
        }

        // have we finished yet in the u or v directions?
        const Eigen::AlignedBox2d& domain = surface.getDomain();
        const bool narrowedDownEnoughU    = std::abs(domain.max().x() - domain.min().x()) < (Hybrid ? mNewtonStartSize : mRequiredDomainSize);
        const bool narrowedDownEnoughV    = std::abs(domain.max().y() - domain.min().y()) < (Hybrid ? mNewtonStartSize : mRequiredDomainSize);
        const bool zeroEpsilonReached     = surface.sample(domain.center()).norm() < mZeroEpsilon;

        // quad small enough? terminate
        if (narrowedDownEnoughU && narrowedDownEnoughV && (Hybrid || zeroEpsilonReached))
        {
            Eigen::Vector2d candidate;
            bool accepted;
            if constexpr (Hybrid)
            {
                // for hybrid solvers we first do a newton refinement from the mean of the current normalized domain.
                const TBezierSurface dU = surface.generatePartialTensorProduct(PartialUV::du);
                const TBezierSurface dV = surface.generatePartialTensorProduct(PartialUV::dv);
                candidate               = newtonRefinement(surface, dU, dV);
                // result must be in normalized domain
                accepted = candidate.x() >= 0 && candidate.x() <= 1 && candidate.y() >= 0 && candidate.y() <= 1;
                if (accepted)
                {
                    // then we translate that to a point to the actual domain and check if it is too close to existing roots.
                    candidate.x() = (domain.max().x() - domain.min().x()) * candidate.x() + domain.min().x();
                    candidate.y() = (domain.max().y() - domain.min().y()) * candidate.y() + domain.min().y();
                    accepted      = !domainCoordinateContained(candidate, uvs);
                }
            }
            else
            {
                // for non-hybrid we use the center of the current domain
                candidate = domain.center();
                // and acceptance is determined by not being too close to existing roots
                accepted = !domainCoordinateContained(candidate, uvs);
            }

            if (accepted)
            {
                uvs.push_back(candidate);
            }
            return;
        }

        // exit when the recursion limit is reached
        if (depth == 0)
        {
            return;
        }

        // there is a crossing -> subdivide...
        bool subdivideU = (narrowedDownEnoughV && !narrowedDownEnoughU) || (depth % 2 == 0 && (!narrowedDownEnoughU || narrowedDownEnoughV));
        if (subdivideU)
        {
            TBezierSurface left;
            TBezierSurface right;
            surface.subdivideU(left, right, domain.center().x());
            solveRecursiveBisection<Hybrid>(left, depth - 1, uvs);
            solveRecursiveBisection<Hybrid>(right, depth - 1, uvs);
        }
        else
        {
            TBezierSurface left;
            TBezierSurface right;
            surface.subdivideV(left, right, domain.center().y());
            solveRecursiveBisection<Hybrid>(left, depth - 1, uvs);
            solveRecursiveBisection<Hybrid>(right, depth - 1, uvs);
        }
    }

    /**
     * @brief Computes and returns the pseudo inverse of given matrix.
     *
     * @param a
     * @param epsilon
     * @return
     */
    template <int RowsIn, int ColsIn>
    Eigen::Matrix<double, ColsIn, RowsIn> pseudoInverse(const Eigen::Matrix<double, RowsIn, ColsIn>& a,
                                                        const double& zeroEpsilon = 1e-10,
                                                        const double& epsilon     = std::numeric_limits<double>::epsilon())
    {
        auto ata = a.transpose() * a;
        if (ata.determinant() > zeroEpsilon)
        {
            return ata.ldlt().solve(a.transpose());
        }
        using MatrixType        = Eigen::Matrix<double, RowsIn, ColsIn>;
        using PseudoInverseType = Eigen::Matrix<double, ColsIn, RowsIn>;

        Eigen::JacobiSVD<MatrixType> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::VectorXd singularValues = svd.singularValues();

        const double tolerance = epsilon * std::max(a.cols(), a.rows()) *
                                 singularValues.maxCoeff();

        PseudoInverseType invertedSingularValues;
        invertedSingularValues.setZero();

        for (Eigen::Index index = 0; index < std::min(ColsIn, RowsIn); ++index)
        {
            double value = singularValues(index);
            if (value > tolerance)
            {
                invertedSingularValues(index, index) = 1.0 / value;
            }
        }

        return (svd.matrixV() * invertedSingularValues * svd.matrixU().adjoint()).eval();
    }

    /**
     * @brief Performs newton iteration steps to refine the given candidate towards the root.
     *
     * @return
     */
    DomainCoord newtonRefinement(TBezierSurface bezierSurface, const TBezierSurface& dU, const TBezierSurface& dV, const DomainCoord& startPoint = DomainCoord(0.5, 0.5))
    {
        // we always start at the center of the normalized domain
        DomainCoord x(0.5, 0.5);
        bezierSurface.setDomain(Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)));
        const Eigen::AlignedBox2d& domain = bezierSurface.getDomain();

        int attempts        = 0;
        double stepSize     = 0.5;
        Value functionValue = bezierSurface.sample(x);
        Eigen::Matrix<double, Components, 2> jacobian;
        jacobian.col(0)                                     = dU.sample(x);
        jacobian.col(1)                                     = dV.sample(x);
        Eigen::Matrix<double, 2, Components> pseudoInverseJ = (jacobian.transpose() * jacobian).ldlt().solve(jacobian.transpose());
        double error                                        = functionValue.norm();

        double newError;
        DomainCoord newX;
        Value newFunctionValue;

        while (error > mZeroEpsilon && attempts < mMaxNewtonIterations)
        {

            // if progress is far too small, terminate since this might be a local minimum
            // if (jacobian.squaredNorm() < 1E-26 || stepSize < 1E-5)
            // {
            //     return DomainCoord(-1, -1);
            // }

            if constexpr (Components == 2)
            {
                newX = x - jacobian.fullPivHouseholderQr().solve(functionValue) * stepSize;
            }
            else
            {
                newX = x - (pseudoInverseJ * functionValue) * stepSize;
            }

            newFunctionValue = bezierSurface.sample(newX);
            newError         = newFunctionValue.norm();

            if (newError > error || newX.x() < domain.min().x() || newX.x() > domain.max().x() || newX.y() < domain.min().y() || newX.y() > domain.max().y())
            {
                stepSize *= 0.5;
            }
            else
            {
                stepSize        = std::max(stepSize * 2.0, 1.0);
                x               = newX;
                error           = newError;
                functionValue   = newFunctionValue;
                jacobian.col(0) = dU.sample(x);
                jacobian.col(1) = dV.sample(x);
                pseudoInverseJ  = (jacobian.transpose() * jacobian).ldlt().solve(jacobian.transpose());
            }
            attempts++;
        }

        if (error < mZeroEpsilon)
        {
            return x;
        }
        return DomainCoord(-1, -1);
    }
};
