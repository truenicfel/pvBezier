#pragma once

#include "bezier_volume.hpp"
#include "binomial_table.hpp"

/**
 * @brief Contains static functions to compute the acceleration volume of a bezier volume using the Jacobian of the volume.
 */
template <typename TBezierVolume>
class BezierAcceleration
{
public:
    // make sure that the bezier volumes have the same degree in all dimensions
    static_assert(TBezierVolume::N == TBezierVolume::M, "N, M and O have to be equal!");
    static_assert(TBezierVolume::M == TBezierVolume::O, "N, M and O have to be equal!");

    // pull in your degrees
    static constexpr int Degree  = TBezierVolume::N;
    static constexpr int Degree2 = 2 * Degree;

    using Value        = typename TBezierVolume::Value;
    using ResultVolume = BezierVolume<2 * Degree, 2 * Degree, 2 * Degree, Value>;

    static constexpr int D  = TBezierVolume::N;
    static constexpr int D2 = 2 * D;
    using InputVol          = BezierVolume<D, D, D, Value>;
    using ResultVol         = BezierVolume<D2, D2, D2, Value>;

    BezierAcceleration()
    {
        buildWeightTables(); // fill weightU_,weightV_,weightW_ as before
    }

    /**
     * @brief Computes the acceleration bezier volume for the given velocity bezier volume and the bezier volumes du, dv and dw (partials of the velocity bezier volume).
     *
     * @tparam TBezierVolume with degree n1, n2, n3 (N, M, O) this will determine the degree of the
     * parameters (n1, n2, n3) and the degree of the resulting acceleration volume (2 * n1, 2 * n2, 2 * n3).
     *
     * @param v The velocity bezier volume.
     * @param du The u partial derivative of the velocity bezier volume.
     * @param dv The v partial derivative of the velocity bezier volume.
     * @param dw The w partial derivative of the velocity bezier volume.
     *
     * @return The acceleration as a bezier volume.
     */
    ResultVolume compute(
        const TBezierVolume& v,
        const TBezierVolume& du,
        const TBezierVolume& dv,
        const TBezierVolume& dw)
    {
        // grab references to input control grids
        auto const& V  = v.getControlPoints();
        auto const& DU = du.getControlPoints();
        auto const& DV = dv.getControlPoints();
        auto const& DW = dw.getControlPoints();

        // Precompute w-direction diagonal weights once (k fixed, K = k + r)
        std::array<std::array<double, Degree + 1>, Degree + 1> wW_diag;
        for (int k = 0; k <= Degree; ++k)
            for (int r = 0; r <= Degree; ++r)
                wW_diag[k][r] = weightW_[k][r][k + r];

        // Prepare output (unchanged)
        ResultVolume result;
        auto& C = result.getControlPoints();
        for (int I = 0; I <= Degree2; ++I)
            for (int J = 0; J <= Degree2; ++J)
                for (int K = 0; K <= Degree2; ++K)
                    C[I][J][K].setZero();

        // Compute acceleration control points (Eq. 27 and Additional Material Algorithm 3)
        // Six-deep convolution with micro-opts and diagonal W cache
        for (int i = 0; i <= Degree; ++i)
            for (int j = 0; j <= Degree; ++j)
                for (int k = 0; k <= Degree; ++k)
                {

                    // Cache Jacobian at (i,j,k)
                    const auto& du_ijk = DU[i][j][k];
                    const auto& dv_ijk = DV[i][j][k];
                    const auto& dw_ijk = DW[i][j][k];

                    // Alias the diagonal weight row for this k (contiguous in r)
                    const double* wWk = wW_diag[k].data();

                    for (int p = 0; p <= Degree; ++p)
                    {
                        const int I     = i + p;
                        const double wU = weightU_[i][p][I];

                        for (int q = 0; q <= Degree; ++q)
                        {
                            const int J      = j + q;
                            const double wUV = wU * weightV_[j][q][J];

                            // Pointer to the contiguous K-line at output
                            auto* rowC = &C[I][J][0];

                            // Inner loop over r rewritten as a contiguous loop over K
                            // r = 0..Degree, K = k + r (runs from k..k+Degree)
                            for (int r = 0; r <= Degree; ++r)
                            {
                                const int K = k + r;

                                // Read velocity
                                const auto& vel = V[p][q][r];

                                // Final weight with diagonal cache
                                const double w = wUV * wWk[r];

                                rowC[K] += w * (du_ijk * vel.x() + dv_ijk * vel.y() + dw_ijk * vel.z());
                            }
                        }
                    }
                }

        return result; // Eq. 26
    }

private:
    std::array<std::array<std::array<double, Degree2 + 1>, Degree + 1>, Degree + 1> weightU_;
    std::array<std::array<std::array<double, Degree2 + 1>, Degree + 1>, Degree + 1> weightV_;
    std::array<std::array<std::array<double, Degree2 + 1>, Degree + 1>, Degree + 1> weightW_;

    void buildWeightTables()
    {
        // Additional Material Algorithm 3 ("Precompute weights")
        auto const& bin = BinomialTable<10>::table;

        // U-direction
        for (int i = 0; i <= Degree; ++i)
            for (int p = 0; p <= Degree; ++p)
            {
                double numUP = bin[Degree][i] * bin[Degree][p];
                for (int I = i; I <= i + Degree; ++I)
                {
                    weightU_[i][p][I] = numUP / bin[Degree2][I];
                }
            }

        // V-direction
        for (int j = 0; j <= Degree; ++j)
            for (int q = 0; q <= Degree; ++q)
            {
                double numVQ = bin[Degree][j] * bin[Degree][q];
                for (int Jv = j; Jv <= j + Degree; ++Jv)
                {
                    weightV_[j][q][Jv] = numVQ / bin[Degree2][Jv];
                }
            }

        // W-direction
        for (int k = 0; k <= Degree; ++k)
            for (int r = 0; r <= Degree; ++r)
            {
                double numKR = bin[Degree][k] * bin[Degree][r];
                for (int K = k; K <= k + Degree; ++K)
                {
                    weightW_[k][r][K] = numKR / bin[Degree2][K];
                }
            }
    }
};
