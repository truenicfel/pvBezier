#pragma once

#include "bezier_surface.hpp"
#include "binomial_table.hpp"

    /**
     * @brief Contains static functions to compute the bezier surface from the cross product of two bezier surfaces.
     *
     * @details This implements equation 8 and 9.
     */
    template <typename TBezierSurfaceV, typename TBezierSurfaceW>
    class BezierCrossQuad
    {
    public:
        static constexpr int NV = TBezierSurfaceV::N;
        static constexpr int MV = TBezierSurfaceV::M;
        static constexpr int NW = TBezierSurfaceW::N;
        static constexpr int MW = TBezierSurfaceW::M;
        static constexpr int NC = NV + NW;
        static constexpr int MC = MV + MW;
        using Value             = typename TBezierSurfaceV::Value;

        using SurfaceOut = BezierSurface<NC, MC, Value>;
        using SurfaceInV = TBezierSurfaceV;
        using SurfaceInW = TBezierSurfaceW;

        BezierCrossQuad()
        {
            buildWeightTable();
        }

        SurfaceOut compute(
            const BezierSurface<NV, MV, Value>& v,
            const BezierSurface<NW, MW, Value>& w)
        {
            // inside compute(v,w):
            auto const& V = v.getControlPoints(); // [0..NV][0..MV]
            auto const& W = w.getControlPoints(); // [0..NW][0..MW]

            // 1) zero the output
            std::array<std::array<Value, MC + 1>, NC + 1> C;
            for (int i = 0; i <= NC; ++i)
                for (int j = 0; j <= MC; ++j)
                    C[i][j].setZero();

            // 2) the four‐deep loop
            for (int k = 0; k <= NV; ++k)
                for (int l = 0; l <= MV; ++l)
                {
                    // compute each control point c_ij (Eq. 9) (here: c_kl)
                    // Additional Material Algorithm 1 ("Compute control points")
                    const Value& Vkl = V[k][l];
                    for (int iW = 0; iW <= NW; ++iW)
                        for (int jW = 0; jW <= MW; ++jW)
                        {
                            int iC     = k + iW;
                            int jC     = l + jW;
                            double wgt = _weight[k][l][iW][jW];
                            // one single cross() and one vector‐scale & add:
                            Value cr = Vkl.cross(W[iW][jW]);
                            C[iC][jC] += wgt * cr;
                        }
                }

            // 3) wrap and return
            return BezierSurface<NC, MC, Value>(C);
        }

    private:
        std::array<std::array<std::array<std::array<double, MW + 1>, NW + 1>, MV + 1>, NV + 1>
            _weight;

        void buildWeightTable()
        {
            // Additional Material Algorithm 1 ("Precompute weights")
            auto const& bin = BinomialTable<10>::table;

            for (int k = 0; k <= NV; ++k)
                for (int l = 0; l <= MV; ++l)
                {
                    for (int iW = 0; iW <= NW; ++iW)
                        for (int jW = 0; jW <= MW; ++jW)
                        {
                            int iC = k + iW;
                            int jC = l + jW;
                            double numer =
                                bin[NV][k] * bin[NW][iW] * bin[MV][l] * bin[MW][jW];
                            double denom =
                                bin[NC][iC] * bin[MC][jC];
                            _weight[k][l][iW][jW] = numer / denom;
                        }
                }
        }
    };
