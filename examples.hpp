#pragma once

#include "pvBezier/bezier_cross_quad.hpp"
#include "pvBezier/bezier_surface.hpp"
#include "pvBezier/bezier_surface_root_finding.hpp"
#include "pvBezier/bezier_volume.hpp"
#include "pvBezier/catmull_rom_interpolant.hpp"
#include "pvBezier/bezier_acceleration.hpp"

#include <chrono>

inline void expect_eq(int a, int b, const char* exprA, const char* exprB, const char* file, int line)
{
    if (a != b)
    {
        std::cerr << file << ":" << line
                  << " EXPECT_EQ failed: "
                  << exprA << " (" << a << ") != "
                  << exprB << " (" << b << ")"
                  << std::endl;
    }
}

inline void expect_near(double a, double b, double epsilon,
                 const char* exprA, const char* exprB, const char* exprEps,
                 const char* file, int line)
{
    if (std::fabs(a - b) > epsilon)
    {
        std::cerr << file << ":" << line
                  << " EXPECT_NEAR failed: |"
                  << exprA << " (" << a << ") - "
                  << exprB << " (" << b << ")| > "
                  << exprEps << " (" << epsilon << ")"
                  << std::endl;
    }
}

inline void expect_true(bool condition,
                 const char* expr,
                 const char* file,
                 int line)
{
    if (!condition)
    {
        std::cerr << file << ":" << line
                  << " EXPECT_TRUE failed: "
                  << expr << " is false"
                  << std::endl;
    }
}

#define EXPECT_EQ(a, b) \
    expect_eq((a), (b), #a, #b, __FILE__, __LINE__)

#define EXPECT_NEAR(a, b, eps) \
    expect_near((a), (b), (eps), #a, #b, #eps, __FILE__, __LINE__)

#define EXPECT_TRUE(expr) \
    expect_true((expr), #expr, __FILE__, __LINE__)

inline void test_solvers()
{
    BiquadraticBezierSurface3d bezierSurface;
    BiquadraticBezierSurface3d::ControlPointsType& controlPoints = bezierSurface.getControlPoints();

    controlPoints[0][0] = Eigen::Vector3d(-1378.8105372499039731337688863277435302734375, -610.680763660446473295451141893863677978515625, 673.19925631300293389358557760715484619140625);
    controlPoints[0][1] = Eigen::Vector3d(-222.715217705232134903781116008758544921875, 142.674505273780596326105296611785888671875, 573.9129224513955023212474770843982696533203125);
    controlPoints[0][2] = Eigen::Vector3d(547.4650942466035985489725135266780853271484375, 1503.077316870347203803248703479766845703125, 345.51448186154794939284329302608966827392578125);
    controlPoints[1][0] = Eigen::Vector3d(-1510.3109463750670329318381845951080322265625, -353.06525097172607274842448532581329345703125, 829.691005278495367747382260859012603759765625);
    controlPoints[1][1] = Eigen::Vector3d(-20.7937261363476437736608204431831836700439453125, -292.9328639675362637717626057565212249755859375, 575.04286343058765851310454308986663818359375);
    controlPoints[1][2] = Eigen::Vector3d(444.26892914429339498383342288434505462646484375, 762.4615891636275364362518303096294403076171875, 52.916172788771206114688538946211338043212890625);
    controlPoints[2][0] = Eigen::Vector3d(-1592.87714667968020876287482678890228271484375, -108.4932429889022387214936316013336181640625, 1009.2418775208938086507259868085384368896484375);
    controlPoints[2][1] = Eigen::Vector3d(227.77883934492382422831724397838115692138671875, -785.4896864954997681707027368247509002685546875, 569.940535021738241994171403348445892333984375);
    controlPoints[2][2] = Eigen::Vector3d(-24.2196107642677844751233351416885852813720703125, 359.67972678065581249029492028057575225830078125, -317.663892151538675534538924694061279296875);

    const double epsilon = 1E-5;
    const Eigen::Vector2d uvExpect(0.9951660204957602, 0.81358986532549682);

    BezierSurfaceRootFinding<BiquadraticBezierSurface3d> bezierRootFinding(1e-10, 1e-10, 1e-3, 300, 1e-3, 1024);

    std::vector<Eigen::Vector2d> uvs;
    // Fig. 4, Sec. 3.4: "Component-wise Clipping"
    bezierRootFinding.clippingSolver<ComponentWise>(bezierSurface, uvs);

    EXPECT_EQ(uvs.size(), 1);

    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs[0];
        if (uvs.size() == 1)
        {
            EXPECT_NEAR(uvExpect.x(), uvHave.x(), epsilon);
            EXPECT_NEAR(uvExpect.y(), uvHave.y(), epsilon);
        }
    }

    uvs.clear();
    // Fig. 5, Sec. 3.4: "Projection-based Clipping"
    bezierRootFinding.clippingSolver<ProjectionBased>(bezierSurface, uvs);

    EXPECT_EQ(uvs.size(), 1);

    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs[0];
        if (uvs.size() == 1)
        {
            EXPECT_NEAR(uvExpect.x(), uvHave.x(), epsilon);
            EXPECT_NEAR(uvExpect.y(), uvHave.y(), epsilon);
        }
    }

    uvs.clear();
    // Bezier Bisection (Fig. 2, Sec. 2.1 "Bezier Bisection")
    bezierRootFinding.bisectionSolver(bezierSurface, uvs);

    EXPECT_EQ(uvs.size(), 1);

    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs[0];
        if (uvs.size() == 1)
        {
            EXPECT_NEAR(uvExpect.x(), uvHave.x(), epsilon);
            EXPECT_NEAR(uvExpect.y(), uvHave.y(), epsilon);
        }
    }

    uvs.clear();
    // Newton (Sec. 2.2 "Local Techniques")
    bezierRootFinding.newtonSolver(bezierSurface, uvs, Eigen::Vector2d(0.5, 0.5));

    EXPECT_EQ(uvs.size(), 1);

    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs[0];
        if (uvs.size() == 1)
        {
            EXPECT_NEAR(uvExpect.x(), uvHave.x(), epsilon);
            EXPECT_NEAR(uvExpect.y(), uvHave.y(), epsilon);
        }
    }

    uvs.clear();
    // Hybrid Subdivision (Sec 4.3)
    bezierRootFinding.hybridBisectionSolver(bezierSurface, uvs);

    EXPECT_EQ(uvs.size(), 1);

    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs[0];
        if (uvs.size() == 1)
        {
            EXPECT_NEAR(uvExpect.x(), uvHave.x(), epsilon);
            EXPECT_NEAR(uvExpect.y(), uvHave.y(), epsilon);
        }
    }

    uvs.clear();
    // Hybrid Clipping (Sec 4.3)
    bezierRootFinding.hybridClippingSolver<ProjectionBased>(bezierSurface, uvs);

    EXPECT_EQ(uvs.size(), 1);

    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs[0];
        if (uvs.size() == 1)
        {
            EXPECT_NEAR(uvExpect.x(), uvHave.x(), epsilon);
            EXPECT_NEAR(uvExpect.y(), uvHave.y(), epsilon);
        }
    }
}

inline void test_solver_performance()
{
    BiquadraticBezierSurface3d bezierSurface;
    BiquadraticBezierSurface3d::ControlPointsType& controlPoints = bezierSurface.getControlPoints();

    controlPoints[0][0] = Eigen::Vector3d(-1378.8105372499039731337688863277435302734375, -610.680763660446473295451141893863677978515625, 673.19925631300293389358557760715484619140625);
    controlPoints[0][1] = Eigen::Vector3d(-222.715217705232134903781116008758544921875, 142.674505273780596326105296611785888671875, 573.9129224513955023212474770843982696533203125);
    controlPoints[0][2] = Eigen::Vector3d(547.4650942466035985489725135266780853271484375, 1503.077316870347203803248703479766845703125, 345.51448186154794939284329302608966827392578125);
    controlPoints[1][0] = Eigen::Vector3d(-1510.3109463750670329318381845951080322265625, -353.06525097172607274842448532581329345703125, 829.691005278495367747382260859012603759765625);
    controlPoints[1][1] = Eigen::Vector3d(-20.7937261363476437736608204431831836700439453125, -292.9328639675362637717626057565212249755859375, 575.04286343058765851310454308986663818359375);
    controlPoints[1][2] = Eigen::Vector3d(444.26892914429339498383342288434505462646484375, 762.4615891636275364362518303096294403076171875, 52.916172788771206114688538946211338043212890625);
    controlPoints[2][0] = Eigen::Vector3d(-1592.87714667968020876287482678890228271484375, -108.4932429889022387214936316013336181640625, 1009.2418775208938086507259868085384368896484375);
    controlPoints[2][1] = Eigen::Vector3d(227.77883934492382422831724397838115692138671875, -785.4896864954997681707027368247509002685546875, 569.940535021738241994171403348445892333984375);
    controlPoints[2][2] = Eigen::Vector3d(-24.2196107642677844751233351416885852813720703125, 359.67972678065581249029492028057575225830078125, -317.663892151538675534538924694061279296875);

    BezierSurfaceRootFinding<BiquadraticBezierSurface3d> bezierRootFinding(1e-10, 1e-10, 1e-3, 300, 1e-3, 1024);

    std::vector<Eigen::Vector2d> uvs;
    // Clipping Component Wise
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10000; i++)
    {
        bezierRootFinding.clippingSolver<ComponentWise>(bezierSurface, uvs);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const long clippingTimeComponentWise      = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Clipping (Component Wise) Time = " << clippingTimeComponentWise << "[ms]" << std::endl;

    // Clipping Distance Based
    uvs.clear();
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10000; i++)
    {
        bezierRootFinding.clippingSolver<ProjectionBased>(bezierSurface, uvs);
    }
    end                            = std::chrono::steady_clock::now();
    long clippingTimeDistanceBased = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Clipping (Distance Based) Time = " << clippingTimeDistanceBased << "[ms]" << std::endl;

    // Subdivision
    uvs.clear();
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10000; i++)
    {
        bezierRootFinding.bisectionSolver(bezierSurface, uvs);
    }
    end                = std::chrono::steady_clock::now();
    long bisectionTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Bisection Time = " << bisectionTime << "[ms]" << std::endl;

    // Newton
    uvs.clear();
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10000; i++)
    {
        bezierRootFinding.bisectionSolver(bezierSurface, uvs);
    }
    end             = std::chrono::steady_clock::now();
    long newtonTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Newton Time = " << newtonTime << "[ms]" << std::endl;

    // Hybrid Clipping
    uvs.clear();
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10000; i++)
    {
        bezierRootFinding.hybridClippingSolver(bezierSurface, uvs);
    }
    end                     = std::chrono::steady_clock::now();
    long hybridClippingTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Hybrid Clipping Time = " << hybridClippingTime << "[ms]" << std::endl;

    // Hybrid Subdivision
    uvs.clear();
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10000; i++)
    {
        bezierRootFinding.hybridBisectionSolver(bezierSurface, uvs);
    }
    end                      = std::chrono::steady_clock::now();
    long hybridBisectionTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Hybrid Bisection Time = " << hybridBisectionTime << "[ms]" << std::endl;

    std::cout << "Distance-based Clipping speedup over..." << std::endl;
    std::cout << "--- Component Wise: " << static_cast<double>(clippingTimeComponentWise) / static_cast<double>(clippingTimeDistanceBased) << std::endl;
    std::cout << "--- Bisection: " << static_cast<double>(bisectionTime) / static_cast<double>(clippingTimeDistanceBased) << std::endl;
    std::cout << "--- Newton: " << static_cast<double>(newtonTime) / static_cast<double>(clippingTimeDistanceBased) << std::endl;
    std::cout << "--- Hybrid Bisection: " << static_cast<double>(hybridBisectionTime) / static_cast<double>(clippingTimeDistanceBased) << std::endl;
    std::cout << "--- Hybrid Clipping: " << static_cast<double>(hybridClippingTime) / static_cast<double>(clippingTimeDistanceBased) << std::endl;
}

inline void test_voxel_trilinear()
{
    // First, the bezier volumes must be initialized. The trilinear interpolation method allows us to
    // use the grid points as they are as control points to our Bezier Volumes. (Eq. 6-7)

    TrilinearBezierVolume3d v;
    TrilinearBezierVolume3d::ControlPointsType& vControlPoints = v.getControlPoints();

    vControlPoints[0][0][0] = Eigen::Vector3d(0.716664, -0.465993, -0.219746);
    vControlPoints[1][0][0] = Eigen::Vector3d(0.811002, -0.36816, -0.267679);
    vControlPoints[0][1][0] = Eigen::Vector3d(0.913761, -0.20279, -0.196658);
    vControlPoints[1][1][0] = Eigen::Vector3d(0.955215, -0.114034, -0.0584682);
    vControlPoints[0][0][1] = Eigen::Vector3d(0.740329, -0.309348, -0.154075);
    vControlPoints[1][0][1] = Eigen::Vector3d(0.781149, -0.427953, -0.359684);
    vControlPoints[0][1][1] = Eigen::Vector3d(0.818845, -0.521896, -0.238994);
    vControlPoints[1][1][1] = Eigen::Vector3d(0.777773, -0.57749, -0.00832442);

    TrilinearBezierVolume3d w;
    TrilinearBezierVolume3d::ControlPointsType& wControlPoints = w.getControlPoints();

    wControlPoints[0][0][0] = Eigen::Vector3d(-454.503, -373.697, 599.975);
    wControlPoints[1][0][0] = Eigen::Vector3d(-210.043, 195.978, 1069.07);
    wControlPoints[0][1][0] = Eigen::Vector3d(-1.24326, 625.098, 837.562);
    wControlPoints[1][1][0] = Eigen::Vector3d(-33.0282, 405.66, 880.974);
    wControlPoints[0][0][1] = Eigen::Vector3d(6.8094, -418.767, -883.431);
    wControlPoints[1][0][1] = Eigen::Vector3d(179.194, 84.3214, -128.959);
    wControlPoints[0][1][1] = Eigen::Vector3d(35.6774, 788.756, -16.9819);
    wControlPoints[1][1][1] = Eigen::Vector3d(-87.8368, 512.57, -452.363);

    // We now have Eq. 6 & 7

    // Test setup:
    std::vector<Eigen::Vector2d> uvs;
    const double epsilon = 1e-4;

    // Cross product and root finding setup:
    BezierCrossQuad<BilinearBezierSurface3d, BilinearBezierSurface3d> crossQuad;
    BezierSurfaceRootFinding<BiquadraticBezierSurface3d> bezierRootFinding(1e-10, 1e-10, 1e-3, 300, 1e-3, 1024);

    // Bottom Face
    BilinearBezierSurface3d vSurface;
    v.getBoundarySurfaceVMin(vSurface);
    BilinearBezierSurface3d wSurface;
    w.getBoundarySurfaceVMin(wSurface);
    // (Eq. 8-9)
    BiquadraticBezierSurface3d crossProduct = crossQuad.compute(vSurface, wSurface);
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }
    // expect one root on this face
    Eigen::Vector2d expect = Eigen::Vector2d(0.748314, 0.239626);
    EXPECT_EQ(uvs.size(), 1);
    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs.at(0);
        EXPECT_NEAR(expect.x(), uvHave.x(), epsilon);
        EXPECT_NEAR(expect.y(), uvHave.y(), epsilon);
    }
    // cleanup
    uvs.clear();

    // Left Face
    v.getBoundarySurfaceUMin(vSurface);
    w.getBoundarySurfaceUMin(wSurface);
    // (Eq. 8-9)
    crossProduct = crossQuad.compute(vSurface, wSurface);
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }

    // expect one root on this face
    expect = Eigen::Vector2d(0.455695, 0.407094);
    EXPECT_EQ(uvs.size(), 1);
    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs.at(0);
        EXPECT_NEAR(expect.x(), uvHave.x(), epsilon);
        EXPECT_NEAR(expect.y(), uvHave.y(), epsilon);
    }
    // cleanup
    uvs.clear();

    // Front Face
    v.getBoundarySurfaceWMin(vSurface);
    w.getBoundarySurfaceWMin(wSurface);
    // (Eq. 8-9)
    crossProduct = crossQuad.compute(vSurface, wSurface);
    crossProduct.recomputeBoundingBox();
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }

    // expect no root on this face
    EXPECT_EQ(uvs.size(), 0);
}

inline void test_voxel_tricubic()
{
    // Here, it is not possible to use the grid points directly as the control points because
    // a tricubic interpolation is used. The grid points are stored in vCells and wCells.

    CatmullRomInterpolant::Cells<Eigen::Vector3d> vCells;

    vCells[0][0][0] = Eigen::Vector3d(0.384905, -0.331021, 0.393825);
    vCells[1][0][0] = Eigen::Vector3d(0.348277, -0.413362, 0.365962);
    vCells[2][0][0] = Eigen::Vector3d(0.326964, -0.4534, 0.255795);
    vCells[3][0][0] = Eigen::Vector3d(0.341653, -0.42509, 0.155603);
    vCells[0][1][0] = Eigen::Vector3d(0.329087, -0.450383, 0.230793);
    vCells[1][1][0] = Eigen::Vector3d(0.340945, -0.420839, 0.149053);
    vCells[2][1][0] = Eigen::Vector3d(0.360521, -0.382515, 0.0886838);
    vCells[3][1][0] = Eigen::Vector3d(0.382397, -0.350426, 0.0549903);
    vCells[0][2][0] = Eigen::Vector3d(0.360575, -0.377413, 0.0833482);
    vCells[1][2][0] = Eigen::Vector3d(0.383187, -0.342312, 0.0525134);
    vCells[2][2][0] = Eigen::Vector3d(0.194541, -0.194253, 0.0180661);
    vCells[3][2][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000172922);
    vCells[0][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000228108);
    vCells[1][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000209707);
    vCells[2][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000191317);
    vCells[3][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000172921);
    vCells[0][0][1] = Eigen::Vector3d(0.398916, -0.288387, 0.648467);
    vCells[1][0][1] = Eigen::Vector3d(0.407706, -0.336648, 0.525273);
    vCells[2][0][1] = Eigen::Vector3d(0.464419, -0.342912, 0.205989);
    vCells[3][0][1] = Eigen::Vector3d(0.522456, -0.35406, -0.131747);
    vCells[0][1][1] = Eigen::Vector3d(0.476781, -0.359284, 0.164483);
    vCells[1][1][1] = Eigen::Vector3d(0.533955, -0.347191, -0.163723);
    vCells[2][1][1] = Eigen::Vector3d(0.604242, -0.2743, -0.199436);
    vCells[3][1][1] = Eigen::Vector3d(0.67171, -0.167851, -0.159466);
    vCells[0][2][1] = Eigen::Vector3d(0.627536, -0.248802, -0.207454);
    vCells[1][2][1] = Eigen::Vector3d(0.680804, -0.15109, -0.146521);
    vCells[2][2][1] = Eigen::Vector3d(0.711689, -0.0849616, -0.0435621);
    vCells[3][2][1] = Eigen::Vector3d(0.703087, -0.0777879, 0.0506719);
    vCells[0][3][1] = Eigen::Vector3d(0.72153, -0.0838632, 0.00731354);
    vCells[1][3][1] = Eigen::Vector3d(0.697467, -0.106539, 0.0691244);
    vCells[2][3][1] = Eigen::Vector3d(0.646464, -0.168068, 0.101593);
    vCells[3][3][1] = Eigen::Vector3d(0.5826, -0.255839, 0.108278);
    vCells[0][0][2] = Eigen::Vector3d(0.484601, -0.0884814, 0.685874);
    vCells[1][0][2] = Eigen::Vector3d(0.506102, -0.10573, 0.552217);
    vCells[2][0][2] = Eigen::Vector3d(0.518997, -0.155266, 0.220147);
    vCells[3][0][2] = Eigen::Vector3d(0.531333, -0.233848, -0.102992);
    vCells[0][1][2] = Eigen::Vector3d(0.532968, -0.148716, 0.188138);
    vCells[1][1][2] = Eigen::Vector3d(0.551587, -0.230482, -0.114794);
    vCells[2][1][2] = Eigen::Vector3d(0.582001, -0.318849, -0.267985);
    vCells[3][1][2] = Eigen::Vector3d(0.60904, -0.358736, -0.205413);
    vCells[0][2][2] = Eigen::Vector3d(0.597367, -0.332034, -0.264587);
    vCells[1][2][2] = Eigen::Vector3d(0.610086, -0.388842, -0.178064);
    vCells[2][2][2] = Eigen::Vector3d(0.579485, -0.430263, -0.00620217);
    vCells[3][2][2] = Eigen::Vector3d(0.542953, -0.474194, 0.163135);
    vCells[0][3][2] = Eigen::Vector3d(0.572444, -0.468671, 0.0727383);
    vCells[1][3][2] = Eigen::Vector3d(0.544172, -0.536679, 0.229456);
    vCells[2][3][2] = Eigen::Vector3d(0.519358, -0.588139, 0.334183);
    vCells[3][3][2] = Eigen::Vector3d(0.510477, -0.631277, 0.377793);
    vCells[0][0][3] = Eigen::Vector3d(0.551293, 0.0301241, 0.642897);
    vCells[1][0][3] = Eigen::Vector3d(0.551752, 0.0384976, 0.584381);
    vCells[2][0][3] = Eigen::Vector3d(0.538636, -0.00932551, 0.385826);
    vCells[3][0][3] = Eigen::Vector3d(0.535415, -0.11538, 0.113912);
    vCells[0][1][3] = Eigen::Vector3d(0.583807, 0.00105604, 0.377048);
    vCells[1][1][3] = Eigen::Vector3d(0.55558, -0.120834, 0.109389);
    vCells[2][1][3] = Eigen::Vector3d(0.56956, -0.270999, -0.0659052);
    vCells[3][1][3] = Eigen::Vector3d(0.604788, -0.389835, -0.0395852);
    vCells[0][2][3] = Eigen::Vector3d(0.591889, -0.274596, -0.00640475);
    vCells[1][2][3] = Eigen::Vector3d(0.608694, -0.4124, 0.00948705);
    vCells[2][2][3] = Eigen::Vector3d(0.599275, -0.496561, 0.174348);
    vCells[3][2][3] = Eigen::Vector3d(0.566664, -0.546236, 0.322577);
    vCells[0][3][3] = Eigen::Vector3d(0.596525, -0.500423, 0.270982);
    vCells[1][3][3] = Eigen::Vector3d(0.57696, -0.559376, 0.394342);
    vCells[2][3][3] = Eigen::Vector3d(0.562996, -0.612187, 0.471696);
    vCells[3][3][3] = Eigen::Vector3d(0.555753, -0.650018, 0.518281);

    CatmullRomInterpolant::Cells<Eigen::Vector3d> wCells;

    wCells[0][0][0] = Eigen::Vector3d(0.0134557, 0.0364571, 0.574216);
    wCells[1][0][0] = Eigen::Vector3d(0.0406675, 0.0130815, 0.488785);
    wCells[2][0][0] = Eigen::Vector3d(0.0836333, 0.033461, 0.203229);
    wCells[3][0][0] = Eigen::Vector3d(0.0512807, -0.0192547, 0.0804576);
    wCells[0][1][0] = Eigen::Vector3d(0.0813105, 0.00729747, 0.199809);
    wCells[1][1][0] = Eigen::Vector3d(0.0520183, -0.0378754, 0.096699);
    wCells[2][1][0] = Eigen::Vector3d(0.26171, -0.221512, 0.0804091);
    wCells[3][1][0] = Eigen::Vector3d(0.312962, -0.2873, 0.0492301);
    wCells[0][2][0] = Eigen::Vector3d(0.49172, -0.435381, 0.0843253);
    wCells[1][2][0] = Eigen::Vector3d(0.207959, -0.205871, 0.0371094);
    wCells[2][2][0] = Eigen::Vector3d(7.35215e-05, -0.0111327, 0.0142914);
    wCells[3][2][0] = Eigen::Vector3d(0.000135638, 0.000111358, -8.48327e-06);
    wCells[0][3][0] = Eigen::Vector3d(0.000185823, 0.000157328, -5.8268e-05);
    wCells[1][3][0] = Eigen::Vector3d(0.000179823, 9.17055e-05, 4.21921e-06);
    wCells[2][3][0] = Eigen::Vector3d(0.000130017, 3.86277e-05, 3.29982e-05);
    wCells[3][3][0] = Eigen::Vector3d(7.98224e-05, -2.00181e-05, 3.66595e-05);
    wCells[0][0][1] = Eigen::Vector3d(0.155001, 0.514089, 0.843615);
    wCells[1][0][1] = Eigen::Vector3d(0.24003, 0.562655, 0.496099);
    wCells[2][0][1] = Eigen::Vector3d(0.0780004, 0.1486, -0.0966623);
    wCells[3][0][1] = Eigen::Vector3d(-0.187626, -0.164309, 0.25562);
    wCells[0][1][1] = Eigen::Vector3d(0.0567082, 0.109739, -0.097436);
    wCells[1][1][1] = Eigen::Vector3d(-0.18867, -0.155126, 0.249057);
    wCells[2][1][1] = Eigen::Vector3d(-0.0871918, 0.081353, 0.443784);
    wCells[3][1][1] = Eigen::Vector3d(0.00247988, 0.241447, 0.304571);
    wCells[0][2][1] = Eigen::Vector3d(-0.0845409, 0.125036, 0.437893);
    wCells[1][2][1] = Eigen::Vector3d(-0.000516096, 0.259487, 0.347683);
    wCells[2][2][1] = Eigen::Vector3d(-0.0137104, 0.168395, 0.365704);
    wCells[3][2][1] = Eigen::Vector3d(-0.00162429, -0.137887, 0.254013);
    wCells[0][3][1] = Eigen::Vector3d(0.00649781, 0.0374063, 0.331684);
    wCells[1][3][1] = Eigen::Vector3d(-0.0116833, -0.244238, 0.165298);
    wCells[2][3][1] = Eigen::Vector3d(0.0227272, -0.297158, 0.0896801);
    wCells[3][3][1] = Eigen::Vector3d(0.12524, -0.167133, 0.0716439);
    wCells[0][0][2] = Eigen::Vector3d(0.337757, 0.685898, -0.0604283);
    wCells[1][0][2] = Eigen::Vector3d(0.28338, 0.585777, -0.394459);
    wCells[2][0][2] = Eigen::Vector3d(0.0484554, 0.126152, -0.487364);
    wCells[3][0][2] = Eigen::Vector3d(0.00396854, -0.140528, -0.302435);
    wCells[0][1][2] = Eigen::Vector3d(0.0426426, 0.0815117, -0.545703);
    wCells[1][1][2] = Eigen::Vector3d(0.00282667, -0.173836, -0.366724);
    wCells[2][1][2] = Eigen::Vector3d(0.0743859, 0.0350029, -0.0535327);
    wCells[3][1][2] = Eigen::Vector3d(0.0262071, 0.267642, 0.0255088);
    wCells[0][2][2] = Eigen::Vector3d(0.0625161, 0.0697091, -0.0937734);
    wCells[1][2][2] = Eigen::Vector3d(0.0148102, 0.327423, -0.00704943);
    wCells[2][2][2] = Eigen::Vector3d(-0.0364622, 0.212775, -0.187782);
    wCells[3][2][2] = Eigen::Vector3d(-0.0272192, 0.00410804, -0.229815);
    wCells[0][3][2] = Eigen::Vector3d(-0.0174583, 0.126733, -0.31148);
    wCells[1][3][2] = Eigen::Vector3d(-0.028702, -0.004114, -0.285703);
    wCells[2][3][2] = Eigen::Vector3d(-0.00655798, -0.0219703, -0.144016);
    wCells[3][3][2] = Eigen::Vector3d(0.0267034, -0.0159486, 0.162332);
    wCells[0][0][3] = Eigen::Vector3d(0.210906, 0.394395, -0.240989);
    wCells[1][0][3] = Eigen::Vector3d(0.124439, 0.293293, -0.410238);
    wCells[2][0][3] = Eigen::Vector3d(0.0256846, 0.0255014, -0.426495);
    wCells[3][0][3] = Eigen::Vector3d(0.0235816, -0.154364, -0.363658);
    wCells[0][1][3] = Eigen::Vector3d(0.0325416, -0.0534158, -0.490571);
    wCells[1][1][3] = Eigen::Vector3d(-0.0364524, -0.227525, -0.391313);
    wCells[2][1][3] = Eigen::Vector3d(0.0377604, -0.112829, -0.188707);
    wCells[3][1][3] = Eigen::Vector3d(0.00749213, 0.131531, 0.0154943);
    wCells[0][2][3] = Eigen::Vector3d(0.0836453, -0.10032, -0.0830725);
    wCells[1][2][3] = Eigen::Vector3d(-0.0128917, 0.140953, -0.00762958);
    wCells[2][2][3] = Eigen::Vector3d(-0.0312827, 0.316187, -0.0578429);
    wCells[3][2][3] = Eigen::Vector3d(0.0888737, 0.332096, -0.207193);
    wCells[0][3][3] = Eigen::Vector3d(0.0329841, 0.391755, -0.034668);
    wCells[1][3][3] = Eigen::Vector3d(0.156468, 0.43512, -0.222207);
    wCells[2][3][3] = Eigen::Vector3d(0.23602, 0.492642, -0.123861);
    wCells[3][3][3] = Eigen::Vector3d(0.225855, 0.565122, 0.0768213);

    // Initialize the Bezier volumes:
    TricubicBezierVolume3d v;
    TricubicBezierVolume3d::ControlPointsType& vControlPoints = v.getControlPoints();
    TricubicBezierVolume3d w;
    TricubicBezierVolume3d::ControlPointsType& wControlPoints = w.getControlPoints();

    // Convert grid points using the Catmull Rom interpolant:
    // Eq. 10-17
    CatmullRomInterpolant::convert(vCells, vControlPoints);
    CatmullRomInterpolant::convert(wCells, wControlPoints);
    // this gives us Eq. 6 & 7

    // Test setup:
    std::vector<Eigen::Vector2d> uvs;
    const double epsilon = 1e-4;

    // Cross product and root finding setup:
    BezierCrossQuad<BicubicBezierSurface3d, BicubicBezierSurface3d> crossQuad;
    BezierSurfaceRootFinding<BezierSurface<6, 6, Eigen::Vector3d>> bezierRootFinding(1e-10, 1e-10, 1e-3, 300, 1e-3, 1024);

    // Bottom Face
    BicubicBezierSurface3d vSurface;
    v.getBoundarySurfaceVMin(vSurface);
    BicubicBezierSurface3d wSurface;
    w.getBoundarySurfaceVMin(wSurface);
    // Eq. 8-9
    BezierSurface<6, 6, Eigen::Vector3d> crossProduct = crossQuad.compute(vSurface, wSurface);
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }
    // expect one root on this face
    Eigen::Vector2d expect = Eigen::Vector2d(0.737744, 0.225455);
    EXPECT_EQ(uvs.size(), 1);
    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs.at(0);
        EXPECT_NEAR(expect.x(), uvHave.x(), epsilon);
        EXPECT_NEAR(expect.y(), uvHave.y(), epsilon);
    }
    // cleanup
    uvs.clear();

    // Left Face
    v.getBoundarySurfaceUMin(vSurface);
    w.getBoundarySurfaceUMin(wSurface);
    // Eq. 8-9
    crossProduct = crossQuad.compute(vSurface, wSurface);
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }

    // expect one root on this face
    expect = Eigen::Vector2d(0.440337, 0.460312);
    EXPECT_EQ(uvs.size(), 1);
    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs.at(0);
        EXPECT_NEAR(expect.x(), uvHave.x(), epsilon);
        EXPECT_NEAR(expect.y(), uvHave.y(), epsilon);
    }
    // cleanup
    uvs.clear();

    // Front Face
    v.getBoundarySurfaceWMin(vSurface);
    w.getBoundarySurfaceWMin(wSurface);
    // Eq. 8-9
    crossProduct = crossQuad.compute(vSurface, wSurface);
    crossProduct.recomputeBoundingBox();
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }

    // expect no root on this face
    EXPECT_EQ(uvs.size(), 0);
}

inline void test_voxel_tricubic_derived_acceleration()
{
    // Here, it is not possible to use the grid points directly as the control points because
    // a tricubic interpolation is used. The grid points are stored in vCells and wCells.
    // Additionally, the acceleration is derived from the Bezier volume for v directly.

    CatmullRomInterpolant::Cells<Eigen::Vector3d> vCells;

    vCells[0][0][0] = Eigen::Vector3d(0.384905, -0.331021, 0.393825);
    vCells[1][0][0] = Eigen::Vector3d(0.348277, -0.413362, 0.365962);
    vCells[2][0][0] = Eigen::Vector3d(0.326964, -0.4534, 0.255795);
    vCells[3][0][0] = Eigen::Vector3d(0.341653, -0.42509, 0.155603);
    vCells[0][1][0] = Eigen::Vector3d(0.329087, -0.450383, 0.230793);
    vCells[1][1][0] = Eigen::Vector3d(0.340945, -0.420839, 0.149053);
    vCells[2][1][0] = Eigen::Vector3d(0.360521, -0.382515, 0.0886838);
    vCells[3][1][0] = Eigen::Vector3d(0.382397, -0.350426, 0.0549903);
    vCells[0][2][0] = Eigen::Vector3d(0.360575, -0.377413, 0.0833482);
    vCells[1][2][0] = Eigen::Vector3d(0.383187, -0.342312, 0.0525134);
    vCells[2][2][0] = Eigen::Vector3d(0.194541, -0.194253, 0.0180661);
    vCells[3][2][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000172922);
    vCells[0][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000228108);
    vCells[1][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000209707);
    vCells[2][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000191317);
    vCells[3][3][0] = Eigen::Vector3d(2.75934e-05, 0, 0.000172921);
    vCells[0][0][1] = Eigen::Vector3d(0.398916, -0.288387, 0.648467);
    vCells[1][0][1] = Eigen::Vector3d(0.407706, -0.336648, 0.525273);
    vCells[2][0][1] = Eigen::Vector3d(0.464419, -0.342912, 0.205989);
    vCells[3][0][1] = Eigen::Vector3d(0.522456, -0.35406, -0.131747);
    vCells[0][1][1] = Eigen::Vector3d(0.476781, -0.359284, 0.164483);
    vCells[1][1][1] = Eigen::Vector3d(0.533955, -0.347191, -0.163723);
    vCells[2][1][1] = Eigen::Vector3d(0.604242, -0.2743, -0.199436);
    vCells[3][1][1] = Eigen::Vector3d(0.67171, -0.167851, -0.159466);
    vCells[0][2][1] = Eigen::Vector3d(0.627536, -0.248802, -0.207454);
    vCells[1][2][1] = Eigen::Vector3d(0.680804, -0.15109, -0.146521);
    vCells[2][2][1] = Eigen::Vector3d(0.711689, -0.0849616, -0.0435621);
    vCells[3][2][1] = Eigen::Vector3d(0.703087, -0.0777879, 0.0506719);
    vCells[0][3][1] = Eigen::Vector3d(0.72153, -0.0838632, 0.00731354);
    vCells[1][3][1] = Eigen::Vector3d(0.697467, -0.106539, 0.0691244);
    vCells[2][3][1] = Eigen::Vector3d(0.646464, -0.168068, 0.101593);
    vCells[3][3][1] = Eigen::Vector3d(0.5826, -0.255839, 0.108278);
    vCells[0][0][2] = Eigen::Vector3d(0.484601, -0.0884814, 0.685874);
    vCells[1][0][2] = Eigen::Vector3d(0.506102, -0.10573, 0.552217);
    vCells[2][0][2] = Eigen::Vector3d(0.518997, -0.155266, 0.220147);
    vCells[3][0][2] = Eigen::Vector3d(0.531333, -0.233848, -0.102992);
    vCells[0][1][2] = Eigen::Vector3d(0.532968, -0.148716, 0.188138);
    vCells[1][1][2] = Eigen::Vector3d(0.551587, -0.230482, -0.114794);
    vCells[2][1][2] = Eigen::Vector3d(0.582001, -0.318849, -0.267985);
    vCells[3][1][2] = Eigen::Vector3d(0.60904, -0.358736, -0.205413);
    vCells[0][2][2] = Eigen::Vector3d(0.597367, -0.332034, -0.264587);
    vCells[1][2][2] = Eigen::Vector3d(0.610086, -0.388842, -0.178064);
    vCells[2][2][2] = Eigen::Vector3d(0.579485, -0.430263, -0.00620217);
    vCells[3][2][2] = Eigen::Vector3d(0.542953, -0.474194, 0.163135);
    vCells[0][3][2] = Eigen::Vector3d(0.572444, -0.468671, 0.0727383);
    vCells[1][3][2] = Eigen::Vector3d(0.544172, -0.536679, 0.229456);
    vCells[2][3][2] = Eigen::Vector3d(0.519358, -0.588139, 0.334183);
    vCells[3][3][2] = Eigen::Vector3d(0.510477, -0.631277, 0.377793);
    vCells[0][0][3] = Eigen::Vector3d(0.551293, 0.0301241, 0.642897);
    vCells[1][0][3] = Eigen::Vector3d(0.551752, 0.0384976, 0.584381);
    vCells[2][0][3] = Eigen::Vector3d(0.538636, -0.00932551, 0.385826);
    vCells[3][0][3] = Eigen::Vector3d(0.535415, -0.11538, 0.113912);
    vCells[0][1][3] = Eigen::Vector3d(0.583807, 0.00105604, 0.377048);
    vCells[1][1][3] = Eigen::Vector3d(0.55558, -0.120834, 0.109389);
    vCells[2][1][3] = Eigen::Vector3d(0.56956, -0.270999, -0.0659052);
    vCells[3][1][3] = Eigen::Vector3d(0.604788, -0.389835, -0.0395852);
    vCells[0][2][3] = Eigen::Vector3d(0.591889, -0.274596, -0.00640475);
    vCells[1][2][3] = Eigen::Vector3d(0.608694, -0.4124, 0.00948705);
    vCells[2][2][3] = Eigen::Vector3d(0.599275, -0.496561, 0.174348);
    vCells[3][2][3] = Eigen::Vector3d(0.566664, -0.546236, 0.322577);
    vCells[0][3][3] = Eigen::Vector3d(0.596525, -0.500423, 0.270982);
    vCells[1][3][3] = Eigen::Vector3d(0.57696, -0.559376, 0.394342);
    vCells[2][3][3] = Eigen::Vector3d(0.562996, -0.612187, 0.471696);
    vCells[3][3][3] = Eigen::Vector3d(0.555753, -0.650018, 0.518281);

    // Initialize the Bezier volume...
    TricubicBezierVolume3d v;
    TricubicBezierVolume3d::ControlPointsType& vControlPoints = v.getControlPoints();

    // Convert grid points using the Catmull Rom interpolant:
    // Eq. 10-17
    CatmullRomInterpolant::convert(vCells, vControlPoints);
    // this gives us Eq. 18

    // Compute the acceleration
    BezierAcceleration<TricubicBezierVolume3d> bezierAcceleration;
    // Eq. 19
    TricubicBezierVolume3d du           = v.generatePartialTensorProduct(PartialUVW(1, 0, 0)); // computed with  Eq. 22-25
    // Eq. 20
    TricubicBezierVolume3d dv           = v.generatePartialTensorProduct(PartialUVW(0, 1, 0)); // computed with  Eq. 22-25
    // Eq. 21
    TricubicBezierVolume3d dw           = v.generatePartialTensorProduct(PartialUVW(0, 0, 1)); // computed with  Eq. 22-25
    // Eq. 26
    BezierVolume<6, 6, 6, Eigen::Vector3d> w = bezierAcceleration.compute(v, du, dv, dw); // computed with Eq. 27

    // Test setup:
    std::vector<Eigen::Vector2d> uvs;
    const double epsilon = 1e-4;

    BezierCrossQuad<BicubicBezierSurface3d, BezierVolume<6, 6, 6, Eigen::Vector3d>> crossQuad;
    BezierSurfaceRootFinding<BezierSurface<9, 9, Eigen::Vector3d>> bezierRootFinding(1e-10, 1e-10, 1e-3, 300, 1e-3, 1024);

    // Bottom Face
    BicubicBezierSurface3d vSurface;
    v.getBoundarySurfaceVMin(vSurface);
    BezierSurface<6, 6, Eigen::Vector3d> wSurface;
    w.getBoundarySurfaceVMin(wSurface);
    // Eq. 8-9
    BezierSurface<9, 9, Eigen::Vector3d> crossProduct = crossQuad.compute(vSurface, wSurface);
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }
    // expect one root on this face
    Eigen::Vector2d expect = Eigen::Vector2d(0.630879, 0.234697);
    EXPECT_EQ(uvs.size(), 1);
    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs.at(0);
        EXPECT_NEAR(expect.x(), uvHave.x(), epsilon);
        EXPECT_NEAR(expect.y(), uvHave.y(), epsilon);
    }
    // cleanup
    uvs.clear();

    // Left Face
    v.getBoundarySurfaceUMin(vSurface);
    w.getBoundarySurfaceUMin(wSurface);
    // Eq. 8-9
    crossProduct = crossQuad.compute(vSurface, wSurface);
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }

    // expect one root on this face
    expect = Eigen::Vector2d(0.955859, 0.664451);
    EXPECT_EQ(uvs.size(), 1);
    if (uvs.size() == 1)
    {
        Eigen::Vector2d uvHave = uvs.at(0);
        EXPECT_NEAR(expect.x(), uvHave.x(), epsilon);
        EXPECT_NEAR(expect.y(), uvHave.y(), epsilon);
    }
    // cleanup
    uvs.clear();

    // Front Face
    v.getBoundarySurfaceWMin(vSurface);
    w.getBoundarySurfaceWMin(wSurface);
    // Eq. 8-9
    crossProduct = crossQuad.compute(vSurface, wSurface);
    crossProduct.recomputeBoundingBox();
    // quickly check if the bounding box even contains zero
    crossProduct.recomputeBoundingBox();
    if (crossProduct.getBoundingBox().contains(Eigen::Vector3d::Zero()))
    {
        // Fig. 5 (Sec. 3.4 "Projection-based Clipping")
        bezierRootFinding.clippingSolver<ProjectionBased>(crossProduct, uvs);
    }

    // expect no root on this face
    EXPECT_EQ(uvs.size(), 0);
}