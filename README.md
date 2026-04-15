# Reference Implementation for Parallel Vectors Extraction using Bézier Clipping

<img width="3806" height="1243" alt="teaser" src="https://github.com/user-attachments/assets/c9412218-bf1a-4734-b6e9-9c12f6a22a5e" />

## About

Reference Implementation for the EuroVis 2026 paper "Parallel Vectors "Parallel Vectors Extraction using Bézier Clipping". 
Important files:

- `examples.hpp` Contains functions with examples of how to solve for roots on Bézier surfaces. The _test_solvers_ and _test_solver_performance_ functions compare all available solvers and their performance on a toy example. Please refer to the paper for complete benchmark results. Further, the file contains _test_voxel_trilinear|tricubic|tricubic_derived_acceleration_ which show how a trilinear, tricubic and tricubic with derived acceleration Bézier volume/surface is built and used for root finding.
- `bezier_surface_root_finding.hpp` Contains all available solvers (Clipping (Component-wise, Projection-based), Bisection, Newton).
- `bezier_cross_quad` Contains the cross product computation from two Bézier surfaces.
- `catmull_rom_interpolant` Contains the cubic interpolant int Bézier form.
- `bezier_acceleration` Contains the derived acceleration computation.
- `bezier_volume` Contains the Bézier volume representation and the partial tensor product computation used for the derived acceleration computation.

Below is a short guide on how to get, build and run the code for Linux-based systems.


### Get the Code

```
git clone https://github.com/fau-vc/pvBezier.git
```

### Build

```
mkdir build
cd build
cmake ..
cmake --build .
cd ..
```

### Run

```
./build/pvBezierMain
```
