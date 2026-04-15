#include <iostream>

#include "examples.hpp"

// TODO: insert as many equation references as possible

int main()
{
    std::cout << "Test if different solvers agree on sample biquadratic Bezier surface..." << std::endl;

    test_solvers();

    std::cout << "Done." << std::endl;

    std::cout << "Test solver runtime by repeating previous experiment 10.000x..." << std::endl;

    test_solver_performance();

    std::cout << "Done." << std::endl;

    std::cout << "Test different interpolators..." << std::endl;

    test_voxel_trilinear();
    test_voxel_tricubic();
    test_voxel_tricubic_derived_acceleration();

    std::cout << "Done." << std::endl;

}
