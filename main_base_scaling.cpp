#include "types.hpp"
#include "physics.hpp"
#include "solver.hpp"
#include "init.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    // fixed parameters
    const int nx_base = 500;
    const int ny_base = 197;
    const int ng = 2;
    const double Lx = 0.225;
    const double Ly_base = 0.089;
    const double cfl = 0.4;
    const double t_end = 0.0011741;

    // variable parameters for scaling
    int nx = nx_base;
    int copies_y = 1;

    // Usage:
    //   ./serial_scaling.exe [nx] [copies_y]
    //   ./omp_scaling.exe    [nx] [copies_y]

    if (argc >= 2) {
        nx = std::atoi(argv[1]);
    }
    if (argc >= 3) {
        copies_y = std::atoi(argv[2]);
    }

    if (nx <= 0 || copies_y <= 0) {
        std::cerr << "Error: nx and copies_y must be positive.\n";
        return 1;
    }

    const int ny = ny_base * copies_y;
    const double Ly = Ly_base * copies_y;

    Grid grid;
    grid.init(nx, ny, ng, Lx, Ly);
    initialize_shock_bubble_weak(grid, copies_y);

    int step = 0;
    double t = 0.0;

    auto t0 = std::chrono::steady_clock::now();

    double dt = compute_dt(grid, cfl);

    while (t < t_end) {
        if (t + dt > t_end) {
            dt = t_end - t;
        }

        advance_one_step(grid, dt);

        t += dt;
        ++step;

        dt = compute_dt(grid, cfl);
    }

    auto t1 = std::chrono::steady_clock::now();
    const double wall =
        std::chrono::duration<double>(t1 - t0).count();

#ifdef _OPENMP
    const int p = omp_get_max_threads();
    std::cout << "mode=omp";
#else
    const int p = 1;
    std::cout << "mode=serial";
#endif

    std::cout << " p=" << p
              << " nx=" << nx
              << " ny=" << ny
              << " copies_y=" << copies_y
              << " steps=" << step
              << " wall_seconds=" << wall
              << std::endl;

    return 0;
}