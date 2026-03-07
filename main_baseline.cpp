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
    int nx = 500;
    int ny = 197;
    const int ng = 2;
    const double Lx = 0.225;
    const double Ly = 0.089;
    const double cfl = 0.4;
    double t_end = 0.0011741;

    if (argc >= 3) {
        nx = std::atoi(argv[1]);
        ny = std::atoi(argv[2]);
    }
    if (argc >= 4) {
        t_end = std::atof(argv[3]);
    }

    Grid grid;
    grid.init(nx, ny, ng, Lx, Ly);
    initialize_shock_bubble(grid);

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
              << " steps=" << step
              << " wall_seconds=" << wall
              << std::endl;

    return 0;
}