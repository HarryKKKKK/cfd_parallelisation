#include "types.hpp"
#include "solver.hpp"
#include "init.hpp"

#include <chrono>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
    // fixed problem setup
    const int nx = 500;
    const int ny = 197;
    const int ng = 2;
    const double Lx = 0.225;
    const double Ly = 0.089;
    const double cfl = 0.4;
    const double t_end = 0.0011741;

    Grid grid;
    grid.init(nx, ny, ng, Lx, Ly);
    initialize_shock_bubble(grid);

    std::cout << std::setprecision(12);
    std::cout << "[INIT] nx=" << nx
              << " ny=" << ny
              << " ng=" << ng
              << " Lx=" << Lx
              << " Ly=" << Ly
              << " dx=" << grid.dx
              << " dy=" << grid.dy
              << " cfl=" << cfl
              << " t_end=" << t_end
              << std::endl;

    int step = 0;
    double t = 0.0;

    auto t0 = std::chrono::steady_clock::now();

    double dt = compute_dt(grid, cfl);

    std::cout << "[INIT] initial dt=" << dt << std::endl;

    while (t < t_end) {
        if (t + dt > t_end) {
            dt = t_end - t;
        }

        if (step % 500 == 0) {
            auto now = std::chrono::steady_clock::now();
            double wall_so_far =
                std::chrono::duration<double>(now - t0).count();

            std::cout << "[LOOP] step=" << step
                      << " t=" << t
                      << " dt=" << dt
                      << " wall_so_far=" << wall_so_far
                      << std::endl;
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