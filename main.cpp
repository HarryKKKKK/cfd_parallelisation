#include "types.hpp"
#include "physics.hpp"
#include "solver.hpp"
#include "init.hpp"
#include "utils.hpp"

#include <filesystem>
#include <iostream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// =============================
//   OUTPUT CONTROL SWITCH
// =============================

constexpr bool ENABLE_OUTPUT = false;

int main() {

#ifdef _OPENMP
    const std::string out_dir = "res/omp";
    #pragma omp parallel
    {
        #pragma omp single
        {
        std::cout << "OpenMP max threads = " << omp_get_max_threads() << "\n";
        std::cout << "OpenMP num threads = " << omp_get_num_threads() << "\n";
        }
    }
#else
    const std::string out_dir = "res/serial";
#endif

    if (ENABLE_OUTPUT) {
        std::filesystem::create_directories(out_dir);
    }

    // ----------------- Init grid -----------------
    Grid grid;
    grid.init(500, 197, 2, 0.225, 0.089);

    initialize_shock_bubble(grid);

    int step = 0;
    double t = 0.0;

    if (ENABLE_OUTPUT) {
        write_grid_csv(grid, make_filename(out_dir, step, t));
    }

    // ----------------- time loop settings -----------------
    const double t_end = 0.0011741;
    const double cfl   = 0.4;

    const int n_outputs = 10;
    const double output_dt = t_end / (n_outputs - 1);
    double next_output_t = output_dt;

    // ----------------- Main loop -----------------
    double dt = compute_dt(grid, cfl);
    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        advance_one_step(grid, dt);

        t += dt;
        step++;

        if (ENABLE_OUTPUT) {
            if (t + 1e-15 >= next_output_t || t >= t_end) {
                write_grid_csv(grid, make_filename(out_dir, step, t));
                next_output_t += output_dt;
            }
        }
    }

    std::cout << "Finished. Steps = " << step << std::endl;

    return 0;
}
