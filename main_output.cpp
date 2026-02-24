#include "types.hpp"
#include "physics.hpp"
#include "solver.hpp"
#include "init.hpp"
#include "utils.hpp"

#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <algorithm>

constexpr bool ENABLE_OUTPUT = true;

int main() {
    const std::string out_dir = "res/serial/hllc";
    if (ENABLE_OUTPUT) {
        std::filesystem::create_directories(out_dir);
    }

    // ----------------- Init grid -----------------
    Grid grid;
    grid.init(500, 197, 2, 0.225, 0.089);

    initialize_shock_bubble(grid);

    int step = 0;
    double t = 0.0;

    // Initial checks
    check_symmetry(grid, t, true);

    if (ENABLE_OUTPUT) {
        write_grid_csv(grid, make_filename(out_dir, step, t));
    }

    // ----------------- time settings -----------------
    const double t_end = 0.0011741;  // physical time corresponding to t*=19
    const double cfl   = 0.4;

    const double tstar_targets[] = {0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 9.4, 12.6, 17.4, 19.0};
    const int n_targets = sizeof(tstar_targets)/sizeof(tstar_targets[0]);

    const double t_ref = t_end / 19.0;
    double targets_phys[n_targets];
    for (int k = 0; k < n_targets; ++k) targets_phys[k] = tstar_targets[k] * t_ref;

    int next_k = 0;

    // ----------------- Main loop -----------------
    double dt = compute_dt(grid, cfl);

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        if (ENABLE_OUTPUT && next_k < n_targets) {
            const double t_out = targets_phys[next_k];
            if (t < t_out && t + dt > t_out) {
                dt = t_out - t;
            }
        }

        advance_one_step(grid, dt);

        t += dt;
        step++;

        if (ENABLE_OUTPUT && next_k < n_targets) {
            const double t_out = targets_phys[next_k];
            if (t + 1e-15 >= t_out) {
                check_symmetry(grid, t, true);
                write_grid_csv(grid, make_filename(out_dir, step, t));
                next_k++;
            }
        }

        dt = compute_dt(grid, cfl);
    }

    std::cout << "Finished. Steps = " << step << std::endl;
    return 0;
}