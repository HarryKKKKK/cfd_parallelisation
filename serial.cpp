#include "types.hpp"
#include "physics.hpp"
#include "solver.hpp"
#include "init.hpp"
#include "utils.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

int main() {
    const std::string out_dir = "res/serial";
    std::filesystem::create_directories(out_dir);

    // ----------------- Init grid -----------------
    Grid grid;
    grid.init(500, 197, 2, 0.225, 0.089);

    initialize_shock_bubble(grid);

    int step = 0;
    double t = 0.0;

    // Write t=0
    write_grid_csv(grid, make_filename(out_dir, step, t));
    std::cout << "Wrote: " << make_filename(out_dir, step, t) << "\n";

    // ----------------- time loop settings -----------------
    const double t_end = 0.0011741;
    const double cfl   = 0.5;

    // output frequency related
    const int n_outputs = 10;
    const double output_dt = t_end / (n_outputs-1);
    double next_output_t = output_dt;

    // ----------------- Main loop -----------------
    while (t < t_end) {
        double dt = compute_dt(grid, cfl);
        if (t + dt > t_end) dt = t_end - t;

        advance_one_step(grid, dt);
        t += dt;
        step++;

        // Output 
        if (t + 1e-15 >= next_output_t || t >= t_end) {
            const std::string fname = make_filename(out_dir, step, t);
            write_grid_csv(grid, fname);
            std::cout << "Wrote: " << fname << "\n";
            next_output_t += output_dt;
        }
    }

    return 0;
}
