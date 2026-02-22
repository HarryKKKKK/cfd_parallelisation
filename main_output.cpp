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

// -----------------------------
// Diagnostic helpers
// -----------------------------
static void debug_dump_rho(const Grid& grid, double t, const char* tag) {
    const int nx = grid.nx;
    const int ny = grid.ny;
    const int ng = grid.ng;

    auto rho_at = [&](int i_phys, int j_phys) -> double {
        // i_phys: 0..nx-1, j_phys: 0..ny-1
        const int I = i_phys + ng;
        const int J = j_phys + ng;
        const Conserved& U = grid.U[grid.idx(I, J)];
        const Primitive W = conserved_to_primitive(U);
        return W.rho;
    };

    // 1) min/max over physical domain
    double rho_min = 1e300, rho_max = -1e300;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const double r = rho_at(i, j);
            rho_min = std::min(rho_min, r);
            rho_max = std::max(rho_max, r);
        }
    }

    // 2) boundary-vs-interior differences (L/R/T/B)
    // compare first physical line vs second physical line
    double max_L = 0.0, max_R = 0.0, max_B = 0.0, max_T = 0.0;
    for (int j = 0; j < ny; ++j) {
        max_L = std::max(max_L, std::abs(rho_at(0, j) - rho_at(1, j)));
        max_R = std::max(max_R, std::abs(rho_at(nx-1, j) - rho_at(nx-2, j)));
    }
    for (int i = 0; i < nx; ++i) {
        max_B = std::max(max_B, std::abs(rho_at(i, 0) - rho_at(i, 1)));
        max_T = std::max(max_T, std::abs(rho_at(i, ny-1) - rho_at(i, ny-2)));
    }

    // 3) plateau / repeated-values detection:
    // bin rho by rounding to 12 decimals and count the top-5 most frequent.
    std::unordered_map<long long, long long> cnt;
    cnt.reserve(static_cast<size_t>(nx) * static_cast<size_t>(ny) / 4);

    auto key12 = [](double x) -> long long {
        // map to integer key of rounded 1e-12
        // (safe for typical rho magnitudes here)
        return static_cast<long long>(std::llround(x * 1e12));
    };

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            cnt[key12(rho_at(i, j))] += 1;
        }
    }

    // extract top-5
    std::vector<std::pair<long long, long long>> top(cnt.begin(), cnt.end());
    std::partial_sort(
        top.begin(),
        top.begin() + std::min<size_t>(5, top.size()),
        top.end(),
        [](auto& a, auto& b){ return a.second > b.second; }
    );
    const int K = static_cast<int>(std::min<size_t>(5, top.size()));

    std::cout << "\n[DEBUG] " << tag << " at t=" << t << "\n";
    std::cout << "  rho_min=" << rho_min << " rho_max=" << rho_max << "\n";
    std::cout << "  max|rho(x0)-rho(x1)| over y = " << max_L << "\n";
    std::cout << "  max|rho(xR)-rho(xR-1)| over y = " << max_R << "\n";
    std::cout << "  max|rho(y0)-rho(y1)| over x = " << max_B << "\n";
    std::cout << "  max|rho(yT)-rho(yT-1)| over x = " << max_T << "\n";

    const long long total = static_cast<long long>(nx) * static_cast<long long>(ny);
    std::cout << "  Top repeated rho (rounded to 1e-12):\n";
    for (int k = 0; k < K; ++k) {
        const double rho_val = top[k].first / 1e12;
        const double frac = static_cast<double>(top[k].second) / static_cast<double>(total);
        std::cout << "    rho≈" << rho_val << "  count=" << top[k].second
                  << "  frac=" << frac << "\n";
    }

    // 4) left boundary column summary (helps if the black band is near x=0)
    // compute min/max along i=0 column (physical)
    double col_min = 1e300, col_max = -1e300;
    for (int j = 0; j < ny; ++j) {
        const double r = rho_at(0, j);
        col_min = std::min(col_min, r);
        col_max = std::max(col_max, r);
    }
    std::cout << "  Left boundary column (i=0) rho_min=" << col_min
              << " rho_max=" << col_max << "\n";
}

int main() {
    const std::string out_dir = "res/serial";
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
    debug_dump_rho(grid, t, "INIT");

    if (ENABLE_OUTPUT) {
        write_grid_csv(grid, make_filename(out_dir, step, t));
    }

    // ----------------- time settings -----------------
    const double t_end = 0.0011741;  // physical time corresponding to t*=19
    const double cfl   = 0.4;

    const double tstar_targets[] = {0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 12.6, 19.0};
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
                debug_dump_rho(grid, t, "OUTPUT");
                write_grid_csv(grid, make_filename(out_dir, step, t));
                next_k++;
            }
        }

        dt = compute_dt(grid, cfl);
    }

    std::cout << "Finished. Steps = " << step << std::endl;
    return 0;
}