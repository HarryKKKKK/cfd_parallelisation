#include "solver.hpp"
#include "physics.hpp"
#include "constants.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

// ----------------- small helpers on Conserved -----------------
static inline Conserved add(const Conserved& a, const Conserved& b) {
    return {a.rho + b.rho, a.rhou + b.rhou, a.rhov + b.rhov, a.E + b.E};
}
static inline Conserved sub(const Conserved& a, const Conserved& b) {
    return {a.rho - b.rho, a.rhou - b.rhou, a.rhov - b.rhov, a.E - b.E};
}
static inline Conserved mul(const Conserved& a, double s) {
    return {a.rho * s, a.rhou * s, a.rhov * s, a.E * s};
}

// ----------------- limiter -----------------
double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

static inline Conserved minmod_vec(const Conserved& a, const Conserved& b) {
    return {
        minmod(a.rho, b.rho),
        minmod(a.rhou, b.rhou),
        minmod(a.rhov, b.rhov),
        minmod(a.E, b.E)
    };
}

// ----------------- boundary conditions -----------------
void apply_boundary_conditions(Grid& grid) {
    DBG_PRINT("apply boundary conditions");

    const int nx_tot = grid.nx + 2 * grid.ng;
    const int ny_tot = grid.ny + 2 * grid.ng;
    const int ng = grid.ng;

    for (int j = 0; j < ny_tot; ++j) {
        for (int g = 0; g < ng; ++g) {
            grid.U[grid.idx(g, j)] = grid.U[grid.idx(ng, j)];
            grid.U[grid.idx(nx_tot-1-g, j)] = grid.U[grid.idx(nx_tot-1-ng, j)];
        }
    }

    for (int i = 0; i < nx_tot; ++i) {
        for (int g = 0; g < ng; ++g) {
            grid.U[grid.idx(i, g)] = grid.U[grid.idx(i, ng)];
            grid.U[grid.idx(i, ny_tot-1-g)] = grid.U[grid.idx(i, ny_tot-1-ng)];
        }
    }
}

// ----------------- CFL timestep -----------------
double compute_dt(const Grid& grid, double cfl) {
    DBG_PRINT("computing dt");

    const int ng = grid.ng;
    const int i0 = ng, i1 = ng + grid.nx - 1;
    const int j0 = ng, j1 = ng + grid.ny - 1;

    double max_lambda = 0.0;

    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            const Primitive W = conserved_to_primitive(grid.U[grid.idx(i, j)]);
            const double a = std::sqrt(phys::gamma * W.p / W.rho);
            max_lambda = std::max(
                max_lambda,
                (std::abs(W.u) + a) / grid.dx +
                (std::abs(W.v) + a) / grid.dy
            );
        }
    }

    if (max_lambda <= 0.0) return std::numeric_limits<double>::infinity();
    return cfl / max_lambda;
}

// ----------------- MUSCL slopes -----------------
static inline Conserved slope_x_safe(const Grid& grid, int i, int j) {
    const int nx_tot = grid.nx + 2 * grid.ng;
    if (i <= 0 || i >= nx_tot - 1) return {0,0,0,0};

    const Conserved& Um = grid.U[grid.idx(i - 1, j)];
    const Conserved& U0 = grid.U[grid.idx(i,     j)];
    const Conserved& Up = grid.U[grid.idx(i + 1, j)];
    return minmod_vec(sub(U0, Um), sub(Up, U0));
}

static inline Conserved slope_y_safe(const Grid& grid, int i, int j) {
    const int ny_tot = grid.ny + 2 * grid.ng;
    if (j <= 0 || j >= ny_tot - 1) return {0,0,0,0};

    const Conserved& Um = grid.U[grid.idx(i, j - 1)];
    const Conserved& U0 = grid.U[grid.idx(i, j    )];
    const Conserved& Up = grid.U[grid.idx(i, j + 1)];
    return minmod_vec(sub(U0, Um), sub(Up, U0));
}

// ----------------- RHS computation -----------------
void compute_rhs(const Grid& grid, std::vector<Conserved>& rhs) {
    DBG_PRINT("computing rhs");

    const int nx_tot = grid.nx + 2 * grid.ng;
    const int ny_tot = grid.ny + 2 * grid.ng;
    const int ng     = grid.ng;

    rhs.assign(nx_tot * ny_tot, Conserved{0,0,0,0});

    std::vector<Conserved> Fx((nx_tot + 1) * ny_tot);
    std::vector<Conserved> Gy(nx_tot * (ny_tot + 1));

    auto Fx_idx = [&](int iface, int j) { return iface + (nx_tot + 1) * j; };
    auto Gy_idx = [&](int i, int jf) { return i + nx_tot * jf; };

    const int iface0 = ng;
    const int iface1 = ng + grid.nx;
    const int jface0 = ng;
    const int jface1 = ng + grid.ny;

    for (int j = 0; j < ny_tot; ++j) {
        for (int iface = iface0; iface <= iface1; ++iface) {
            int iL = iface - 1;
            int iR = iface;

            const Conserved ULc = grid.U[grid.idx(iL, j)];
            const Conserved URc = grid.U[grid.idx(iR, j)];

            const Conserved UL = add(ULc, mul(slope_x_safe(grid, iL, j), 0.5));
            const Conserved UR = sub(URc, mul(slope_x_safe(grid, iR, j), 0.5));

            Fx[Fx_idx(iface, j)] = hll_flux_x(UL, UR);
        }
    }

    for (int i = 0; i < nx_tot; ++i) {
        for (int jface = jface0; jface <= jface1; ++jface) {
            int jB = jface - 1;
            int jT = jface;

            const Conserved UBc = grid.U[grid.idx(i, jB)];
            const Conserved UTc = grid.U[grid.idx(i, jT)];

            const Conserved UL = add(UBc, mul(slope_y_safe(grid, i, jB), 0.5));
            const Conserved UR = sub(UTc, mul(slope_y_safe(grid, i, jT), 0.5));

            Gy[Gy_idx(i, jface)] = hll_flux_y(UL, UR);
        }
    }

    for (int j = ng; j < ng + grid.ny; ++j) {
        for (int i = ng; i < ng + grid.nx; ++i) {
            rhs[grid.idx(i,j)] =
                sub(
                    sub({0,0,0,0}, mul(sub(Fx[Fx_idx(i+1,j)], Fx[Fx_idx(i,j)]), 1.0/grid.dx)),
                    mul(sub(Gy[Gy_idx(i,j+1)], Gy[Gy_idx(i,j)]), 1.0/grid.dy)
                );
        }
    }
}

// ----------------- RK2 time stepping -----------------
void advance_one_step(Grid& grid, double dt) {
    DBG_PRINT("advance one step");

    apply_boundary_conditions(grid);

    std::vector<Conserved> k1;
    compute_rhs(grid, k1);

    Grid tmp = grid;
    for (int j = grid.ng; j < grid.ng + grid.ny; ++j)
        for (int i = grid.ng; i < grid.ng + grid.nx; ++i)
            tmp.U[tmp.idx(i,j)] = add(grid.U[grid.idx(i,j)], mul(k1[grid.idx(i,j)], dt));

    apply_boundary_conditions(tmp);

    std::vector<Conserved> k2;
    compute_rhs(tmp, k2);

    for (int j = grid.ng; j < grid.ng + grid.ny; ++j)
        for (int i = grid.ng; i < grid.ng + grid.nx; ++i)
            grid.U[grid.idx(i,j)] =
                add(grid.U[grid.idx(i,j)], mul(add(k1[grid.idx(i,j)], k2[grid.idx(i,j)]), 0.5*dt));

    apply_boundary_conditions(grid);
}
