#include "init.hpp"
#include "constants.hpp"
#include "physics.hpp"

#include <cmath>

static void fill_transmissive_ghosts(Grid& g) {
    const int ng = g.ng;
    const int nx = g.nx;
    const int ny = g.ny;

    const int nx_tot = nx + 2 * ng;
    const int ny_tot = ny + 2 * ng;

    const int iL = ng;
    const int iR = ng + nx - 1;
    const int jB = ng;
    const int jT = ng + ny - 1;

    // Left/Right ghost columns (copy nearest interior column)
    for (int J = 0; J < ny_tot; ++J) {
        for (int gcol = 0; gcol < ng; ++gcol) {
            g.U[g.idx(gcol, J)] = g.U[g.idx(iL, J)];
            g.U[g.idx(nx_tot - 1 - gcol, J)] = g.U[g.idx(iR, J)];
        }
    }

    // Bottom/Top ghost rows (copy nearest interior row)
    for (int I = 0; I < nx_tot; ++I) {
        for (int grow = 0; grow < ng; ++grow) {
            g.U[g.idx(I, grow)] = g.U[g.idx(I, jB)];
            g.U[g.idx(I, ny_tot - 1 - grow)]  = g.U[g.idx(I, jT)];
        }
    }
}

void initialize_shock_bubble(Grid& grid) {
    const double a0 = std::sqrt(phys::gamma * phys::p0 / phys::rho_air);

    const double M2 = phys::mach * phys::mach;

    const double rho2 = phys::rho_air * ((phys::gamma + 1.0) * M2) /
                        ((phys::gamma - 1.0) * M2 + 2.0);

    const double p2 = phys::p0 * (1.0 + 2.0 * phys::gamma / (phys::gamma + 1.0) *
                                  (M2 - 1.0));

    const double u2 = 2.0 * a0 / (phys::gamma + 1.0) *
                      (phys::mach - 1.0 / phys::mach);

    const double R2 = phys::bubble_r * phys::bubble_r;

    // ------------------------------------------------------------
    // ONLY fill INTERIOR cells (do NOT touch ghost cells)
    // This avoids subtle asymmetry / inconsistencies from ghost init.
    //
    // Interior indices in storage are:
    //   I = ng .. ng+nx-1
    //   J = ng .. ng+ny-1
    // Physical cell-centre coords:
    //   x = x0 + (i + 0.5)*dx   where i = 0..nx-1
    //   y = y0 + (j + 0.5)*dy   where j = 0..ny-1
    // ------------------------------------------------------------
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {

            const int I = i + grid.ng;
            const int J = j + grid.ng;

            const double x = grid.x0 + (i + 0.5) * grid.dx;
            const double y = grid.y0 + (j + 0.5) * grid.dy;

            Primitive W;

            // Post-shock air (left of shock)
            if (x < phys::shock_x) {
                W.rho = rho2;
                W.u   = u2;
                W.v   = 0.0;
                W.p   = p2;
            } else {
                // Pre-shock air
                W.rho = phys::rho_air;
                W.u   = 0.0;
                W.v   = 0.0;
                W.p   = phys::p0;
            }

            // Helium bubble overrides density only
            const double rx = x - phys::bubble_x;
            const double ry = y - phys::bubble_y;
            const double r2 = rx * rx + ry * ry;

            // Use <= and squared distance to avoid tiny floating asymmetries
            if (r2 <= R2) {
                W.rho = phys::rho_helium;
            }

            grid.U[grid.idx(I, J)] = primitive_to_conserved(W);
        }
    }

    fill_transmissive_ghosts(grid);
}

void initialize_shock_bubble_weak(Grid& grid, int copies_y) {
    const double a0 = std::sqrt(phys::gamma * phys::p0 / phys::rho_air);

    const double M2 = phys::mach * phys::mach;

    const double rho2 = phys::rho_air * ((phys::gamma + 1.0) * M2) /
                        ((phys::gamma - 1.0) * M2 + 2.0);

    const double p2 = phys::p0 * (1.0 + 2.0 * phys::gamma / (phys::gamma + 1.0) *
                                  (M2 - 1.0));

    const double u2 = 2.0 * a0 / (phys::gamma + 1.0) *
                      (phys::mach - 1.0 / phys::mach);

    const double R2 = phys::bubble_r * phys::bubble_r;

    // base cell height of the original single-bubble problem
    const double Ly_cell = 0.089;

    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {

            const int I = i + grid.ng;
            const int J = j + grid.ng;

            const double x = grid.x0 + (i + 0.5) * grid.dx;
            const double y = grid.y0 + (j + 0.5) * grid.dy;

            // Map global y back into one repeated cell [0, Ly_cell)
            double y_cell = std::fmod(y, Ly_cell);
            if (y_cell < 0.0) {
                y_cell += Ly_cell;
            }

            Primitive W;

            // Post-shock air (left of shock)
            if (x < phys::shock_x) {
                W.rho = rho2;
                W.u   = u2;
                W.v   = 0.0;
                W.p   = p2;
            } else {
                // Pre-shock air
                W.rho = phys::rho_air;
                W.u   = 0.0;
                W.v   = 0.0;
                W.p   = phys::p0;
            }

            // Helium bubble overrides density only
            // IMPORTANT:
            // use y_cell instead of y, so the bubble is duplicated in y-direction
            const double rx = x - phys::bubble_x;
            const double ry = y_cell - phys::bubble_y;
            const double r2 = rx * rx + ry * ry;

            if (r2 <= R2) {
                W.rho = phys::rho_helium;
            }

            grid.U[grid.idx(I, J)] = primitive_to_conserved(W);
        }
    }

    fill_transmissive_ghosts(grid);
}