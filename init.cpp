#include "init.hpp"
#include "constants.hpp"
#include "physics.hpp"

#include <cmath>

void initialize_shock_bubble(Grid& grid) {
    const double a0 = std::sqrt(phys::gamma * phys::p0 / phys::rho_air);

    const double rho2 = phys::rho_air * ((phys::gamma + 1.0) * phys::mach * phys::mach) /
                        ((phys::gamma - 1.0) * phys::mach * phys::mach + 2.0);

    const double p2 = phys::p0 * (1.0 + 2.0 * phys::gamma / (phys::gamma + 1.0) *
                                  (phys::mach * phys::mach - 1.0));

    const double u2 = 2.0 * a0 / (phys::gamma + 1.0) *
                      (phys::mach - 1.0 / phys::mach);

    // Fill cells
    for (int j = 0; j < grid.ny + 2 * grid.ng; ++j) {
        for (int i = 0; i < grid.nx + 2 * grid.ng; ++i) {

            const double x = grid.x0 + (i - grid.ng + 0.5) * grid.dx;
            const double y = grid.y0 + (j - grid.ng + 0.5) * grid.dy;

            Primitive W;

            // Post-shock air
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

            // Helium bubble overrides density only (pressure/velocity unchanged)
            const double rx = x - phys::bubble_x;
            const double ry = y - phys::bubble_y;
            const double r  = std::sqrt(rx * rx + ry * ry);

            if (r < phys::bubble_r) {
                W.rho = phys::rho_helium;
            }

            grid.U[grid.idx(i, j)] = primitive_to_conserved(W);
        }
    }
}
