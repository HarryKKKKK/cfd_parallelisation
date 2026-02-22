#include "utils.hpp"
#include "physics.hpp"   
#include <fstream>
#include <iomanip>
#include <iostream>

void write_grid_csv(const Grid& grid,
                    const std::string& filename)
{
    std::ofstream file(filename);
    file << std::setprecision(10);

    // CSV header
    file << "x,y,rho,u,v,p\n";

    for (int j = grid.ng; j < grid.ny + grid.ng; ++j) {
        for (int i = grid.ng; i < grid.nx + grid.ng; ++i) {

            const double x = grid.x0 + (i - grid.ng + 0.5) * grid.dx;
            const double y = grid.y0 + (j - grid.ng + 0.5) * grid.dy;

            const Conserved& U = grid.U[grid.idx(i, j)];
            const Primitive W = conserved_to_primitive(U);

            file << x << ","
                 << y << ","
                 << W.rho << ","
                 << W.u   << ","
                 << W.v   << ","
                 << W.p   << "\n";
        }
    }

    file.close();
    std::cout << "Wrote: " << filename << "\n";
}

std::string make_filename(const std::string& dir, int step, double t) {
    std::ostringstream oss;
    oss << dir << "/step_" << std::setw(6) << std::setfill('0') << step
        << "_t_" << std::fixed << std::setprecision(6) << t << ".csv";
    return oss.str();
}

double check_symmetry(const Grid& grid, double t, bool verbose = true)
{
    const int nx = grid.nx;
    const int ny = grid.ny;
    const int ng = grid.ng;

    double max_err = 0.0;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            const int I1 = i + ng;
            const int J1 = j + ng;

            const int I2 = i + ng;
            const int J2 = (ny - 1 - j) + ng;

            const double rho1 = grid.U[grid.idx(I1, J1)].rho;
            const double rho2 = grid.U[grid.idx(I2, J2)].rho;

            const double err = std::abs(rho1 - rho2);

            if (err > max_err)
                max_err = err;
        }
    }

    if (verbose) {
        std::cout << "Symmetry check at t = "
                  << t
                  << " : max |rho - rho_flip| = "
                  << max_err << std::endl;
    }

    return max_err;
}