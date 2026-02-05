#include "utils.hpp"
#include "physics.hpp"   
#include <fstream>
#include <iomanip>

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
}

std::string make_filename(const std::string& dir, int step, double t) {
    std::ostringstream oss;
    oss << dir << "/step_" << std::setw(6) << std::setfill('0') << step
        << "_t_" << std::fixed << std::setprecision(6) << t << ".csv";
    return oss.str();
}