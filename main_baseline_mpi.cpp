#include "types.hpp"
#include "physics.hpp"
#include "init.hpp"
#include "mpi_solver.hpp"

#include <mpi.h>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int nx_global = 500;
    int ny_global = 197;
    const int ng = 2;
    const double Lx = 0.225;
    const double Ly = 0.089;
    const double cfl = 0.4;
    double t_end = 0.0011741;

    // optional command line:
    // mpirun -np 4 ./mpi_scaling.exe [nx] [ny] [t_end]
    if (argc >= 3) {
        nx_global = std::atoi(argv[1]);
        ny_global = std::atoi(argv[2]);
    }
    if (argc >= 4) {
        t_end = std::atof(argv[3]);
    }

    MpiDomain mp = make_mpi_domain_y_slab(nx_global, ny_global, ng, MPI_COMM_WORLD);

    const double dx = Lx / nx_global;
    const double dy = Ly / ny_global;
    const double y0_local = mp.y0_global * dy;

    Grid grid;
    grid.init(nx_global, mp.ny_local, ng, Lx, Ly, 0.0, y0_local);
    grid.dx = dx;
    grid.dy = dy;

    initialize_shock_bubble(grid);

    exchange_halo_y_mpi(grid, mp);
    apply_boundary_conditions_mpi(grid, mp);

    int step = 0;
    double t = 0.0;

    MPI_Barrier(mp.comm);
    const double t0 = MPI_Wtime();

    double dt = compute_dt_mpi(grid, mp, cfl);

    while (t < t_end) {
        if (t + dt > t_end) {
            dt = t_end - t;
        }

        advance_one_step_mpi(grid, mp, dt);

        t += dt;
        ++step;

        dt = compute_dt_mpi(grid, mp, cfl);
    }

    MPI_Barrier(mp.comm);
    const double t1 = MPI_Wtime();
    const double wall_local = t1 - t0;

    double wall_max = 0.0;
    MPI_Reduce(&wall_local, &wall_max, 1, MPI_DOUBLE, MPI_MAX, 0, mp.comm);

    if (mp.rank == 0) {
        std::cout << "mode=mpi"
                  << " p=" << mp.size
                  << " nx=" << nx_global
                  << " ny=" << ny_global
                  << " steps=" << step
                  << " wall_seconds=" << wall_max
                  << std::endl;
    }

    MPI_Finalize();
    return 0;
}