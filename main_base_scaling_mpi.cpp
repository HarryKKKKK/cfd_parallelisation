#include "types.hpp"
#include "physics.hpp"
#include "init.hpp"
#include "mpi_solver.hpp"

#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    const int nx_global = 500;
    const int ny_base = 197;
    const int ng = 2;
    const double Lx = 0.225;
    const double Ly_base = 0.089;
    const double cfl = 0.4;
    const double t_end = 0.0011741;

    int copies_y = world_size; 
    if (argc >= 2) {
        copies_y = std::atoi(argv[1]);
    }

    if (copies_y <= 0) {
        if (world_rank == 0) {
            std::cerr << "Error: copies_y must be positive.\n";
        }
        MPI_Finalize();
        return 1;
    }

    const int ny_global = ny_base * copies_y;
    const double Ly = Ly_base * copies_y;

    MpiDomain mp = make_mpi_domain_y_slab(nx_global, ny_global, ng, MPI_COMM_WORLD);

    const double dx = Lx / nx_global;
    const double dy = Ly / ny_global;
    const double y0_local = mp.y0_global * dy;

    Grid grid;
    grid.init(nx_global, mp.ny_local, ng, Lx, Ly, 0.0, y0_local);
    grid.dx = dx;
    grid.dy = dy;

    initialize_shock_bubble_weak(grid, copies_y);

    exchange_halo_y_mpi(grid, mp);
    apply_boundary_conditions_mpi(grid, mp);

    if (mp.rank == 0) {
        std::cout << std::setprecision(12);
        std::cout << "[INIT] nx=" << nx_global
                  << " ny=" << ny_global
                  << " copies_y=" << copies_y
                  << " ng=" << ng
                  << " Lx=" << Lx
                  << " Ly=" << Ly
                  << " dx=" << grid.dx
                  << " dy=" << grid.dy
                  << " cfl=" << cfl
                  << " t_end=" << t_end
                  << std::endl;
    }

    int step = 0;
    double t = 0.0;

    MPI_Barrier(mp.comm);
    const double t0 = MPI_Wtime();

    double dt = compute_dt_mpi(grid, mp, cfl);

    if (mp.rank == 0) {
        std::cout << "[INIT] initial dt=" << dt << std::endl;
    }

    while (t < t_end) {
        if (t + dt > t_end) {
            dt = t_end - t;
        }

        // 仅由 Rank 0 负责打印 LOOP 进度条，使用 MPI_Wtime
        if (step % 500 == 0 && mp.rank == 0) {
            double wall_so_far = MPI_Wtime() - t0;
            std::cout << "[LOOP] step=" << step
                      << " t=" << t
                      << " dt=" << dt
                      << " wall_so_far=" << wall_so_far
                      << std::endl;
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
                  << " copies_y=" << copies_y
                  << " steps=" << step
                  << " wall_seconds=" << wall_max
                  << std::endl;
    }

    MPI_Finalize();
    return 0;
}