#ifndef MPI_SOLVER_HPP
#define MPI_SOLVER_HPP

#include "types.hpp"
#include <mpi.h>

// y-slab decomposition (full x, chunked y)
struct MpiDomain {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0, size = 1;

  int nx_global = 0;
  int ny_global = 0;
  int ng = 2;

  int ny_local  = 0;   // interior rows on this rank
  int y0_global = 0;   // starting global y (interior) index

  int nbr_down = MPI_PROC_NULL; // rank-1
  int nbr_up   = MPI_PROC_NULL; // rank+1
};

MpiDomain make_mpi_domain_y_slab(int nx_global, int ny_global, int ng,
                                 MPI_Comm comm = MPI_COMM_WORLD);

// CFL timestep using GLOBAL max(|v|+cs)
double compute_dt_mpi(const Grid& grid_local, const MpiDomain& mp, double cfl);

// Halo exchange (y only, exchange ng rows including x-ghost columns)
void exchange_halo_y_mpi(Grid& grid_local, const MpiDomain& mp);

// Transmissive BC on physical boundaries:
// - x boundaries: always (since x not decomposed)
// - y boundaries: only on global bottom/top ranks
void apply_boundary_conditions_mpi(Grid& grid_local, const MpiDomain& mp);

// Limiter (EXACT same as serial signature)
double minmod(double dL, double dR);

// One step: ghost-cell MUSCL-Hancock + HLL flux (loop structure same as serial)
void advance_one_step_mpi(Grid& grid_local, const MpiDomain& mp, double dt);

#endif // MPI_SOLVER_HPP