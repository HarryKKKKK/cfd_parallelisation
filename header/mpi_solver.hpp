#ifndef MPI_SOLVER_HPP
#define MPI_SOLVER_HPP

#include "types.hpp"
#include <mpi.h>

// ============================================================
// y-slab decomposition: full x-range on each rank, chunked in y
// ============================================================
struct MpiDomain {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0;
  int size = 1;

  int nx_global = 0;
  int ny_global = 0;
  int ng = 2;

  int ny_local  = 0;   // local interior rows
  int y0_global = 0;   // global start row of local interior

  // neighbour ranks in y-direction
  int nbr_down = MPI_PROC_NULL; // smaller global y
  int nbr_up   = MPI_PROC_NULL; // larger global y
};

// Build y-slab decomposition
MpiDomain make_mpi_domain_y_slab(int nx_global, int ny_global, int ng,
                                 MPI_Comm comm = MPI_COMM_WORLD);

// CFL timestep using GLOBAL max wave speed
double compute_dt_mpi(const Grid& grid, const MpiDomain& mp, double cfl);

// Physical transmissive BC
// - x boundaries: always on every rank
// - y boundaries: only on global bottom/top ranks
void apply_boundary_conditions_mpi(Grid& grid, const MpiDomain& mp);

// Exchange y halo rows between neighbouring ranks
void exchange_halo_y_mpi(Grid& grid, const MpiDomain& mp);

// Limiter (same as serial)
double minmod(double a, double b);

// One-step advance using the SAME MUSCL-Hancock + HLLC algorithm
void advance_one_step_mpi(Grid& grid, const MpiDomain& mp, double dt);

#endif // MPI_SOLVER_HPP