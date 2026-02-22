#pragma once
#include "types.hpp"
#include <mpi.h>

struct MpiDecomp {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0;
  int size = 1;

  int nbr_down = MPI_PROC_NULL; // rank-1
  int nbr_up   = MPI_PROC_NULL; // rank+1

  int nx_global = 0;
  int ny_global = 0;

  int ny_local  = 0;   // interior rows on this rank
  int y0_global = 0;   // starting global interior row index for this rank

  int ng = 2;
};

MpiDecomp make_y_slab_decomp(int nx_global, int ny_global, int ng,
                            MPI_Comm comm = MPI_COMM_WORLD);

// Halo exchange for ng rows (y direction)
void exchange_halo_y(Grid& grid, const MpiDecomp& mp);

// Apply transmissive BC on global physical boundaries only
void apply_physical_bc_mpi(Grid& grid, const MpiDecomp& mp);

// Local max(|v|+c) for dt
double compute_amax_local(const Grid& grid);

// Global dt using MPI_Allreduce (MAX)
double compute_dt_mpi(const Grid& grid, const MpiDecomp& mp, double cfl);

// One step of MUSCL-Hancock + HLL (MPI version, uses ghost cells)
void advance_one_step_mpi(Grid& grid, const MpiDecomp& mp, double dt);
