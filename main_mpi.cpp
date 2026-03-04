#include "types.hpp"
#include "physics.hpp"
#include "mpi_solver.hpp"
#include "init.hpp"
#include "utils.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

constexpr bool ENABLE_OUTPUT = true;

// Local interior pack: [ny_local][nx_global] contiguous
static std::vector<Conserved> pack_local_interior(const Grid& grid) {
  const int ng = grid.ng;
  const int nx = grid.nx;
  const int ny = grid.ny;

  std::vector<Conserved> buf(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));

  for (int j = 0; j < ny; ++j) {
    const int J = ng + j;
    const Conserved* src = &grid.U[ grid.idx(ng, J) ];
    std::copy(src, src + nx, &buf[static_cast<std::size_t>(j) * nx]);
  }
  return buf;
}

// MPI datatype for Conserved (same approach as mpi_solver.cpp, but local to main)
static MPI_Datatype mpi_conserved_type_main() {
  static MPI_Datatype T = MPI_DATATYPE_NULL;
  static bool inited = false;
  if (!inited) {
    Conserved dummy{};
    int blocklen[4] = {1,1,1,1};
    MPI_Aint disp[4], base;
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.rho,  &disp[0]);
    MPI_Get_address(&dummy.rhou, &disp[1]);
    MPI_Get_address(&dummy.rhov, &disp[2]);
    MPI_Get_address(&dummy.E,    &disp[3]);
    for (int k=0;k<4;++k) disp[k] -= base;
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(4, blocklen, disp, types, &T);
    MPI_Type_commit(&T);
    inited = true;
  }
  return T;
}

// On rank0: unpack global interior into full Grid interior
static void unpack_global_interior(Grid& g_full, const std::vector<Conserved>& all) {
  const int ng = g_full.ng;
  const int nx = g_full.nx;
  const int ny = g_full.ny;

  for (int j = 0; j < ny; ++j) {
    const Conserved* src = &all[static_cast<std::size_t>(j) * nx];
    Conserved* dst = &g_full.U[g_full.idx(ng, ng + j)];
    std::copy(src, src + nx, dst);
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const std::string out_dir = "res/mpi/hllc";
  if (ENABLE_OUTPUT && rank == 0) {
    std::filesystem::create_directories(out_dir);
  }

  // ----------------- Global settings (same as serial) -----------------
  const int NX_GLOBAL = 500;
  const int NY_GLOBAL = 197;
  const int NG        = 2;
  const double X_LEN  = 0.225;
  const double Y_LEN  = 0.089;

  // ----------------- Init MPI domain + local grid -----------------
  MpiDomain mp = make_mpi_domain_y_slab(NX_GLOBAL, NY_GLOBAL, NG, MPI_COMM_WORLD);

  Grid grid; // local grid on each rank
  grid.init(NX_GLOBAL, mp.ny_local, NG, X_LEN, Y_LEN);

  // ----------------- To be EXACTLY serial-identical at t=0: rank0 init full grid then scatter interior -----------------
  MPI_Datatype T = mpi_conserved_type_main();

  // compute sendcounts/displs (Conserved units) on rank0
  std::vector<int> sendcounts, displs;
  if (rank == 0) {
    sendcounts.resize(size);
    displs.resize(size);
    for (int r = 0; r < size; ++r) {
      const int base = NY_GLOBAL / size;
      const int rem  = NY_GLOBAL % size;
      const int ny_r = base + (r < rem ? 1 : 0);
      const int y0_r = r * base + std::min(r, rem);
      sendcounts[r] = ny_r * NX_GLOBAL;
      displs[r]     = y0_r * NX_GLOBAL;
    }
  }

  // rank0 builds full interior buffer and scatters
  std::vector<Conserved> local_interior(static_cast<std::size_t>(NX_GLOBAL) * mp.ny_local);

  if (size == 1) {
    initialize_shock_bubble(grid);
  } else {
    if (rank == 0) {
      Grid g_full;
      g_full.init(NX_GLOBAL, NY_GLOBAL, NG, X_LEN, Y_LEN);
      initialize_shock_bubble(g_full);

      std::vector<Conserved> full_interior(static_cast<std::size_t>(NX_GLOBAL) * NY_GLOBAL);
      for (int j = 0; j < NY_GLOBAL; ++j) {
        const Conserved* src = &g_full.U[g_full.idx(NG, NG + j)];
        Conserved* dst = &full_interior[static_cast<std::size_t>(j) * NX_GLOBAL];
        std::copy(src, src + NX_GLOBAL, dst);
      }

      MPI_Scatterv(full_interior.data(),
                   sendcounts.data(), displs.data(), T,
                   local_interior.data(), static_cast<int>(local_interior.size()), T,
                   0, MPI_COMM_WORLD);
    } else {
      MPI_Scatterv(nullptr, nullptr, nullptr, T,
                   local_interior.data(), static_cast<int>(local_interior.size()), T,
                   0, MPI_COMM_WORLD);
    }

    // unpack local interior into grid (leave ghosts; solver will fill)
    for (int j = 0; j < mp.ny_local; ++j) {
      Conserved* dst = &grid.U[grid.idx(NG, NG + j)];
      const Conserved* src = &local_interior[static_cast<std::size_t>(j) * NX_GLOBAL];
      std::copy(src, src + NX_GLOBAL, dst);
    }
  }

  int step = 0;
  double t = 0.0;

  // ----------------- output helper: gather interior to rank0 and write with existing utils -----------------
  auto write_global_snapshot = [&](int step_now, double t_now) {
    if (!ENABLE_OUTPUT) return;

    std::vector<Conserved> local = pack_local_interior(grid);

    std::vector<int> recvcounts, rdispls;
    std::vector<Conserved> all;

    if (rank == 0) {
      recvcounts.resize(size);
      rdispls.resize(size);
      for (int r = 0; r < size; ++r) {
        const int base = NY_GLOBAL / size;
        const int rem  = NY_GLOBAL % size;
        const int ny_r = base + (r < rem ? 1 : 0);
        const int y0_r = r * base + std::min(r, rem);
        recvcounts[r] = ny_r * NX_GLOBAL;
        rdispls[r]    = y0_r * NX_GLOBAL;
      }
      all.resize(static_cast<std::size_t>(NX_GLOBAL) * NY_GLOBAL);
    }

    MPI_Gatherv(local.data(), static_cast<int>(local.size()), T,
                rank == 0 ? all.data() : nullptr,
                rank == 0 ? recvcounts.data() : nullptr,
                rank == 0 ? rdispls.data() : nullptr,
                T, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      Grid g_full;
      g_full.init(NX_GLOBAL, NY_GLOBAL, NG, X_LEN, Y_LEN);
      unpack_global_interior(g_full, all);

      // We do NOT need serial apply_boundary_conditions for output if write_grid_csv writes interior.
      // If your write_grid_csv uses ghosts, then you should refactor BC into a shared bc.cpp.
      write_grid_csv(g_full, make_filename(out_dir, step_now, t_now));
    }
  };

  // Write t=0 (same behavior as serial)
  if (ENABLE_OUTPUT) {
    write_global_snapshot(step, t);
  }

  // ----------------- Physical time settings (same as serial) -----------------
  const double cfl = 0.4;

  const double T_collision = get_collision_time();
  const double T_ref       = get_time_ref();

  const double t_end = T_collision + 19.0 * T_ref;

  const double tstar_targets[] =
      {0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 9.4, 12.6, 17.4, 19.0};

  const int n_targets = static_cast<int>(sizeof(tstar_targets) / sizeof(tstar_targets[0]));

  double targets_phys[n_targets];
  for (int k = 0; k < n_targets; ++k) {
    targets_phys[k] = T_collision + tstar_targets[k] * T_ref;
  }

  int next_k = 0;

  // ----------------- Main loop (same structure as serial) -----------------
  double dt = compute_dt_mpi(grid, mp, cfl);

  while (t < t_end) {
    if (t + dt > t_end)
      dt = t_end - t;

    // force exact hit of next output time
    if (ENABLE_OUTPUT && next_k < n_targets) {
      const double t_out = targets_phys[next_k];
      if (t < t_out && t + dt > t_out) {
        dt = t_out - t;
      }
    }

    advance_one_step_mpi(grid, mp, dt);

    t += dt;
    step++;

    if (ENABLE_OUTPUT && next_k < n_targets) {
      const double t_out = targets_phys[next_k];

      if (t + 1e-15 >= t_out) {
        double t_star = (t - T_collision) / T_ref;

        if (rank == 0) {
          std::cout << "Writing t* = "
                    << t_star
                    << "  (physical t = "
                    << t << ")\n";
        }

        write_global_snapshot(step, t);

        next_k++;
      }
    }

    dt = compute_dt_mpi(grid, mp, cfl);
  }

  if (rank == 0) {
    std::cout << "Finished. Steps = " << step << std::endl;
  }

  MPI_Finalize();
  return 0;
}