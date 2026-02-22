#include "types.hpp"
#include "init.hpp"
#include "utils.hpp"
#include "mpi_solver.hpp"

#include <mpi.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

// Pack local interior (no ghosts) into a contiguous buffer [ny_local][nx_global]
static std::vector<Conserved> pack_local_interior(const Grid& grid) {
  const int ng = grid.ng;
  const int nx = grid.nx; // interior nx (should equal nx_global)
  const int ny = grid.ny; // local interior ny

  std::vector<Conserved> buf(static_cast<size_t>(nx) * static_cast<size_t>(ny));

  for (int j = 0; j < ny; ++j) {
    // local interior row j corresponds to grid index J = ng + j
    const int J = ng + j;
    const Conserved* row_ptr = &grid.U[grid.idx(ng, J)];
    std::copy(row_ptr, row_ptr + nx, &buf[static_cast<size_t>(j) * nx]);
  }
  return buf;
}

// Build a global Grid on rank0 and write a single CSV using your existing writer
static void write_global_csv_rank0(const std::vector<Conserved>& global_interior,
                                  int nx_global, int ny_global,
                                  int ng,
                                  double Lx, double Ly,
                                  double dx, double dy,
                                  const std::string& filename) {
  Grid g;
  g.init(nx_global, ny_global, ng, Lx, Ly); // will set dx/dy internally, we override next
  g.dx = dx;
  g.dy = dy;
  g.x0 = 0.0;
  g.y0 = 0.0;

  // Fill interior (no ghosts) into g.U
  for (int j = 0; j < ny_global; ++j) {
    const int J = ng + j;
    for (int i = 0; i < nx_global; ++i) {
      const int I = ng + i;
      g.U[g.idx(I, J)] = global_interior[static_cast<size_t>(j) * nx_global + i];
    }
  }

  // Apply transmissive BC on all four sides for completeness
  // (If your writer only writes interior, this is not strictly necessary.)
  // You can reuse your serial BC if it is in solver.cpp, but we keep it simple:
  // left/right
  const int nx_tot = nx_global + 2 * ng;
  const int ny_tot = ny_global + 2 * ng;
  for (int J = 0; J < ny_tot; ++J) {
    for (int gcol = 0; gcol < ng; ++gcol) {
      g.U[g.idx(gcol, J)]              = g.U[g.idx(ng, J)];
      g.U[g.idx(nx_tot - 1 - gcol, J)] = g.U[g.idx(nx_tot - 1 - ng, J)];
    }
  }
  // bottom/top
  for (int I = 0; I < nx_tot; ++I) {
    for (int grow = 0; grow < ng; ++grow) {
      g.U[g.idx(I, grow)]              = g.U[g.idx(I, ng)];
      g.U[g.idx(I, ny_tot - 1 - grow)] = g.U[g.idx(I, ny_tot - 1 - ng)];
    }
  }

  write_grid_csv(g, filename);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // ----------------- Global problem settings (match your serial main) -----------------
  const int nx_global = 500;
  const int ny_global = 197;
  const int ng        = 2;

  const double Lx = 0.225;
  const double Ly = 0.089;

  const double t_end = 0.0011741;
  const double cfl   = 0.5;

  const int n_outputs = 10;
  const std::string out_dir = "res/mpi";

  // ----------------- MPI decomposition -----------------
  MpiDecomp mp = make_y_slab_decomp(nx_global, ny_global, ng, MPI_COMM_WORLD);

  if (mp.rank == 0) {
    std::filesystem::create_directories(out_dir);
    std::cout << "MPI size = " << mp.size << "\n";
    std::cout << "Writing output to: " << out_dir << "\n";
  }

  // Global spacing (IMPORTANT: use global dx/dy on all ranks)
  const double dx = Lx / nx_global;
  const double dy = Ly / ny_global;

  // Local y0 in physical coordinates (so init uses correct y)
  const double y0_local = mp.y0_global * dy;

  // ----------------- Init local grid -----------------
  Grid grid;
  // NOTE: grid.init computes dx/dy from Lx/Ly and local ny; we override to global next.
  grid.init(nx_global, mp.ny_local, ng, Lx, Ly, 0.0, y0_local);
  grid.dx = dx;
  grid.dy = dy;

  // Initialise shock-bubble on local subdomain using physical coords
  initialize_shock_bubble(grid);

  // Ensure ghosts are valid before first output/step
  exchange_halo_y(grid, mp);
  apply_physical_bc_mpi(grid, mp);

  // ----------------- Output schedule -----------------
  int step = 0;
  double t = 0.0;
  const int output_every = std::max(1, static_cast<int>(std::ceil((t_end / (n_outputs)) / 1e-12))); // placeholder; we'll do time-based below

  // We'll output when crossing evenly spaced times
  std::vector<double> out_times;
  out_times.reserve(n_outputs + 1);
  for (int k = 0; k <= n_outputs; ++k) {
    out_times.push_back(t_end * (static_cast<double>(k) / n_outputs));
  }
  int next_out_idx = 0;

  auto do_output = [&](int step_id, double time_now) {
    // Pack local interior
    std::vector<Conserved> sendbuf = pack_local_interior(grid);

    // Prepare recv on rank0
    std::vector<int> recvcounts, displs;
    std::vector<Conserved> recvbuf;

    if (mp.rank == 0) {
      recvcounts.resize(mp.size);
      displs.resize(mp.size);

      // Recompute the same ny_local per rank to form recvcounts/displs
      const int base = ny_global / mp.size;
      const int rem  = ny_global % mp.size;

      int disp = 0;
      for (int r = 0; r < mp.size; ++r) {
        const int ny_r = base + (r < rem ? 1 : 0);
        recvcounts[r] = nx_global * ny_r; // in Conserved units
        displs[r]     = disp;
        disp         += recvcounts[r];
      }

      recvbuf.resize(static_cast<size_t>(nx_global) * static_cast<size_t>(ny_global));
    }

    // Gather to rank0
    MPI_Gatherv(sendbuf.data(),
                static_cast<int>(sendbuf.size()),
                mpi_conserved_type(),
                (mp.rank == 0 ? recvbuf.data() : nullptr),
                (mp.rank == 0 ? recvcounts.data() : nullptr),
                (mp.rank == 0 ? displs.data() : nullptr),
                mpi_conserved_type(),
                0,
                mp.comm);

    // Rank0 writes single CSV
    if (mp.rank == 0) {
      const std::string fname = make_filename(out_dir, step_id, time_now);
      write_global_csv_rank0(recvbuf, nx_global, ny_global, ng, Lx, Ly, dx, dy, fname);
      std::cout << "Wrote: " << fname << "\n";
    }
  };

  // Output at t=0
  do_output(step, t);
  next_out_idx = 1;

  // ----------------- Time loop -----------------
  while (t < t_end) {
    // global dt
    const double dt = compute_dt_mpi(grid, mp, cfl);

    // clamp last step
    const double dt_eff = (t + dt > t_end) ? (t_end - t) : dt;

    advance_one_step_mpi(grid, mp, dt_eff);

    t += dt_eff;
    step++;

    // Output when crossing scheduled times
    while (next_out_idx < static_cast<int>(out_times.size()) && t >= out_times[next_out_idx] - 1e-15) {
      do_output(step, t);
      next_out_idx++;
    }
  }

  if (mp.rank == 0) {
    std::cout << "Finished. Steps = " << step << "\n";
  }

  MPI_Finalize();
  return 0;
}
