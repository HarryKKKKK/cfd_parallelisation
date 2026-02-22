#include "mpi_solver.hpp"
#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// ------------------------------
// Small helpers
// ------------------------------
static inline Conserved addC(const Conserved& a, const Conserved& b) {
  return {a.rho + b.rho, a.rhou + b.rhou, a.rhov + b.rhov, a.E + b.E};
}
static inline Conserved subC(const Conserved& a, const Conserved& b) {
  return {a.rho - b.rho, a.rhou - b.rhou, a.rhov - b.rhov, a.E - b.E};
}
static inline Conserved mulC(const Conserved& a, double s) {
  return {a.rho * s, a.rhou * s, a.rhov * s, a.E * s};
}

static inline double minmod(double a, double b) {
  if (a * b <= 0.0) return 0.0;
  return (std::abs(a) < std::abs(b)) ? a : b;
}
static inline Conserved minmod_vec(const Conserved& a, const Conserved& b) {
  return {
    minmod(a.rho,  b.rho),
    minmod(a.rhou, b.rhou),
    minmod(a.rhov, b.rhov),
    minmod(a.E,    b.E)
  };
}

// MPI datatype for Conserved = 4 doubles
static MPI_Datatype mpi_conserved_type() {
  static MPI_Datatype T = MPI_DATATYPE_NULL;
  static bool inited = false;
  if (!inited) {
    MPI_Type_contiguous(4, MPI_DOUBLE, &T);
    MPI_Type_commit(&T);
    inited = true;
  }
  return T;
}

MpiDecomp make_y_slab_decomp(int nx_global, int ny_global, int ng, MPI_Comm comm) {
  MpiDecomp mp;
  mp.comm = comm;
  MPI_Comm_rank(comm, &mp.rank);
  MPI_Comm_size(comm, &mp.size);

  mp.nx_global = nx_global;
  mp.ny_global = ny_global;
  mp.ng = ng;

  const int base = ny_global / mp.size;
  const int rem  = ny_global % mp.size;

  mp.ny_local  = base + (mp.rank < rem ? 1 : 0);
  mp.y0_global = mp.rank * base + std::min(mp.rank, rem);

  mp.nbr_down = (mp.rank == 0) ? MPI_PROC_NULL : (mp.rank - 1);
  mp.nbr_up   = (mp.rank == mp.size - 1) ? MPI_PROC_NULL : (mp.rank + 1);

  return mp;
}

void exchange_halo_y(Grid& grid, const MpiDecomp& mp) {
  const int ng = grid.ng;
  if (ng != mp.ng) throw std::runtime_error("exchange_halo_y: grid.ng != mp.ng");
  if (grid.ny != mp.ny_local) throw std::runtime_error("exchange_halo_y: grid.ny != mp.ny_local");

  const int nx_tot = grid.nx + 2 * ng;
  const int count  = ng * nx_tot; // in Conserved units
  MPI_Datatype T   = mpi_conserved_type();

  // interior bottom slab starts at j = ng
  Conserved* send_down = &grid.U[grid.idx(0, ng)];
  // interior top slab starts at j = ng + ny - ng
  Conserved* send_up   = &grid.U[grid.idx(0, ng + grid.ny - ng)];

  // ghost bottom slab starts at j = 0
  Conserved* recv_down = &grid.U[grid.idx(0, 0)];
  // ghost top slab starts at j = ng + ny
  Conserved* recv_up   = &grid.U[grid.idx(0, ng + grid.ny)];

  MPI_Sendrecv(send_down, count, T, mp.nbr_down, 100,
               recv_down, count, T, mp.nbr_down, 101,
               mp.comm, MPI_STATUS_IGNORE);

  MPI_Sendrecv(send_up, count, T, mp.nbr_up, 101,
               recv_up, count, T, mp.nbr_up, 100,
               mp.comm, MPI_STATUS_IGNORE);
}

void apply_physical_bc_mpi(Grid& grid, const MpiDecomp& mp) {
  const int ng = grid.ng;
  const int nx_tot = grid.nx + 2 * ng;
  const int ny_tot = grid.ny + 2 * ng;

  // x transmissive (all ranks)
  for (int j = 0; j < ny_tot; ++j) {
    for (int g = 0; g < ng; ++g) {
      grid.U[grid.idx(g, j)]              = grid.U[grid.idx(ng, j)];
      grid.U[grid.idx(nx_tot - 1 - g, j)] = grid.U[grid.idx(nx_tot - 1 - ng, j)];
    }
  }

  // y transmissive only at global boundaries
  if (mp.nbr_down == MPI_PROC_NULL) {
    for (int i = 0; i < nx_tot; ++i) {
      for (int g = 0; g < ng; ++g) {
        grid.U[grid.idx(i, g)] = grid.U[grid.idx(i, ng)];
      }
    }
  }
  if (mp.nbr_up == MPI_PROC_NULL) {
    for (int i = 0; i < nx_tot; ++i) {
      for (int g = 0; g < ng; ++g) {
        grid.U[grid.idx(i, ny_tot - 1 - g)] = grid.U[grid.idx(i, ny_tot - 1 - ng)];
      }
    }
  }
}

double compute_amax_local(const Grid& grid) {
  const int ng = grid.ng;
  const int i0 = ng, i1 = ng + grid.nx - 1;
  const int j0 = ng, j1 = ng + grid.ny - 1;

  double amax = 0.0;
  for (int j = j0; j <= j1; ++j) {
    for (int i = i0; i <= i1; ++i) {
      const Conserved& U = grid.U[grid.idx(i, j)];
      const Primitive W  = conserved_to_primitive(U);
      const double vmag  = std::sqrt(W.u*W.u + W.v*W.v);
      const double cs    = sound_speed(W);
      amax = std::max(amax, vmag + cs);
    }
  }
  return amax;
}

double compute_dt_mpi(const Grid& grid, const MpiDecomp& mp, double cfl) {
  const double h = std::min(grid.dx, grid.dy);
  const double a_local = compute_amax_local(grid);

  double a_global = 0.0;
  MPI_Allreduce(&a_local, &a_global, 1, MPI_DOUBLE, MPI_MAX, mp.comm);

  if (a_global <= 0.0) return 1e-12;
  return cfl * h / a_global;
}

// ------------------------------
// MUSCL-Hancock + HLL with ghost cells
// ------------------------------
namespace {
struct WorkspaceMPI {
  int nx=-1, ny=-1, ng=-1;
  std::size_t nTot=0;

  std::vector<Conserved> Umx, Upx, Umy, Upy;
  std::vector<Conserved> Uhat_mx, Uhat_px, Uhat_my, Uhat_py;
  std::vector<Conserved> Fx, Gy;

  void ensure(const Grid& grid) {
    const int nx_tot = grid.nx + 2*grid.ng;
    const int ny_tot = grid.ny + 2*grid.ng;
    const std::size_t n = static_cast<std::size_t>(nx_tot) * static_cast<std::size_t>(ny_tot);
    if (grid.nx==nx && grid.ny==ny && grid.ng==ng && n==nTot) return;

    nx = grid.nx; ny = grid.ny; ng = grid.ng; nTot = n;
    Umx.resize(nTot); Upx.resize(nTot); Umy.resize(nTot); Upy.resize(nTot);
    Uhat_mx.resize(nTot); Uhat_px.resize(nTot); Uhat_my.resize(nTot); Uhat_py.resize(nTot);

    // Interfaces on interior domain only:
    // x-interfaces: (nx+1)*ny over interior rows
    Fx.resize(static_cast<std::size_t>(nx+1) * static_cast<std::size_t>(ny));
    // y-interfaces: nx*(ny+1)
    Gy.resize(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny+1));
  }
};

static WorkspaceMPI ws;

static inline std::size_t tot_id(const Grid& grid, int I, int J) {
  // I,J are in [0..nx+2ng-1], [0..ny+2ng-1]
  return static_cast<std::size_t>(grid.idx(I, J));
}
static inline std::size_t fx_id(int nx, int iface_i, int j) {
  // iface_i = 0..nx, j = 0..ny-1 (interior indexing)
  return static_cast<std::size_t>(iface_i + (nx+1)*j);
}
static inline std::size_t gy_id(int nx, int i, int iface_j) {
  // i = 0..nx-1, iface_j = 0..ny
  return static_cast<std::size_t>(i + nx*iface_j);
}
} // namespace

void advance_one_step_mpi(Grid& grid, const MpiDecomp& mp, double dt) {
  // Update halos & physical BC before using neighbours
  exchange_halo_y(grid, mp);
  apply_physical_bc_mpi(grid, mp);

  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  ws.ensure(grid);

  const double half = 0.5;
  const double cx   = dt / (2.0 * grid.dx);
  const double cy   = dt / (2.0 * grid.dy);

  // (I)+(II) reconstruction + Hancock predictor on interior cells only
  for (int j = 0; j < ny; ++j) {
    const int J = ng + j;
    for (int i = 0; i < nx; ++i) {
      const int I = ng + i;
      const auto id = tot_id(grid, I, J);

      const Conserved& Uc = grid.U[grid.idx(I, J)];
      const Conserved& Ul = grid.U[grid.idx(I-1, J)];
      const Conserved& Ur = grid.U[grid.idx(I+1, J)];
      const Conserved& Ub = grid.U[grid.idx(I, J-1)];
      const Conserved& Ut = grid.U[grid.idx(I, J+1)];

      const Conserved sx = minmod_vec(subC(Uc, Ul), subC(Ur, Uc));
      const Conserved sy = minmod_vec(subC(Uc, Ub), subC(Ut, Uc));

      const Conserved Umx = subC(Uc, mulC(sx, half));
      const Conserved Upx = addC(Uc, mulC(sx, half));
      const Conserved Umy = subC(Uc, mulC(sy, half));
      const Conserved Upy = addC(Uc, mulC(sy, half));

      ws.Umx[id] = Umx; ws.Upx[id] = Upx;
      ws.Umy[id] = Umy; ws.Upy[id] = Upy;

      const Conserved dF   = subC(flux_x(Umx), flux_x(Upx));
      const Conserved dG   = subC(flux_y(Umy), flux_y(Upy));
      const Conserved corr = addC(mulC(dF, cx), mulC(dG, cy));

      ws.Uhat_mx[id] = addC(Umx, corr);
      ws.Uhat_px[id] = addC(Upx, corr);
      ws.Uhat_my[id] = addC(Umy, corr);
      ws.Uhat_py[id] = addC(Upy, corr);
    }
  }

  // (IIIa) x-interfaces fluxes (iface=0..nx, j=0..ny-1)
  for (int j = 0; j < ny; ++j) {
    const int J = ng + j;
    for (int iface = 0; iface <= nx; ++iface) {
      const Conserved *ULp, *URp;

      if (iface == 0) {
        const auto id0 = tot_id(grid, ng + 0, J);
        ULp = &ws.Uhat_px[id0];
        URp = &ws.Uhat_mx[id0];
      } else if (iface == nx) {
        const auto idN = tot_id(grid, ng + (nx - 1), J);
        ULp = &ws.Uhat_px[idN];
        URp = &ws.Uhat_mx[idN];
      } else {
        const auto idL = tot_id(grid, ng + (iface - 1), J);
        const auto idR = tot_id(grid, ng + iface,       J);
        ULp = &ws.Uhat_px[idL];
        URp = &ws.Uhat_mx[idR];
      }

      ws.Fx[fx_id(nx, iface, j)] = hll_flux_x(*ULp, *URp);
    }
  }

  // (IIIb) y-interfaces fluxes (iface=0..ny, i=0..nx-1)
  for (int iface = 0; iface <= ny; ++iface) {
    const int Jiface = ng + iface; // interface between Jiface-1 and Jiface
    for (int i = 0; i < nx; ++i) {
      const int I = ng + i;

      const Conserved *ULp, *URp;

      if (iface == 0) {
        // bottom boundary interface: between ghost row (ng-1) and first interior (ng)
        const auto id0 = tot_id(grid, I, ng + 0);
        // Using transmissive BC, this is effectively interior vs itself (already copied into ghosts)
        ULp = &ws.Uhat_py[id0];
        URp = &ws.Uhat_my[id0];
      } else if (iface == ny) {
        const auto idN = tot_id(grid, I, ng + (ny - 1));
        ULp = &ws.Uhat_py[idN];
        URp = &ws.Uhat_my[idN];
      } else {
        const auto idB = tot_id(grid, I, ng + (iface - 1));
        const auto idT = tot_id(grid, I, ng + iface);
        ULp = &ws.Uhat_py[idB];
        URp = &ws.Uhat_my[idT];
      }

      ws.Gy[gy_id(nx, i, iface)] = hll_flux_y(*ULp, *URp);
    }
  }

  // (IV) FV update (interior only)
  const double dtdx = dt / grid.dx;
  const double dtdy = dt / grid.dy;

  for (int j = 0; j < ny; ++j) {
    const int J = ng + j;
    for (int i = 0; i < nx; ++i) {
      const int I = ng + i;

      const Conserved& F_im = ws.Fx[fx_id(nx, i,   j)];
      const Conserved& F_ip = ws.Fx[fx_id(nx, i+1, j)];

      const Conserved& G_jm = ws.Gy[gy_id(nx, i, j)];
      const Conserved& G_jp = ws.Gy[gy_id(nx, i, j+1)];

      const Conserved update =
        addC(mulC(subC(F_ip, F_im), dtdx),
             mulC(subC(G_jp, G_jm), dtdy));

      grid.U[grid.idx(I, J)] = subC(grid.U[grid.idx(I, J)], update);
    }
  }

  // refresh halos for next step consistency
  exchange_halo_y(grid, mp);
  apply_physical_bc_mpi(grid, mp);
}
