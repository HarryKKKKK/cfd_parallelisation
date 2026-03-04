#include "mpi_solver.hpp"
#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

// ------------------------------
// MPI datatype for Conserved (safe w/ padding)
// ------------------------------
static MPI_Datatype mpi_conserved_type() {
  static MPI_Datatype T = MPI_DATATYPE_NULL;
  static bool inited = false;
  if (!inited) {
    Conserved dummy{};
    int blocklen[4] = {1, 1, 1, 1};
    MPI_Aint disp[4], base;
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.rho,  &disp[0]);
    MPI_Get_address(&dummy.rhou, &disp[1]);
    MPI_Get_address(&dummy.rhov, &disp[2]);
    MPI_Get_address(&dummy.E,    &disp[3]);
    for (int k = 0; k < 4; ++k) disp[k] -= base;
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(4, blocklen, disp, types, &T);
    MPI_Type_commit(&T);
    inited = true;
  }
  return T;
}

// ------------------------------
// Decomposition
// ------------------------------
MpiDomain make_mpi_domain_y_slab(int nx_global, int ny_global, int ng, MPI_Comm comm) {
  MpiDomain mp;
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

// ------------------------------
// Halo exchange in y (exchange ng rows including x ghost columns)
// ------------------------------
void exchange_halo_y_mpi(Grid& g, const MpiDomain& mp) {
  const int ng = g.ng;
  if (ng != mp.ng) throw std::runtime_error("exchange_halo_y_mpi: grid.ng != mp.ng");
  if (g.nx != mp.nx_global) throw std::runtime_error("exchange_halo_y_mpi: grid.nx != mp.nx_global");
  if (g.ny != mp.ny_local) throw std::runtime_error("exchange_halo_y_mpi: grid.ny != mp.ny_local");

  const int nx_tot = g.nx + 2 * ng;
  const int count  = nx_tot * ng; // Conserved entries in one slab

  MPI_Datatype T = mpi_conserved_type();

  Conserved* send_down = &g.U[g.idx(0, ng)];
  Conserved* recv_down = &g.U[g.idx(0, 0)];

  Conserved* send_up   = &g.U[g.idx(0, ng + g.ny - ng)];
  Conserved* recv_up   = &g.U[g.idx(0, ng + g.ny)];

  MPI_Sendrecv(send_down, count, T, mp.nbr_down, 100,
               recv_down, count, T, mp.nbr_down, 101,
               mp.comm, MPI_STATUS_IGNORE);

  MPI_Sendrecv(send_up,   count, T, mp.nbr_up,   101,
               recv_up,   count, T, mp.nbr_up,   100,
               mp.comm, MPI_STATUS_IGNORE);
}

// ------------------------------
// Transmissive BC on physical boundaries only
// ------------------------------
void apply_boundary_conditions_mpi(Grid& g, const MpiDomain& mp) {
  const int ng = g.ng;
  if (ng <= 0) return;

  const int nx = g.nx;
  const int ny = g.ny;
  const int nx_tot = nx + 2 * ng;
  const int ny_tot = ny + 2 * ng;

  const int iL = ng;
  const int iR = ng + nx - 1;
  const int jB = ng;
  const int jT = ng + ny - 1;

  // Left/right ghost columns: always (x not decomposed)
  for (int J = 0; J < ny_tot; ++J) {
    for (int gc = 0; gc < ng; ++gc) {
      g.U[g.idx(gc, J)]              = g.U[g.idx(iL, J)];
      g.U[g.idx(nx_tot - 1 - gc, J)] = g.U[g.idx(iR, J)];
    }
  }

  // Bottom ghost rows: only on global bottom
  if (mp.nbr_down == MPI_PROC_NULL) {
    for (int I = 0; I < nx_tot; ++I) {
      for (int gr = 0; gr < ng; ++gr) {
        g.U[g.idx(I, gr)] = g.U[g.idx(I, jB)];
      }
    }
  }

  // Top ghost rows: only on global top
  if (mp.nbr_up == MPI_PROC_NULL) {
    for (int I = 0; I < nx_tot; ++I) {
      for (int gr = 0; gr < ng; ++gr) {
        g.U[g.idx(I, ny_tot - 1 - gr)] = g.U[g.idx(I, jT)];
      }
    }
  }
}

// ------------------------------
// minmod limiter (EXACT same as serial)
// ------------------------------
static inline double minmod2(double a, double b) {
  if (a * b <= 0.0) return 0.0;
  return (std::abs(a) < std::abs(b)) ? a : b;
}

double minmod(double dL, double dR) {
  const double a = 0.5 * (dL + dR);
  const double b = 2.0 * dL;
  const double c = 2.0 * dR;
  return minmod2(a, minmod2(b, c));
}

static inline Conserved minmod_vec(const Conserved& a, const Conserved& b) {
  return {
    minmod(a.rho,  b.rho),
    minmod(a.rhou, b.rhou),
    minmod(a.rhov, b.rhov),
    minmod(a.E,    b.E)
  };
}

// ------------------------------
// dt = cfl*h / global_max(|v|+cs)
// ------------------------------
double compute_dt_mpi(const Grid& grid, const MpiDomain& mp, double cfl) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  const double h = std::min(grid.dx, grid.dy);
  double amax_local = 0.0;

  for (int j = ng; j < ny + ng; ++j) {
    for (int i = ng; i < nx + ng; ++i) {
      const Conserved& U = grid.U[ grid.idx(i, j) ];
      const Primitive  W = conserved_to_primitive(U);

      const double vmag = std::sqrt(W.u * W.u + W.v * W.v);
      const double cs   = sound_speed(W);

      amax_local = std::max(amax_local, vmag + cs);
    }
  }

  double amax = 0.0;
  MPI_Allreduce(&amax_local, &amax, 1, MPI_DOUBLE, MPI_MAX, mp.comm);

  if (amax <= 0.0) return 1e-12;
  return cfl * h / amax;
}

// ------------------------------
// Workspace (same as serial)
// ------------------------------
namespace {

struct Workspace {
  int nx = -1, ny = -1, ng = -1;
  std::size_t nU = 0;

  std::vector<Conserved> Umx, Upx, Umy, Upy;
  std::vector<Conserved> Uhat_mx, Uhat_px, Uhat_my, Uhat_py;

  std::vector<Conserved> Fx, Gy;

  void ensure(const Grid& grid) {
    if (grid.nx == nx && grid.ny == ny && grid.ng == ng && grid.U.size() == nU) return;
    nx = grid.nx; ny = grid.ny; ng = grid.ng; nU = grid.U.size();

    Umx.resize(nU); Upx.resize(nU); Umy.resize(nU); Upy.resize(nU);
    Uhat_mx.resize(nU); Uhat_px.resize(nU); Uhat_my.resize(nU); Uhat_py.resize(nU);

    Fx.resize(static_cast<std::size_t>(nx + 1) * static_cast<std::size_t>(ny));
    Gy.resize(static_cast<std::size_t>(nx)     * static_cast<std::size_t>(ny + 1));
  }
};

static Workspace ws;

static inline std::size_t fx_id(int nx, int iface_i, int j) {
  return static_cast<std::size_t>(iface_i + (nx + 1) * j);
}
static inline std::size_t gy_id(int nx, int i, int iface_j) {
  return static_cast<std::size_t>(i + nx * iface_j);
}

} // namespace

// ------------------------------
// One step: identical algorithmic structure to serial
// ------------------------------
void advance_one_step_mpi(Grid& grid, const MpiDomain& mp, double dt) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  if (ng < 2) {
    throw std::runtime_error("advance_one_step_mpi requires ng >= 2 (same as serial).");
  }

  // 1) fill ghost cells: (a) halo exchange, (b) transmissive BC on physical boundaries
  exchange_halo_y_mpi(grid, mp);
  apply_boundary_conditions_mpi(grid, mp);

  ws.ensure(grid);

  const double cx = dt / (2.0 * grid.dx);
  const double cy = dt / (2.0 * grid.dy);

  // predictor states on: I = ng-1 .. ng+nx, J = ng-1 .. ng+ny (inclusive)
  const int Imin = ng - 1;
  const int Imax = ng + nx;
  const int Jmin = ng - 1;
  const int Jmax = ng + ny;

  // (I) Reconstruction + (II) Hancock predictor
  for (int J = Jmin; J <= Jmax; ++J) {
    for (int I = Imin; I <= Imax; ++I) {
      const auto id = grid.idx(I, J);

      const Conserved& Uc = grid.U[id];

      const Conserved& Ul = grid.U[ grid.idx(I - 1, J) ];
      const Conserved& Ur = grid.U[ grid.idx(I + 1, J) ];
      const Conserved& Ub = grid.U[ grid.idx(I, J - 1) ];
      const Conserved& Ut = grid.U[ grid.idx(I, J + 1) ];

      const Conserved sx = minmod_vec(Uc - Ul, Ur - Uc);
      const Conserved sy = minmod_vec(Uc - Ub, Ut - Uc);

      const Conserved Umx = Uc - (sx * 0.5);
      const Conserved Upx = Uc + (sx * 0.5);
      const Conserved Umy = Uc - (sy * 0.5);
      const Conserved Upy = Uc + (sy * 0.5);

      ws.Umx[id] = Umx; ws.Upx[id] = Upx;
      ws.Umy[id] = Umy; ws.Upy[id] = Upy;

      const Conserved dF = flux_x(Umx) - flux_x(Upx);
      const Conserved dG = flux_y(Umy) - flux_y(Upy);

      const Conserved corr = (dF * cx) + (dG * cy);

      ws.Uhat_mx[id] = Umx + corr;
      ws.Uhat_px[id] = Upx + corr;
      ws.Uhat_my[id] = Umy + corr;
      ws.Uhat_py[id] = Upy + corr;
    }
  }

  // (IIIa) x-interfaces HLL fluxes: iface = 0..nx, j = 0..ny-1
  for (int j = 0; j < ny; ++j) {
    const int J = j + ng;
    for (int iface = 0; iface <= nx; ++iface) {
      const int IL = ng + iface - 1;
      const int IR = ng + iface;

      const auto idL = grid.idx(IL, J);
      const auto idR = grid.idx(IR, J);

      const Conserved& UL = ws.Uhat_px[idL];
      const Conserved& UR = ws.Uhat_mx[idR];

      ws.Fx[ fx_id(nx, iface, j) ] = hll_flux_x(UL, UR);
    }
  }

  // (IIIb) y-interfaces HLL fluxes: iface = 0..ny, i = 0..nx-1
  for (int iface = 0; iface <= ny; ++iface) {
    const int JB = ng + iface - 1;
    const int JT = ng + iface;

    for (int i = 0; i < nx; ++i) {
      const int I = i + ng;

      const auto idB = grid.idx(I, JB);
      const auto idT = grid.idx(I, JT);

      const Conserved& UL = ws.Uhat_py[idB];
      const Conserved& UR = ws.Uhat_my[idT];

      ws.Gy[ gy_id(nx, i, iface) ] = hll_flux_y(UL, UR);
    }
  }

  // (IV) FV update on PHYSICAL cells only
  const double dtdx = dt / grid.dx;
  const double dtdy = dt / grid.dy;

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const Conserved& F_im = ws.Fx[ fx_id(nx, i,   j) ];
      const Conserved& F_ip = ws.Fx[ fx_id(nx, i+1, j) ];

      const Conserved& G_jm = ws.Gy[ gy_id(nx, i, j) ];
      const Conserved& G_jp = ws.Gy[ gy_id(nx, i, j+1) ];

      const Conserved update = ((F_ip - F_im) * dtdx) + ((G_jp - G_jm) * dtdy);

      const int I = i + ng;
      const int J = j + ng;
      const auto id = grid.idx(I, J);

      grid.U[id] = grid.U[id] - update;
    }
  }
}