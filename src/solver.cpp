// solver.cpp (no-negative-index MUSCL-Hancock)
#include "solver.hpp"
#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


// ------------------------------
// Small helpers
// ------------------------------
static inline Conserved add(const Conserved& a, const Conserved& b) {
  return {a.rho + b.rho, a.rhou + b.rhou, a.rhov + b.rhov, a.E + b.E};
}
static inline Conserved sub(const Conserved& a, const Conserved& b) {
  return {a.rho - b.rho, a.rhou - b.rhou, a.rhov - b.rhov, a.E - b.E};
}
static inline Conserved mul(const Conserved& a, double s) {
  return {a.rho * s, a.rhou * s, a.rhov * s, a.E * s};
}

double minmod(double a, double b) {
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

// ------------------------------
// CFL timestep (basic: dt = cfl*min(dx,dy)/max(|v|+cs))
// ------------------------------
double compute_dt(const Grid& grid, double cfl) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const double h = std::min(grid.dx, grid.dy);
  double amax = 0.0;

  #pragma omp parallel for collapse(2) reduction(max:amax) schedule(static)
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const Conserved& U = grid.U[ grid.idx(i,j) ];
      const Primitive W  = conserved_to_primitive(U);
      const double vmag  = std::sqrt(W.u*W.u + W.v*W.v);
      const double cs    = sound_speed(W);
      amax = std::max(amax, vmag + cs);
    }
  }
  if (amax <= 0.0) return 1e-12;
  return cfl * h / amax;
}

// ------------------------------
// Boundary conditions
// NOTE: With this solver, BCs are enforced WITHOUT ghost cells.
// So this function can be empty or kept for other parts.
// We'll keep it as a no-op to avoid negative index usage.
// ------------------------------
void apply_boundary_conditions(Grid&) {
  // no-op: BC handled in flux assembly below
}

// ------------------------------
// Optional RHS not used here
// ------------------------------
void compute_rhs(const Grid& grid, std::vector<Conserved>& rhs) {
  rhs.assign(grid.U.size(), Conserved{0,0,0,0});
}

// ------------------------------
// MUSCL-Hancock (Toro-style) without ghost / negative indexing
// ------------------------------
namespace {
struct Workspace {
  int nx = -1, ny = -1;
  std::size_t nU = 0;

  std::vector<Conserved> Umx, Upx, Umy, Upy;
  std::vector<Conserved> Uhat_mx, Uhat_px, Uhat_my, Uhat_py;

  // Interface fluxes on domain boundaries included:
  // Fx: (nx+1)*ny corresponds to i=0..nx (left boundary to right boundary)
  // Gy: nx*(ny+1) corresponds to j=0..ny
  std::vector<Conserved> Fx, Gy;

  void ensure(const Grid& grid) {
    if (grid.nx == nx && grid.ny == ny && grid.U.size() == nU) return;
    nx = grid.nx; ny = grid.ny; nU = grid.U.size();
    Umx.resize(nU); Upx.resize(nU); Umy.resize(nU); Upy.resize(nU);
    Uhat_mx.resize(nU); Uhat_px.resize(nU); Uhat_my.resize(nU); Uhat_py.resize(nU);
    Fx.resize((nx+1) * ny);
    Gy.resize(nx * (ny+1));
  }
};

static Workspace ws;

static inline std::size_t fx_id(int nx, int iface_i, int j) {
  return static_cast<std::size_t>(iface_i + (nx+1)*j);
}
static inline std::size_t gy_id(int nx, int i, int iface_j) {
  return static_cast<std::size_t>(i + nx*iface_j);
}

// Safe accessors for neighbor cell indices WITHOUT negative indexing.
// For transmissive BC, out-of-range neighbor uses nearest in-range cell.
static inline int clamp_i(int i, int nx) { return std::max(0, std::min(nx-1, i)); }
static inline int clamp_j(int j, int ny) { return std::max(0, std::min(ny-1, j)); }

} // namespace

void advance_one_step(Grid& grid, double dt) {
  const int nx = grid.nx;
  const int ny = grid.ny;

  ws.ensure(grid);

  const double half = 0.5;
  const double cx   = dt / (2.0 * grid.dx);
  const double cy   = dt / (2.0 * grid.dy);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // ------------------------------------------------------------
    // (I) Reconstruction + (II) Hancock predictor
    // ------------------------------------------------------------
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const auto id = grid.idx(i, j);
        const Conserved& Uc = grid.U[id];

        const int il = clamp_i(i - 1, nx);
        const int ir = clamp_i(i + 1, nx);
        const int jb = clamp_j(j - 1, ny);
        const int jt = clamp_j(j + 1, ny);

        const Conserved& Ul = grid.U[ grid.idx(il, j) ];
        const Conserved& Ur = grid.U[ grid.idx(ir, j) ];
        const Conserved& Ub = grid.U[ grid.idx(i, jb) ];
        const Conserved& Ut = grid.U[ grid.idx(i, jt) ];

        const Conserved sx = minmod_vec(sub(Uc, Ul), sub(Ur, Uc));
        const Conserved sy = minmod_vec(sub(Uc, Ub), sub(Ut, Uc));

        const Conserved Umx = sub(Uc, mul(sx, half));
        const Conserved Upx = add(Uc, mul(sx, half));
        const Conserved Umy = sub(Uc, mul(sy, half));
        const Conserved Upy = add(Uc, mul(sy, half));

        ws.Umx[id] = Umx; ws.Upx[id] = Upx;
        ws.Umy[id] = Umy; ws.Upy[id] = Upy;

        const Conserved dF   = sub(flux_x(Umx), flux_x(Upx));
        const Conserved dG   = sub(flux_y(Umy), flux_y(Upy));
        const Conserved corr = add(mul(dF, cx), mul(dG, cy));

        ws.Uhat_mx[id] = add(Umx, corr);
        ws.Uhat_px[id] = add(Upx, corr);
        ws.Uhat_my[id] = add(Umy, corr);
        ws.Uhat_py[id] = add(Upy, corr);
      }
    }

    // ------------------------------------------------------------
    // (IIIa) x-interfaces HLL fluxes: iface = 0..nx, j=0..ny-1
    // ------------------------------------------------------------
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int j = 0; j < ny; ++j) {
      for (int iface = 0; iface <= nx; ++iface) {
        const Conserved *ULp, *URp;

        if (iface == 0) {
          const auto id0 = grid.idx(0, j);
          ULp = &ws.Uhat_px[id0];
          URp = &ws.Uhat_mx[id0];
        } else if (iface == nx) {
          const auto idN = grid.idx(nx - 1, j);
          ULp = &ws.Uhat_px[idN];
          URp = &ws.Uhat_mx[idN];
        } else {
          const auto idL = grid.idx(iface - 1, j);
          const auto idR = grid.idx(iface,     j);
          ULp = &ws.Uhat_px[idL];
          URp = &ws.Uhat_mx[idR];
        }

        ws.Fx[ fx_id(nx, iface, j) ] = hll_flux_x(*ULp, *URp);
      }
    }

    // ------------------------------------------------------------
    // (IIIb) y-interfaces HLL fluxes: iface = 0..ny, i=0..nx-1
    // ------------------------------------------------------------
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int iface = 0; iface <= ny; ++iface) {
      for (int i = 0; i < nx; ++i) {
        const Conserved *ULp, *URp;

        if (iface == 0) {
          const auto id0 = grid.idx(i, 0);
          ULp = &ws.Uhat_py[id0];
          URp = &ws.Uhat_my[id0];
        } else if (iface == ny) {
          const auto idN = grid.idx(i, ny - 1);
          ULp = &ws.Uhat_py[idN];
          URp = &ws.Uhat_my[idN];
        } else {
          const auto idB = grid.idx(i, iface - 1);
          const auto idT = grid.idx(i, iface);
          ULp = &ws.Uhat_py[idB];
          URp = &ws.Uhat_my[idT];
        }

        ws.Gy[ gy_id(nx, i, iface) ] = hll_flux_y(*ULp, *URp);
      }
    }

    // ------------------------------------------------------------
    // (IV) FV update (in-place)
    // ------------------------------------------------------------
    const double dtdx = dt / grid.dx;
    const double dtdy = dt / grid.dy;

#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const Conserved& F_im = ws.Fx[ fx_id(nx, i,   j) ];
        const Conserved& F_ip = ws.Fx[ fx_id(nx, i+1, j) ];

        const Conserved& G_jm = ws.Gy[ gy_id(nx, i, j) ];
        const Conserved& G_jp = ws.Gy[ gy_id(nx, i, j+1) ];

        const Conserved update =
          add(mul(sub(F_ip, F_im), dtdx),
              mul(sub(G_jp, G_jm), dtdy));

        const auto id = grid.idx(i, j);
        grid.U[id] = sub(grid.U[id], update);
      }
    }
  } // end parallel region
}
