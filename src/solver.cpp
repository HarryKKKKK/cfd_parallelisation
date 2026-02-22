// solver.cpp (ghost-cell MUSCL-Hancock, transmissive BC, HLL flux)
#include "solver.hpp"
#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
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
// Transmissive (zero-gradient) ghost fill
// For ng >= 1. We will REQUIRE ng >= 2 in advance_one_step()
// because we compute predictor on (physical + 1 ghost layer) and need (I±1,J±1).
// ------------------------------
void apply_boundary_conditions(Grid& g) {
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

  // Left / Right ghost columns: copy nearest interior column
  for (int J = 0; J < ny_tot; ++J) {
    for (int gc = 0; gc < ng; ++gc) {
      g.U[g.idx(gc, J)]              = g.U[g.idx(iL, J)];
      g.U[g.idx(nx_tot - 1 - gc, J)] = g.U[g.idx(iR, J)];
    }
  }

  // Bottom / Top ghost rows: copy nearest interior row
  for (int I = 0; I < nx_tot; ++I) {
    for (int gr = 0; gr < ng; ++gr) {
      g.U[g.idx(I, gr)]              = g.U[g.idx(I, jB)];
      g.U[g.idx(I, ny_tot - 1 - gr)] = g.U[g.idx(I, jT)];
    }
  }
}

// ------------------------------
// CFL timestep (scan ONLY physical cells)
// dt = cfl*min(dx,dy)/max(|v|+cs)
// ------------------------------
double compute_dt(const Grid& grid, double cfl) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  const double h = std::min(grid.dx, grid.dy);
  double amax = 0.0;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(max:amax) schedule(static)
#endif
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const int I = i + ng;
      const int J = j + ng;

      const Conserved& U = grid.U[ grid.idx(I, J) ];
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
// Optional RHS not used here
// ------------------------------
void compute_rhs(const Grid& grid, std::vector<Conserved>& rhs) {
  rhs.assign(grid.U.size(), Conserved{0,0,0,0});
}

// ------------------------------
// MUSCL-Hancock with ghost cells + HLL flux
// - predictor computed on (physical + 1 ghost layer)
// - flux computed on all physical interfaces
// - update only physical cells
// REQUIRE ng >= 2
// ------------------------------
namespace {

struct Workspace {
  int nx = -1, ny = -1, ng = -1;
  std::size_t nU = 0;

  std::vector<Conserved> Umx, Upx, Umy, Upy;
  std::vector<Conserved> Uhat_mx, Uhat_px, Uhat_my, Uhat_py;

  // Interface fluxes over PHYSICAL interfaces:
  // Fx: (nx+1)*ny for x-faces iface=0..nx, j=0..ny-1
  // Gy: nx*(ny+1) for y-faces i=0..nx-1, iface=0..ny
  std::vector<Conserved> Fx, Gy;

  void ensure(const Grid& grid) {
    if (grid.nx == nx && grid.ny == ny && grid.ng == ng && grid.U.size() == nU) return;
    nx = grid.nx; ny = grid.ny; ng = grid.ng; nU = grid.U.size();

    Umx.resize(nU); Upx.resize(nU); Umy.resize(nU); Upy.resize(nU);
    Uhat_mx.resize(nU); Uhat_px.resize(nU); Uhat_my.resize(nU); Uhat_py.resize(nU);

    Fx.resize((nx + 1) * ny);
    Gy.resize(nx * (ny + 1));
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

void advance_one_step(Grid& grid, double dt) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  if (ng < 2) {
    throw std::runtime_error("advance_one_step requires ng >= 2 for ghost-cell MUSCL-Hancock predictor.");
  }

  // 1) fill ghost cells (transmissive)
  apply_boundary_conditions(grid);

  ws.ensure(grid);

  const double cx = dt / (2.0 * grid.dx);
  const double cy = dt / (2.0 * grid.dy);

  // We compute predictor states on: I = ng-1 .. ng+nx, J = ng-1 .. ng+ny
  // This covers one ghost layer around the physical domain so interface fluxes
  // at iface=0 or iface=nx (and iface=0/ny in y) can read ws.Uhat_* safely.
  const int Imin = ng - 1;
  const int Imax = ng + nx;     // inclusive
  const int Jmin = ng - 1;
  const int Jmax = ng + ny;     // inclusive

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // ------------------------------------------------------------
    // (I) Reconstruction + (II) Hancock predictor
    // on (physical + 1 ghost layer)
    // ------------------------------------------------------------
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int J = Jmin; J <= Jmax; ++J) {
      for (int I = Imin; I <= Imax; ++I) {
        const auto id = grid.idx(I, J);
        const Conserved& Uc = grid.U[id];

        const Conserved& Ul = grid.U[ grid.idx(I - 1, J) ];
        const Conserved& Ur = grid.U[ grid.idx(I + 1, J) ];
        const Conserved& Ub = grid.U[ grid.idx(I, J - 1) ];
        const Conserved& Ut = grid.U[ grid.idx(I, J + 1) ];

        const Conserved sx = minmod_vec(sub(Uc, Ul), sub(Ur, Uc));
        const Conserved sy = minmod_vec(sub(Uc, Ub), sub(Ut, Uc));

        const Conserved Umx = sub(Uc, mul(sx, 0.5));
        const Conserved Upx = add(Uc, mul(sx, 0.5));
        const Conserved Umy = sub(Uc, mul(sy, 0.5));
        const Conserved Upy = add(Uc, mul(sy, 0.5));

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
    // (IIIa) x-interfaces HLL fluxes: iface = 0..nx, j = 0..ny-1
    // Use ghost cells so NO special casing needed.
    // Interface between left cell (I=ng+iface-1) and right cell (I=ng+iface)
    // ------------------------------------------------------------
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int j = 0; j < ny; ++j) {
      const int J = j + ng;
      for (int iface = 0; iface <= nx; ++iface) {
        const int IL = ng + iface - 1;
        const int IR = ng + iface;

        const auto idL = grid.idx(IL, J);
        const auto idR = grid.idx(IR, J);

        const Conserved& UL = ws.Uhat_px[idL]; // +x face of left cell
        const Conserved& UR = ws.Uhat_mx[idR]; // -x face of right cell

        ws.Fx[ fx_id(nx, iface, j) ] = hll_flux_x(UL, UR);
      }
    }

    // ------------------------------------------------------------
    // (IIIb) y-interfaces HLL fluxes: iface = 0..ny, i = 0..nx-1
    // Interface between bottom cell (J=ng+iface-1) and top cell (J=ng+iface)
    // ------------------------------------------------------------
#ifdef _OPENMP
#pragma omp for collapse(2) schedule(static)
#endif
    for (int iface = 0; iface <= ny; ++iface) {
      const int JB = ng + iface - 1;
      const int JT = ng + iface;

      for (int i = 0; i < nx; ++i) {
        const int I = i + ng;

        const auto idB = grid.idx(I, JB);
        const auto idT = grid.idx(I, JT);

        const Conserved& UL = ws.Uhat_py[idB]; // +y face of bottom cell
        const Conserved& UR = ws.Uhat_my[idT]; // -y face of top cell

        ws.Gy[ gy_id(nx, i, iface) ] = hll_flux_y(UL, UR);
      }
    }

    // ------------------------------------------------------------
    // (IV) FV update (in-place) on PHYSICAL cells only
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

        const int I = i + ng;
        const int J = j + ng;
        const auto id = grid.idx(I, J);

        grid.U[id] = sub(grid.U[id], update);
      }
    }
  } // end parallel region
}