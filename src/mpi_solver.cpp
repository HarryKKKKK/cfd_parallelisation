#include "mpi_solver.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

struct Decoded {
  double rho, u, v, p, a;
};

static inline Decoded decode_state(const Conserved& U) {
  Decoded d;
  d.rho = U.rho;

#ifndef NDEBUG
  if (d.rho <= 0.0) {
    throw std::runtime_error("Non-positive density encountered.");
  }
#endif

  const double inv_rho = 1.0 / d.rho;
  d.u = U.rhou * inv_rho;
  d.v = U.rhov * inv_rho;

  const double kinetic = 0.5 * d.rho * (d.u * d.u + d.v * d.v);
  d.p = (phys::gamma - 1.0) * (U.E - kinetic);

#ifndef NDEBUG
  if (d.p <= 0.0) {
    throw std::runtime_error("Non-positive pressure encountered.");
  }
#endif

  d.a = std::sqrt(phys::gamma * d.p * inv_rho);
  return d;
}

static inline Primitive cons_to_prim_local(const Conserved& U) {
  const auto d = decode_state(U);
  Primitive W;
  W.rho = d.rho;
  W.u   = d.u;
  W.v   = d.v;
  W.p   = d.p;
  return W;
}

static inline Conserved flux_x_decoded(const Conserved& U, double u, double p) {
  Conserved F;
  F.rho  = U.rhou;
  F.rhou = U.rhou * u + p;
  F.rhov = U.rhov * u;
  F.E    = (U.E + p) * u;
  return F;
}

static inline Conserved flux_y_decoded(const Conserved& U, double v, double p) {
  Conserved G;
  G.rho  = U.rhov;
  G.rhou = U.rhou * v;
  G.rhov = U.rhov * v + p;
  G.E    = (U.E + p) * v;
  return G;
}

static inline Conserved hllc_flux_x(
    const Conserved& UL, const Primitive& WL,
    const Conserved& UR, const Primitive& WR) {

  const double aL = std::sqrt(phys::gamma * WL.p / WL.rho);
  const double aR = std::sqrt(phys::gamma * WR.p / WR.rho);

  const double SL = std::min(WL.u - aL, WR.u - aR);
  const double SR = std::max(WL.u + aL, WR.u + aR);

  const Conserved FL = flux_x_decoded(UL, WL.u, WL.p);
  const Conserved FR = flux_x_decoded(UR, WR.u, WR.p);

  if (SL >= 0.0) return FL;
  if (SR <= 0.0) return FR;

  const double denom =
      WL.rho * (SL - WL.u) - WR.rho * (SR - WR.u);

  const double eps = 1e-14;
  const double inv_denom =
      1.0 / (std::abs(denom) < eps ? (denom >= 0.0 ? eps : -eps) : denom);

  const double Sstar =
      (WR.p - WL.p
       + WL.rho * WL.u * (SL - WL.u)
       - WR.rho * WR.u * (SR - WR.u)) * inv_denom;

  const double rhoL_star =
      WL.rho * (SL - WL.u) / (SL - Sstar);

  Conserved UL_star;
  UL_star.rho  = rhoL_star;
  UL_star.rhou = rhoL_star * Sstar;
  UL_star.rhov = rhoL_star * WL.v;
  UL_star.E =
      rhoL_star *
      ( UL.E / WL.rho
        + (Sstar - WL.u) *
          (Sstar + WL.p / (WL.rho * (SL - WL.u))) );

  const double rhoR_star =
      WR.rho * (SR - WR.u) / (SR - Sstar);

  Conserved UR_star;
  UR_star.rho  = rhoR_star;
  UR_star.rhou = rhoR_star * Sstar;
  UR_star.rhov = rhoR_star * WR.v;
  UR_star.E =
      rhoR_star *
      ( UR.E / WR.rho
        + (Sstar - WR.u) *
          (Sstar + WR.p / (WR.rho * (SR - WR.u))) );

  if (Sstar >= 0.0) {
    return FL + (UL_star - UL) * SL;
  } else {
    return FR + (UR_star - UR) * SR;
  }
}

static inline Conserved hllc_flux_y(
    const Conserved& UL, const Primitive& WL,
    const Conserved& UR, const Primitive& WR) {

  const double aL = std::sqrt(phys::gamma * WL.p / WL.rho);
  const double aR = std::sqrt(phys::gamma * WR.p / WR.rho);

  const double SL = std::min(WL.v - aL, WR.v - aR);
  const double SR = std::max(WL.v + aL, WR.v + aR);

  const Conserved GL = flux_y_decoded(UL, WL.v, WL.p);
  const Conserved GR = flux_y_decoded(UR, WR.v, WR.p);

  if (SL >= 0.0) return GL;
  if (SR <= 0.0) return GR;

  const double denom =
      WL.rho * (SL - WL.v) - WR.rho * (SR - WR.v);

  const double eps = 1e-14;
  const double inv_denom =
      1.0 / (std::abs(denom) < eps ? (denom >= 0.0 ? eps : -eps) : denom);

  const double Sstar =
      (WR.p - WL.p
       + WL.rho * WL.v * (SL - WL.v)
       - WR.rho * WR.v * (SR - WR.v)) * inv_denom;

  const double rhoL_star =
      WL.rho * (SL - WL.v) / (SL - Sstar);

  Conserved UL_star;
  UL_star.rho  = rhoL_star;
  UL_star.rhou = rhoL_star * WL.u;
  UL_star.rhov = rhoL_star * Sstar;
  UL_star.E =
      rhoL_star *
      ( UL.E / WL.rho
        + (Sstar - WL.v) *
          (Sstar + WL.p / (WL.rho * (SL - WL.v))) );

  const double rhoR_star =
      WR.rho * (SR - WR.v) / (SR - Sstar);

  Conserved UR_star;
  UR_star.rho  = rhoR_star;
  UR_star.rhou = rhoR_star * WR.u;
  UR_star.rhov = rhoR_star * Sstar;
  UR_star.E =
      rhoR_star *
      ( UR.E / WR.rho
        + (Sstar - WR.v) *
          (Sstar + WR.p / (WR.rho * (SR - WR.v))) );

  if (Sstar >= 0.0) {
    return GL + (UL_star - UL) * SL;
  } else {
    return GR + (UR_star - UR) * SR;
  }
}

// ============================================================
// Row workspace
// ============================================================

struct RowData {
  std::vector<Conserved> Lc, Rc, Dc, Uc;
  std::vector<Primitive> Lp, Rp, Dp, Up;

  void resize(int nx) {
    Lc.resize(nx); Rc.resize(nx); Dc.resize(nx); Uc.resize(nx);
    Lp.resize(nx); Rp.resize(nx); Dp.resize(nx); Up.resize(nx);
  }
};

static inline int gid(const Grid& g, int I, int J) {
  return g.idx(I, J);
}

static inline void build_row_predictor(
    const Grid& g,
    int J,
    double cx,
    double cy,
    RowData& row) {

  const int nx = g.nx;
  const int ng = g.ng;

  row.resize(nx);

  for (int i = 0; i < nx; ++i) {
    const int I = i + ng;

    const Conserved& Uc = g.U[gid(g, I,     J)];
    const Conserved& Ul = g.U[gid(g, I - 1, J)];
    const Conserved& Ur = g.U[gid(g, I + 1, J)];
    const Conserved& Ub = g.U[gid(g, I, J - 1)];
    const Conserved& Ut = g.U[gid(g, I, J + 1)];

    const Conserved sx = {
      minmod(Uc.rho  - Ul.rho,  Ur.rho  - Uc.rho),
      minmod(Uc.rhou - Ul.rhou, Ur.rhou - Uc.rhou),
      minmod(Uc.rhov - Ul.rhov, Ur.rhov - Uc.rhov),
      minmod(Uc.E    - Ul.E,    Ur.E    - Uc.E)
    };

    const Conserved sy = {
      minmod(Uc.rho  - Ub.rho,  Ut.rho  - Uc.rho),
      minmod(Uc.rhou - Ub.rhou, Ut.rhou - Uc.rhou),
      minmod(Uc.rhov - Ub.rhov, Ut.rhov - Uc.rhov),
      minmod(Uc.E    - Ub.E,    Ut.E    - Uc.E)
    };

    const Conserved L = Uc - sx * 0.5;
    const Conserved R = Uc + sx * 0.5;
    const Conserved D = Uc - sy * 0.5;
    const Conserved U = Uc + sy * 0.5;

    const auto dL = decode_state(L);
    const auto dR = decode_state(R);
    const auto dD = decode_state(D);
    const auto dU = decode_state(U);

    const Conserved FxL = flux_x_decoded(L, dL.u, dL.p);
    const Conserved FxR = flux_x_decoded(R, dR.u, dR.p);
    const Conserved GyD = flux_y_decoded(D, dD.v, dD.p);
    const Conserved GyU = flux_y_decoded(U, dU.v, dU.p);

    const Conserved corr = (FxL - FxR) * cx + (GyD - GyU) * cy;

    row.Lc[i] = L + corr;
    row.Rc[i] = R + corr;
    row.Dc[i] = D + corr;
    row.Uc[i] = U + corr;

    row.Lp[i] = cons_to_prim_local(row.Lc[i]);
    row.Rp[i] = cons_to_prim_local(row.Rc[i]);
    row.Dp[i] = cons_to_prim_local(row.Dc[i]);
    row.Up[i] = cons_to_prim_local(row.Uc[i]);
  }
}

static inline void compute_flux_x_row(
    const RowData& row,
    std::vector<Conserved>& Fx) {

  const int nx = static_cast<int>(row.Lc.size());
  Fx.resize(nx + 1);

  for (int iface = 0; iface <= nx; ++iface) {
    const int iL = iface - 1;
    const int iR = iface;

    Fx[iface] = hllc_flux_x(
        row.Rc[iL < 0 ? 0 : iL], row.Rp[iL < 0 ? 0 : iL],
        row.Lc[iR >= nx ? nx - 1 : iR], row.Lp[iR >= nx ? nx - 1 : iR]);
  }
}

static inline void compute_flux_y_between_rows(
    const RowData& rowB,
    const RowData& rowT,
    std::vector<Conserved>& Gy) {

  const int nx = static_cast<int>(rowB.Lc.size());
  Gy.resize(nx);

  for (int i = 0; i < nx; ++i) {
    Gy[i] = hllc_flux_y(
        rowB.Uc[i], rowB.Up[i],
        rowT.Dc[i], rowT.Dp[i]);
  }
}

static MPI_Datatype make_conserved_type() {
  MPI_Datatype type;
  MPI_Type_contiguous(4, MPI_DOUBLE, &type);
  MPI_Type_commit(&type);
  return type;
}

} // namespace

// ============================================================
// Domain decomposition
// ============================================================

MpiDomain make_mpi_domain_y_slab(int nx_global, int ny_global, int ng, MPI_Comm comm) {
  MpiDomain mp;
  mp.comm = comm;
  mp.nx_global = nx_global;
  mp.ny_global = ny_global;
  mp.ng = ng;

  MPI_Comm_rank(comm, &mp.rank);
  MPI_Comm_size(comm, &mp.size);

  const int base = ny_global / mp.size;
  const int rem  = ny_global % mp.size;

  mp.ny_local  = base + (mp.rank < rem ? 1 : 0);
  mp.y0_global = mp.rank * base + std::min(mp.rank, rem);

  mp.nbr_down = (mp.rank == 0) ? MPI_PROC_NULL : mp.rank - 1;
  mp.nbr_up   = (mp.rank == mp.size - 1) ? MPI_PROC_NULL : mp.rank + 1;

  return mp;
}

// ============================================================
// Limiter
// ============================================================

double minmod(double a, double b) {
  if (a * b <= 0.0) return 0.0;
  return (std::abs(a) < std::abs(b)) ? a : b;
}

// ============================================================
// Physical BC only (same transmissive semantics as serial)
// ============================================================

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

  for (int J = 0; J < ny_tot; ++J) {
    for (int gc = 0; gc < ng; ++gc) {
      g.U[g.idx(gc, J)]              = g.U[g.idx(iL, J)];
      g.U[g.idx(nx_tot - 1 - gc, J)] = g.U[g.idx(iR, J)];
    }
  }

  if (mp.nbr_down == MPI_PROC_NULL) {
    for (int I = 0; I < nx_tot; ++I) {
      for (int gr = 0; gr < ng; ++gr) {
        g.U[g.idx(I, gr)] = g.U[g.idx(I, jB)];
      }
    }
  }

  if (mp.nbr_up == MPI_PROC_NULL) {
    for (int I = 0; I < nx_tot; ++I) {
      for (int gr = 0; gr < ng; ++gr) {
        g.U[g.idx(I, ny_tot - 1 - gr)] = g.U[g.idx(I, jT)];
      }
    }
  }
}

// ============================================================
// Halo exchange in y only
// ============================================================

void exchange_halo_y_mpi(Grid& g, const MpiDomain& mp) {
  const int ng = g.ng;
  if (ng <= 0) return;

  const int nx = g.nx;
  const int ny = g.ny;
  const int nx_tot = nx + 2 * ng;
  const int jB = ng;
  const int jT = ng + ny - 1;

  MPI_Datatype cons_type = make_conserved_type();

  std::vector<Conserved> send_down(static_cast<std::size_t>(ng) * nx_tot);
  std::vector<Conserved> recv_down(static_cast<std::size_t>(ng) * nx_tot);
  std::vector<Conserved> send_up(static_cast<std::size_t>(ng) * nx_tot);
  std::vector<Conserved> recv_up(static_cast<std::size_t>(ng) * nx_tot);

  for (int gr = 0; gr < ng; ++gr) {
    const int Jsrc = jB + gr;
    for (int I = 0; I < nx_tot; ++I) {
      send_down[static_cast<std::size_t>(gr) * nx_tot + I] = g.U[g.idx(I, Jsrc)];
    }
  }

  for (int gr = 0; gr < ng; ++gr) {
    const int Jsrc = jT - ng + 1 + gr;
    for (int I = 0; I < nx_tot; ++I) {
      send_up[static_cast<std::size_t>(gr) * nx_tot + I] = g.U[g.idx(I, Jsrc)];
    }
  }

  MPI_Request reqs[4];
  int nreq = 0;

  if (mp.nbr_down != MPI_PROC_NULL) {
    MPI_Irecv(recv_down.data(), static_cast<int>(recv_down.size()), cons_type,
              mp.nbr_down, 101, mp.comm, &reqs[nreq++]);
    MPI_Isend(send_down.data(), static_cast<int>(send_down.size()), cons_type,
              mp.nbr_down, 202, mp.comm, &reqs[nreq++]);
  }

  if (mp.nbr_up != MPI_PROC_NULL) {
    MPI_Irecv(recv_up.data(), static_cast<int>(recv_up.size()), cons_type,
              mp.nbr_up, 202, mp.comm, &reqs[nreq++]);
    MPI_Isend(send_up.data(), static_cast<int>(send_up.size()), cons_type,
              mp.nbr_up, 101, mp.comm, &reqs[nreq++]);
  }

  if (nreq > 0) {
    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
  }

  if (mp.nbr_down != MPI_PROC_NULL) {
    for (int gr = 0; gr < ng; ++gr) {
      for (int I = 0; I < nx_tot; ++I) {
        g.U[g.idx(I, gr)] = recv_down[static_cast<std::size_t>(gr) * nx_tot + I];
      }
    }
  }

  if (mp.nbr_up != MPI_PROC_NULL) {
    for (int gr = 0; gr < ng; ++gr) {
      for (int I = 0; I < nx_tot; ++I) {
        g.U[g.idx(I, jT + 1 + gr)] = recv_up[static_cast<std::size_t>(gr) * nx_tot + I];
      }
    }
  }

  MPI_Type_free(&cons_type);
}

// ============================================================
// CFL dt: same formula as serial, with global MPI reduction
// ============================================================

double compute_dt_mpi(const Grid& grid, const MpiDomain& mp, double cfl) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  const double h = std::min(grid.dx, grid.dy);
  double amax_local = 0.0;

  for (int j = ng; j < ny + ng; ++j) {
    for (int i = ng; i < nx + ng; ++i) {
      const Conserved& U = grid.U[grid.idx(i, j)];
      const auto d = decode_state(U);

      const double vmag = std::sqrt(d.u * d.u + d.v * d.v);
      amax_local = std::max(amax_local, vmag + d.a);
    }
  }

  double amax_global = 0.0;
  MPI_Allreduce(&amax_local, &amax_global, 1, MPI_DOUBLE, MPI_MAX, mp.comm);

  if (amax_global <= 0.0) return 1e-12;
  return cfl * h / amax_global;
}

// ============================================================
// One-step advance: same algorithmic kernel as serial
// ============================================================
void advance_one_step_mpi(Grid& grid, const MpiDomain& mp, double dt) {
  const int nx = grid.nx;
  const int ny = grid.ny;
  const int ng = grid.ng;

  apply_boundary_conditions_mpi(grid, mp);
  exchange_halo_y_mpi(grid, mp);
  apply_boundary_conditions_mpi(grid, mp); 

  const double cx = dt / (2.0 * grid.dx);
  const double cy = dt / (2.0 * grid.dy);
  const double dtdx = dt / grid.dx;
  const double dtdy = dt / grid.dy;

  RowData row_B, row_C, row_T;
  std::vector<Conserved> Fx, Gy_B, Gy_T;

  for (int j = 0; j < ny; ++j) {
    const int J = j + ng;

    build_row_predictor(grid, J - 1, cx, cy, row_B);
    build_row_predictor(grid, J,     cx, cy, row_C);
    build_row_predictor(grid, J + 1, cx, cy, row_T);

    compute_flux_x_row(row_C, Fx);
    compute_flux_y_between_rows(row_B, row_C, Gy_B);
    compute_flux_y_between_rows(row_C, row_T, Gy_T);

    for (int i = 0; i < nx; ++i) {
      const Conserved update =
          (Fx[i + 1] - Fx[i]) * dtdx
        + (Gy_T[i]    - Gy_B[i]) * dtdy;

      const int I = i + ng;
      const int id = grid.idx(I, J);
      
      grid.U_new[id] = grid.U[id] - update; 
    }
  }

  std::swap(grid.U, grid.U_new);
}