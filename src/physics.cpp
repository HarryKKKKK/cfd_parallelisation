#include "physics.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// A small POD to hold decoded state (no allocations)
struct Decoded {
    double rho, u, v, p, a;
};

// Decode from Conserved once: rho,u,v,p,a
static inline Decoded decode_state(const Conserved& U) {
    Decoded d;
    d.rho = U.rho;

    // If you want max speed, avoid throws inside hot loops.
    // For debugging, keep them behind a macro.
#ifndef NDEBUG
    if (d.rho <= 0.0) throw std::runtime_error("Non-positive density encountered.");
#endif

    const double inv_rho = 1.0 / d.rho;
    d.u = U.rhou * inv_rho;
    d.v = U.rhov * inv_rho;

    const double kinetic = 0.5 * d.rho * (d.u*d.u + d.v*d.v);
    d.p = (phys::gamma - 1.0) * (U.E - kinetic);

#ifndef NDEBUG
    if (d.p <= 0.0) throw std::runtime_error("Non-positive pressure encountered.");
#endif

    d.a = std::sqrt(phys::gamma * d.p * inv_rho);
    return d;
}

// Fluxes using decoded values (no repeated primitive conversion)
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


double pressure_from_conserved(const Conserved& U) {
    // p = (gamma - 1) * (E - 0.5*rho*(u^2+v^2))
    const double rho = U.rho;
    if (rho <= 0.0) throw std::runtime_error("Non-positive density encountered.");

    const double u = U.rhou / rho;
    const double v = U.rhov / rho;

    const double kinetic = 0.5 * rho * (u*u + v*v);
    const double p = (phys::gamma - 1.0) * (U.E - kinetic);

    if (p <= 0.0) throw std::runtime_error("Non-positive pressure encountered.");
    return p;
}

double sound_speed_from_conserved(const Conserved& U) {
    const double p = pressure_from_conserved(U);
    const double rho = U.rho;
    return std::sqrt(phys::gamma * p / rho);
}

Primitive conserved_to_primitive(const Conserved& U) {
    Primitive W;
    W.rho = U.rho;

    if (W.rho <= 0.0) throw std::runtime_error("Non-positive density in conserved_to_primitive.");

    W.u = U.rhou / W.rho;
    W.v = U.rhov / W.rho;
    W.p = pressure_from_conserved(U);
    return W;
}

Conserved primitive_to_conserved(const Primitive& W) {
    if (W.rho <= 0.0) throw std::runtime_error("Non-positive density in primitive_to_conserved.");
    if (W.p   <= 0.0) throw std::runtime_error("Non-positive pressure in primitive_to_conserved.");

    Conserved U;
    U.rho  = W.rho;
    U.rhou = W.rho * W.u;
    U.rhov = W.rho * W.v;

    // Total energy: E = p/(gamma-1) + 0.5*rho*(u^2+v^2)
    const double internal = W.p / (phys::gamma - 1.0);
    const double kinetic  = 0.5 * W.rho * (W.u*W.u + W.v*W.v);
    U.E = internal + kinetic;

    return U;
}

Conserved flux_x(const Conserved& U) {
    const auto d = decode_state(U);
    return flux_x_decoded(U, d.u, d.p);
}

Conserved flux_y(const Conserved& U) {
    const auto d = decode_state(U);
    return flux_y_decoded(U, d.v, d.p);
}


double max_wave_speed_x(const Conserved& U) {
    const auto d = decode_state(U);
    return std::abs(d.u) + d.a;
}
double max_wave_speed_y(const Conserved& U) {
    const auto d = decode_state(U);
    return std::abs(d.v) + d.a;
}


Conserved hll_flux_x(const Conserved& UL, const Conserved& UR) {
    // // HLL flux
    // const auto L = decode_state(UL);
    // const auto R = decode_state(UR);

    // const double SL = std::min(L.u - L.a, R.u - R.a);
    // const double SR = std::max(L.u + L.a, R.u + R.a);

    // const Conserved FL = flux_x_decoded(UL, L.u, L.p);
    // const Conserved FR = flux_x_decoded(UR, R.u, R.p);

    // if (SL >= 0.0) return FL;
    // if (SR <= 0.0) return FR;

    // const double inv = 1.0 / (SR - SL);

    // Conserved FH;
    // FH.rho  = (SR*FL.rho  - SL*FR.rho  + SL*SR*(UR.rho  - UL.rho )) * inv;
    // FH.rhou = (SR*FL.rhou - SL*FR.rhou + SL*SR*(UR.rhou - UL.rhou)) * inv;
    // FH.rhov = (SR*FL.rhov - SL*FR.rhov + SL*SR*(UR.rhov - UL.rhov)) * inv;
    // FH.E    = (SR*FL.E    - SL*FR.E    + SL*SR*(UR.E    - UL.E   )) * inv;
    // return FH;

    // HLLC flux
    const auto L = decode_state(UL);
    const auto R = decode_state(UR);

    const double SL = std::min(L.u - L.a, R.u - R.a);
    const double SR = std::max(L.u + L.a, R.u + R.a);

    const Conserved FL = flux_x_decoded(UL, L.u, L.p);
    const Conserved FR = flux_x_decoded(UR, R.u, R.p);

    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    const double denom =
        L.rho * (SL - L.u) - R.rho * (SR - R.u);

    const double eps = 1e-14;
    const double inv_denom =
        1.0 / (std::abs(denom) < eps ? (denom >= 0 ? eps : -eps) : denom);

    const double Sstar =
        (R.p - L.p
         + L.rho * L.u * (SL - L.u)
         - R.rho * R.u * (SR - R.u)) * inv_denom;

    // Left star state
    const double rhoL_star =
        L.rho * (SL - L.u) / (SL - Sstar);

    Conserved UL_star;
    UL_star.rho  = rhoL_star;
    UL_star.rhou = rhoL_star * Sstar;
    UL_star.rhov = rhoL_star * L.v;
    UL_star.E =
        rhoL_star *
        ( UL.E / L.rho
          + (Sstar - L.u) *
            (Sstar + L.p / (L.rho * (SL - L.u))) );

    // Right star state
    const double rhoR_star =
        R.rho * (SR - R.u) / (SR - Sstar);

    Conserved UR_star;
    UR_star.rho  = rhoR_star;
    UR_star.rhou = rhoR_star * Sstar;
    UR_star.rhov = rhoR_star * R.v;
    UR_star.E =
        rhoR_star *
        ( UR.E / R.rho
          + (Sstar - R.u) *
            (Sstar + R.p / (R.rho * (SR - R.u))) );

    if (Sstar >= 0.0)
        return FL + (UL_star - UL) * SL;
    else
        return FR + (UR_star - UR) * SR;
}

Conserved hll_flux_y(const Conserved& UL, const Conserved& UR) {
    // // HLL flux
    // const auto L = decode_state(UL);
    // const auto R = decode_state(UR);

    // const double SL = std::min(L.v - L.a, R.v - R.a);
    // const double SR = std::max(L.v + L.a, R.v + R.a);

    // const Conserved GL = flux_y_decoded(UL, L.v, L.p);
    // const Conserved GR = flux_y_decoded(UR, R.v, R.p);

    // if (SL >= 0.0) return GL;
    // if (SR <= 0.0) return GR;

    // const double inv = 1.0 / (SR - SL);

    // Conserved GH;
    // GH.rho  = (SR*GL.rho  - SL*GR.rho  + SL*SR*(UR.rho  - UL.rho )) * inv;
    // GH.rhou = (SR*GL.rhou - SL*GR.rhou + SL*SR*(UR.rhou - UL.rhou)) * inv;
    // GH.rhov = (SR*GL.rhov - SL*GR.rhov + SL*SR*(UR.rhov - UL.rhov)) * inv;
    // GH.E    = (SR*GL.E    - SL*GR.E    + SL*SR*(UR.E    - UL.E   )) * inv;
    // return GH;

    // HLLC flux
    const auto L = decode_state(UL);
    const auto R = decode_state(UR);

    const double SL = std::min(L.v - L.a, R.v - R.a);
    const double SR = std::max(L.v + L.a, R.v + R.a);

    const Conserved GL = flux_y_decoded(UL, L.v, L.p);
    const Conserved GR = flux_y_decoded(UR, R.v, R.p);

    if (SL >= 0.0) return GL;
    if (SR <= 0.0) return GR;

    const double denom =
        L.rho * (SL - L.v) - R.rho * (SR - R.v);

    const double eps = 1e-14;
    const double inv_denom =
        1.0 / (std::abs(denom) < eps ? (denom >= 0 ? eps : -eps) : denom);

    const double Sstar =
        (R.p - L.p
         + L.rho * L.v * (SL - L.v)
         - R.rho * R.v * (SR - R.v)) * inv_denom;

    const double rhoL_star =
        L.rho * (SL - L.v) / (SL - Sstar);

    Conserved UL_star;
    UL_star.rho  = rhoL_star;
    UL_star.rhou = rhoL_star * L.u;
    UL_star.rhov = rhoL_star * Sstar;
    UL_star.E =
        rhoL_star *
        ( UL.E / L.rho
          + (Sstar - L.v) *
            (Sstar + L.p / (L.rho * (SL - L.v))) );

    const double rhoR_star =
        R.rho * (SR - R.v) / (SR - Sstar);

    Conserved UR_star;
    UR_star.rho  = rhoR_star;
    UR_star.rhou = rhoR_star * R.u;
    UR_star.rhov = rhoR_star * Sstar;
    UR_star.E =
        rhoR_star *
        ( UR.E / R.rho
          + (Sstar - R.v) *
            (Sstar + R.p / (R.rho * (SR - R.v))) );

    if (Sstar >= 0.0)
        return GL + (UL_star - UL) * SL;
    else
        return GR + (UR_star - UR) * SR;
}

double sound_speed(const Primitive& W)
{
    const double rho = std::max(W.rho, 1e-14);
    const double p   = std::max(W.p,   1e-14);

    return std::sqrt(phys::gamma * p / rho);
}
