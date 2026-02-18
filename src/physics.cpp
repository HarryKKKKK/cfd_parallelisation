#include "physics.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

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
    const Primitive W = conserved_to_primitive(U);

    Conserved F;
    F.rho  = U.rhou;
    F.rhou = U.rhou * W.u + W.p;
    F.rhov = U.rhov * W.u;
    F.E    = (U.E + W.p) * W.u;
    return F;
}

Conserved flux_y(const Conserved& U) {
    const Primitive W = conserved_to_primitive(U);

    Conserved G;
    G.rho  = U.rhov;
    G.rhou = U.rhou * W.v;
    G.rhov = U.rhov * W.v + W.p;
    G.E    = (U.E + W.p) * W.v;
    return G;
}

double max_wave_speed_x(const Conserved& U) {
    const Primitive W = conserved_to_primitive(U);
    const double a = std::sqrt(phys::gamma * W.p / W.rho);
    return std::abs(W.u) + a;
}

double max_wave_speed_y(const Conserved& U) {
    const Primitive W = conserved_to_primitive(U);
    const double a = std::sqrt(phys::gamma * W.p / W.rho);
    return std::abs(W.v) + a;
}

Conserved hll_flux_x(const Conserved& UL, const Conserved& UR) {
    const Primitive WL = conserved_to_primitive(UL);
    const Primitive WR = conserved_to_primitive(UR);

    const double aL = std::sqrt(phys::gamma * WL.p / WL.rho);
    const double aR = std::sqrt(phys::gamma * WR.p / WR.rho);

    const double SL = std::min(WL.u - aL, WR.u - aR);
    const double SR = std::max(WL.u + aL, WR.u + aR);

    const Conserved FL = flux_x(UL);
    const Conserved FR = flux_x(UR);

    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    // (SR*FL - SL*FR + SL*SR*(UR-UL)) / (SR-SL)
    const double inv = 1.0 / (SR - SL);

    Conserved FH;
    FH.rho  = (SR*FL.rho  - SL*FR.rho  + SL*SR*(UR.rho  - UL.rho )) * inv;
    FH.rhou = (SR*FL.rhou - SL*FR.rhou + SL*SR*(UR.rhou - UL.rhou)) * inv;
    FH.rhov = (SR*FL.rhov - SL*FR.rhov + SL*SR*(UR.rhov - UL.rhov)) * inv;
    FH.E    = (SR*FL.E    - SL*FR.E    + SL*SR*(UR.E    - UL.E   )) * inv;

    return FH;
}

Conserved hll_flux_y(const Conserved& UL, const Conserved& UR) {
    const Primitive WL = conserved_to_primitive(UL);
    const Primitive WR = conserved_to_primitive(UR);

    const double aL = std::sqrt(phys::gamma * WL.p / WL.rho);
    const double aR = std::sqrt(phys::gamma * WR.p / WR.rho);

    const double SL = std::min(WL.v - aL, WR.v - aR);
    const double SR = std::max(WL.v + aL, WR.v + aR);

    const Conserved GL = flux_y(UL);
    const Conserved GR = flux_y(UR);

    if (SL >= 0.0) return GL;
    if (SR <= 0.0) return GR;

    const double inv = 1.0 / (SR - SL);

    Conserved GH;
    GH.rho  = (SR*GL.rho  - SL*GR.rho  + SL*SR*(UR.rho  - UL.rho )) * inv;
    GH.rhou = (SR*GL.rhou - SL*GR.rhou + SL*SR*(UR.rhou - UL.rhou)) * inv;
    GH.rhov = (SR*GL.rhov - SL*GR.rhov + SL*SR*(UR.rhov - UL.rhov)) * inv;
    GH.E    = (SR*GL.E    - SL*GR.E    + SL*SR*(UR.E    - UL.E   )) * inv;

    return GH;
}
