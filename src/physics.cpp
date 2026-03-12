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

// double max_wave_speed_x(const Conserved& U) {
//     const auto d = decode_state(U);
//     return std::abs(d.u) + d.a;
// }
// double max_wave_speed_y(const Conserved& U) {
//     const auto d = decode_state(U);
//     return std::abs(d.v) + d.a;
// }

double sound_speed(const Primitive& W)
{
    const double rho = std::max(W.rho, 1e-14);
    const double p   = std::max(W.p,   1e-14);

    return std::sqrt(phys::gamma * p / rho);
}
