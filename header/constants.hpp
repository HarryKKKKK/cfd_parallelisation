#pragma once

namespace phys {
constexpr double gamma = 1.4;

// Initial conditions
constexpr double rho_air    = 1.29;
constexpr double rho_helium = 0.214;
constexpr double p0         = 1.01325e5;
constexpr double mach       = 1.22;

// Shock-bubble
constexpr double bubble_r = 0.025;
constexpr double bubble_x = 0.035;
constexpr double bubble_y = 0.0445;

constexpr double shock_x  = 0.005;

}
