#pragma once
#include "types.hpp"
#include "constants.hpp"

#include <iostream>
#include <string>
#include <cmath>

#ifdef DEBUG
    #define DBG_PRINT(x) \
        do { std::cerr << "[DBG] " << __FILE__ << ":" << __LINE__ \
                       << " | " << x << std::endl; } while (0)
#else
    #define DBG_PRINT(x) do {} while (0)
#endif

void write_grid_csv(const Grid& grid, const std::string& filename);
std::string make_filename(const std::string& dir, int step, double t);
double check_symmetry(const Grid& grid, double t, bool verbose);

/// Ambient sound speed:
/// s = sqrt(gamma * p0 / rho_air)
inline double get_sound_speed()
{
    return std::sqrt(
        phys::gamma * phys::p0 / phys::rho_air
    );
}

/// Shock speed:
/// U_s = M_s * s
inline double get_shock_speed()
{
    return phys::mach * get_sound_speed();
}

/// Bubble left boundary:
/// x_bubble_left = bubble_x - bubble_r
inline double get_bubble_left()
{
    return phys::bubble_x - phys::bubble_r;
}

/// Collision time:
///
/// T_collision
/// = (x_bubble_left - shock_x) / U_s
///
/// = ( (bubble_x - bubble_r) - shock_x )
///   / ( mach * sqrt(gamma * p0 / rho_air) )
///
inline double get_collision_time()
{
    const double Us = get_shock_speed();
    return (get_bubble_left() - phys::shock_x) / Us;
}

/// Reference time (used in paper):
///
/// T_ref = r / U_s
///
/// = bubble_r / (mach * sqrt(gamma * p0 / rho_air))
///
inline double get_time_ref()
{
    return phys::bubble_r / get_shock_speed();
}