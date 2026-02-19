#pragma once
#include "types.hpp"
#include <vector>

// CFL timestep
double compute_dt(const Grid& grid, double cfl);

// Boundary conditions (transmissive)
void apply_boundary_conditions(Grid& grid);

// Limiter
double minmod(double a, double b);

// Time advance (RK2)
void advance_one_step(Grid& grid, double dt);
