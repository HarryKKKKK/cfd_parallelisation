#pragma once
#include "types.hpp"

// --- Conversion between conserved and primitive ---
Primitive conserved_to_primitive(const Conserved& U);
Conserved primitive_to_conserved(const Primitive& W);

// --- Euler fluxes ---
Conserved flux_x(const Conserved& U);
Conserved flux_y(const Conserved& U);
Conserved hll_flux_x(const Conserved& UL, const Conserved& UR);
Conserved hll_flux_y(const Conserved& UL, const Conserved& UR);

// --- Wave-speed helpers ---
double pressure_from_conserved(const Conserved& U);
double sound_speed_from_conserved(const Conserved& U);
double max_wave_speed_x(const Conserved& U);
double max_wave_speed_y(const Conserved& U);
double sound_speed(const Primitive& W);