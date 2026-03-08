#pragma once
#include "types.hpp"

// --- Conversion between conserved and primitive ---
Primitive conserved_to_primitive(const Conserved& U);
Conserved primitive_to_conserved(const Primitive& W);

// --- Wave-speed helpers ---
double pressure_from_conserved(const Conserved& U);
double sound_speed_from_conserved(const Conserved& U);
double sound_speed(const Primitive& W);