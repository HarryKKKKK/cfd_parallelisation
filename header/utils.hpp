#pragma once
#include "types.hpp"

#include <iostream>
#include <string>

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