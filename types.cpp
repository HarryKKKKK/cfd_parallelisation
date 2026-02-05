#include "types.hpp"

void Grid::init(int nx_, int ny_, int ng_,
                double Lx, double Ly,
                double x0_, double y0_) {
    nx = nx_;
    ny = ny_;
    ng = ng_;
    x0 = x0_;
    y0 = y0_;
    dx = Lx / nx;
    dy = Ly / ny;
    U.resize((nx + 2*ng) * (ny + 2*ng));
}
