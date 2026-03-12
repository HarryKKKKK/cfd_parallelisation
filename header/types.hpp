#pragma once
#include <vector>

struct Conserved {
    double rho;  
    double rhou;  
    double rhov; 
    double E;   

    Conserved operator+(const Conserved& other) const {
        return {rho + other.rho,
                rhou + other.rhou,
                rhov + other.rhov,
                E + other.E};
    }

    Conserved operator-(const Conserved& other) const {
        return {rho - other.rho,
                rhou - other.rhou,
                rhov - other.rhov,
                E - other.E};
    }

    Conserved operator*(double s) const {
        return {rho * s,
                rhou * s,
                rhov * s,
                E * s};
    }
};

struct Primitive {
    double rho;
    double u;
    double v;
    double p;
};

struct Grid {
    int nx, ny;
    int ng;
    double dx, dy;
    double x0, y0;

    std::vector<Conserved> U;
    std::vector<Conserved> U_new;

    void init(int nx_, int ny_, int ng_, double Lx, double Ly, double x0_ = 0.0, double y0_ = 0.0);

    inline int idx(int i, int j) const {
        return i + (nx + 2*ng) * j;
    }
};