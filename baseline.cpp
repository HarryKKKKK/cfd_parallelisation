#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <array>
#include <chrono>
#include <limits>
#include <string>

// ========== 1. Initial Data ==========

const double GAMMA = 1.4;
const double Mach = 1.22;
const double rho_air = 1.29;
const double rho_helium = 0.214;
const double p_atm = 101325.0;

const double bubble_radius = 0.025;
const double bubble_x = 0.035;
const double bubble_y = 0.0445;
const double shock_x = 0.005;

const double x_start = 0.0;
const double x_end = 0.225;
const double y_start = 0.0;
double y_end = 0.089;

const double y_end_base = 0.089; 
int num_bubbles = 1;

const int nGhosts = 2;

const double CFL = 0.45;
const double T_end = 0.0011741;

int nx = 500;
int ny = 197;

double dx;
double dy;
int total_x;
int total_y;
int n_grid;

inline int index(int i, int j) {
    return j * total_x + i;
}

// ========== 2. Data Structure ==========

struct Primitive {
    double rho;
    double u;
    double v;
    double p;
};

struct Conserved {
    double rho;
    double mom_u;
    double mom_v;
    double E;
    
    Conserved operator+(const Conserved& other) const {
        return {rho + other.rho, mom_u + other.mom_u, mom_v + other.mom_v, E + other.E};
    }
    
    Conserved operator-(const Conserved& other) const {
        return {rho - other.rho, mom_u - other.mom_u, mom_v - other.mom_v, E - other.E};
    }

    Conserved operator*(double scalar) const {
        return {rho * scalar, mom_u * scalar, mom_v * scalar, E * scalar};
    }
};

struct Flux {
    double rho;
    double mom_u;
    double mom_v;
    double energy;
};

// ========== 3. Physical Function ==========

inline double computeSoundSpeed(const Primitive& prim) {
    return std::sqrt(GAMMA * prim.p / prim.rho);
}

inline double computeTotalEnergy(const Primitive &prim) {
    return prim.p / (GAMMA - 1.0) + 0.5 * prim.rho * (prim.u * prim.u + prim.v * prim.v);
}

inline Conserved primToCons(const Primitive& prim) {
    Conserved cons;
    cons.rho    = prim.rho;
    cons.mom_u  = prim.rho * prim.u;
    cons.mom_v  = prim.rho * prim.v;
    cons.E      = computeTotalEnergy(prim);
    return cons;
}

inline Primitive consToPrim(const Conserved& cons){
    Primitive prim;
    prim.rho = cons.rho;
    double inv_rho = 1.0 / prim.rho;
    prim.u      = cons.mom_u * inv_rho;
    prim.v      = cons.mom_v * inv_rho;
    double kinetic = 0.5 * inv_rho * (cons.mom_u * cons.mom_u + cons.mom_v * cons.mom_v);
    prim.p      = (cons.E - kinetic) * (GAMMA - 1.0);
    return prim;
}

void updatePrimitivesGlobal(const std::vector<Conserved>& cons_grid, std::vector<Primitive>& prim_grid) {
    for (int i = 0; i < n_grid; ++i) {
        prim_grid[i] = consToPrim(cons_grid[i]);
    }
}

// ========== 4. Flux Functions ==========

inline Flux computeFluxFromCons(const Conserved& U) {
    Flux flux;
    double inv_rho = 1.0 / U.rho;
    double u = U.mom_u * inv_rho;
    double v = U.mom_v * inv_rho;
    double kinetic = 0.5 * inv_rho * (U.mom_u * U.mom_u + U.mom_v * U.mom_v);
    double p = (GAMMA - 1.0) * (U.E - kinetic);

    flux.rho    = U.mom_u;
    flux.mom_u  = U.mom_u * u + p;
    flux.mom_v  = U.mom_u * v;
    flux.energy = u * (U.E + p);
    return flux;
}

inline Flux computeFluxFromBoth(const Primitive& p, const Conserved& U) {
    Flux f;
    f.rho    = U.mom_u;
    f.mom_u  = U.mom_u * p.u + p.p;
    f.mom_v  = U.mom_u * p.v;
    f.energy = p.u * (U.E + p.p);
    return f;
}

inline Flux computeFluxYFromCons(const Conserved& U) {
    Flux flux;
    double inv_rho = 1.0 / U.rho;
    double u = U.mom_u * inv_rho;
    double v = U.mom_v * inv_rho;
    double kinetic = 0.5 * inv_rho * (U.mom_u * U.mom_u + U.mom_v * U.mom_v);
    double p = (GAMMA - 1.0) * (U.E - kinetic);

    flux.rho    = U.mom_v;
    flux.mom_u  = U.mom_v * u;
    flux.mom_v  = U.mom_v * v + p;
    flux.energy = v * (U.E + p);
    return flux;
}

inline Flux computeFluxYFromBoth(const Primitive& p, const Conserved& U) {
    Flux f;
    f.rho    = U.mom_v;
    f.mom_u  = U.mom_v * p.u;
    f.mom_v  = U.mom_v * p.v + p.p;
    f.energy = p.v * (U.E + p.p);
    return f;
}

// ========== 5. Slope Limiter ==========

inline double getVanLeerLimiter(double deltaL, double deltaR){
    if (std::abs(deltaR) < 1e-12) return 0.0;
    double r = deltaL / deltaR;
    double xiL = 2.0 * r / (1.0 + r);
    double xiR = 2.0 / (1.0 + r);
    if (r <= 0.0) return 0.0;
    else return xiL < xiR ? xiL : xiR; 
}

// ========== 6. MUSCL Data Reconstruction ==========

inline double getSlopeDelta(double uL, double u_i, double uR){
    double deltaL = u_i - uL;
    double deltaR = uR - u_i;
    double delta_i = 0.5 * (deltaL + deltaR);
    double xi = getVanLeerLimiter(deltaL, deltaR);
    return xi * delta_i;
}

inline Conserved getSlopeDeltaVector(const Conserved& uL, const Conserved& u_i, const Conserved& uR){
    Conserved delta_vec;
    delta_vec.rho   = getSlopeDelta(uL.rho,   u_i.rho,   uR.rho);
    delta_vec.mom_u = getSlopeDelta(uL.mom_u, u_i.mom_u, uR.mom_u);
    delta_vec.mom_v = getSlopeDelta(uL.mom_v, u_i.mom_v, uR.mom_v);
    delta_vec.E     = getSlopeDelta(uL.E,     u_i.E,     uR.E);
    return delta_vec;
}

inline std::array<Conserved,2> dataReconstruction(const Conserved& uL, const Conserved& u_i, const Conserved& uR){
    Conserved delta = getSlopeDeltaVector(uL, u_i, uR);
    Conserved U_minus = u_i - delta * 0.5;
    Conserved U_plus  = u_i + delta * 0.5;
    return {U_minus, U_plus};
}

// ========== 7. HLLC Approximation ==========

inline Conserved getHLLCState(const Primitive& prim, const Conserved& U, double S_K, double S_Star){
    Conserved U_HLLC_K;
    double rho_star_K = prim.rho * (S_K - prim.u) / (S_K - S_Star);
    U_HLLC_K.rho   = rho_star_K;
    U_HLLC_K.mom_u = rho_star_K * S_Star;
    U_HLLC_K.mom_v = rho_star_K * prim.v;
    U_HLLC_K.E     = rho_star_K * (U.E / prim.rho + (S_Star - prim.u) * (S_Star + prim.p / (prim.rho * (S_K - prim.u))));
    return U_HLLC_K;
}

inline Conserved getHLLCStateY(const Primitive& prim, const Conserved& U, double S_K, double S_Star){
    Conserved U_HLLC_K;
    double rho_star_K = prim.rho * (S_K - prim.v) / (S_K - S_Star);
    U_HLLC_K.rho   = rho_star_K;
    U_HLLC_K.mom_u = rho_star_K * prim.u; 
    U_HLLC_K.mom_v = rho_star_K * S_Star; 
    U_HLLC_K.E     = rho_star_K * (U.E / prim.rho + (S_Star - prim.v) * (S_Star + prim.p / (prim.rho * (S_K - prim.v))));
    return U_HLLC_K;
}

inline Flux fluxHLLC(const Primitive& state_L, const Conserved& U_L, 
                     const Primitive& state_R, const Conserved& U_R){
    double rho_L = state_L.rho, u_L = state_L.u, p_L = state_L.p;
    double rho_R = state_R.rho, u_R = state_R.u, p_R = state_R.p;

    double a_L = computeSoundSpeed(state_L);
    double a_R = computeSoundSpeed(state_R);
    double rho_bar = 0.5 * (rho_L + rho_R);
    double a_bar = 0.5 * (a_L + a_R);

    double p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * (rho_bar * a_bar);
    double p_star = p_pvrs > 0.0 ? p_pvrs : 0.0;

    double q_L = (p_star <= p_L) ? 1.0 : std::sqrt(1.0 + (GAMMA + 1.0) / (2.0 * GAMMA) * (p_star / p_L - 1));
    double S_L = u_L - a_L * q_L;
    if (S_L >= 0.0) return computeFluxFromBoth(state_L, U_L);

    double q_R = (p_star <= p_R) ? 1.0 : std::sqrt(1.0 + (GAMMA + 1.0) / (2.0 * GAMMA) * (p_star / p_R - 1));
    double S_R = u_R + a_R * q_R;
    if (S_R <= 0.0) return computeFluxFromBoth(state_R, U_R);

    double S_Star = (p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R)) / (rho_L*(S_L - u_L) - rho_R*(S_R - u_R));

    Flux flux;
    if (S_L < 0.0 && S_Star >= 0.0) {
        Flux fL = computeFluxFromBoth(state_L, U_L);
        Conserved U_HLLC_L = getHLLCState(state_L, U_L, S_L, S_Star);
        flux.rho    = fL.rho    + S_L * (U_HLLC_L.rho - U_L.rho);
        flux.mom_u  = fL.mom_u  + S_L * (U_HLLC_L.mom_u - U_L.mom_u);
        flux.mom_v  = fL.mom_v  + S_L * (U_HLLC_L.mom_v - U_L.mom_v);
        flux.energy = fL.energy + S_L * (U_HLLC_L.E - U_L.E);
    } else { 
        Flux fR = computeFluxFromBoth(state_R, U_R);
        Conserved U_HLLC_R = getHLLCState(state_R, U_R, S_R, S_Star);
        flux.rho    = fR.rho    + S_R * (U_HLLC_R.rho - U_R.rho);
        flux.mom_u  = fR.mom_u  + S_R * (U_HLLC_R.mom_u - U_R.mom_u);
        flux.mom_v  = fR.mom_v  + S_R * (U_HLLC_R.mom_v - U_R.mom_v);
        flux.energy = fR.energy + S_R * (U_HLLC_R.E - U_R.E);
    } 
    return flux;
}

inline Flux fluxHLLC_Y(const Primitive& state_L, const Conserved& U_L, 
                       const Primitive& state_R, const Conserved& U_R){
    double rho_L = state_L.rho, v_L = state_L.v, p_L = state_L.p;
    double rho_R = state_R.rho, v_R = state_R.v, p_R = state_R.p;

    double a_L = computeSoundSpeed(state_L);
    double a_R = computeSoundSpeed(state_R);
    double rho_bar = 0.5 * (rho_L + rho_R);
    double a_bar = 0.5 * (a_L + a_R);

    double p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (v_R - v_L) * (rho_bar * a_bar);
    double p_star = p_pvrs > 0.0 ? p_pvrs : 0.0;

    double q_L = (p_star <= p_L) ? 1.0 : std::sqrt(1.0 + (GAMMA + 1.0) / (2.0 * GAMMA) * (p_star / p_L - 1));
    double S_L = v_L - a_L * q_L;
    if (S_L >= 0.0) return computeFluxYFromBoth(state_L, U_L);

    double q_R = (p_star <= p_R) ? 1.0 : std::sqrt(1.0 + (GAMMA + 1.0) / (2.0 * GAMMA) * (p_star / p_R - 1));
    double S_R = v_R + a_R * q_R;
    if (S_R <= 0.0) return computeFluxYFromBoth(state_R, U_R);

    double S_Star = (p_R - p_L + rho_L*v_L*(S_L - v_L) - rho_R*v_R*(S_R - v_R)) / (rho_L*(S_L - v_L) - rho_R*(S_R - v_R));

    Flux flux;
    if (S_L < 0.0 && S_Star >= 0.0) {
        Flux fL = computeFluxYFromBoth(state_L, U_L);
        Conserved U_HLLC_L = getHLLCStateY(state_L, U_L, S_L, S_Star);
        flux.rho    = fL.rho    + S_L * (U_HLLC_L.rho - U_L.rho);
        flux.mom_u  = fL.mom_u  + S_L * (U_HLLC_L.mom_u - U_L.mom_u);
        flux.mom_v  = fL.mom_v  + S_L * (U_HLLC_L.mom_v - U_L.mom_v);
        flux.energy = fL.energy + S_L * (U_HLLC_L.E - U_L.E);
    } else { 
        Flux fR = computeFluxYFromBoth(state_R, U_R);
        Conserved U_HLLC_R = getHLLCStateY(state_R, U_R, S_R, S_Star);
        flux.rho    = fR.rho    + S_R * (U_HLLC_R.rho - U_R.rho);
        flux.mom_u  = fR.mom_u  + S_R * (U_HLLC_R.mom_u - U_R.mom_u);
        flux.mom_v  = fR.mom_v  + S_R * (U_HLLC_R.mom_v - U_R.mom_v);
        flux.energy = fR.energy + S_R * (U_HLLC_R.E - U_R.E);
    } 
    return flux;
}

// ========== 8. Time Step Calculation ==========

double computeTimeStep(const std::vector<Primitive> &prim_grid) {
    double max_wave_speed = 0.0;

    for (int j = nGhosts; j < total_y - nGhosts; ++j) {
        for (int i = nGhosts; i < total_x - nGhosts; ++i) {
            int idx = index(i, j);
            const Primitive &p = prim_grid[idx];
            double cs = std::sqrt(GAMMA * p.p / p.rho);
            double curr_wave_speed = std::sqrt(p.u * p.u + p.v * p.v) + cs;
            if (curr_wave_speed > max_wave_speed) max_wave_speed = curr_wave_speed;
        }
    }

    if (max_wave_speed <= 1e-12) return std::numeric_limits<double>::infinity();
    double min_dx_dy = dx < dy ? dx : dy;
    return CFL * min_dx_dy / max_wave_speed;
}

// ========== 9. Boundary Condition ==========

void applyBoundaryConditions(std::vector<Conserved> &cons_grid) {
    for (int j = 0; j < total_y; j++) {
        for (int i = 0; i < nGhosts; i++) {
            cons_grid[index(i, j)] = cons_grid[index(nGhosts, j)];
            cons_grid[index(total_x - 1 - i, j)] = cons_grid[index(total_x - 1 - nGhosts, j)];
        }
    }
    for (int j = 0; j < nGhosts; j++) {
        for (int i = 0; i < total_x; i++) {
            cons_grid[index(i, j)] = cons_grid[index(i, nGhosts)];
            cons_grid[index(i, total_y - 1 - j)] = cons_grid[index(i, total_y - 1 - nGhosts)];
        }
    }
}

// ========== 10. Initial Conditions ==========

void initializeShockBubble(std::vector<Conserved> &cons_grid) {
    double rho1 = rho_air;
    double p1 = p_atm;
    double c_air = std::sqrt(GAMMA * p1 / rho1);
    double rho_post = rho1 * ((GAMMA + 1.0) * Mach * Mach) / ((GAMMA - 1.0) * Mach * Mach + 2.0);
    double p_post = p1 * (1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (Mach * Mach - 1.0));
    double u_post = (2.0 * c_air / (GAMMA + 1.0)) * (Mach - 1.0 / Mach);

    std::cout << "Initialize Shock Bubble Interaction (" << num_bubbles << " Bubbles)..." << std::endl;

    for (int j = 0; j < total_y; j++) {
        for (int i = 0; i < total_x; i++) {
            double y = y_start + (j - nGhosts + 0.5) * dy;
            double x = x_start + (i - nGhosts + 0.5) * dx;
            int idx = index(i, j);

            Primitive temp_prim;
            if (x < shock_x) temp_prim = {rho_post, u_post, 0.0, p_post};
            else temp_prim = {rho1, 0.0, 0.0, p1};

            for (int b = 0; b < num_bubbles; b++) {
                double current_bubble_y = bubble_y + b * y_end_base;
                if (std::sqrt((x - bubble_x)*(x - bubble_x) +(y - current_bubble_y)*(y - current_bubble_y)) < bubble_radius) {
                    temp_prim.rho = rho_helium;
                    break;
                }
            }
            cons_grid[idx] = primToCons(temp_prim);
        }
    } 
}

// ========== 11. Output Function ==========

void dataOutput(const std::vector<Primitive>& prim_grid, int step, double time) {
    std::string filename = "output_" + std::to_string(step) + ".dat";
    std::ofstream output(filename);
    output<<std::scientific<<std::setprecision(18);
    output << "Variables= x y rho u v p" << "\n";
    output << "nx = " << nx << " ny =" << ny << "\n";

    for (int j = nGhosts; j < ny + nGhosts; j++) {
        double y = y_start + (j - nGhosts + 0.5) * dy;
        for (int i = nGhosts; i < nx + nGhosts; i++) {
            double x = x_start + (i - nGhosts + 0.5) * dx;  
            const Primitive& p = prim_grid[index(i, j)];
            output << x << " " << y << " " << p.rho << " " << p.u << " " << p.v << " " << p.p << "\n";
        }
    }
    output.close();
    std::cout << "Data saved to " << filename << " (Time: " << time << ")" << std::endl;
}

// ========== 12. Update Solution (Serial Cache-Blocked) ==========

// 定义行数据缓存区（仅包含单行数据，锁定在 CPU L1/L2 Cache）
struct RowData {
    std::vector<Primitive> L, R, D, U;
    std::vector<Conserved> L_cons, R_cons, D_cons, U_cons;
    RowData(int size) : L(size), R(size), D(size), U(size), L_cons(size), R_cons(size), D_cons(size), U_cons(size) {}
};

// 定义私有工作区
struct Workspace {
    RowData row_curr;
    RowData row_next;
    std::vector<Flux> flux_x;
    std::vector<Flux> flux_y_curr;
    std::vector<Flux> flux_y_prev;

    Workspace(int size) : row_curr(size), row_next(size), flux_x(size), flux_y_curr(size), flux_y_prev(size) {}
};

// 单行重构函数
inline void compute_U_bar_row(const std::vector<Conserved>& cons_grid, int total_x, int nx, int nGhosts,
                              int j, double rx_half, double ry_half, RowData& row) {
    for (int i = nGhosts - 1; i <= nx + nGhosts; i++) {
        int idx = j * total_x + i;
        const Conserved& U_ij = cons_grid[idx];

        std::array<Conserved, 2> rec_x = dataReconstruction(cons_grid[idx - 1], U_ij, cons_grid[idx + 1]);
        std::array<Conserved, 2> rec_y = dataReconstruction(cons_grid[idx - total_x], U_ij, cons_grid[idx + total_x]);

        Flux F_L = computeFluxFromCons(rec_x[0]);
        Flux F_R = computeFluxFromCons(rec_x[1]);
        Flux G_D = computeFluxYFromCons(rec_y[0]);
        Flux G_U = computeFluxYFromCons(rec_y[1]);

        Conserved dU;
        dU.rho   = rx_half * (F_L.rho    - F_R.rho)    + ry_half * (G_D.rho    - G_U.rho);
        dU.mom_u = rx_half * (F_L.mom_u  - F_R.mom_u)  + ry_half * (G_D.mom_u  - G_U.mom_u);
        dU.mom_v = rx_half * (F_L.mom_v  - F_R.mom_v)  + ry_half * (G_D.mom_v  - G_U.mom_v);
        dU.E     = rx_half * (F_L.energy - F_R.energy) + ry_half * (G_D.energy - G_U.energy);

        row.L_cons[i] = rec_x[0] + dU; row.L[i] = consToPrim(row.L_cons[i]);
        row.R_cons[i] = rec_x[1] + dU; row.R[i] = consToPrim(row.R_cons[i]);
        row.D_cons[i] = rec_y[0] + dU; row.D[i] = consToPrim(row.D_cons[i]);
        row.U_cons[i] = rec_y[1] + dU; row.U[i] = consToPrim(row.U_cons[i]);
    }
}

// 融合后的主更新函数
void updateSolution(const std::vector<Conserved>& cons_grid, 
                    std::vector<Conserved>& new_grid, double ratio_x, double ratio_y) 
{   
    const double rx_half = 0.5 * ratio_x; 
    const double ry_half = 0.5 * ratio_y;

    Workspace work(total_x); 
    
    int j_start = nGhosts;
    int j_end = ny + nGhosts - 1;

    // -- 1. 初始化滑动窗口 --
    compute_U_bar_row(cons_grid, total_x, nx, nGhosts, j_start - 1, rx_half, ry_half, work.row_curr);
    compute_U_bar_row(cons_grid, total_x, nx, nGhosts, j_start, rx_half, ry_half, work.row_next);
    
    for (int i = nGhosts; i <= nx + nGhosts - 1; i++) {
        work.flux_y_prev[i] = fluxHLLC_Y(work.row_curr.U[i], work.row_curr.U_cons[i], 
                                         work.row_next.D[i], work.row_next.D_cons[i]);
    }
    std::swap(work.row_curr, work.row_next); 

    // -- 2. 核心迭代 --
    for (int j = j_start; j <= j_end; j++) {
        
        compute_U_bar_row(cons_grid, total_x, nx, nGhosts, j + 1, rx_half, ry_half, work.row_next);

        for (int i = nGhosts - 1; i <= nx + nGhosts - 1; i++) { 
            work.flux_x[i] = fluxHLLC(work.row_curr.R[i], work.row_curr.R_cons[i], 
                                      work.row_curr.L[i + 1], work.row_curr.L_cons[i + 1]);
        }

        for (int i = nGhosts; i <= nx + nGhosts - 1; i++) {
            work.flux_y_curr[i] = fluxHLLC_Y(work.row_curr.U[i], work.row_curr.U_cons[i], 
                                             work.row_next.D[i], work.row_next.D_cons[i]); 
        }

        for (int i = nGhosts; i <= nx + nGhosts - 1; i++) {
            int idx = j * total_x + i;
            const Conserved& U_old = cons_grid[idx]; 
            
            new_grid[idx].rho   = U_old.rho   - ratio_x * (work.flux_x[i].rho   - work.flux_x[i - 1].rho)   - ratio_y * (work.flux_y_curr[i].rho   - work.flux_y_prev[i].rho);
            new_grid[idx].mom_u = U_old.mom_u - ratio_x * (work.flux_x[i].mom_u - work.flux_x[i - 1].mom_u) - ratio_y * (work.flux_y_curr[i].mom_u - work.flux_y_prev[i].mom_u);
            new_grid[idx].mom_v = U_old.mom_v - ratio_x * (work.flux_x[i].mom_v - work.flux_x[i - 1].mom_v) - ratio_y * (work.flux_y_curr[i].mom_v - work.flux_y_prev[i].mom_v);
            new_grid[idx].E     = U_old.E     - ratio_x * (work.flux_x[i].energy- work.flux_x[i - 1].energy)- ratio_y * (work.flux_y_curr[i].energy- work.flux_y_prev[i].energy);
        }

        std::swap(work.row_curr, work.row_next); 
        std::swap(work.flux_y_prev, work.flux_y_curr); 
    }
}

// ========== 13. Main Function ==========

int main(int argc, char* argv[]) {

    // 默认值
    nx = 500;
    ny = 197;
    num_bubbles = 1;

    // 解析终端参数，完美解绑网格与气泡数
    if (argc > 1) nx = std::atoi(argv[1]); 
    if (argc > 2) ny = std::atoi(argv[2]); 
    if (argc > 3) num_bubbles = std::atoi(argv[3]); 

    y_end = num_bubbles * y_end_base;

    dx = (x_end - x_start) / nx;
    dy = (y_end - y_start) / ny;

    total_x = nx + 2 * nGhosts;
    total_y = ny + 2 * nGhosts;
    n_grid = total_x * total_y;

    std::cout << "========================================" << std::endl;
    std::cout << "  2D Euler Solver (MUSCL-Hancock HLLC) " << std::endl;
    std::cout << "  Resolution: " << nx << " x " << ny << std::endl;
    std::cout << "  Bubbles Count: " << num_bubbles << std::endl;
    std::cout << "  Execution: Pure Serial" << std::endl;
    std::cout << "  Maths: Cache Blocking (Loop Fusion) Opt" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<Conserved> new_grid(n_grid);
    std::vector<Conserved> cons_grid(n_grid);
    std::vector<Primitive> prim_grid(n_grid);

    initializeShockBubble(cons_grid);
    applyBoundaryConditions(cons_grid);
    updatePrimitivesGlobal(cons_grid, prim_grid); 

    double t = 0.0;
    int step = 0;

    std::cout << "Starting simulation..." << std::endl;

    double total_compute_time = 0.0;

    while (t < T_end) {

        auto step_start = std::chrono::high_resolution_clock::now();

        double dt = computeTimeStep(prim_grid);
        if (t + dt > T_end) { dt = T_end - t; }

        updateSolution(cons_grid, new_grid, dt/dx, dt/dy);
        
        std::swap(cons_grid, new_grid);
        applyBoundaryConditions(cons_grid);
        updatePrimitivesGlobal(cons_grid, prim_grid);

        t += dt;
        step++;

        auto step_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> step_diff = step_end - step_start;
        total_compute_time += step_diff.count();
    }
    
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Simulation Finished." << std::endl;
    std::cout << "Total Computation Time (No I/O): " << std::fixed << std::setprecision(4) << total_compute_time << "s." << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;

    return 0;
}