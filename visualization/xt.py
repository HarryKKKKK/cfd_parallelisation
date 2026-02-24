import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ========= config (match your setup) =========
NX, NY = 500, 197
X_LEN, Y_LEN = 0.225, 0.089

T_END_PHYS = 0.0011741
T_END_DIM  = 19.0

PATTERN  = "step_*_t_*.csv"
OUT_NAME = "xt.png"

# bubble setting
R_BUBBLE  = 0.025
BUBBLE_CX = 0.035
BUBBLE_CY = 0.0445

# FV cell-center coords
dx = X_LEN / NX
dy = Y_LEN / NY
xs = (np.arange(NX) + 0.5) * dx
ys = (np.arange(NY) + 0.5) * dy

# pick y index closest to bubble centerline
j_mid = int(np.argmin(np.abs(ys - BUBBLE_CY)))

# --- interface threshold: choose mid value between bubble rho and ambient rho
# You can change these two if you know your initial values:
RHO_BUBBLE0  = 0.214   # from your earlier debug: rho_min ~0.214
RHO_AMBIENT0 = 1.29    # from your earlier debug: top repeated ~1.29
RHO_THR = 0.5 * (RHO_BUBBLE0 + RHO_AMBIENT0)

# --- tracking window settings ---
# initial bubble x-range
X_UP0   = BUBBLE_CX - R_BUBBLE
X_DN0   = BUBBLE_CX + R_BUBBLE
# search half-width around previous interface position
WIN_HALF = 0.02   # meters; tune 0.015~0.03 if needed

# jet center: search between interfaces and a bit downstream
JET_PAD_LEFT  = 0.001
JET_PAD_RIGHT = 0.01
# ============================================

def parse_time(fname: str) -> float:
    base = os.path.basename(fname)
    t_str = base.split("_t_")[-1].replace(".csv", "")
    return float(t_str)

def load_rho(df: pd.DataFrame) -> np.ndarray:
    # rho stored row-major (NY, NX) in your scripts
    return df["rho"].to_numpy(dtype=float).reshape((NY, NX))

def find_crossings(x, y, thr):
    """Return all x where y crosses thr (linear interpolation)."""
    s = y - thr
    idx = np.where(s[:-1] * s[1:] <= 0)[0]
    xs_cross = []
    for k in idx:
        s0, s1 = s[k], s[k+1]
        if s0 == s1:
            xs_cross.append(x[k])
        else:
            alpha = s0 / (s0 - s1)
            xs_cross.append(x[k] + alpha * (x[k+1] - x[k]))
    return xs_cross

def pick_nearest(candidates, x_prev):
    if not candidates:
        return np.nan
    candidates = np.array(candidates, dtype=float)
    return candidates[np.argmin(np.abs(candidates - x_prev))]

def clip_window(x_center, half_width):
    xL = max(0.0, x_center - half_width)
    xR = min(xs[-1], x_center + half_width)
    iL = int(np.searchsorted(xs, xL, side="left"))
    iR = int(np.searchsorted(xs, xR, side="right"))
    return iL, iR

def main(indir="res/serial"):
    files = sorted(glob.glob(os.path.join(indir, PATTERN)))
    if not files:
        raise SystemExit("No csv files found.")

    t_ref = T_END_PHYS / T_END_DIM
    t_phys_all = np.array([parse_time(f) for f in files], dtype=float)
    t_dim_all  = t_phys_all / t_ref

    # sort by time just in case
    order = np.argsort(t_dim_all)
    t_dim_all = t_dim_all[order]
    files = [files[i] for i in order]

    # trajectories
    ts, x_up, x_dn, x_jet = [], [], [], []

    # previous positions (start from initial bubble)
    xup_prev = X_UP0
    xdn_prev = X_DN0

    for f, tdim in zip(files, t_dim_all):
        df = pd.read_csv(f)
        rho = load_rho(df)
        prof = rho[j_mid, :]

        # --- upstream: search crossings near previous xup ---
        iL, iR = clip_window(xup_prev, WIN_HALF)
        crosses_up = find_crossings(xs[iL:iR], prof[iL:iR], RHO_THR)
        xup = pick_nearest(crosses_up, xup_prev)

        # --- downstream: search crossings near previous xdn ---
        iL, iR = clip_window(xdn_prev, WIN_HALF)
        crosses_dn = find_crossings(xs[iL:iR], prof[iL:iR], RHO_THR)
        xdn = pick_nearest(crosses_dn, xdn_prev)

        # if one side failed, skip this frame (keeps curves clean)
        if not np.isfinite(xup) or not np.isfinite(xdn):
            continue

        # enforce ordering (upstream < downstream); if violated, skip
        if xup >= xdn:
            continue

        # --- jet center: look for minimum rho between interfaces (+ small pads) ---
        jetL = max(0.0, xup + JET_PAD_LEFT)
        jetR = min(xs[-1], xdn + JET_PAD_RIGHT)
        iJL = int(np.searchsorted(xs, jetL, side="left"))
        iJR = int(np.searchsorted(xs, jetR, side="right"))
        if iJR <= iJL + 2:
            xj = np.nan
        else:
            k = int(np.argmin(prof[iJL:iJR]))
            xj = xs[iJL + k]

        ts.append(tdim)
        x_up.append(xup)
        x_dn.append(xdn)
        x_jet.append(xj)

        # update prev for tracking
        xup_prev = xup
        xdn_prev = xdn

    ts = np.array(ts)
    x_up = np.array(x_up)
    x_dn = np.array(x_dn)
    x_jet = np.array(x_jet)

    # ---- plot ----
    plt.figure(figsize=(7,5))
    plt.plot(ts, x_up,  label=f"Upstream interface (rho={RHO_THR:.3g})")
    plt.plot(ts, x_dn,  label=f"Downstream interface (rho={RHO_THR:.3g})")
    plt.plot(ts, x_jet, label="Jet center (min rho between interfaces)")

    plt.xlabel("t (dimensionless)")
    plt.ylabel("x (m)")
    plt.xlim(0, T_END_DIM)
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    out_path = os.path.join(indir, OUT_NAME)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print("Saved:", out_path)
    print("Notes:")
    print(f"  Using y ≈ {ys[j_mid]:.6g} (closest to bubble y_c={BUBBLE_CY})")
    print(f"  RHO_THR = {RHO_THR:.6g} (mid of {RHO_BUBBLE0} and {RHO_AMBIENT0})")
    print(f"  WIN_HALF = {WIN_HALF} m (increase if tracking breaks)")

if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv)>1 else "res/serial")