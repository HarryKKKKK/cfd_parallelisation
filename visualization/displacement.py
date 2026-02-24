import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NX, NY = 500, 197
X_LEN, Y_LEN = 0.225, 0.089
T_END_PHYS = 0.0011741
T_END_DIM  = 19.0
PATTERN = "step_*_t_*.csv"
OUT_NAME = "displacement.png"

# bubble setting
BUBBLE_CY = 0.0445

# FV coords
dx = X_LEN / NX
dy = Y_LEN / NY
xs = (np.arange(NX) + 0.5) * dx
ys = (np.arange(NY) + 0.5) * dy
j_mid = int(np.argmin(np.abs(ys - BUBBLE_CY)))

# tune this once
RHO_THR = 1.0
X_WINDOW = (0.0, 0.12)

def parse_time(fname: str) -> float:
    base = os.path.basename(fname)
    t_str = base.split("_t_")[-1].replace(".csv", "")
    return float(t_str)

def find_crossings(x, y, thr):
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

def main(indir="res/serial"):
    files = sorted(glob.glob(os.path.join(indir, PATTERN)))
    if not files:
        raise SystemExit("No csv files found.")

    t_ref = T_END_PHYS / T_END_DIM
    t_phys = np.array([parse_time(f) for f in files], dtype=float)
    t_dim  = t_phys / t_ref

    xL, xR = X_WINDOW
    iL = int(np.searchsorted(xs, xL, side="left"))
    iR = int(np.searchsorted(xs, xR, side="right"))

    ts, up, down = [], [], []

    for f, t in zip(files, t_dim):
        df = pd.read_csv(f)
        rho = df["rho"].to_numpy(float).reshape((NY, NX))
        prof = rho[j_mid, :]
        crosses = find_crossings(xs[iL:iR], prof[iL:iR], RHO_THR)
        if len(crosses) >= 2:
            ts.append(t); up.append(crosses[0]); down.append(crosses[-1])

    ts = np.array(ts)
    up = np.array(up)
    down = np.array(down)

    # displacement relative to first available point
    up0, down0 = up[0], down[0]
    dup = up - up0
    ddown = down - down0

    plt.figure(figsize=(6,4.5))
    plt.plot(ts, dup, label="Upstream displacement")
    plt.plot(ts, ddown, label="Downstream displacement")
    plt.xlabel("t (dimensionless)")
    plt.ylabel("Δx (m)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    out_path = os.path.join(indir, OUT_NAME)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print("Saved:", out_path)
    print("NOTE: tune RHO_THR if displacement jumps/noisy.")

if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv)>1 else "res/serial")