import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import argparse
import os

def parse_time_from_filename(fname: str) -> float:
    t_str = fname.split("_t_")[-1].replace(".csv", "")
    return float(t_str)

def load_frame(fname):
    return pd.read_csv(fname)

def find_nearest_files(files, t_dim_all, targets):
    chosen = []
    chosen_t = []
    for tt in targets:
        k = int(np.argmin(np.abs(t_dim_all - tt)))
        chosen.append(files[k])
        chosen_t.append(t_dim_all[k])
    return chosen, chosen_t

def compute_uv(df):
    cols = set(df.columns)
    if "u" in cols and "v" in cols:
        return df["u"].values, df["v"].values
    # fallback to conserved
    if "rho" in cols and "rhou" in cols and "rhov" in cols:
        rho = df["rho"].values
        u = df["rhou"].values / rho
        v = df["rhov"].values / rho
        return u, v
    raise ValueError("CSV must contain either (u,v) or (rho,rhou,rhov).")

# -----------------------------
# CLI
# -----------------------------
parser = argparse.ArgumentParser(description="Fig9-style vorticity contour panels from CSV outputs")
parser.add_argument("mode", choices=["serial", "omp", "mpi"], help="res/<mode>")
parser.add_argument("--targets", type=float, nargs="*", default=[0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 12.6, 19.0],
                    help="Dimensionless target times t* (default: paper-like list)")
parser.add_argument("--T_end_phys", type=float, default=0.0011741, help="Physical end time in seconds")
parser.add_argument("--T_end_dim",  type=float, default=19.0, help="Dimensionless end time (default 19)")
# vorticity contour control
parser.add_argument("--omega_max", type=float, default=3000.0,
                    help="Max |omega| for contour levels (default 3000). Tune to match paper look.")
parser.add_argument("--n_levels", type=int, default=41, help="Number of vorticity contour levels (odd recommended)")
parser.add_argument("--linewidth", type=float, default=0.5)
args = parser.parse_args()

csv_dir = f"res/{args.mode}"
files = sorted(glob.glob(f"{csv_dir}/step_*.csv"))
if len(files) == 0:
    raise ValueError(f"No files: {csv_dir}/step_*.csv")

t_phys = np.array([parse_time_from_filename(f) for f in files], dtype=float)
t_ref  = args.T_end_phys / args.T_end_dim
t_dim  = t_phys / t_ref

targets = np.array(args.targets, dtype=float)
chosen_files, chosen_t_dim = find_nearest_files(files, t_dim, targets)

# grid info
df0 = load_frame(chosen_files[0])
xs = np.unique(df0["x"].values)
ys = np.unique(df0["y"].values)
nx, ny = len(xs), len(ys)
extent = [xs.min(), xs.max(), ys.min(), ys.max()]
dx = xs[1] - xs[0]
dy = ys[1] - ys[0]

# symmetric levels around 0
omega_max = float(args.omega_max)
levels = np.linspace(-omega_max, omega_max, args.n_levels)

# -----------------------------
# plotting
# -----------------------------
N = len(chosen_files)
fig, axes = plt.subplots(N, 1, figsize=(6, 2.0 * N))
if N == 1:
    axes = [axes]

for ax, fname, tstar in zip(axes, chosen_files, chosen_t_dim):
    df = load_frame(fname)
    u_flat, v_flat = compute_uv(df)

    u = u_flat.reshape(ny, nx)
    v = v_flat.reshape(ny, nx)

    # vorticity omega = dv/dx - du/dy
    dv_dx = np.gradient(v, dx, axis=1)
    du_dy = np.gradient(u, dy, axis=0)
    omega = dv_dx - du_dy

    ax.contour(omega, levels=levels, colors="k", linewidths=args.linewidth, extent=extent)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"t={tstar:.1f}", fontsize=12, pad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.tight_layout(h_pad=0.2)

out_path = os.path.join(csv_dir, f"fig9_vorticity_contours_{args.mode}.png")
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")