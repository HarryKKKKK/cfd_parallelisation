import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import math
import argparse
import os

# =============================
# command line argument
# =============================
parser = argparse.ArgumentParser(
    description="Plot schlieren evolution (log10(|grad rho|)) from CFD output CSVs"
)
parser.add_argument(
    "mode",
    choices=["serial", "omp", "mpi"],
    help="Select which result folder to use (res/<mode>)"
)
args = parser.parse_args()

csv_dir = f"res/{args.mode}"

# =============================
# settings
# =============================
n_panels = 10
ncols = 2
cmap = "gray"

# =============================
# load csv files
# =============================
files = sorted(glob.glob(f"{csv_dir}/step_*.csv"))
if len(files) < n_panels:
    raise ValueError(f"Not enough CSV files in {csv_dir}: found {len(files)}, need {n_panels}")

# pick evenly spaced snapshots
indices = np.linspace(0, len(files) - 1, n_panels, dtype=int)
files = [files[i] for i in indices]

# =============================
# grid info from first file
# =============================
df0 = pd.read_csv(files[0])
xs = np.unique(df0["x"].values)
ys = np.unique(df0["y"].values)

nx = len(xs)
ny = len(ys)

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]

# =============================
# plotting
# =============================
nrows = math.ceil(n_panels / ncols)

fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(10, 14),
    sharex=True,
    sharey=True,
    constrained_layout=True
)

axes = axes.flatten()

for k, (ax, fname) in enumerate(zip(axes, files)):
    df = pd.read_csv(fname)
    rho = df["rho"].values.reshape(ny, nx)

    # -------- schlieren --------
    drho_dx = np.gradient(rho, dx, axis=1)
    drho_dy = np.gradient(rho, dy, axis=0)
    schlieren = np.sqrt(drho_dx**2 + drho_dy**2)

    im = ax.imshow(
        np.log10(schlieren + 1e-12),
        origin="lower",
        extent=[xs.min(), xs.max(), ys.min(), ys.max()],
        aspect="equal",
        cmap=cmap
    )

    # panel label
    label = chr(97 + k)  # a, b, c, ...
    t_str = fname.split("_t_")[-1].replace(".csv", "")
    ax.set_title(f"({label})  t = {t_str}", fontsize=10)

    ax.set_xticks([])
    ax.set_yticks([])

# hide unused axes
for ax in axes[len(files):]:
    ax.axis("off")

# shared colorbar
cbar = fig.colorbar(im, ax=axes, shrink=0.9)
cbar.set_label(r"$\log_{10}(|\nabla \rho|)$")

fig.suptitle(
    f"Shock–Bubble Interaction ({args.mode.upper()}): Schlieren Evolution",
    fontsize=16
)

# =============================
# save
# =============================
out_path = os.path.join(csv_dir, f"schlieren_{args.mode}.png")
plt.savefig(out_path, dpi=300)
print(f"Saved plot to {out_path}")
