import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import math
import argparse
import os

# -----------------------------
# command line argument
# -----------------------------
parser = argparse.ArgumentParser(
    description="Plot density evolution from CFD output"
)
parser.add_argument(
    "mode",
    choices=["serial", "omp", "mpi"],
    help="Select which result folder to use"
)

args = parser.parse_args()

csv_dir = f"res/{args.mode}"
n_panels = 10
ncols = 2
cmap = "viridis"

# -----------------------------
# load files
# -----------------------------
files = sorted(glob.glob(f"{csv_dir}/step_*.csv"))

if len(files) < n_panels:
    raise ValueError("Not enough CSV files in folder")

# pick evenly spaced snapshots
indices = np.linspace(0, len(files)-1, n_panels, dtype=int)
files = [files[i] for i in indices]

# -----------------------------
# grid info
# -----------------------------
df0 = pd.read_csv(files[0])
xs = np.unique(df0["x"].values)
ys = np.unique(df0["y"].values)
nx, ny = len(xs), len(ys)

# -----------------------------
# plotting
# -----------------------------
nrows = math.ceil(n_panels / ncols)

fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(10, 14),
    sharex=True,
    sharey=True,
    constrained_layout=True
)

axes = axes.flatten()

for ax, fname in zip(axes, files):
    df = pd.read_csv(fname)
    rho = df["rho"].values.reshape(ny, nx)

    im = ax.imshow(
        rho,
        origin="lower",
        extent=[xs.min(), xs.max(), ys.min(), ys.max()],
        aspect="equal",
        cmap=cmap
    )

    # extract time from filename
    t_str = fname.split("_t_")[-1].replace(".csv", "")
    ax.set_title(f"t = {t_str}", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])

# hide unused axes
for ax in axes[len(files):]:
    ax.axis("off")

# shared colorbar
cbar = fig.colorbar(im, ax=axes, shrink=0.9)
cbar.set_label(r"$\rho$")

fig.suptitle(
    f"Shock–Bubble Interaction ({args.mode.upper()}): Density Evolution",
    fontsize=16
)

output_file = os.path.join(csv_dir, f"density_plot_{args.mode}.png")
plt.savefig(output_file, dpi=300)
print(f"Saved plot to {output_file}")
