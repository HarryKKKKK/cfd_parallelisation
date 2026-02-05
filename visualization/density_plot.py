import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import math

# -----------------------------
# settings
# -----------------------------
csv_dir = "res/serial"
n_panels = 10
ncols = 2
cmap = "viridis"

# -----------------------------
# load files
# -----------------------------
files = sorted(glob.glob(f"{csv_dir}/step_*.csv"))
assert len(files) >= n_panels, "Not enough CSV files"

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
    figsize=(10, 14),   # tall figure
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
    "Shock–Bubble Interaction: Density Evolution",
    fontsize=16
)

plt.savefig(f"{csv_dir}/density_plot.png")
