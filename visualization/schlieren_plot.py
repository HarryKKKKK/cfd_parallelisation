import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import math

# =============================
# settings
# =============================
csv_dir = "res/serial"
n_panels = 10
ncols = 2
cmap = "gray"

# =============================
# load csv files
# =============================
files = sorted(glob.glob(f"{csv_dir}/step_*.csv"))
assert len(files) >= n_panels, "Not enough CSV files"

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
    "Shock–Bubble Interaction: Schlieren Evolution",
    fontsize=16
)

# =============================
# save / show
# =============================
# plt.savefig("schlieren_2x5.pdf", dpi=300)
plt.savefig(f"{csv_dir}/schlieren.png", dpi=300)
