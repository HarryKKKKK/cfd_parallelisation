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
C = 200.0          # scaling constant in exp(-|grad rho| / C)
cmap = "gray"

# =============================
# load csv files
# =============================
files = sorted(glob.glob(f"{csv_dir}/step_*.csv"))
indices = np.linspace(0, len(files) - 1, n_panels, dtype=int)
files = [files[i] for i in indices]

# =============================
# grid info
# =============================
df0 = pd.read_csv(files[0])
xs = np.unique(df0["x"].values)
ys = np.unique(df0["y"].values)

nx, ny = len(xs), len(ys)
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

    # gradient magnitude
    drho_dx = np.gradient(rho, dx, axis=1)
    drho_dy = np.gradient(rho, dy, axis=0)
    grad_rho = np.sqrt(drho_dx**2 + drho_dy**2)

    # exponential indicator
    M = np.exp(-grad_rho / C)

    im = ax.imshow(
        M,
        origin="lower",
        extent=[xs.min(), xs.max(), ys.min(), ys.max()],
        aspect="equal",
        cmap=cmap,
        vmin=0.0,
        vmax=1.0
    )

    label = chr(97 + k)
    t_str = fname.split("_t_")[-1].replace(".csv", "")
    ax.set_title(f"({label})  t = {t_str}", fontsize=10)

    ax.set_xticks([])
    ax.set_yticks([])

for ax in axes[len(files):]:
    ax.axis("off")

cbar = fig.colorbar(im, ax=axes, shrink=0.9)
cbar.set_label(r"$M = \exp(-|\nabla \rho| / 200)$")

fig.suptitle(
    "Shock–Bubble Interaction: Exponential Density-Gradient Indicator",
    fontsize=16
)

plt.savefig(f"{csv_dir}/m_indicator.png", dpi=300)
