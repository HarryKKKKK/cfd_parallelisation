#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ================= 固定配置 =================
TARGETS = [0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 9.4, 12.6, 17.4]

NX, NY = 500, 197
X_LEN, Y_LEN = 0.225, 0.089

LEVELS = np.linspace(0.1, 2.8, 45)

# 物理时间终点（和 main 一致）
T_END_PHYS = 0.0011741
T_END_DIM  = 19.0

PATTERN = "step_*_t_*.csv"
OUT_NAME = "density_landscape.png"

# ---- 虚线小球（初始 bubble） ----
R_BUBBLE = 0.025
BUBBLE_CX = 0.035
BUBBLE_CY = 0.0445
DRAW_BUBBLE = True
BUBBLE_LINEWIDTH = 0.8

# ---- domain origin ----
X0, Y0 = 0.0, 0.0

# ---- 布局设置 ----
N_COLS = 2          # 可改成 3
TITLE_FONTSIZE = 12
CONTOUR_LINEWIDTH = 0.5
SPINE_LINEWIDTH = 1.0
# ===============================


def parse_time(fname: str) -> float:
    base = os.path.basename(fname)
    t_str = base.split("_t_")[-1].replace(".csv", "")
    return float(t_str)


def main(indir: str = "res/serial/hllc"):
    out_path = os.path.join(indir, OUT_NAME)

    files = sorted(glob.glob(os.path.join(indir, PATTERN)))
    if not files:
        raise SystemExit(f"No files found: {os.path.join(indir, PATTERN)}")

    # 映射：t_dim = t_phys / t_ref
    t_ref = T_END_PHYS / T_END_DIM
    t_phys_all = np.array([parse_time(f) for f in files], dtype=float)
    t_dim_all = t_phys_all / t_ref

    def load_rho_by_dim_t(t_dim_target: float):
        idx = int(np.argmin(np.abs(t_dim_all - t_dim_target)))
        df = pd.read_csv(files[idx])

        # 如果你的输出 reshape 方向不对，可试：
        # rho = df["rho"].to_numpy(dtype=float).reshape((NX, NY)).T
        rho = df["rho"].to_numpy(dtype=float).reshape((NY, NX))
        return rho, files[idx], t_phys_all[idx], t_dim_all[idx]

    # ---------------------------------------------------------
    # meshgrid: FV cell-centers
    # ---------------------------------------------------------
    dx = X_LEN / NX
    dy = Y_LEN / NY
    xs = X0 + (np.arange(NX) + 0.5) * dx
    ys = Y0 + (np.arange(NY) + 0.5) * dy
    X, Y = np.meshgrid(xs, ys)

    # 横向布局
    n_plots = len(TARGETS)
    n_cols = N_COLS
    n_rows = math.ceil(n_plots / n_cols)

    # 2列时比较适合 A4 横板/论文横向图
    fig_width = 12 if n_cols == 2 else 15
    fig_height = 2.4 * n_rows
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))

    # 统一转成一维，方便索引
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    else:
        axes = np.array(axes).reshape(-1)

    # 预计算 bubble 虚线圆
    if DRAW_BUBBLE:
        theta = np.linspace(0, 2 * np.pi, 400)
        xb = BUBBLE_CX + R_BUBBLE * np.cos(theta)
        yb = BUBBLE_CY + R_BUBBLE * np.sin(theta)

    for i, t_dim_target in enumerate(TARGETS):
        ax = axes[i]
        rho, chosen_file, chosen_t_phys, chosen_t_dim = load_rho_by_dim_t(t_dim_target)

        ax.contour(X, Y, rho, levels=LEVELS, colors="k", linewidths=CONTOUR_LINEWIDTH)

        if DRAW_BUBBLE:
            ax.plot(xb, yb, "k--", linewidth=BUBBLE_LINEWIDTH)

        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"t={t_dim_target:g}", fontsize=TITLE_FONTSIZE, pad=3)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(SPINE_LINEWIDTH)
            spine.set_color("black")

        print(
            f"[{i}] target t_dim={t_dim_target:g} -> chosen t_dim={chosen_t_dim:g}, "
            f"t_phys={chosen_t_phys:.6g}, file={os.path.basename(chosen_file)}"
        )

    # 删除多余空子图
    for j in range(n_plots, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout(pad=0.6, w_pad=0.8, h_pad=0.8)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"图片已生成：{out_path}")


if __name__ == "__main__":
    import sys
    indir = sys.argv[1] if len(sys.argv) > 1 else "res/serial/hllc"
    main(indir)