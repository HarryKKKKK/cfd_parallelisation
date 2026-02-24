import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ================= 固定配置 =================
TARGETS = [0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 9.4, 12.6, 17.4, 19.0]

NX, NY = 500, 197
X_LEN, Y_LEN = 0.225, 0.089

# 时间映射（和 density 一样）
T_END_PHYS = 0.0011741
T_END_DIM  = 19.0

PATTERN = "step_*_t_*.csv"
OUT_NAME = "pressure.png"

# ---- bubble（你的 setting） ----
R_BUBBLE  = 0.025
BUBBLE_CX = 0.035
BUBBLE_CY = 0.0445
DRAW_BUBBLE = True
BUBBLE_LINEWIDTH = 0.8

# =====================================================

def parse_time(fname: str) -> float:
    base = os.path.basename(fname)
    t_str = base.split("_t_")[-1].replace(".csv", "")
    return float(t_str)

def main(indir: str = "res/serial"):
    out_path = os.path.join(indir, OUT_NAME)

    files = sorted(glob.glob(os.path.join(indir, PATTERN)))
    if not files:
        raise SystemExit(f"No files found: {os.path.join(indir, PATTERN)}")

    # t_dim 映射
    t_ref = T_END_PHYS / T_END_DIM
    t_phys_all = np.array([parse_time(f) for f in files], dtype=float)
    t_dim_all  = t_phys_all / t_ref

    def load_p_by_dim_t(t_dim_target: float):
        idx = int(np.argmin(np.abs(t_dim_all - t_dim_target)))
        df = pd.read_csv(files[idx])
        p = df["p"].to_numpy(dtype=float).reshape((NY, NX))
        return p, files[idx], t_phys_all[idx], t_dim_all[idx]

    # ---- FV cell-center 坐标 ----
    dx = X_LEN / NX
    dy = Y_LEN / NY
    xs = (np.arange(NX) + 0.5) * dx
    ys = (np.arange(NY) + 0.5) * dy
    X, Y = np.meshgrid(xs, ys)

    # ---- 自动生成合理 contour levels ----
    # 先扫一遍所有帧确定压力范围
    Ps = []
    metas = []
    for t in TARGETS:
        p, f, tphys, tdim = load_p_by_dim_t(t)
        Ps.append(p)
        metas.append((t, f, tphys, tdim))

    pall = np.concatenate([p.ravel() for p in Ps])
    pmin = np.percentile(pall, 1)
    pmax = np.percentile(pall, 99)

    LEVELS = np.linspace(pmin, pmax, 30)

    # ---- 预计算 bubble 圆 ----
    if DRAW_BUBBLE:
        theta = np.linspace(0, 2*np.pi, 400)
        xb = BUBBLE_CX + R_BUBBLE * np.cos(theta)
        yb = BUBBLE_CY + R_BUBBLE * np.sin(theta)

    fig, axes = plt.subplots(len(TARGETS), 1, figsize=(6, 18))
    if len(TARGETS) == 1:
        axes = [axes]

    for i, (p, meta) in enumerate(zip(Ps, metas)):
        t_target, chosen_file, chosen_t_phys, chosen_t_dim = meta
        ax = axes[i]

        ax.contour(X, Y, p, levels=LEVELS, colors="k", linewidths=0.5)

        if DRAW_BUBBLE:
            ax.plot(xb, yb, "k--", linewidth=BUBBLE_LINEWIDTH)

        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"t={t_target:g}", fontsize=12, pad=2)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
            spine.set_color("black")

        print(f"[{i}] target t={t_target:g} -> chosen t={chosen_t_dim:g}, "
              f"T={chosen_t_phys:.6g}, file={os.path.basename(chosen_file)}")

    plt.tight_layout(h_pad=0.2)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"图片已生成：{out_path}")

if __name__ == "__main__":
    import sys
    indir = sys.argv[1] if len(sys.argv) > 1 else "res/serial"
    main(indir)