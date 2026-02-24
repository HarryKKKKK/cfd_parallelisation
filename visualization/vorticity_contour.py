import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ================= 固定配置 =================
TARGETS = [0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 9.4, 12.6, 17.4, 19.0]

NX, NY = 500, 197
X_LEN, Y_LEN = 0.225, 0.089

# 物理时间终点（和你 main 一致）
T_END_PHYS = 0.0011741
T_END_DIM  = 19.0

PATTERN  = "step_*_t_*.csv"
OUT_NAME = "vorticity.png"

# ---- 虚线小球（初始 bubble） ----
# setting: Bubble radius 0.025m, centre (0.035m, 0.0445m)
R_BUBBLE  = 0.025
BUBBLE_CX = 0.035
BUBBLE_CY = 0.0445
DRAW_BUBBLE = True
BUBBLE_LINEWIDTH = 0.8

# ---- Vorticity contour config ----
N_LEVELS_EACH_SIDE = 8
ROBUST_PCTL = 99.5
LINEWIDTH = 0.6
# ===============================

def parse_time(fname: str) -> float:
    base = os.path.basename(fname)
    t_str = base.split("_t_")[-1].replace(".csv", "")
    return float(t_str)

def pick_column(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def load_uv(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    u_col = pick_column(df, ["u", "U", "velx", "vx"])
    v_col = pick_column(df, ["v", "V", "vely", "vy"])

    if u_col is not None and v_col is not None:
        u = df[u_col].to_numpy(dtype=float).reshape((NY, NX))
        v = df[v_col].to_numpy(dtype=float).reshape((NY, NX))
        return u, v

    rho_col  = pick_column(df, ["rho", "density"])
    rhou_col = pick_column(df, ["rhou", "rho_u", "mx", "momx"])
    rhov_col = pick_column(df, ["rhov", "rho_v", "my", "momy"])
    if rho_col is None or rhou_col is None or rhov_col is None:
        raise KeyError(
            "Need columns either (u,v) or (rho,rhou,rhov).\n"
            f"Found: {list(df.columns)}"
        )

    rho  = df[rho_col].to_numpy(dtype=float).reshape((NY, NX))
    rhou = df[rhou_col].to_numpy(dtype=float).reshape((NY, NX))
    rhov = df[rhov_col].to_numpy(dtype=float).reshape((NY, NX))

    rho_safe = np.maximum(rho, 1e-14)
    u = rhou / rho_safe
    v = rhov / rho_safe
    return u, v

def compute_vorticity(u: np.ndarray, v: np.ndarray, dx: float, dy: float) -> np.ndarray:
    dv_dy, dv_dx = np.gradient(v, dy, dx, edge_order=1)
    du_dy, du_dx = np.gradient(u, dy, dx, edge_order=1)
    return dv_dx - du_dy

def main(indir: str = "res/serial"):
    out_path = os.path.join(indir, OUT_NAME)

    files = sorted(glob.glob(os.path.join(indir, PATTERN)))
    if not files:
        raise SystemExit(f"No files found: {os.path.join(indir, PATTERN)}")

    # 映射：t_dim = t_phys / t_ref
    t_ref = T_END_PHYS / T_END_DIM
    t_phys_all = np.array([parse_time(f) for f in files], dtype=float)
    t_dim_all  = t_phys_all / t_ref

    # meshgrid（保持和你的 density 脚本一致）
    xs = np.linspace(0, X_LEN, NX)
    ys = np.linspace(0, Y_LEN, NY)
    X, Y = np.meshgrid(xs, ys)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    def load_w_by_dim_t(t_dim_target: float):
        idx = int(np.argmin(np.abs(t_dim_all - t_dim_target)))
        df = pd.read_csv(files[idx])
        u, v = load_uv(df)
        w = compute_vorticity(u, v, dx, dy)
        return w, files[idx], t_phys_all[idx], t_dim_all[idx]

    if DRAW_BUBBLE:
        theta = np.linspace(0, 2*np.pi, 400)
        xb = BUBBLE_CX + R_BUBBLE * np.cos(theta)
        yb = BUBBLE_CY + R_BUBBLE * np.sin(theta)

    # 先算所有帧，用于统一 levels
    Ws, metas = [], []
    for t_dim_target in TARGETS:
        w, f, tphys, tdim = load_w_by_dim_t(t_dim_target)
        Ws.append(w)
        metas.append((t_dim_target, f, tphys, tdim))

    w_all = np.concatenate([w.ravel() for w in Ws])
    wmax = np.percentile(np.abs(w_all), ROBUST_PCTL)
    if not np.isfinite(wmax) or wmax <= 0:
        wmax = np.max(np.abs(w_all)) + 1e-12

    levels_pos = np.linspace(wmax / N_LEVELS_EACH_SIDE, wmax, N_LEVELS_EACH_SIDE)
    levels_neg = -levels_pos[::-1]

    fig, axes = plt.subplots(len(TARGETS), 1, figsize=(6, 18))
    if len(TARGETS) == 1:
        axes = [axes]

    for i, (w, meta) in enumerate(zip(Ws, metas)):
        t_dim_target, chosen_file, chosen_t_phys, chosen_t_dim = meta
        ax = axes[i]

        ax.contour(X, Y, w, levels=levels_pos, colors="k",
                   linewidths=LINEWIDTH, linestyles="solid")
        ax.contour(X, Y, w, levels=levels_neg, colors="k",
                   linewidths=LINEWIDTH, linestyles="dashed")

        if DRAW_BUBBLE:
            ax.plot(xb, yb, "k--", linewidth=BUBBLE_LINEWIDTH)

        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"t={t_dim_target:g}", fontsize=12, pad=2)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
            spine.set_color("black")

        print(f"[{i}] target t_dim={t_dim_target:g} -> chosen t_dim={chosen_t_dim:g}, "
              f"t_phys={chosen_t_phys:.6g}, file={os.path.basename(chosen_file)}")

    plt.tight_layout(h_pad=0.2)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"图片已生成：{out_path}")

if __name__ == "__main__":
    import sys
    indir = sys.argv[1] if len(sys.argv) > 1 else "res/serial"
    main(indir)