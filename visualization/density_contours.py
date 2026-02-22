import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

# ================= 配置 =================
targets = [0.6, 1.2, 1.8, 3.0, 4.6, 6.2, 7.8, 12.6, 19.0]

nx, ny = 500, 197
x_len, y_len = 0.225, 0.089

levels = np.linspace(0.1, 2.8, 45)

# 物理时间终点（和你 main 一致）
T_end_phys = 0.0011741
T_end_dim  = 19.0
t_ref = T_end_phys / T_end_dim
# ========================================


# ================= 读取 CSV =================
files = sorted(glob.glob("res/serial/step_*.csv"))

def parse_time(fname):
    t_str = fname.split("_t_")[-1].replace(".csv", "")
    return float(t_str)

t_phys_all = np.array([parse_time(f) for f in files])
t_dim_all  = t_phys_all / t_ref

def load_data(t_dim_target):
    # 找到最接近该无量纲时间的文件
    idx = np.argmin(np.abs(t_dim_all - t_dim_target))
    df = pd.read_csv(files[idx])
    rho = df["rho"].values.reshape((ny, nx))
    return rho
# ============================================


# 创建 meshgrid（⚠️ 这是关键修正）
xs = np.linspace(0, x_len, nx)
ys = np.linspace(0, y_len, ny)
X, Y = np.meshgrid(xs, ys)


fig, axes = plt.subplots(len(targets), 1, figsize=(6, 18))

for i, t in enumerate(targets):
    ax = axes[i]
    rho = load_data(t)

    # ✅ 使用 X, Y 传给 contour（避免 extent 贴边问题）
    ax.contour(X, Y, rho,
               levels=levels,
               colors='k',
               linewidths=0.5)

    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"t={t}", fontsize=12, pad=2)

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_visible(False)

plt.tight_layout(h_pad=0.2)
plt.savefig("Fig4_MySolver_Verification.png",
            dpi=300,
            bbox_inches='tight')

print("图片已生成：Fig4_MySolver_Verification.png")