#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_strong.py <path_to_strong_csv>")
        sys.exit(1)

    csv_path = sys.argv[1]
    if not os.path.exists(csv_path):
        print(f"Error: File '{csv_path}' not found.")
        sys.exit(1)

    # 读取并清洗数据
    df = pd.read_csv(csv_path)
    df = df[df['status'] == 'ok']
    
    # 按 mode 和 p 求平均执行时间
    agg = df.groupby(['mode', 'p'])['wall_seconds'].mean().reset_index()

    # 获取纯串行基准时间 T_serial
    try:
        t_serial = agg[(agg['mode'] == 'serial') & (agg['p'] == 1)]['wall_seconds'].values[0]
    except IndexError:
        print("Error: Missing 'serial' mode data at p=1.")
        sys.exit(1)

    # 【核心修改】：强制将所有模式在 p=1 时的运行时间替换为 serial 的时间
    # 这保证了在 p=1 时，Speedup 绝对为 1.0，Efficiency 绝对为 100%
    agg.loc[agg['p'] == 1, 'wall_seconds'] = t_serial

    # 计算加速比和效率
    agg['speedup'] = t_serial / agg['wall_seconds']
    agg['efficiency'] = (agg['speedup'] / agg['p']) * 100

    cores = sorted(agg['p'].unique())

    # ================== 图 1: Speedup (加速比) ==================
    plt.figure(figsize=(8, 6))
    for m, marker, color in zip(['mpi', 'omp'], ['o', 's'], ['#1f77b4', '#ff7f0e']):
        data = agg[agg['mode'] == m]
        if not data.empty:
            plt.plot(data['p'], data['speedup'], marker=marker, label=f'{m.upper()} Speedup', color=color, linewidth=2)
    
    plt.plot(cores, cores, 'k--', label='Ideal Speedup', alpha=0.7)
    plt.title('Strong Scaling: Speedup vs Cores', fontsize=14)
    plt.xlabel('Number of Cores (p)', fontsize=12)
    plt.ylabel('Speedup (S)', fontsize=12)
    plt.xticks(cores)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('strong_speedup.png', dpi=300)
    print("Saved -> strong_speedup.png")

    # ================== 图 2: Efficiency (并行效率) ==================
    plt.figure(figsize=(8, 6))
    for m, marker, color in zip(['mpi', 'omp'], ['o', 's'], ['#1f77b4', '#ff7f0e']):
        data = agg[agg['mode'] == m]
        if not data.empty:
            plt.plot(data['p'], data['efficiency'], marker=marker, label=f'{m.upper()} Efficiency', color=color, linewidth=2)
    
    plt.axhline(y=100, color='k', linestyle='--', label='Ideal Efficiency (100%)', alpha=0.7)
    plt.title('Strong Scaling: Parallel Efficiency vs Cores', fontsize=14)
    plt.xlabel('Number of Cores (p)', fontsize=12)
    plt.ylabel('Efficiency (%)', fontsize=12)
    plt.ylim(0, 110)
    plt.xticks(cores)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('strong_efficiency.png', dpi=300)
    print("Saved -> strong_efficiency.png")

if __name__ == '__main__':
    main()