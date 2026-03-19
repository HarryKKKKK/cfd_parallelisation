#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_weak.py <path_to_weak_csv>")
        sys.exit(1)

    csv_path = sys.argv[1]
    if not os.path.exists(csv_path):
        print(f"Error: File '{csv_path}' not found.")
        sys.exit(1)

    # 读取并清洗数据
    df = pd.read_csv(csv_path)
    df = df[df['status'] == 'ok']
    
    agg = df.groupby(['mode', 'p'])['wall_seconds'].mean().reset_index()

    # 获取纯串行基准时间 T_serial
    try:
        t_serial = agg[(agg['mode'] == 'serial') & (agg['p'] == 1)]['wall_seconds'].values[0]
    except IndexError:
        print("Error: Missing 'serial' mode data at p=1.")
        sys.exit(1)

    # 【核心修改】：强制将所有模式在 p=1 时的运行时间替换为 serial 的时间
    # 保证弱扩展性的基准线完美平齐
    agg.loc[agg['p'] == 1, 'wall_seconds'] = t_serial

    # 计算弱扩展性效率: E_weak = T1 / Tp * 100%
    agg['efficiency'] = (t_serial / agg['wall_seconds']) * 100

    cores = sorted(agg['p'].unique())

    # ================== 图 1: Wall-clock Time (执行时间) ==================
    plt.figure(figsize=(8, 6))
    for m, marker, color in zip(['mpi', 'omp'], ['o', 's'], ['#1f77b4', '#ff7f0e']):
        data = agg[agg['mode'] == m]
        if not data.empty:
            plt.plot(data['p'], data['wall_seconds'], marker=marker, label=f'{m.upper()} Time', color=color, linewidth=2)
    
    plt.axhline(y=t_serial, color='k', linestyle='--', label=f'Ideal Time ({t_serial:.1f}s)', alpha=0.7)
    plt.title('Weak Scaling: Execution Time vs Cores', fontsize=14)
    plt.xlabel('Number of Cores (p) & Problem Size Multiplier', fontsize=12)
    plt.ylabel('Wall-clock Time (seconds)', fontsize=12)
    plt.xticks(cores)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('weak_walltime.png', dpi=300)
    print("Saved -> weak_walltime.png")

    # ================== 图 2: Weak Efficiency (弱扩展性效率) ==================
    plt.figure(figsize=(8, 6))
    for m, marker, color in zip(['mpi', 'omp'], ['o', 's'], ['#1f77b4', '#ff7f0e']):
        data = agg[agg['mode'] == m]
        if not data.empty:
            plt.plot(data['p'], data['efficiency'], marker=marker, label=f'{m.upper()} Efficiency', color=color, linewidth=2)
    
    plt.axhline(y=100, color='k', linestyle='--', label='Ideal Efficiency (100%)', alpha=0.7)
    plt.title('Weak Scaling: Parallel Efficiency vs Cores', fontsize=14)
    plt.xlabel('Number of Cores (p) & Problem Size Multiplier', fontsize=12)
    plt.ylabel('Weak Efficiency (%)', fontsize=12)
    plt.ylim(0, 110)
    plt.xticks(cores)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('weak_efficiency.png', dpi=300)
    print("Saved -> weak_efficiency.png")

if __name__ == '__main__':
    main()