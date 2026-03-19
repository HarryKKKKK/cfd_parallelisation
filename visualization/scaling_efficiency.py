#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
    resolutions = [
        (500, 197),
        (1000, 394),
        (1500, 591),
        (2000, 788)
    ]

    all_data = []

    print("loading data...")
    for nx, ny in resolutions:
        filename = f"results/resolution/res_{nx}_raw.csv"
        if not os.path.exists(filename):
            print(f"Error: cannot find file {filename}")
            continue

        df = pd.read_csv(filename)
        df = df[df['status'] == 'ok'].copy()

        df['res_str'] = f"{nx}x{ny}"

        serial_data = df[(df['mode'] == 'serial') & (df['p'] == 1)]
        if not serial_data.empty:
            t_serial = serial_data['wall_seconds'].mean()
        else:
            t_serial = df[df['p'] == 1]['wall_seconds'].mean()
            print(f"Warning: {filename} lacks serial data, fallback to p=1 average")

        df['efficiency'] = t_serial / (df['p'] * df['wall_seconds'])
        all_data.append(df)

    if not all_data:
        print("No data found")
        return

    master_df = pd.concat(all_data, ignore_index=True)
    res_labels = [f"{nx}x{ny}" for nx, ny in resolutions]
    target_cores = [4, 8, 16]

    color_mpi = '#377eb8'
    color_omp = '#ff7f00'

    for p in target_cores:
        fig, ax = plt.subplots(figsize=(8, 4.5))

        mpi_df = master_df[(master_df['mode'] == 'mpi') & (master_df['p'] == p)]
        mpi_means = mpi_df.groupby('res_str')['efficiency'].mean().reindex(res_labels)

        omp_df = master_df[(master_df['mode'] == 'omp') & (master_df['p'] == p)]
        omp_means = omp_df.groupby('res_str')['efficiency'].mean().reindex(res_labels)

        ax.plot(
            res_labels, mpi_means,
            marker='o', label='MPI',
            color=color_mpi, linewidth=1.5, markersize=6
        )
        ax.plot(
            res_labels, omp_means,
            marker='s', label='OpenMP',
            color=color_omp, linewidth=1.5, markersize=6
        )

        ax.set_title(f"Parallel efficiency vs resolution ({p} cores)", fontsize=14)
        ax.set_xlabel("Resolution", fontsize=12)
        ax.set_ylabel("Parallel efficiency", fontsize=12)
        ax.set_ylim(0, 1.05)
        ax.grid(True, linestyle='-', alpha=0.3)
        ax.legend(loc='lower left', fontsize=11)

        plt.tight_layout()

        output_img = f'efficiency_vs_resolution_{p}cores.png'
        plt.savefig(output_img, dpi=300, bbox_inches='tight')
        print(f"plot saved: {output_img}")
        plt.close(fig)

if __name__ == '__main__':
    main()