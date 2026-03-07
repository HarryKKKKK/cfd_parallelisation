#!/usr/bin/env python3
"""
Visualize strong scaling results from a CSV file.

Usage:
    python plot_strong_scaling.py /path/to/results.csv

Input CSV columns expected:
    mode,p,run_id,nx,ny,t_end,wall_seconds

Outputs (saved in the same directory as the CSV):
    strong_scaling_time.png
    strong_scaling_speedup.png
    strong_scaling_efficiency.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt


def load_and_validate(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    required_cols = {
        "mode", "p", "run_id", "nx", "ny", "t_end", "wall_seconds"
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Clean types
    df["mode"] = df["mode"].astype(str).str.strip().str.lower()
    df["p"] = pd.to_numeric(df["p"], errors="raise").astype(int)
    df["wall_seconds"] = pd.to_numeric(df["wall_seconds"], errors="raise")

    # Keep only positive timing values
    df = df[df["wall_seconds"] > 0].copy()
    if df.empty:
        raise ValueError("No valid positive wall_seconds values found.")

    return df


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["mode", "p"], as_index=False)
        .agg(
            mean_wall=("wall_seconds", "mean"),
            std_wall=("wall_seconds", "std"),
            n_runs=("wall_seconds", "count"),
        )
        .sort_values(["mode", "p"])
    )

    # Replace NaN std (e.g. only one run) with 0
    summary["std_wall"] = summary["std_wall"].fillna(0.0)
    return summary


def get_serial_baseline(summary: pd.DataFrame) -> float:
    serial_rows = summary[(summary["mode"] == "serial") & (summary["p"] == 1)]
    if serial_rows.empty:
        raise ValueError("Could not find baseline row: mode='serial', p=1")
    return float(serial_rows["mean_wall"].iloc[0])


def add_scaling_metrics(summary: pd.DataFrame, serial_baseline: float) -> pd.DataFrame:
    out = summary.copy()
    out["speedup"] = serial_baseline / out["mean_wall"]
    out["efficiency"] = out["speedup"] / out["p"]
    return out


def plot_time(summary: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 5))

    for mode, group in summary.groupby("mode"):
        group = group.sort_values("p")
        plt.errorbar(
            group["p"],
            group["mean_wall"],
            yerr=group["std_wall"],
            marker="o",
            capsize=4,
            label=mode,
        )

    plt.xlabel("Processes / Threads (p)")
    plt.ylabel("Wall time (s)")
    plt.title("Strong Scaling: Wall Time")
    plt.xticks(sorted(summary["p"].unique()))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_speedup(summary: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 5))

    p_values = sorted(summary["p"].unique())
    plt.plot(p_values, p_values, linestyle="--", label="Ideal")

    for mode, group in summary.groupby("mode"):
        group = group.sort_values("p")
        plt.plot(
            group["p"],
            group["speedup"],
            marker="o",
            label=mode,
        )

    plt.xlabel("Processes / Threads (p)")
    plt.ylabel("Speedup")
    plt.title("Strong Scaling: Speedup")
    plt.xticks(p_values)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_efficiency(summary: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 5))

    for mode, group in summary.groupby("mode"):
        group = group.sort_values("p")
        plt.plot(
            group["p"],
            group["efficiency"],
            marker="o",
            label=mode,
        )

    plt.axhline(1.0, linestyle="--", label="Ideal")
    plt.xlabel("Processes / Threads (p)")
    plt.ylabel("Parallel efficiency")
    plt.title("Strong Scaling: Efficiency")
    plt.xticks(sorted(summary["p"].unique()))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot strong scaling results from a CSV file."
    )
    parser.add_argument(
        "csv_path",
        type=Path,
        help="Path to results CSV file",
    )
    args = parser.parse_args()

    try:
        csv_path = args.csv_path.resolve()
        out_dir = csv_path.parent

        df = load_and_validate(csv_path)
        summary = summarize(df)
        serial_baseline = get_serial_baseline(summary)
        summary = add_scaling_metrics(summary, serial_baseline)

        # Save summary table too, useful for debugging/reporting
        summary_csv = out_dir / "strong_scaling_summary.csv"
        summary.to_csv(summary_csv, index=False)

        plot_time(summary, out_dir / "strong_scaling_time.png")
        plot_speedup(summary, out_dir / "strong_scaling_speedup.png")
        plot_efficiency(summary, out_dir / "strong_scaling_efficiency.png")

        print(f"Loaded: {csv_path}")
        print(f"Baseline serial time (p=1): {serial_baseline:.6f} s")
        print(f"Saved summary: {summary_csv}")
        print(f"Saved figure: {out_dir / 'strong_scaling_time.png'}")
        print(f"Saved figure: {out_dir / 'strong_scaling_speedup.png'}")
        print(f"Saved figure: {out_dir / 'strong_scaling_efficiency.png'}")
        return 0

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())