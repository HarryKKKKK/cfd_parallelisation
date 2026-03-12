import os
import sys
import pandas as pd
import matplotlib.pyplot as plt


def summarise_strong(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["mode", "p"], as_index=False)
          .agg(
              mean_wall_seconds=("wall_seconds", "mean"),
              std_wall_seconds=("wall_seconds", "std"),
              min_wall_seconds=("wall_seconds", "min"),
              max_wall_seconds=("wall_seconds", "max"),
              runs=("wall_seconds", "count"),
          )
          .sort_values(["mode", "p"])
    )

    baseline = (
        summary[summary["p"] == 1][["mode", "mean_wall_seconds"]]
        .rename(columns={"mean_wall_seconds": "t1"})
    )

    summary = summary.merge(baseline, on="mode", how="left")
    summary["speedup"] = summary["t1"] / summary["mean_wall_seconds"]
    summary["efficiency"] = summary["speedup"] / summary["p"]

    return summary


def plot_runtime(summary: pd.DataFrame, outdir: str) -> None:
    plt.figure(figsize=(7, 5))
    for mode in summary["mode"].unique():
        d = summary[summary["mode"] == mode]
        plt.errorbar(
            d["p"], d["mean_wall_seconds"], yerr=d["std_wall_seconds"],
            marker="o", capsize=4, label=mode
        )
    plt.xlabel("p")
    plt.ylabel("Runtime (s)")
    plt.title("Strong scaling: runtime vs p")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "strong_runtime_vs_p.png"), dpi=200)
    plt.close()


def plot_speedup(summary: pd.DataFrame, outdir: str) -> None:
    plt.figure(figsize=(7, 5))

    p_all = sorted(summary["p"].unique())
    plt.plot(p_all, p_all, linestyle="--", label="ideal")

    for mode in summary["mode"].unique():
        d = summary[summary["mode"] == mode]
        plt.plot(d["p"], d["speedup"], marker="o", label=mode)

    plt.xlabel("p")
    plt.ylabel("Speedup")
    plt.title("Strong scaling: speedup vs p")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "strong_speedup_vs_p.png"), dpi=200)
    plt.close()


def plot_efficiency(summary: pd.DataFrame, outdir: str) -> None:
    plt.figure(figsize=(7, 5))

    for mode in summary["mode"].unique():
        d = summary[summary["mode"] == mode]
        plt.plot(d["p"], d["efficiency"], marker="o", label=mode)

    plt.axhline(1.0, linestyle="--", label="ideal")
    plt.xlabel("p")
    plt.ylabel("Parallel efficiency")
    plt.title("Strong scaling: efficiency vs p")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "strong_efficiency_vs_p.png"), dpi=200)
    plt.close()


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 eval_strong_scaling.py <csv_folder>")
        sys.exit(1)

    folder = sys.argv[1]
    infile = os.path.join(folder, "strong_scaling_raw.csv")

    if not os.path.exists(infile):
        print(f"Error: file not found: {infile}")
        sys.exit(1)

    df = pd.read_csv(infile)
    summary = summarise_strong(df)

    summary_file = os.path.join(folder, "strong_scaling_summary.csv")
    summary.to_csv(summary_file, index=False)

    plot_runtime(summary, folder)
    plot_speedup(summary, folder)
    plot_efficiency(summary, folder)

    print(f"[OK] Summary written to: {summary_file}")
    print(f"[OK] Plots written to: {folder}")


if __name__ == "__main__":
    main()