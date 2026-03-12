#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import glob
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd


def list_step_files(d: str, pattern: str = "step_*.csv") -> List[str]:
    return sorted(glob.glob(os.path.join(d, pattern)))


def basename_set(files: List[str]) -> set:
    return set(os.path.basename(p) for p in files)


def get_numeric_columns(df: pd.DataFrame, exclude: List[str]) -> List[str]:
    cols = []
    for c in df.columns:
        if c in exclude:
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            cols.append(c)
    return cols


def file_mse(
    f1: str,
    f2: str,
    columns: Optional[List[str]],
    exclude_cols: List[str],
    require_same_shape: bool,
) -> Tuple[float, int]:
    """
    Returns: (mse, n_elements_used)
    """
    df1 = pd.read_csv(f1)
    df2 = pd.read_csv(f2)

    if require_same_shape and df1.shape != df2.shape:
        raise ValueError(f"Shape mismatch: {os.path.basename(f1)}: {df1.shape} vs {df2.shape}")

    # Determine columns to compare
    if columns is None:
        c1 = get_numeric_columns(df1, exclude_cols)
        c2 = get_numeric_columns(df2, exclude_cols)
        cols = [c for c in c1 if c in c2]
    else:
        missing1 = [c for c in columns if c not in df1.columns]
        missing2 = [c for c in columns if c not in df2.columns]
        if missing1 or missing2:
            raise ValueError(
                f"Missing columns in {os.path.basename(f1)}: {missing1}; "
                f"in {os.path.basename(f2)}: {missing2}"
            )
        cols = columns

    if not cols:
        raise ValueError(f"No comparable numeric columns found in {os.path.basename(f1)}")

    # Align by common columns; ensure same length (if shapes differ but require_same_shape=False)
    n = min(len(df1), len(df2))
    a = df1.loc[: n - 1, cols].to_numpy(dtype=np.float64, copy=False)
    b = df2.loc[: n - 1, cols].to_numpy(dtype=np.float64, copy=False)

    diff = a - b
    mse = float(np.mean(diff * diff))
    return mse, diff.size


def main():
    ap = argparse.ArgumentParser(
        description="Compute total MSE across matching step_*.csv files in two directories."
    )
    ap.add_argument("dir1", help="First directory (e.g., res/serial/hllc)")
    ap.add_argument("dir2", help="Second directory (e.g., res/mpi/hllc)")
    ap.add_argument("--pattern", default="step_*.csv", help="Glob pattern (default: step_*.csv)")
    ap.add_argument(
        "--exclude", default="x,y",
        help="Comma-separated columns to exclude when auto-picking numeric columns (default: x,y)"
    )
    ap.add_argument(
        "--cols", default=None,
        help="Comma-separated list of columns to compare explicitly (e.g., rho,rhou,rhov,E). "
             "If not set, compares all numeric cols except excluded."
    )
    ap.add_argument(
        "--require-same-shape", action="store_true",
        help="If set, error out when csv shapes differ. Otherwise compares up to min(nrows)."
    )
    ap.add_argument(
        "--quiet", action="store_true",
        help="If set, only prints totals (no per-file lines)."
    )

    args = ap.parse_args()

    d1 = args.dir1
    d2 = args.dir2

    exclude_cols = [c.strip() for c in args.exclude.split(",") if c.strip()]
    cols = None
    if args.cols is not None:
        cols = [c.strip() for c in args.cols.split(",") if c.strip()]

    files1 = list_step_files(d1, args.pattern)
    files2 = list_step_files(d2, args.pattern)

    b1 = basename_set(files1)
    b2 = basename_set(files2)

    common = sorted(b1 & b2)
    only1 = sorted(b1 - b2)
    only2 = sorted(b2 - b1)

    if not common:
        print("No matching files found.")
        print(f"dir1 matched: {len(files1)} files")
        print(f"dir2 matched: {len(files2)} files")
        return

    total_mse_sum = 0.0
    total_elements = 0
    per_file = []

    for name in common:
        f1 = os.path.join(d1, name)
        f2 = os.path.join(d2, name)
        mse, n_elem = file_mse(
            f1, f2,
            columns=cols,
            exclude_cols=exclude_cols,
            require_same_shape=args.require_same_shape
        )
        per_file.append((name, mse, n_elem))
        total_mse_sum += mse
        total_elements += n_elem

    if not args.quiet:
        print(f"Compared {len(common)} matching files.")
        if only1:
            print(f"[WARN] Only in dir1 ({len(only1)}): {only1[:10]}{' ...' if len(only1)>10 else ''}")
        if only2:
            print(f"[WARN] Only in dir2 ({len(only2)}): {only2[:10]}{' ...' if len(only2)>10 else ''}")
        print()

        for name, mse, n_elem in per_file:
            print(f"{name}: MSE={mse:.6e}  (elements={n_elem})")

        print("\n--- Totals ---")

    print(f"TOTAL_MSE_SUM = {total_mse_sum:.12e}")
    print(f"TOTAL_ELEMENTS_COMPARED = {total_elements}")
    print(f"MEAN_OF_FILE_MSE = {total_mse_sum/len(common):.12e}")


if __name__ == "__main__":
    main()