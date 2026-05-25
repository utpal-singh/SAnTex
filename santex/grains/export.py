"""
Export grain analysis results to various formats.

Formats
-------
  CSV          : summary DataFrame → CSV
  Neper .ori   : grain orientations for Neper FEM tessellation software
  Neper .tess  : simplified tessellation (centroid-based Voronoi; requires Neper)
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import os


def to_csv(grains, filepath: str, extra_cols: list[str] | None = None) -> None:
    """
    Export grain summary statistics to a CSV file.

    Parameters
    ----------
    grains     : Grain2D object
    filepath   : output CSV path
    extra_cols : additional columns from grains.data to include (per-pixel, averaged per grain)
    """
    df = grains.summary_df()

    if extra_cols:
        for col in extra_cols:
            if col in grains.data.columns:
                vals = []
                for gid in grains._valid_ids:
                    px = grains._pixel_idx[int(gid)]
                    vals.append(grains.data[col].iloc[px].mean())
                df[col + "_mean"] = vals

    df.to_csv(filepath, index=False)
    print(f"Exported {len(df)} grains → {filepath}")


def to_neper_ori(grains, filepath: str) -> None:
    """
    Write a Neper .ori file (grain orientations in Bunge Euler angles, degrees).

    Format per line: phi1  PHI  phi2  (space-delimited, degrees)
    Neper grain numbering starts at 1 and must match tessellation.
    """
    me = grains.mean_euler   # (n_grains, 3) degrees
    lines = []
    for i in range(grains.n_grains):
        lines.append(f"{me[i, 0]:.4f}  {me[i, 1]:.4f}  {me[i, 2]:.4f}")
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Neper .ori written ({grains.n_grains} grains) → {filepath}")


def to_neper_tess(grains, filepath: str) -> None:
    """
    Write a minimal Neper .tess file describing a 2-D tessellation.

    This is a *centroid-based* approximation — it correctly encodes grain
    count, phase assignment, and orientation, but the cell geometry is derived
    from pixel-based centroids rather than exact polygonal boundaries.
    Sufficient for Neper to assign properties without re-tessellating.
    """
    n = grains.n_grains
    centroids = grains.centroid      # (n, 2)
    areas     = grains.area          # (n,)
    phase_ids = grains.phase_id      # (n,)
    me        = grains.mean_euler    # (n, 3)

    x_min = grains.data['X'].min()
    x_max = grains.data['X'].max()
    y_min = grains.data['Y'].min()
    y_max = grains.data['Y'].max()

    lines = ["***tess",
             " **format",
             " 2.0",
             " **general",
             f" 2  {n}",
             " **cell",
             f" {n}",
             " *id"]
    lines.append(" ".join(str(i + 1) for i in range(n)))
    lines += [" *ori",
              " euler-bunge"]
    for i in range(n):
        lines.append(f"  {me[i, 0]:.4f}  {me[i, 1]:.4f}  {me[i, 2]:.4f}")
    lines += [" *phase"]
    lines.append(" ".join(str(p) for p in phase_ids))
    lines.append("***end")

    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Neper .tess written ({n} grains) → {filepath}")


def to_ctf(grains, filepath: str) -> None:
    """
    Re-export EBSD data with grain_id column appended, in CTF-style CSV format.
    Useful for downstream tools that need grain labels attached to pixel data.
    """
    df = grains.data.copy()
    df["GrainID"] = grains.grain_ids
    df.to_csv(filepath, index=False)
    print(f"EBSD + GrainID exported → {filepath}")
