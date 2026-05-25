"""
ODF export functions.

Supported formats
-----------------
- POPLA (.epf): Fortran-style format read by PopLA and BEARTEX
- MTEX (.txt):  Tab-delimited phi1, PHI, phi2, ODF columns
- VPSC:         Euler angles + weights for VPSC / MTEX import
"""

from __future__ import annotations
import numpy as np


def to_popla(odf, filepath: str, resolution_deg: float = 5.0) -> None:
    """
    Export ODF to POPLA .epf format.

    Parameters
    ----------
    odf       : ODF instance
    filepath  : output file path (should end in .epf)
    resolution_deg : grid resolution in degrees
    """
    phi1 = np.arange(0, 360, resolution_deg)
    PHI  = np.arange(0, 91,  resolution_deg)
    phi2 = np.arange(0, 360, resolution_deg)

    lines = []
    lines.append("(2I5,F10.4)")    # format line
    lines.append(f" {len(phi2):4d} {len(PHI):4d}")

    for p2 in phi2:
        for P in PHI:
            row_vals = []
            for p1 in phi1:
                grid = np.array([[p1, P, p2]])
                val = odf.evaluate(grid)[0]
                row_vals.append(val)
            # 8 values per line, 10.4f format
            for i in range(0, len(row_vals), 8):
                line_vals = row_vals[i:i + 8]
                lines.append("".join(f"{v:10.4f}" for v in line_vals))

    with open(filepath, "w") as f:
        f.write("\n".join(lines))


def to_mtex_txt(odf, filepath: str, resolution_deg: float = 5.0) -> None:
    """
    Export ODF as a tab-delimited text file (phi1, PHI, phi2, f) readable by MTEX.

    Parameters
    ----------
    odf        : ODF instance
    filepath   : output file path
    resolution_deg : grid resolution in degrees
    """
    phi1 = np.arange(0, 360, resolution_deg)
    PHI  = np.arange(0, 91,  resolution_deg)
    phi2 = np.arange(0, 360, resolution_deg)
    p1_g, P_g, p2_g = np.meshgrid(phi1, PHI, phi2)
    grid = np.stack([p1_g.ravel(), P_g.ravel(), p2_g.ravel()], axis=1)

    f_vals = odf.evaluate(grid, batch_size=500)

    with open(filepath, "w") as fh:
        fh.write("phi1\tPHI\tphi2\tODF\n")
        for (p1v, Pv, p2v), fv in zip(grid, f_vals):
            fh.write(f"{p1v:.2f}\t{Pv:.2f}\t{p2v:.2f}\t{fv:.6f}\n")


def to_vpsc(odf, filepath: str, n_orientations: int = 1000,
            rng=None) -> None:
    """
    Export ODF as a set of discrete orientations in VPSC format.

    Parameters
    ----------
    odf           : ODF instance
    filepath      : output file path
    n_orientations: number of orientations to sample
    rng           : numpy random Generator (optional)
    """
    euler_deg = odf.random_sampling(n_orientations, rng=rng)
    weight    = 1.0 / n_orientations

    with open(filepath, "w") as f:
        f.write("Discrete orientations from santex ODF\n")
        f.write("B  <- Bunge Euler angles\n")
        f.write(f"  {n_orientations}\n")
        for e in euler_deg:
            f.write(f"  {e[0]:10.4f}  {e[1]:10.4f}  {e[2]:10.4f}  {weight:10.8f}\n")


def pole_figure_to_csv(vectors: np.ndarray, filepath: str,
                       weights: np.ndarray | None = None) -> None:
    """Export pole figure vectors to CSV."""
    import pandas as pd
    df = pd.DataFrame(vectors, columns=["x", "y", "z"])
    if weights is not None:
        df["weight"] = weights
    df.to_csv(filepath, index=False)
