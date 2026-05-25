"""
Shared stereonet figure builder — MTEX-style solid fills.

All anisotropy stereonets in SAnTex call ``make_stereonet_figure``.
Three visual styles are supported:

  "Smooth fill"     go.Heatmap(zsmooth='best')            — seamless colour wash
  "Filled contours" go.Contour(coloring='heatmap')        — MTEX-default look
  "Scatter (dots)"  go.Scattergl with markers             — fast scatter

**Projection**: Lambert equal-area (Schmidt net), matching MTEX default.
  r = √2 · sin(θ/2),  where θ is the polar angle from the upper pole.
  The equator maps to r = 1; areas are correctly preserved.

Every style gets:
  • Solid circular border
  • Dashed inner rings at 30 ° and 60 ° from pole
  • Tick marks every 30 ° on the boundary
  • N / E / S / W compass labels
"""

from __future__ import annotations
import numpy as np
import plotly.graph_objects as go

# Exported for combo-box population
STEREONET_STYLES: list[str] = [
    "Smooth fill",
    "Filled contours",
    "Scatter (dots)",
]

# ── helpers ──────────────────────────────────────────────────────────────────

def _circle_xy(n: int = 300):
    t = np.linspace(0, 2 * np.pi, n)
    return np.cos(t), np.sin(t)


def _ring_r(theta_deg: float) -> float:
    """Lambert equal-area radius for a small-circle at *theta_deg* from the pole.

    r = √2 · sin(θ/2)   — matches MTEX's default Schmidt-net projection.
    """
    return np.sqrt(2.0) * np.sin(np.deg2rad(theta_deg) / 2.0)


def _add_stereonet_decorations(fig: go.Figure) -> None:
    """Add compass labels, grid rings, and tick marks to *fig* in-place."""

    # ── compass labels ────────────────────────────────────────────────────
    for label, px, py in [("N", 0.0, 1.13), ("E", 1.13, 0.0),
                           ("S", 0.0, -1.13), ("W", -1.13, 0.0)]:
        fig.add_annotation(
            x=px, y=py, text=f"<b>{label}</b>",
            showarrow=False,
            font=dict(size=12, color="black"),
            xref="x", yref="y",
        )

    # ── inner dashed rings ────────────────────────────────────────────────
    for deg in (30, 60):
        r = _ring_r(deg)
        fig.add_shape(
            type="circle", xref="x", yref="y",
            x0=-r, y0=-r, x1=r, y1=r,
            line=dict(color="rgba(0,0,0,0.30)", width=0.8, dash="dot"),
            layer="above",
        )

    # ── outer boundary circle ─────────────────────────────────────────────
    cx, cy = _circle_xy()
    fig.add_trace(go.Scatter(
        x=cx, y=cy,
        mode="lines",
        line=dict(color="black", width=2),
        hoverinfo="skip",
        showlegend=False,
        name="__border__",
    ))

    # ── tick marks every 30 ° ─────────────────────────────────────────────
    for deg in range(0, 360, 30):
        a = np.deg2rad(deg)
        fig.add_shape(
            type="line", xref="x", yref="y",
            x0=0.93 * np.cos(a), y0=0.93 * np.sin(a),
            x1=1.00 * np.cos(a), y1=1.00 * np.sin(a),
            line=dict(color="black", width=1),
            layer="above",
        )

    # ── cross-hairs (faint) ───────────────────────────────────────────────
    for x0, y0, x1, y1 in [(-1, 0, 1, 0), (0, -1, 0, 1)]:
        fig.add_shape(
            type="line", xref="x", yref="y",
            x0=x0, y0=y0, x1=x1, y1=y1,
            line=dict(color="rgba(0,0,0,0.18)", width=0.8, dash="dot"),
            layer="below",
        )

    # ── projection label ─────────────────────────────────────────────────
    fig.add_annotation(
        x=-1.18, y=-1.18,
        text="<i>equal-area</i>",
        showarrow=False,
        font=dict(size=9, color="grey"),
        xref="x", yref="y",
    )


def _base_layout(title: str) -> dict:
    return dict(
        title=dict(text=title, font=dict(size=13), x=0.5, xanchor="center"),
        plot_bgcolor="white",
        margin=dict(l=40, r=40, t=50, b=40),
        xaxis=dict(
            range=[-1.25, 1.25], showgrid=False, zeroline=False,
            showticklabels=False, fixedrange=True,
            scaleanchor="y",
        ),
        yaxis=dict(
            range=[-1.25, 1.25], showgrid=False, zeroline=False,
            showticklabels=False, fixedrange=True,
        ),
    )


# ── public API ────────────────────────────────────────────────────────────────

def make_stereonet_figure(
    data: dict,
    scalar: str,
    style: str = "Smooth fill",
    colorscale: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    show_colorbar: bool = True,
    n_contours: int = 10,
    title: str = "",
    # scatter-mode only
    pt_size: int = 3,
) -> go.Figure:
    """
    Build a Plotly stereonet figure.

    Parameters
    ----------
    data : dict
        From ``AnisotropyBackend.compute_stereonet_grid`` (keys xi, yi,
        plus velocity arrays) for "Smooth fill" / "Filled contours", **or**
        from ``compute_stereonet_data`` (keys x, y, velocity arrays) for
        "Scatter (dots)".
    scalar : str
        Key inside *data* to plot (e.g. "vp", "avs").
    style : str
        One of :data:`STEREONET_STYLES`.
    """

    unit  = "km/s" if scalar in ("vp", "vs1", "vs2") else (
            "%" if scalar == "avs" else "")
    vals_raw = data[scalar]
    if unit == "km/s":
        vals_raw = vals_raw / 1000.0

    _vmin = float(np.nanmin(vals_raw)) if vmin is None else vmin
    _vmax = float(np.nanmax(vals_raw)) if vmax is None else vmax

    cb_spec = dict(
        title=dict(text=f"{scalar} ({unit})", font=dict(size=11)),
        thickness=14, len=0.75,
    ) if show_colorbar else None

    fig = go.Figure()

    # ── style-specific trace ──────────────────────────────────────────────
    if style == "Smooth fill":
        fig.add_trace(go.Heatmap(
            x=data["xi"], y=data["yi"], z=vals_raw,
            zsmooth="best",
            zmin=_vmin, zmax=_vmax,
            colorscale=colorscale,
            showscale=show_colorbar,
            colorbar=cb_spec,
            hovertemplate=(
                f"x=%{{x:.3f}}<br>y=%{{y:.3f}}<br>"
                f"{scalar}=%{{z:.3f}} {unit}<extra></extra>"
            ),
            name=scalar,
        ))

    elif style == "Filled contours":
        fig.add_trace(go.Contour(
            x=data["xi"], y=data["yi"], z=vals_raw,
            zmin=_vmin, zmax=_vmax,
            colorscale=colorscale,
            ncontours=n_contours,
            contours=dict(
                coloring="heatmap",
                showlabels=True,
                labelfont=dict(size=9, color="white"),
            ),
            line=dict(width=0.8, color="rgba(0,0,0,0.5)"),
            showscale=show_colorbar,
            colorbar=cb_spec,
            connectgaps=False,
            hovertemplate=(
                f"x=%{{x:.3f}}<br>y=%{{y:.3f}}<br>"
                f"{scalar}=%{{z:.3f}} {unit}<extra></extra>"
            ),
            name=scalar,
        ))

    else:  # "Scatter (dots)" — legacy / raw-data path
        x_sc   = data["x"]
        y_sc   = data["y"]
        vals_sc = vals_raw  # already divided above
        fig.add_trace(go.Scattergl(
            x=x_sc, y=y_sc,
            mode="markers",
            marker=dict(
                color=vals_sc,
                colorscale=colorscale,
                cmin=_vmin, cmax=_vmax,
                size=pt_size,
                colorbar=cb_spec,
                showscale=show_colorbar,
            ),
            hovertemplate=(
                f"x=%{{x:.3f}}<br>y=%{{y:.3f}}<br>"
                f"{scalar}=%{{marker.color:.3f}} {unit}<extra></extra>"
            ),
            name=scalar,
        ))

    # ── decorations & layout ──────────────────────────────────────────────
    _add_stereonet_decorations(fig)
    fig.update_layout(**_base_layout(title))

    return fig
