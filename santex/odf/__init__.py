"""
santex.odf — Orientation Distribution Function and Pole Figure module.

Public API
----------
ODF            : kernel density ODF from EBSD orientations
FibreODF       : model fibre ODF
odf_from_ebsd  : convenience constructor

Pole figures / IPF:
  calc_pole_figure        : PF vectors from orientations
  calc_inverse_pole_figure: IPF crystal directions
  ipf_map_colors          : per-pixel IPF RGB colours for EBSD maps
  stereo_project          : stereographic projection

Components / properties:
  component_table   : volume fractions for standard components
  volume_fraction   : fraction near a single component

ODF reconstruction:
  PoleFigureReconstruction : iterative WIMV from experimental PFs

Export:
  to_popla, to_mtex_txt, to_vpsc, pole_figure_to_csv
"""

from .odf import ODF, FibreODF, odf_from_ebsd
from .pole_figure import (
    calc_pole_figure,
    calc_inverse_pole_figure,
    ipf_map_colors,
    stereo_project,
    sphere_kde,
    pole_figure_grid,
    ipf_color,
)
from .components import component_table, volume_fraction, get_components
from .kernel import (
    de_la_vallee_poussin,
    halfwidth_to_kappa,
    kappa_to_halfwidth,
)
from .reconstruction import PoleFigureReconstruction
from .export import to_popla, to_mtex_txt, to_vpsc, pole_figure_to_csv

__all__ = [
    "ODF", "FibreODF", "odf_from_ebsd",
    "calc_pole_figure", "calc_inverse_pole_figure",
    "ipf_map_colors", "stereo_project", "sphere_kde", "pole_figure_grid",
    "ipf_color",
    "component_table", "volume_fraction", "get_components",
    "de_la_vallee_poussin", "halfwidth_to_kappa", "kappa_to_halfwidth",
    "PoleFigureReconstruction",
    "to_popla", "to_mtex_txt", "to_vpsc", "pole_figure_to_csv",
]
