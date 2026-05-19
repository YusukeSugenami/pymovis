"""
PyMoVis: Molecular Orbital Visualization System

A Python package for visualizing molecular orbitals from quantum chemistry calculations.
Supports Gaussian fchk and cube file formats.

Features:
    - Load and parse molecular data from Gaussian output files
    - Evaluate molecular orbitals on custom grids
    - Visualize MO isosurfaces with 3D rendering using PyVista
    - Customize camera positions and visualization parameters

Example:
    >>> from pymovis import savemo
    >>> savemo("water.fchk", ["HOMO", "LUMO"], iso=0.05)

"""

from .pymovis import cli, savemo, save_mo_plot, build_volume_grid, build_molstruct
from .load_mol import load_inp, MoleculeData
from .camera_utils import (
    infer_atom_axis_camera_fit_from_coords,
    infer_axis_camera_fit_from_coords,
)
from .settings import ELEMENT_DATA, ORBITAL_COLORS

__version__ = "0.1.0"
__author__ = "PyMoVis Contributors"

__all__ = [
    "savemo",
    "cli",
    "save_mo_plot",
    "build_volume_grid",
    "build_molstruct",
    "load_inp",
    "MoleculeData",
    "infer_atom_axis_camera_fit_from_coords",
    "infer_axis_camera_fit_from_coords",
    "ELEMENT_DATA",
    "ORBITAL_COLORS",
]
