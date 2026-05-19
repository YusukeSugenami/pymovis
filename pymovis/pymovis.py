import os
from typing import Any, cast

import numpy as np
import periodictable
import pyvista as pv
import argparse

from .settings import *
from .load_mol import load_inp


##############################################
# PyVista settings
##############################################

pv.OFF_SCREEN = True
# pv.global_theme.off_screen = True
# pv.global_theme.rendering.smooth_shading = True
# pv.global_theme.trame.show = False
# pv.start_xvfb()
pv.global_theme.allow_empty_mesh = True


##############################################
# Visualize with PyVista and save image files
##############################################

def build_molstruct(plotter: pv.Plotter, atomnos: np.ndarray, coords: np.ndarray) -> None:
    """
    Add atoms and bonds to a PyVista plotter.
    Parameters:
        plotter: PyVista Plotter object to which the molecule will be added.
        atomnos: Array of atomic numbers for each atom in the molecule.
        coords: Array of atomic coordinates (shape: [num_atoms, 3]).
    Returns:
        None
    """

    atomnos = np.asarray(atomnos)
    coords = np.asarray(coords, dtype=float)
    num_atoms = len(atomnos)
    symbols = [periodictable.elements[atomnos[i]].symbol for i in range(num_atoms)]

    # Adding atoms
    for atom_index in range(num_atoms):
        position = coords[atom_index]
        symbol = symbols[atom_index]
        data = ELEMENT_DATA.get(symbol, DEFAULT_DATA)
        color = data["color"]
        radius = data["radius"]
        sphere = pv.Sphere(center=position, radius=radius)
        plotter.add_mesh(sphere, color=color, smooth_shading=True)

    # Adding bonds
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            symbol_i = symbols[i]
            symbol_j = symbols[j]
            radius_i = ELEMENT_DATA.get(symbol_i, DEFAULT_DATA)["radius"]
            radius_j = ELEMENT_DATA.get(symbol_j, DEFAULT_DATA)["radius"]
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < 2.0 * (radius_i + radius_j):  # Determine bonds
                start = coords[i]
                end = coords[j]
                center = (start + end) / 2
                height = float(distance)
                direction = end - start
                cyl = pv.Cylinder(
                    center=center,
                    direction=direction,
                    radius=BOND_RADIUS,
                    height=height,
                )
                plotter.add_mesh(cyl, color=BOND_DEFAULT_COLOR, smooth_shading=True)


def build_volume_grid(
    mo_cube: np.ndarray,
    origin: np.ndarray,
    grid_vecs: np.ndarray,
) -> pv.StructuredGrid:
    """
    Build a structured grid for MO isosurface rendering.
    Parameters:
        mo_cube: 3D array of MO values on the grid (shape: [nx, ny, nz]).
        origin: 3D coordinates of the grid origin (shape: [3]).
        grid_vecs: 3x3 array of grid vectors defining the grid spacing and orientation (shape: [3, 3]).
    Returns:
        A PyVista StructuredGrid object containing the grid points and MO values.
    """

    nx, ny, nz = mo_cube.shape
    origin = np.asarray(origin, dtype=float)
    grid_x, grid_y, grid_z = np.meshgrid(
        np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij"
    )
    points = (
        origin[None, None, None, :]
        + grid_x[..., None] * grid_vecs[0]
        + grid_y[..., None] * grid_vecs[1]
        + grid_z[..., None] * grid_vecs[2]
    )
    grid = pv.StructuredGrid(points[..., 0], points[..., 1], points[..., 2])
    grid.point_data["mo"] = mo_cube.flatten(order="F")

    return grid


def save_mo_plot(
    atomnos : np.ndarray,
    coords : np.ndarray,
    mo_cube : np.ndarray,
    origin : np.ndarray,
    grid_vecs : np.ndarray,
    out_name : str,
    iso : float = 0.05,
    camera_pos : list[float] | None = None,
    camera_focal : list[float] | None = None,
    camera_up : list[float] | None = None,
    transparent_background : bool = True,
) -> None:
    """
    Render a MO isosurface and save it to an image file.
    Parameters:
        atomnos: Array of atomic numbers for each atom in the molecule.
        coords: Array of atomic coordinates (shape: [num_atoms, 3]).
        mo_cube: 3D array of MO values on the grid (shape: [nx, ny, nz]).
        origin: 3D coordinates of the grid origin (shape: [3]).
        grid_vecs: 3x3 array of grid vectors defining the grid spacing and orientation (shape: [3, 3]).
        out_name: Output file name for the saved image.
        iso: Isosurface value for rendering the MO (default: 0.05).
        camera_pos: Optional list of 3 floats specifying the camera position (default: None).
        camera_focal: Optional list of 3 floats specifying the camera focal point (default: None).
        camera_up: Optional list of 3 floats specifying the camera up vector (default: None).
        transparent_background: Whether to use a transparent background for the saved image (default: True).
    Returns:
        None
    """

    plotter : Any = pv.Plotter(off_screen=True)
    plotter.set_background(BACKGROUND_COLOR)
    grid = build_volume_grid(mo_cube, origin, grid_vecs)
    positive_iso = cast(pv.PolyData, grid.contour([iso], scalars="mo"))
    negative_iso = cast(pv.PolyData, grid.contour([-iso], scalars="mo"))
    plotter.add_mesh(positive_iso, color=ORBITAL_COLORS["pos"], opacity=ORBITAL_OPACITY)
    plotter.add_mesh(negative_iso, color=ORBITAL_COLORS["neg"], opacity=ORBITAL_OPACITY)
    build_molstruct(plotter, atomnos, coords)

    if camera_pos is None and camera_focal is None and camera_up is None:
        plotter.enable_anti_aliasing(ANTI_ALIASING)
        plotter.screenshot(
            out_name, scale=IMAGE_QUALITY, transparent_background=transparent_background
        )
        print(f"Saved {out_name}\ncamera setting: {plotter.camera_position}")

    elif not camera_pos is None and not camera_focal is None and not camera_up is None:
        plotter.enable_anti_aliasing(ANTI_ALIASING)
        plotter.camera_position = [camera_pos, camera_focal, camera_up]
        plotter.screenshot(
            out_name, scale=IMAGE_QUALITY, transparent_background=transparent_background
        )
        print(f"Saved {out_name}\ncamera setting: {plotter.camera_position}")

    else:
        plotter.screenshot(
            out_name, scale=1, transparent_background=transparent_background
        )
        pos, focus, up = plotter.camera_position

        if camera_pos is not None:
            pos = camera_pos

        if camera_focal is not None:
            focus = camera_focal

        if camera_up is not None:
            up = camera_up

        plotter.enable_anti_aliasing(ANTI_ALIASING)
        plotter.camera_position = [pos, focus, up]
        plotter.screenshot(
            out_name, scale=IMAGE_QUALITY, transparent_background=transparent_background
        )
        print(f"Saved {out_name}\ncamera setting: {plotter.camera_position}")


##############################################
# main
##############################################


def savefig(
    inp_file: str,
    mo_index: list[int | str],
    out_name: list[str] | None = None,
    iso: float = 0.03,
    camera_pos: list[float] | None = None,
    camera_focal: list[float] | None = None,
    camera_up: list[float] | None = None,
    transparent_background: bool = True,
) -> None:
    """
    Load molecular data, evaluate selected orbitals, and save plots.
    Parameters:
        inp_file: Path to the input file (Gaussian fchk or cube) containing molecular data and MO information.
        mo_index: List of MO indices to visualize. Can be integers (starting from 0) or strings (e.g., "HOMO", "LUMO+2").
        out_name: Optional list of output file names for the saved images. If None, names will be generated automatically based on the input file name and MO index.
        iso: Isosurface value for rendering the MO (default: 0.03).
        camera_pos: Optional list of 3 floats specifying the camera position (default: None).
        camera_focal: Optional list of 3 floats specifying the camera focal point (default: None).
        camera_up: Optional list of 3 floats specifying the camera up vector (default: None).
        transparent_background: Whether to use a transparent background for the saved images (default: True).
    Returns:
        None
    """

    # Determine output file names based on input file name and MO index
    basename, ext = os.path.basename(inp_file).rsplit(".", 1)
    # if output name is not provided, generate it based on the input file name and MO index
    if out_name is None: 
        if mo_index[0] == -1:
            out_name = [f"{basename}.{DEFAULT_OUT_SUFFIX}"] 
        else:
            out_name = [f"{basename}_{mo_idx}.{DEFAULT_OUT_SUFFIX}" for mo_idx in mo_index] 
    else:
        if len(mo_index) != len(out_name):
            print("warning: length of mo_index and output_name do not match. set output name automatically")
            out_name = [f"{basename}_{mo_idx}.{DEFAULT_OUT_SUFFIX}" for mo_idx in mo_index]
        else: 
            # if output name is provided and has the same length as mo_index, use the provided names. 
            # If any name does not end with an image suffix, append the default suffix.
            out_name = [f"{name}.{DEFAULT_OUT_SUFFIX}" if not name.endswith(IMAGE_SUFFIXES) else name for name in out_name]
    
    moldata = load_inp(inp_file)
    for mo_idx, image_name in zip(mo_index, out_name):
        moldata.evaluate_mo_on_grid(
            mo_idx,
            padding=PADDING,
            grid_setting=GRID_SETTING,
            grid_size=GRID_SIZE,
            grid_shape=GRID_SHAPE
        )
        if (
            moldata.atomnos is None
            or moldata.coords is None
            or moldata.mo_cube is None
            or moldata.origin is None
            or moldata.grid_vecs is None
        ):
            raise ValueError("molecular grid data is incomplete")
        save_mo_plot(
            moldata.atomnos,
            moldata.coords,
            moldata.mo_cube,
            moldata.origin,
            moldata.grid_vecs,
            image_name,
            iso=iso,
            camera_pos=camera_pos,
            camera_focal=camera_focal,
            camera_up=camera_up,
            transparent_background=transparent_background,
        )
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="MO picture generator from Gaussian fchk or cube file."
    )
    for arg_config in PARSER_ARGS:
        parser.add_argument(*arg_config["args"], **arg_config["kwargs"])
    args = parser.parse_args()
    savefig(
        args.input_file,
        args.mo,
        out_name=args.out,
        iso=args.iso,
        camera_pos=args.camera_pos,
        camera_focal=args.camera_focal,
        camera_up=args.camera_up,
        transparent_background=args.transparent,
    )
