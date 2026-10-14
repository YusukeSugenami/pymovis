"""PyVista scene construction and image rendering."""

from typing import Any, cast

import numpy as np
import periodictable
import pyvista as pv

from .settings import (
    ANTI_ALIASING,
    BACKGROUND_COLOR,
    BOND_DEFAULT_COLOR,
    BOND_RADIUS,
    DEFAULT_DATA,
    ELEMENT_DATA,
    IMAGE_QUALITY,
    ORBITAL_COLORS,
    ORBITAL_OPACITY,
)

pv.OFF_SCREEN = True
pv.global_theme.allow_empty_mesh = True


def build_molstruct(
    plotter: pv.Plotter,
    atomnos: np.ndarray,
    coords: np.ndarray,
) -> None:
    """Add atom spheres and distance-based bonds to ``plotter``."""
    atomnos = np.asarray(atomnos)
    coords = np.asarray(coords, dtype=float)
    symbols = [
        periodictable.elements[atom_number].symbol for atom_number in atomnos
    ]

    for atom_index, symbol in enumerate(symbols):
        data = ELEMENT_DATA.get(symbol, DEFAULT_DATA)
        sphere = pv.Sphere(
            center=coords[atom_index],
            radius=data["radius"],
        )
        plotter.add_mesh(sphere, color=data["color"], smooth_shading=True)

    for i in range(len(atomnos)):
        for j in range(i + 1, len(atomnos)):
            radius_i = ELEMENT_DATA.get(symbols[i], DEFAULT_DATA)["radius"]
            radius_j = ELEMENT_DATA.get(symbols[j], DEFAULT_DATA)["radius"]
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < 2.0 * (radius_i + radius_j):
                direction = coords[j] - coords[i]
                cylinder = pv.Cylinder(
                    center=(coords[i] + coords[j]) / 2,
                    direction=direction,
                    radius=BOND_RADIUS,
                    height=float(distance),
                )
                plotter.add_mesh(
                    cylinder,
                    color=BOND_DEFAULT_COLOR,
                    smooth_shading=True,
                )


def build_volume_grid(
    mo_cube: np.ndarray,
    origin: np.ndarray,
    grid_vecs: np.ndarray,
) -> pv.StructuredGrid:
    """Create a PyVista grid from scalar values and grid basis vectors."""
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
    atomnos: np.ndarray,
    coords: np.ndarray,
    mo_cube: np.ndarray | None,
    origin: np.ndarray | None,
    grid_vecs: np.ndarray | None,
    out_name: str,
    iso: float = 0.05,
    camera_pos: list[float] | np.ndarray | None = None,
    camera_focal: list[float] | np.ndarray | None = None,
    camera_up: list[float] | np.ndarray | None = None,
    transparent_background: bool = True,
    without_mo: bool = False,
) -> None:
    """Render molecular structure and optional orbital isosurfaces to a file."""
    plotter = pv.Plotter(off_screen=True)
    cast(Any, plotter).set_background(BACKGROUND_COLOR)
    if not without_mo and mo_cube is not None and origin is not None and grid_vecs is not None:
        grid = build_volume_grid(mo_cube, origin, grid_vecs)
        positive_iso = cast(pv.PolyData, grid.contour([iso], scalars="mo"))
        negative_iso = cast(pv.PolyData, grid.contour([-iso], scalars="mo"))
        plotter.add_mesh(
            positive_iso,
            color=ORBITAL_COLORS["pos"],
            opacity=ORBITAL_OPACITY,
        )
        plotter.add_mesh(
            negative_iso,
            color=ORBITAL_COLORS["neg"],
            opacity=ORBITAL_OPACITY,
        )
    build_molstruct(plotter, atomnos, coords)

    if camera_pos is None and camera_focal is None and camera_up is None:
        plotter.enable_anti_aliasing(ANTI_ALIASING)
        plotter.screenshot(
            out_name,
            scale=IMAGE_QUALITY,
            transparent_background=transparent_background,
        )
    elif camera_pos is not None and camera_focal is not None and camera_up is not None:
        plotter.enable_anti_aliasing(ANTI_ALIASING)
        cast(Any, plotter).camera_position = [camera_pos, camera_focal, camera_up]
        plotter.screenshot(
            out_name,
            scale=IMAGE_QUALITY,
            transparent_background=transparent_background,
        )
    else:
        plotter.screenshot(
            out_name,
            scale=1,
            transparent_background=transparent_background,
        )
        pos, focus, up = plotter.camera_position
        if camera_pos is not None:
            pos = camera_pos
        if camera_focal is not None:
            focus = camera_focal
        if camera_up is not None:
            up = camera_up
        plotter.enable_anti_aliasing(ANTI_ALIASING)
        cast(Any, plotter).camera_position = [pos, focus, up]
        plotter.screenshot(
            out_name,
            scale=IMAGE_QUALITY,
            transparent_background=transparent_background,
        )
    print(f"Saved {out_name}")


__all__ = ["build_molstruct", "build_volume_grid", "save_mo_plot"]
