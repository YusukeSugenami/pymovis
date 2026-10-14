import os
import argparse

import numpy as np

from .settings import (
    DEFAULT_OUT_SUFFIX,
    GRID_SETTING,
    GRID_SHAPE,
    GRID_SIZE,
    IMAGE_SUFFIXES,
    PADDING,
    PARSER_ARGS,
)
from .load_mol import load_inp, MoleculeData
from .renderer import save_mo_plot
from .camera_utils import (
    resolve_camera_position,
)


def savemo(
    inp_file: str | None = None,
    mo_index: list[int | str] | int | str = -1,
    out_name: list[str] | None = None,
    iso: float = 0.03,
    camera_pos: list[float] | None = None,
    camera_focal: list[float] | None = None,
    camera_up: list[float] | None = None,
    camera_axis: str | None = None,
    camera_focal_atoms: list[int] | int | None = None,
    camera_axis_atoms: list[int] | int | None = None,
    camera_up_atoms: list[int] | int | None = None,
    camera_axis_reference: list[float] | None = None,
    camera_up_reference: list[float] | None = None,
    transparent_background: bool = True,
    without_mo: bool = False,
    atomnos: np.ndarray | list[int] | None = None,
    atom_symbols: list[str] | None = None,
    coords: np.ndarray | list[list[float]] | None = None,
    coords_unit: str = "Bohr",
    basis: dict | list | str | None = None,
    mo_coeff: np.ndarray | None = None,
    mo_coeff_beta: np.ndarray | None = None,
    mo_occ: np.ndarray | None = None,
    mo_occ_beta: np.ndarray | None = None,
    mo_energy: np.ndarray | None = None,
    mo_energy_beta: np.ndarray | None = None,
    charge: int = 0,
    multiplicity: int = 1,
    mo_cube: np.ndarray | None = None,
    origin: np.ndarray | None = None,
    grid_vecs: np.ndarray | None = None,
) -> None:
    """Load molecular data, evaluate selected orbitals, and save images.

    Input may be a supported file path or direct molecular data. Camera
    options are resolved by
    :func:`pymovis.camera_utils.resolve_camera_position`.
    """

    # Determine output file names based on input file name and MO index
    if isinstance(mo_index, (int, str)):
        mo_index = [mo_index]
    if inp_file is not None:
        basename = os.path.splitext(os.path.basename(inp_file))[0]
    else:
        basename = "molecule"
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
    
    if inp_file is not None:
        moldata = load_inp(inp_file)
    else:
        if coords is None or (atomnos is None and atom_symbols is None):
            raise ValueError("direct input requires coords and atomnos or atom_symbols")
        moldata = MoleculeData(
            file_type="cube" if mo_cube is not None else "arrays",
            coords=np.asarray(coords, dtype=float),
            coords_unit=coords_unit,
            atomnos=None if atomnos is None else np.asarray(atomnos, dtype=int),
            atom_symbols=atom_symbols,
            basis=basis,
            mo_coeff=mo_coeff,
            mo_coeff_beta=mo_coeff_beta,
            mo_occ=mo_occ,
            mo_occ_beta=mo_occ_beta,
            mo_energy=mo_energy,
            mo_energy_beta=mo_energy_beta,
            charge=charge,
            multiplicity=multiplicity,
            mo_cube=mo_cube,
            origin=origin,
            grid_vecs=grid_vecs,
        )
    if without_mo:
        moldata.mo_info = False

    resolved_camera_pos, resolved_camera_focal, resolved_camera_up = resolve_camera_position(
        moldata.coords,
        camera_pos=camera_pos,
        camera_focal=camera_focal,
        camera_up=camera_up,
        camera_axis=camera_axis,
        camera_focal_atoms=camera_focal_atoms,
        camera_axis_atoms=camera_axis_atoms,
        camera_up_atoms=camera_up_atoms,
        camera_axis_reference=camera_axis_reference,
        camera_up_reference=camera_up_reference,
        padding=PADDING,
    )
    
    for mo_idx, image_name in zip(mo_index, out_name):
        moldata.evaluate_mo_on_grid(
            mo_idx,
            padding=PADDING,
            grid_setting=GRID_SETTING,
            grid_size=GRID_SIZE,
            grid_shape=GRID_SHAPE
        )

        save_mo_plot(
            moldata.atomnos,
            moldata.coords,
            moldata.mo_cube,
            moldata.origin,
            moldata.grid_vecs,
            image_name,
            iso=iso,
            camera_pos=resolved_camera_pos,
            camera_focal=resolved_camera_focal,
            camera_up=resolved_camera_up,
            transparent_background=transparent_background,
            without_mo=without_mo
        )


def cli() -> None:
    """Parse command-line arguments and invoke :func:`savemo`."""
    parser = argparse.ArgumentParser(
        description="MO picture generator from Gaussian fchk or cube file."
    )
    for arg_config in PARSER_ARGS:
        parser.add_argument(*arg_config["args"], **arg_config["kwargs"])
    args = parser.parse_args()
    savemo(
        args.input_file,
        args.mo,
        out_name=args.out,
        iso=args.iso,
        camera_pos=args.camera_pos,
        camera_focal=args.camera_focal,
        camera_up=args.camera_up,
        camera_axis=args.camera_axis,
        camera_focal_atoms=args.camera_focal_atoms,
        camera_axis_atoms=args.camera_axis_atoms,
        camera_up_atoms=args.camera_up_atoms,
        camera_axis_reference=args.camera_axis_reference,
        camera_up_reference=args.camera_up_reference,
        transparent_background=args.transparent,
        without_mo=args.without_mo
    )

if __name__ == "__main__":
    cli()
