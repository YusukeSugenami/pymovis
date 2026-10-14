# PyMovis

PyMovis is a Python package for visualizing molecular orbitals from quantum chemistry calculations. It reads Gaussian formatted checkpoint files, cube files, and PySCF checkpoint files, then renders 3D orbital isosurfaces with PyVista.

## Installation

You can install PyMovis with pip:

```bash
pip install pymovis
```

or from source:

```bash
git clone https://github.com/YusukeSugenami/pymovis.git
cd pymovis
pip install -e .
```

## Dependencies

PyMovis requires Python 3.12 or later and installs these runtime
dependencies:

- [NumPy](https://numpy.org/) `>=1.20.0`
- [periodictable](https://pypi.org/project/periodictable/) `>=2.1.0`
- [PySCF](https://pyscf.org/) `>=2.13.0`
- [PyVista](https://pyvista.org/) `>=0.48.4`

For development and testing, install the optional development dependencies:

```bash
pip install -e ".[dev]"
```

This additionally installs pytest, pytest-cov, Black, Ruff, and isort.

## Quick Start

### Command line

The package installs a `pymovis` command.

```bash
pymovis test.fchk HOMO
```

Visualize multiple orbitals:

```bash
pymovis test.fchk HOMO LUMO
```

Save to specific output files:

```bash
pymovis test.fchk HOMO LUMO -o homo.png lumo.png
```

### Python API

```python
from pymovis import savemo

savemo("test.fchk", ["HOMO", "LUMO"])
```

## Camera Control

PyMoVis supports three levels of camera control.

### 1. Manual camera settings

Use `camera_pos`, `camera_focal`, and `camera_up` directly:

```bash
pymovis test.fchk HOMO \
  -cp 10 10 10 \
  -cf 0 0 0 \
  -cu 0 1 0
```

### 2. Automatic axis-based camera

Use `--camera_axis` to place the camera so the molecule fits in the frame.
The current renderer defaults are used for this calculation:
- window size: `1024 x 768`
- aspect ratio: `4:3`
- vertical field of view: `30 deg`

```bash
pymovis test.fchk HOMO --camera_axis X
```

This option looks at the molecule from the negative side of the selected axis and uses the corresponding screen orientation:
- `X` view: right is `+Y`, up is `+Z`
- `Y` view: right is `+X`, up is `+Z`
- `Z` view: right is `+X`, up is `+Y`

You can invert the axis direction with prefix '-' (e.g. -X)

### 3. Atom-index based camera control

You can define the camera center and axes from atom indices. Indices are 0-based.

#### Camera center

Use `--camera_focal_atoms`.

- One atom index: the camera center is that atom position
- Two or more atom indices: the camera center is the centroid of those atoms

Example:

```bash
pymovis test.fchk HOMO --camera_focal_atoms 0
pymovis test.fchk HOMO --camera_focal_atoms 0 1 2
```

#### Camera axis and up direction

Use `--camera_axis_atoms` for the camera view axis and `--camera_up_atoms` for the screen up direction.

- Two atom indices: the axis is the line through those atoms
- Three or more atom indices: the axis is the normal vector of the best-fit plane through those atoms

Example:

```bash
pymovis test.fchk HOMO \
  --camera_focal_atoms 0 1 2 \
  --camera_axis_atoms 0 1 \
  --camera_up_atoms 0 1 3
```

#### Fixing axis sign with reference vectors

When using three or more atoms to define an axis (via best-fit plane normal), the sign of the resulting axis is ambiguous. You can fix this using `--camera_axis_reference` and `--camera_up_reference`.

Each reference option accepts:
- **Two atom indices**: `[idx1 idx2]` → the reference vector is `coords[idx2] - coords[idx1]`
- **Three floats**: `[X Y Z]` → the reference vector is given directly

The axis is flipped if its dot product with the reference vector is negative.

Example using atom indices:

```bash
pymovis test.fchk HOMO \
  --camera_axis_atoms 0 1 2 \
  --camera_up_atoms 3 4 5 \
  --camera_axis_reference 0 1 \
  --camera_up_reference 2 5
```

Example using direct coordinates:

```bash
pymovis test.fchk HOMO \
  --camera_axis_atoms 0 1 2 \
  --camera_up_atoms 3 4 5 \
  --camera_axis_reference 0 0 1 \
  --camera_up_reference 0 1 0
```

### Camera precedence

The current behavior is:

- If `--camera_axis_atoms` and `--camera_up_atoms` are provided, those atom-based axes are used
- Otherwise, if `--camera_axis` is provided, the axis-based automatic camera is used
- If `--camera_focal` is provided, it becomes the camera center used by the automatic camera calculation
- If `--camera_up` is provided, it overrides the final up vector
- If none of the automatic options are provided, manual `camera_pos`, `camera_focal`, and `camera_up` are used as-is

## CLI Options

```text
Usage: pymovis input_file [mo ...] [options]

Positional arguments:
  input_file              Gaussian fchk or cube file to visualize
  mo                      MO indices to visualize
                          Use integer indices starting from 0 or labels such as HOMO/LUMO

Options:
  -o, --out               Output file names
  -i, --iso               Isosurface value
  -t, --transparent       Transparent background flag
  -cp, --camera_pos       Manual camera position
  -cf, --camera_focal     Manual camera focal point
  -cu, --camera_up        Manual camera up vector
  -ca, --camera_axis      Automatic camera axis: X / Y / Z. Use '-' prefix (e.g. -X) to invert the axis direction.
  --camera_focal_atoms    Atom indices used to define the camera focal point
  --camera_axis_atoms     Atom indices used to define the camera view axis
  --camera_up_atoms       Atom indices used to define the camera up axis
  --camera_axis_reference Reference vector for camera view axis sign (2 atom indices or 3 coordinates)
  --camera_up_reference   Reference vector for camera up axis sign (2 atom indices or 3 coordinates)
  -womo, --without_mo     Only render molecule structure without MO isosurface
```

## Supported File Formats

| Format | Extension | Description |
| --- | --- | --- |
| Gaussian Formatted Checkpoint | `.fchk` | Gaussian MO and geometry data |
| Cube File | `.cube` | Cube grid data |
| PySCF checkpoint | `.chk` | PySCF SCF checkpoint |

## Python API

### `savemo`

```python
savemo(
    inp_file,
    mo_index,
    out_name=None,
    iso=0.03,
    camera_pos=None,
    camera_focal=None,
    camera_up=None,
    camera_axis=None,
    camera_focal_atoms=None,
    camera_axis_atoms=None,
    camera_up_atoms=None,
    camera_axis_reference=None,
    camera_up_reference=None,
    transparent_background=True,
    without_mo=False,
)
```

- `inp_file`: input file path
- `mo_index`: list of MO indices or labels
- `out_name`: output file names
- `iso`: isosurface value
- `camera_axis`: automatic axis-based view (`X`, `Y`, `Z`)
- `camera_focal_atoms`: atom indices for the camera focal point
- `camera_axis_atoms`: atom indices for the view axis
- `camera_up_atoms`: atom indices for the up axis
- `camera_axis_reference`: reference vector for camera axis sign (2 atom indices or 3D coordinates)
- `camera_up_reference`: reference vector for camera up axis sign (2 atom indices or 3D coordinates)
- `transparent_background`: save transparent background picture if True
- `without_mo`: save only molecular structure if True

### `load_inp`

### Direct molecular data

`savemo` can also be used without an input file.  Provide atomic data, a
PySCF-compatible basis, and alpha MO coefficients directly:

```python
import numpy as np
from pymovis import savemo

savemo(
    mo_index=["HOMO", "LUMO"],
    atom_symbols=["O", "H", "H"],
    coords=np.array([[0, 0, 0], [0, 1.4, 1.1], [0, -1.4, 1.1]]),
    coords_unit="Bohr",
    basis={"H": [[0, (1.0, 1.0)]], "O": [[0, (1.0, 1.0)]]},
    mo_coeff=coefficients,
    camera_axis="Z",
)
```

For a precomputed cube, pass `atomnos` or `atom_symbols`, `coords`,
`mo_cube`, `origin`, and `grid_vecs` instead of `basis` and `mo_coeff`.
`molecule=MoleculeData(...)` is also accepted when the data object is already
constructed.

### Charged and open-shell systems

Formatted checkpoint files now use the alpha and beta electron counts and
read both alpha and beta MO coefficient arrays when present.  `HOMO` and
`LUMO` are resolved independently from the occupied orbital arrays rather
than from the sum of atomic numbers.  Select a spin channel with
`"alpha:HOMO"`, `"beta:HOMO"`, `"A:LUMO"`, or `"B:LUMO"`.  Direct and PySCF
checkpoint inputs accept `charge`, `multiplicity`, `mo_occ`, `mo_occ_beta`,
`mo_coeff_beta`, `mo_energy`, and `mo_energy_beta`.

```python
from pymovis import load_inp
mol = load_inp("test.fchk")
```

Returns a `MoleculeData` object with parsed geometry, basis, and MO data.

## Camera Utilities

If you want to use the camera logic directly, the package exposes helper functions such as:

```python
from pymovis import (
    infer_axis_camera_fit_from_coords,
    infer_atom_axis_camera_fit_from_coords,
    resolve_camera_position,
)
```

`resolve_camera_position` applies the same precedence rules as `savemo` and
returns `(camera_pos, camera_focal, camera_up)`. Atom-defined axes take
precedence over `camera_axis`, while an explicitly supplied `camera_up`
overrides an inferred up vector.

## Package structure

The package separates responsibilities into focused modules:

- `parser/`: file-format-specific readers for `.fchk`, `.cube`, and `.chk`
- `load_mol.py`: normalized `MoleculeData`, orbital selection, and grid evaluation
- `camera_utils.py`: camera geometry, fitting, and option resolution
- `renderer.py`: PyVista scene construction and image output
- `pymovis.py`: public orchestration and command-line entry point

## Configuration

Default rendering and parser settings are defined in [pymovis/settings.py](pymovis/settings.py).

Important values include:

- `PADDING`
- `GRID_SETTING`
- `GRID_SHAPE`
- `IMAGE_QUALITY`
- `BACKGROUND_COLOR`
- `PARSER_ARGS`
