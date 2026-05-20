# PyMovis

PyMovis is a Python package for visualizing molecular orbitals from quantum chemistry calculations. It reads Gaussian formatted checkpoint files, cube files, and PySCF checkpoint files, then renders 3D orbital isosurfaces with PyVista.

## Features

- Visualize molecular orbitals from `fchk`, `cube`, and `chk` files
- Select orbitals by index or by labels such as `HOMO`, `LUMO`, `HOMO-1`, `LUMO+2`
- Render atoms, bonds, and orbital isosurfaces in 3D
- Use automatic camera placement with axis-based views
- Specify camera center and camera axes from atom indices
- Save high-quality images with transparent backgrounds

## Installation

### From source

```bash
git clone https://github.com/YusukeSugenami/pymovis.git
cd pymovis
pip install -e .
```

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
  -b, --basis             Basis set name for manual specification
  -t, --transparent       Transparent background flag
  -cp, --camera_pos       Manual camera position
  -cf, --camera_focal     Manual camera focal point
  -cu, --camera_up        Manual camera up vector
  -ca, --camera_axis      Automatic camera axis: X / Y / Z
  --camera_focal_atoms    Atom indices used to define the camera focal point
  --camera_axis_atoms     Atom indices used to define the camera view axis
  --camera_up_atoms       Atom indices used to define the camera up axis
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
- `transparent_background`: save transparent background picture if True
- `without_mo`: save only molecular structure if True

### `load_inp`

```python
from pymovis import load_inp
mol = load_inp("test.fchk")
```

Returns a `MoleculeData` object with parsed geometry, basis, and MO data.

## Camera Utilities

If you want to use the camera logic directly, the package exposes helper functions such as:

```python
from pymovis import (
    infer_camera_from_coords,
    infer_axis_camera_fit_from_coords,
    infer_atom_axis_camera_fit_from_coords,
)
```

## Configuration

Default rendering and parser settings are defined in [pymovis/settings.py](pymovis/settings.py).

Important values include:

- `PADDING`
- `GRID_SETTING`
- `GRID_SHAPE`
- `IMAGE_QUALITY`
- `BACKGROUND_COLOR`
- `PARSER_ARGS`

## Example Workflows

### Visualize the HOMO from a file

```bash
pymovis water.fchk HOMO
```

### Visualize HOMO and LUMO with automatic axis view

```bash
pymovis water.fchk HOMO LUMO --camera_axis Z
```

### Visualize with atom-defined camera center and axes

```bash
pymovis water.fchk HOMO \
  --camera_focal_atoms 0 1 2 \
  --camera_axis_atoms 0 1 \
  --camera_up_atoms 0 1 3
```

### Save transparent output

```bash
pymovis water.fchk HOMO --transparent True
```

