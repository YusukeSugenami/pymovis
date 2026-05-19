# PyMoVis: Molecular Orbital Visualization System

[![Python Version](https://img.shields.io/badge/Python-3.12+-blue.svg)](https://www.python.org/downloads/)

PyMoVis is a Python package for visualizing molecular orbitals (MOs) from quantum chemistry calculations. It reads molecular data from Gaussian output files and renders high-quality 3D isosurface visualizations using PyVista.

## Features

- **Multiple Input Formats**: Support for Gaussian Formatted Checkpoint (fchk), cube files, and PySCF checkpoint files
- **Flexible MO Selection**: Identify orbitals by index or using keywords (e.g., `HOMO`, `LUMO`, `HOMO-1`, `LUMO+2`)
- **3D Visualization**: Render molecular orbital isosurfaces with PyVista
- **Customizable Rendering**: Control camera positions, isosurface values, background colors, and more
- **Automatic Grid Generation**: Evaluate orbitals on adaptive grids around the molecule
- **Atom and Bond Visualization**: Display molecular structure with CPK coloring and bond detection

## Installation

### From PyPI
```bash
pip install pymovis
```

### From Source
```bash
git clone https://github.com/YusukeSugenami/pymovis.git
cd movis
pip install -e .
```

### Development Installation
```bash
pip install -e ".[dev]"
```

## Quick Start

### Basic Usage

Visualize multiple orbitals:
```python
from pymovis import savefig

savefig("molecule.fchk", ["HOMO", "LUMO", "LUMO+1"], iso=0.05)
```

### Command Line

```bash
pymovis molecule.fchk HOMO -o homo.png -i 0.05 --transparent True
```

## Usage Examples

### Python API

#### Basic MO Visualization
```python
from movis import main

# Visualize HOMO and LUMO with custom isosurface value
main("water.fchk", ["HOMO", "LUMO"], iso=0.03)
```

#### Advanced Rendering with Camera Control
```python
from movis import main

# Custom camera position, focal point, and up vector
main(
    "water.fchk",
    ["HOMO"],
    iso=0.05,
    camera_pos=[10, 10, 10],
    camera_focal=[0, 0, 0],
    camera_up=[0, 1, 0],
    transparent_background=True
)
```

#### Loading and Manipulating Molecular Data
```python
from movis import load_inp

# Load molecular data
mol_data = load_inp("molecule.fchk")

# Access molecular properties
print(f"Atomic numbers: {mol_data.atomnos}")
print(f"Coordinates shape: {mol_data.coords.shape}")
print(f"Basis name: {mol_data.basis_name}")

# Evaluate orbital on custom grid
mol_data.evaluate_mo_on_grid(
    "HOMO",
    padding=3.0,
    grid_setting="shape",
    grid_shape=[150, 150, 150]
)
```

#### Using Camera Utilities
```python
from movis import infer_camera_from_coords
import numpy as np

# Automatically infer camera position based on molecular geometry
coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
camera_pos, focal, up = infer_camera_from_coords(coords)
```

### Command Line Options

```
Usage: movis input_file [mo_indices] [options]

Positional arguments:
  input_file              Gaussian fchk or cube file to visualize
  mo                      MO indices (optional, multiple allowed)
                         Specify by index (0, 1, ...) or keyword (HOMO, LUMO, HOMO-1, etc.)

Optional arguments:
  -o, --out              Output file names (one per MO)
  -i, --iso              Isosurface value (default: 0.05)
  -b, --basis            Basis set name (e.g., 6-31G(d,p))
  -t, --transparent      Transparent background (default: False)
  -cp, --camera_pos      Camera position (X Y Z)
  -cf, --camera_focal    Focal point (X Y Z)
  -cu, --camera_up       Up vector (X Y Z)
```

## Supported File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| Gaussian Formatted Checkpoint | `.fchk` | Standard Gaussian output format |
| Cube File | `.cube` | Electron density or orbital cube file |
| PySCF Checkpoint | `.chk` | PySCF checkpoint file with SCF results |

## API Reference

### Main Functions

#### `main(inp_file, mo_index, out_name=None, iso=0.03, camera_pos=None, camera_focal=None, camera_up=None, transparent_background=True)`

Main function to visualize molecular orbitals.

**Parameters:**
- `inp_file` (str): Path to input file (fchk, cube, or chk)
- `mo_index` (list): List of MO indices or keywords
- `out_name` (list, optional): Output file names
- `iso` (float): Isosurface value (default: 0.03)
- `camera_pos` (list, optional): Camera position [X, Y, Z]
- `camera_focal` (list, optional): Focal point [X, Y, Z]
- `camera_up` (list, optional): Up vector [X, Y, Z]
- `transparent_background` (bool): Use transparent background

#### `load_inp(inp_file)`

Load and parse molecular data from file.

**Parameters:**
- `inp_file` (str): Path to input file

**Returns:**
- `MoleculeData`: Object containing molecular information

### Data Classes

#### `MoleculeData`

Dataclass containing molecular information.

**Attributes:**
- `file_type` (str): Type of input file
- `atoms` (list): Atom symbols and coordinates
- `atomnos` (ndarray): Atomic numbers
- `coords` (ndarray): Atomic coordinates
- `coords_unit` (str): Coordinate unit ('Bohr' or 'Angstrom')
- `mo_coeff` (ndarray): MO coefficients in AO basis
- `basis` (dict): Basis set information
- `cart` (bool): Cartesian (True) or spherical (False) basis functions

## Configuration

### Settings

Edit visualization parameters in `movis/settings.py`:

```python
# Grid settings
GRID_SETTING = "shape"      # "size" or "shape"
GRID_SHAPE = [100, 100, 100]
GRID_SIZE = [0.3333, 0.3333, 0.3333]

# Visualization
ORBITAL_OPACITY = 1.0
ORBITAL_COLORS = {"pos": "blue", "neg": "red"}
IMAGE_QUALITY = 4

# Atom rendering
ELEMENT_DATA = {...}        # CPK colors and radii
BOND_RADIUS = 0.15
```

## Examples

### Example 1: Basic Water Molecule Visualization

```bash
# Download example file (if available)
movis water.fchk HOMO LUMO -o homo.png lumo.png
```

### Example 2: Batch Processing

```python
from pathlib import Path
from movis import main

# Process all fchk files in a directory
for fchk_file in Path(".").glob("*.fchk"):
    print(f"Processing {fchk_file}...")
    main(str(fchk_file), ["HOMO", "LUMO"])
```

## Testing

Run the test suite:

```bash
pytest
```

With coverage:

```bash
pytest --cov=movis
```



## Acknowledgments

- Built with [PyVista](https://docs.pyvista.org/) for 3D visualization
- Uses [PySCF](https://pyscf.org/) for quantum chemistry calculations
- Inspired by molecular visualization tools in the quantum chemistry community


## References

- [Gaussian09/16 User's Reference](https://gaussian.com/)
- [PySCF Documentation](https://pyscf.org/)
- [PyVista Documentation](https://docs.pyvista.org/)

**Last Updated**: May 2026
**Version**: 0.1.0
