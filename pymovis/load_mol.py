"""Normalized molecular data model and orbital grid evaluation."""

import os
import re
from collections.abc import Sequence
from typing import Any

import numpy as np
from pyscf import gto
from .parser import load_cube, load_fchk, load_pyscfchk
import periodictable

class MoleculeData:
    """Normalized molecular geometry and molecular-orbital data.

    Parsing is delegated to :mod:`pymovis.parser`; this class owns unit
    conversion, orbital selection, AO evaluation, and MO grid generation.
    """
    def __init__(
        self, 
        file_type: str, 
        coords: np.ndarray,
        coords_unit: str = "Bohr",
        atom_symbols: Sequence[str] | None = None,
        atomnos: np.ndarray | None = None,
        basis: Any = None,
        shell_type: str | None = "sph",
        ao_shell_types: list[int] | None = None,
        mo_coeff: np.ndarray | None = None,
        mo_coeff_beta: np.ndarray | None = None,
        mo_occ: np.ndarray | None = None,
        mo_occ_beta: np.ndarray | None = None,
        mo_energy: np.ndarray | None = None,
        mo_energy_beta: np.ndarray | None = None,
        charge: int = 0,
        multiplicity: int = 1,
        origin: np.ndarray | None = None,
        grid_vecs: np.ndarray | None = None,
        mo_cube: np.ndarray | None = None
    ) -> None:
        
        # Basic molecular information (required)
        self.file_type = file_type
        self.coords = coords
        self.coords_unit = coords_unit

        if atomnos is None:
            if atom_symbols is None:
                raise ValueError("either atomnos or atom_symbols must be provided")
            atomnos = np.array([periodictable.elements.symbol(s).number for s in atom_symbols], dtype=int)
        else:
            if atom_symbols is None:
                atom_symbols = [periodictable.elements[z].symbol for z in atomnos]

        self.atom_symbols: list[str] = list(atom_symbols)
        self.atomnos : np.ndarray = atomnos
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        self.mo_coeff_alpha = mo_coeff
        self.mo_coeff_beta = mo_coeff_beta
        self.mo_occ_alpha = None if mo_occ is None else np.asarray(mo_occ, dtype=float)
        self.mo_occ_beta = None if mo_occ_beta is None else np.asarray(mo_occ_beta, dtype=float)
        self.mo_energy_alpha = None if mo_energy is None else np.asarray(mo_energy, dtype=float)
        self.mo_energy_beta = None if mo_energy_beta is None else np.asarray(mo_energy_beta, dtype=float)

        # MO information (optional, but required for MO visualization)
        if file_type == "cube":
            if mo_cube is None or origin is None or grid_vecs is None:
                print("Warning: incomplete cube data, MO visualization may not work")
                self.mo_info = False
            else:
                self.mo_info = True
                self.mo_identity = True
                self.origin = origin
                self.grid_vecs = grid_vecs
                self.mo_cube = mo_cube
        else:
            if mo_coeff is None or basis is None:
                print("Warning: MO coefficients or basis information is not provided.")
                self.mo_info = False
            else:
                if shell_type == "mixed":
                    if ao_shell_types is None:
                        print("Warning: shell_type is 'mixed' but ao_shell_types is not provided.")
                        self.mo_info = False
                    else:
                        self.mo_info = True
                        self.mo_identity = False
                        self.basis = basis
                        self.shell_type = shell_type
                        self.ao_shell_types = ao_shell_types
                        self.mo_coeff = mo_coeff
                elif shell_type not in ("cart", "sph"):
                    print(f"Warning: unrecognized shell_type '{shell_type}', expected 'cart', 'sph', or 'mixed'.")
                    self.mo_info = False
                else:
                    self.mo_info = True
                    self.mo_identity = False
                    self.basis = basis
                    self.mo_coeff = mo_coeff
                    self.shell_type = shell_type              
                self._initialize_occupations()

        if self.coords_unit == "Angstrom":
            self._convert_A2B()

    
    def _convert_A2B(self) -> None:
        """Convert stored coordinates from Angstrom to Bohr in place."""

        if self.coords_unit == "Bohr":
            return
        self.coords = self.coords / 0.5291772083
        self.coords_unit = "Bohr"

        return

    def _convert_B2A(self) -> None:
        """Convert stored coordinates from Bohr to Angstrom in place."""

        if self.coords_unit == "Angstrom":
            return
        self.coords = self.coords * 0.5291772083
        self.coords_unit = "Angstrom"

        return

    def _initialize_occupations(self) -> None:
        """Infer missing occupations from electron count and multiplicity."""
        if self.mo_coeff_alpha is None:
            return
        nmo = self.mo_coeff_alpha.shape[1]
        total_electrons = int(round(np.sum(self.atomnos) - self.charge))
        if self.multiplicity < 1:
            raise ValueError("multiplicity must be a positive integer")
        if (total_electrons + self.multiplicity - 1) % 2:
            raise ValueError(
                "charge and multiplicity imply a non-integer number of alpha electrons"
            )
        n_alpha = (total_electrons + self.multiplicity - 1) // 2
        n_beta = total_electrons - n_alpha
        if n_alpha < 0 or n_beta < 0:
            raise ValueError("charge and multiplicity imply a negative electron count")
        if n_alpha > nmo or (self.mo_coeff_beta is not None and n_beta > self.mo_coeff_beta.shape[1]):
            raise ValueError("the supplied MO coefficient arrays have too few orbitals")
        if self.mo_occ_alpha is None:
            self.mo_occ_alpha = np.array([1.0 if i < n_alpha else 0.0 for i in range(nmo)])
        if self.mo_coeff_beta is not None and self.mo_occ_beta is None:
            self.mo_occ_beta = np.array(
                [1.0 if i < n_beta else 0.0 for i in range(self.mo_coeff_beta.shape[1])]
            )

    def _parse_orbitalindex(self, mo_index: int | str) -> tuple[str, int]:
        """Resolve an orbital index or HOMO/LUMO label to spin and index."""
        self._initialize_occupations()
        spin = "alpha"
        label = str(mo_index).strip()
        lower_label = label.lower()
        for prefix in ("alpha:", "a:", "beta:", "b:"):
            if lower_label.startswith(prefix):
                spin = "beta" if prefix in ("beta:", "b:") else "alpha"
                label = label[len(prefix):].strip()
                break
        if lower_label.endswith("_alpha") or lower_label.endswith("-alpha"):
            spin, label = "alpha", label[:-6]
        elif lower_label.endswith("_beta") or lower_label.endswith("-beta"):
            spin, label = "beta", label[:-5]

        try:
            index = int(label)
            if index < 0:
                raise ValueError("mo index must be non-negative")
            return spin, index
        except ValueError as _:
            label = label.upper()
            occupations = (
                self.mo_occ_beta if spin == "beta" and self.mo_occ_beta is not None
                else self.mo_occ_alpha
            )
            coefficients = (
                self.mo_coeff_beta if spin == "beta" and self.mo_coeff_beta is not None
                else self.mo_coeff_alpha
            )
            if occupations is None or coefficients is None:
                raise ValueError(f"{spin} orbitals are not available")
            occupied = np.flatnonzero(occupations > 1e-8)
            homo_index = int(occupied[-1]) if occupied.size else -1
            if label.startswith("HOMO"):
                if label == "HOMO":
                    index = homo_index
                else:
                    m = re.fullmatch(r"HOMO([+-]\d+)", label)
                    if m:
                        index = homo_index + int(m.group(1))
                    else:
                        index = None
            elif label.startswith("LUMO"):
                if label == "LUMO":
                    index = homo_index + 1
                else:
                    m = re.fullmatch(r"LUMO([+-]\d+)", label)
                    if m:
                        index = homo_index + 1 + int(m.group(1))
                    else:
                        index = None
            else:
                index = None
            if index is None:
                raise ValueError(f"failed to parse mo index {mo_index}")
            if index < 0 or index >= coefficients.shape[1]:
                raise ValueError(f"mo index {mo_index} is out of bounds")
            
            return spin, index

    def _build_mol(self, cart: bool | None = None) -> gto.Mole:
        """Build a PySCF molecule from the stored geometry and basis."""

        mol = gto.Mole()
        mol.unit = self.coords_unit
        mol.charge = self.charge
        mol.spin = self.multiplicity - 1

        atom_str = ""
        for symbol, (x, y, z) in zip(self.atom_symbols, self.coords):
            atom_str += f"{symbol} {x} {y} {z};"
        mol.atom = atom_str

        if cart is not None:
            mol.cart = cart
        else:
            mol.cart = self.shell_type == "cart"
        mol.basis = self.basis
        mol.build()

        return mol

    def _check_mo_identity(self, tol: float = 1e-5) -> None:
        """Validate MO orthonormality in the AO overlap metric.

        Formatted checkpoint files contain decimal representations of both
        orbital coefficients and basis data.  Rebuilding the overlap matrix
        with those values can therefore leave a small residual, particularly
        for post-SCF or TD calculations.  ``tol`` is applied to the
        Frobenius-norm residual divided by the number of orbitals.
        """

        if self.mo_identity:
            return
        
        mol = self._build_mol()
        S = mol.intor("int1e_ovlp")
        C = self.mo_coeff
        nmo = C.shape[1]
        I = np.eye(nmo)
        err = np.linalg.norm(C.T @ S @ C - I) / nmo

        if err > tol:
            raise ValueError(f"mo coefficients are not orthogonal with error {err}")
        else:
            self.mo_identity = True
            return

    def _check_mo_identity_on_grid(
        self,
        basis_vals: np.ndarray,
        grid_size: Sequence[float],
        tol: float = 1e-2,
    ) -> None:
        """
        Check MO orthonormality on a real-space grid with uniform volume weights.

        The check is based on Psi^T W Psi ~= I, where W = dV * I for a regular grid.
        Due to finite box size and grid resolution, this is only an approximate test.
        """

        if self.mo_identity:
            return
        psi_vals = basis_vals @ self.mo_coeff # (Ngrid, Nmo)
        dx, dy, dz = grid_size
        dV = abs(dx * dy * dz)
        gram = psi_vals.T @ psi_vals * dV
        nmo = self.mo_coeff.shape[1]
        I = np.eye(nmo)
        err = np.linalg.norm(gram - I) / nmo
        if err > tol:
            raise ValueError(f"MO orthonormality check on grid failed with error {err}")
        else:
            self.mo_identity = True
        return 

    def _evaluate_ao_on_grid(self, grid_points: np.ndarray) -> np.ndarray:
        """Evaluate all stored AO basis functions at Cartesian grid points."""

        if self.shell_type == "mixed":
            mol_cart = self._build_mol(cart=True)
            mol_sph = self._build_mol(cart=False)
            basis_cart = mol_cart.eval_gto("GTOval_cart", grid_points)
            basis_sph = mol_sph.eval_gto("GTOval_sph", grid_points)
            blocks = []
            ic = 0
            isph = 0
            for shell_type in self.ao_shell_types:
                l = abs(shell_type)
                ncart = (l + 1) * (l + 2) // 2
                nsph = 2 * l + 1
                cart_block = basis_cart[:, ic : ic + ncart]
                sph_block = basis_sph[:, isph : isph + nsph]

                # Keep the AO representation defined in fchk for each shell.
                if shell_type > 1:
                    blocks.append(cart_block)
                else:
                    blocks.append(sph_block)

                ic += ncart
                isph += nsph

            basis_vals = np.hstack(blocks)
        else:
            mol = self._build_mol()
            basis_vals = mol.eval_gto("GTOval", grid_points)

        return basis_vals

    def evaluate_mo_on_grid(
        self,
        mo_index: int | str,
        padding: float = 2.0,
        grid_setting: str = "size",
        grid_size: Sequence[float] = (0.3333, 0.3333, 0.3333),
        grid_shape: Sequence[int] = (100, 100, 100),
    ) -> None:
        """Evaluate one selected MO on a regular Cartesian grid."""

        if self.file_type == "cube":
            return
        if not self.mo_info:
            self.origin = None
            self.grid_vecs = None
            self.mo_cube = None
            return

        spin, mo_index = self._parse_orbitalindex(mo_index)
        coefficients = (
            self.mo_coeff_beta if spin == "beta" and self.mo_coeff_beta is not None
            else self.mo_coeff_alpha
        )
        if coefficients is None:
            raise ValueError(f"{spin} MO coefficients are not available")

        if self.shell_type != "mixed":
            self._check_mo_identity()
        self._convert_A2B()

        min_xyz = self.coords.min(axis=0) - padding
        max_xyz = self.coords.max(axis=0) + padding

        if grid_setting == "size":
            dx, dy, dz = grid_size
            grid_x = np.arange(min_xyz[0], max_xyz[0], dx)
            grid_y = np.arange(min_xyz[1], max_xyz[1], dy)
            grid_z = np.arange(min_xyz[2], max_xyz[2], dz)
            nx = len(grid_x)
            ny = len(grid_y)
            nz = len(grid_z)

        elif grid_setting == "shape":
            nx, ny, nz = grid_shape
            grid_x = np.linspace(min_xyz[0], max_xyz[0], nx, retstep=False)
            grid_y = np.linspace(min_xyz[1], max_xyz[1], ny, retstep=False)
            grid_z = np.linspace(min_xyz[2], max_xyz[2], nz, retstep=False)
            dx = grid_x[1] - grid_x[0]
            dy = grid_y[1] - grid_y[0]
            dz = grid_z[1] - grid_z[0]

        else:
            raise ValueError("not supported")

        grid_shape = [nx, ny, nz]
        grid_vecs = np.array([[dx, 0.0, 0.0], [0.0, dy, 0.0], [0.0, 0.0, dz]])
        grid_points = np.array(np.meshgrid(grid_x, grid_y, grid_z, indexing="ij")).reshape(3, -1).T  # (N,3)
        basis_vals = self._evaluate_ao_on_grid(grid_points)

        if self.shell_type == "mixed":
            self._check_mo_identity_on_grid(basis_vals, grid_size=[dx, dy, dz])
        mo = (basis_vals @ coefficients[:, mo_index]).reshape(grid_shape)

        self.origin = min_xyz
        self.grid_vecs = grid_vecs
        self.mo_cube = mo

        return


def load_inp(inp_file: str) -> MoleculeData:
    """Load a supported molecular file as normalized ``MoleculeData``."""

    _, ext = os.path.basename(inp_file).rsplit(".", 1)
    if ext.lower() == "cube":
        cube_info = load_cube(inp_file)
        data = MoleculeData(**cube_info)
    elif ext.lower() == "fchk":
        fchk_info = load_fchk(inp_file)
        data = MoleculeData(**fchk_info)
    elif ext.lower() == 'chk':
        pyscf_info = load_pyscfchk(inp_file)
        data = MoleculeData(**pyscf_info)
    else:
        raise ValueError("not supported yet")

    return data
