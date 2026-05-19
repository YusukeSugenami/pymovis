import os
import re
from dataclasses import dataclass
import numpy as np
from pyscf import gto
from .molfile_parser import load_cube, load_fchk, load_pyscfchk

@dataclass
class MoleculeData:
    file_type: str | None = None
    atoms: list | None = None
    coords: np.ndarray | None = None
    coords_unit: str = "Bohr"
    atomnos: np.ndarray | None = None
    basis: list | None = None
    basis_name: str | None = None
    cart: bool = False
    mixed_shells: bool = False
    ao_shell_types: list[int] | None = None
    mo_coeff: np.ndarray | None = None
    origin: np.ndarray | None = None
    grid_vecs: np.ndarray | None = None
    mo_cube: np.ndarray | None = None
    mo_identity: bool = False

    def convert_A2B(self) -> None:

        if self.coords is None or self.coords_unit == "Bohr":
            return
        self.coords = self.coords / 0.5291772083
        self.coords_unit = "Bohr"

        return

    def convert_B2A(self) -> None:

        if self.coords is None or self.coords_unit == "Angstrom":
            return
        self.coords = self.coords * 0.5291772083
        self.coords_unit = "Angstrom"

        return

    def build_mol(self, cart: bool | None = None) -> gto.Mole:

        mol = gto.Mole()
        mol.unit = self.coords_unit

        if self.atoms is None:
            if self.atomnos is None or self.coords is None:
                raise ValueError("not enough information to build molecule")
            atom_str = ""
            for atom, (x, y, z) in zip(self.atomnos, self.coords):
                atom_str += f"{atom} {x} {y} {z};"
            mol.atom = atom_str
        else:
            mol.atom = self.atoms

        if self.basis is None:
            if self.basis_name is None:
                raise ValueError("basis information is required to build molecule")
            mol.basis = self.basis_name
        else:
            mol.basis = self.basis

        mol.cart = self.cart if cart is None else cart
        mol.build()

        return mol

    def check_mo_identity(self) -> None:

        if self.mo_identity:
            return
        if self.mo_coeff is None:
            raise ValueError("mo_coeff is required to check mo identity")
        
        mol = self.build_mol()
        S = mol.intor("int1e_ovlp")
        C = self.mo_coeff
        nao = S.shape[0]
        I = np.eye(nao)
        err = np.linalg.norm(C.T @ S @ C - I) / nao

        if err > 1e-6:
            raise ValueError(f"mo coefficients are not orthogonal with error {err}")
        else:
            self.mo_identity = True
            return

    def parse_orbitalindex(self, mo_index : int | str) -> int:

        try:
            mo_index = int(mo_index)
            return mo_index
        except ValueError as _:
            if self.atomnos is None:
                raise ValueError("atomnos is required to parse mo index in string format")
            homo_index = int(sum(self.atomnos) / 2 - 1)
            mo_index = str(mo_index).strip().upper()

            if mo_index.startswith("HOMO"):
                if mo_index == "HOMO":
                    return homo_index
                m = re.match(r"HOMO([+-]\d+)", mo_index)
                if m:
                    return homo_index + int(m.group(1))

            if mo_index.startswith("LUMO"):
                if mo_index == "LUMO":
                    return homo_index + 1
                m = re.match(r"LUMO([+-]\d+)", mo_index)
                if m:
                    return homo_index + 1 + int(m.group(1))

            raise ValueError(f"invalid format {mo_index} for mo_index")

    def _evaluate_ao_on_points(self, grid_points: np.ndarray) -> np.ndarray:

        if self.mo_coeff is None:
            raise ValueError("mo_coeff is required to evaluate AO values")

        if self.mixed_shells:
            if self.ao_shell_types is None:
                raise ValueError("ao_shell_types is required for mixed-shell AO evaluation")

            mol_cart = self.build_mol(cart=True)
            mol_sph = self.build_mol(cart=False)
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
            mol = self.build_mol()
            basis_vals = mol.eval_gto("GTOval", grid_points)

        if basis_vals.shape[1] != self.mo_coeff.shape[0]:
            raise ValueError(
                f"AO size mismatch: evaluated {basis_vals.shape[1]} AOs, but mo_coeff has {self.mo_coeff.shape[0]} rows"
            )

        return basis_vals

    def check_mo_identity_on_grid(
        self,
        mo_indices: list[int | str] | None = None,
        padding: float = 2.0,
        grid_setting: str = "size",
        grid_size: list[float] = [0.3333, 0.3333, 0.3333],
        grid_shape: list[int] = [100, 100, 100],
        tol: float = 1e-2,
    ) -> tuple[bool, float, np.ndarray]:
        """
        Check MO orthonormality on a real-space grid with uniform volume weights.

        The check is based on Psi^T W Psi ~= I, where W = dV * I for a regular grid.
        Due to finite box size and grid resolution, this is only an approximate test.
        """

        if self.file_type == "cube":
            raise ValueError("grid identity check is only available when MO coefficients are present")
        if self.coords is None or self.mo_coeff is None:
            raise ValueError("not enough information to evaluate MO identity on grid")

        self.convert_A2B()

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

        grid = np.array(np.meshgrid(grid_x, grid_y, grid_z, indexing="ij"))
        grid_points = grid.reshape(3, -1).T
        basis_vals = self._evaluate_ao_on_points(grid_points)

        if mo_indices is None:
            mo_cols = np.arange(self.mo_coeff.shape[1], dtype=int)
        else:
            mo_cols = np.array([self.parse_orbitalindex(i) for i in mo_indices], dtype=int)

        psi_vals = basis_vals @ self.mo_coeff[:, mo_cols]
        dV = abs(dx * dy * dz)
        gram = psi_vals.T @ psi_vals * dV

        I = np.eye(len(mo_cols))
        err = np.linalg.norm(gram - I) / len(mo_cols)
        return bool(err < tol), float(err), gram

    def evaluate_mo_on_grid(
        self,
        mo_index : int | str,
        padding : float = 2.0,
        grid_setting : str = "size",
        grid_size : list[float] = [0.3333, 0.3333, 0.3333],
        grid_shape : list[int] = [100, 100, 100],
    ) -> None:

        if self.file_type == "cube":
            return
        if self.coords is None or self.mo_coeff is None:
            raise ValueError("not enough information to evaluate MO on grid")

        mo_index = self.parse_orbitalindex(mo_index)
        if mo_index < 0:
            raise ValueError("mo index is invalid")

        if not self.mixed_shells:
            self.check_mo_identity()
        self.convert_A2B()

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
        grid = np.array(np.meshgrid(grid_x, grid_y, grid_z, indexing="ij"))
        grid_points = grid.reshape(3, -1).T  # (N,3)

        basis_vals = self._evaluate_ao_on_points(grid_points)
        mo = basis_vals @ self.mo_coeff[:, mo_index]

        self.origin = min_xyz
        self.grid_vecs = grid_vecs
        self.mo_cube = mo.reshape(grid_shape)

        return


def load_inp(inp_file : str) -> MoleculeData:

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
