import os
import re
from dataclasses import dataclass
import numpy as np
from pyscf import gto
from molfile_parser import load_cube, load_fchk, load_pyscfchk

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

    def build_mol(self) -> gto.Mole:

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

        mol.cart = self.cart
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
            mo_index = mo_index.strip().upper()

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

        grid_shape = (nx, ny, nz)
        grid_vecs = np.array([[dx, 0.0, 0.0], [0.0, dy, 0.0], [0.0, 0.0, dz]])
        grid = np.array(np.meshgrid(grid_x, grid_y, grid_z, indexing="ij"))
        grid_points = grid.reshape(3, -1).T  # (N,3)

        mol = self.build_mol()
        basis_vals = mol.eval_gto("GTOval", grid_points)  # shape = (N, nbas)
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
