import os
import re
import numpy as np
from pyscf import gto
from .molfile_parser import load_cube, load_fchk, load_pyscfchk
import periodictable

class MoleculeData:
    def __init__(
        self, 
        file_type: str, 
        coords: np.ndarray,
        coords_unit: str = "Bohr",
        atom_symbols: list | None = None,
        atomnos: np.ndarray | None = None,
        basis: list | str | None = None,
        shell_type: str | None = None,
        ao_shell_types: list[int] | None = None,
        mo_coeff: np.ndarray | None = None,
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

        self.atom_symbols : list[str] = atom_symbols
        self.atomnos : np.ndarray = atomnos

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

    
    def _convert_A2B(self) -> None:

        if self.coords_unit == "Bohr":
            return
        self.coords = self.coords / 0.5291772083
        self.coords_unit = "Bohr"

        return

    def _convert_B2A(self) -> None:

        if self.coords_unit == "Angstrom":
            return
        self.coords = self.coords * 0.5291772083
        self.coords_unit = "Angstrom"

        return

    def _parse_orbitalindex(self, mo_index : int | str) -> int:

        try:
            mo_index = int(mo_index)
            if mo_index < 0:
                raise ValueError("mo index must be non-negative")
            return mo_index
        except ValueError as _:
            homo_index = int(sum(self.atomnos) / 2 - 1)
            mo_index = str(mo_index).strip().upper()
            if mo_index.startswith("HOMO"):
                if mo_index == "HOMO":
                    index = homo_index
                else:
                    m = re.match(r"HOMO([+-]\d+)", mo_index)
                    if m:
                        index = homo_index + int(m.group(1))
                    else:
                        index = None
            elif mo_index.startswith("LUMO"):
                if mo_index == "LUMO":
                    index = homo_index + 1
                else:
                    m = re.match(r"LUMO([+-]\d+)", mo_index)
                    if m:
                        index = homo_index + 1 + int(m.group(1))
                    else:
                        index = None
            else:
                index = None
            if index is None:
                raise ValueError(f"failed to parse mo index {mo_index}")
            if index < 0 or index >= self.mo_coeff.shape[1]:
                raise ValueError(f"mo index {mo_index} is out of bounds")
            
            return index

    def _build_mol(self, cart: bool | None = None) -> gto.Mole:

        mol = gto.Mole()
        mol.unit = self.coords_unit

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

    def _check_mo_identity(self, tol : float = 1e-6) -> None:

        if self.mo_identity:
            return
        
        mol = self._build_mol()
        S = mol.intor("int1e_ovlp")
        C = self.mo_coeff
        nao = S.shape[0]
        I = np.eye(nao)
        err = np.linalg.norm(C.T @ S @ C - I) / nao

        if err > tol:
            raise ValueError(f"mo coefficients are not orthogonal with error {err}")
        else:
            self.mo_identity = True
            return

    def _check_mo_identity_on_grid(self, basis_vals: np.ndarray, grid_size: list[float], tol: float = 1e-2) -> None:
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
        I = np.eye(self.mo_coeff.shape[1])
        err = np.linalg.norm(gram - I) / self.mo_coeff.shape[1]
        if err > tol:
            raise ValueError(f"MO orthonormality check on grid failed with error {err}")
        else:
            self.mo_identity = True
        return 

    def _evaluate_ao_on_grid(self, grid_points: np.ndarray) -> np.ndarray:

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
        mo_index : int | str,
        padding : float = 2.0,
        grid_setting : str = "size",
        grid_size : list[float] = [0.3333, 0.3333, 0.3333],
        grid_shape : list[int] = [100, 100, 100],
    ) -> None:

        if self.file_type == "cube":
            return
        if not self.mo_info:
            self.origin = None
            self.grid_vecs = None
            self.mo_cube = None
            return

        mo_index = self._parse_orbitalindex(mo_index)

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
        mo = (basis_vals @ self.mo_coeff[:, mo_index]).reshape(grid_shape)

        self.origin = min_xyz
        self.grid_vecs = grid_vecs
        self.mo_cube = mo

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
