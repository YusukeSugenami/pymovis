import numpy as np
import periodictable
from pyscf import lib

###############################################################################
# fchk parser
###############################################################################

def read_fchk_array(lines, key, dtype=float):
    """
    fchk array reader

    Parameters
    ----------
    lines : list[str]
    key : str
    dtype : type

    Returns: ndarray
    """

    for i, line in enumerate(lines):
        if key in line:
            # number of elements
            n = int(line.split()[-1])
            vals = []
            j = i + 1
            while len(vals) < n:
                vals.extend(lines[j].split())
                j += 1
            return np.array(vals[:n], dtype=dtype)

    raise ValueError(f"{key} not found")


def load_fchk(inp_file : str):

    with open(inp_file) as f:
        lines = f.readlines()
    
    # ----------------------------------
    # parse geometries
    # ----------------------------------

    atomnos = read_fchk_array(lines, "Atomic numbers", int) # array of atomic numbers starts from 1
    coords = read_fchk_array(lines, "Current cartesian coordinates").reshape(-1, 3) # bohr
    coords_unit = 'Bohr'

    atoms = []
    for i, (z, r) in enumerate(zip(atomnos, coords)):
        symbol = periodictable.elements[z]
        atoms.append([f'{symbol}{i}', r])

    
    # ----------------------------------
    # parse basis functions
    # ----------------------------------

    shell_types = read_fchk_array(lines, "Shell types", int)
    nprim_per_shell = read_fchk_array(lines, "Number of primitives per shell", int)
    shell_to_atom = read_fchk_array(lines, "Shell to atom map", int)
    primitive_exponents = read_fchk_array(lines, "Primitive exponents") 
    contraction_coeffs = read_fchk_array(lines, "Contraction coefficients")

    # SP shell coefficients
    try:
        p_con_coeffs = read_fchk_array(lines, "P(S=P) Contraction coefficients")
    except:
        p_con_coeffs = None
        
    if p_con_coeffs is not None and len(p_con_coeffs) != len(primitive_exponents):
        raise ValueError("Length of P(S=P) Contraction coefficients does not match number of primitives")
    
    # check orbital type: spherical / cartesian
    has_cartesian = any(x > 1 for x in shell_types)
    has_spherical = any(x < -1 for x in shell_types)
    if has_cartesian and has_spherical:
        raise RuntimeError("Mixed Cartesian/spherical shells detected")
    cart = has_cartesian


    # ----------------------------------
    # construct pyscf basis
    # ----------------------------------

    prim_ptr = 0
    current_ao = 0
    pyscf_basis = {}
    fchk_basis = {}
    
    shell_ao_map = {
         0 : [0],
         1 : [0, 1, 2],
        -2 : [4, 2, 0, 1, 3], #shperical d
         2 : [0, 3, 4, 1, 5, 2], # cartecian d
        -3 : [6, 4, 2, 0, 1, 3, 5], # spherical f
         3 : [0, 4, 5, 3, 9, 6, 1, 8, 7, 2], # cartesian f
        -4 : [8, 6, 4, 2, 0, 1, 3, 5, 7], # spherical g
         4 : [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0] # cartesian g
    }

    for ish, stype in enumerate(shell_types):
        nprim = nprim_per_shell[ish]
        atom_idx = shell_to_atom[ish] - 1 # 0 start index for python
        atom_symbol = atoms[atom_idx][0]
        exps = primitive_exponents[prim_ptr : prim_ptr + nprim]
        coeffs = contraction_coeffs[prim_ptr : prim_ptr + nprim]

        if p_con_coeffs is not None:
            pcoeffs = p_con_coeffs[prim_ptr : prim_ptr + nprim]
        else:
            pcoeffs = None

        prim_ptr += nprim
        pyscf_basis.setdefault(atom_symbol, [])
        fchk_basis.setdefault(atom_symbol, [])

        # SP shell
        if stype == -1:
            if pcoeffs is None:
                raise ValueError("P(S=P) Contraction coefficients are required for SP shells")
            # S shell
            s_block: list = [0]
            for e, cs in zip(exps, coeffs):
                s_block.append([float(e), float(cs)])
            pyscf_basis[atom_symbol].append(s_block)

            fchk_basis[atom_symbol].append((0, [current_ao]))
            current_ao += 1

            # P shell
            p_block: list = [1]
            for e, cp in zip(exps, pcoeffs):
                p_block.append([float(e), float(cp)])
            pyscf_basis[atom_symbol].append(p_block)

            fchk_basis[atom_symbol].append((1, list(range(current_ao, current_ao + 3))))
            current_ao += 3

        # normal shells
        else:
            l = abs(stype)
            block = [l]
            for i, (e, c) in enumerate(zip(exps, coeffs)):
                block.append([float(e), float(c)]) 
            pyscf_basis[atom_symbol].append(block)
            
            ao_map = shell_ao_map[stype]
            ao_map = [m + current_ao for m in ao_map]
            fchk_basis[atom_symbol].append((l, ao_map))
            current_ao += len(ao_map)

    for atom_symbol in pyscf_basis:
        pyscf_basis[atom_symbol].sort(key=lambda x: x[0])

    
    # ----------------------------------
    # parse mo coefficients
    # ----------------------------------
    
    nao = current_ao
    coeffs = read_fchk_array(lines, "Alpha MO coefficients")
    nmo = coeffs.size // nao
    C = coeffs.reshape(nmo, nao).T # [ao, mo]

    ao_ordering = []
    ao_scaling = []
    for atom_symbol in fchk_basis:
        # sort shells by absolute angular momentum (keep S, P before D, etc.)
        fchk_basis[atom_symbol].sort(key=lambda x: abs(x[0]))

        for shell in fchk_basis[atom_symbol]:
            ang = shell[0]
            ao_indices = shell[1]
            ao_ordering.extend(ao_indices)
            if ang == 2 and len(ao_indices) == 6:  # cartesian d
                const = np.sqrt(15 / (4 * np.pi))
                ao_scaling.extend([
                    const / np.sqrt(3), # dxx
                    const,              # dxy
                    const,              # dxz
                    const / np.sqrt(3),  # dyy
                    const,              # dyz
                    const / np.sqrt(3)  # dzz
                ])
            elif ang == 3 and len(ao_indices) == 10: # cartesian f
                const = np.sqrt(105 / (4 * np.pi))
                ao_scaling.extend([
                    const / np.sqrt(15), #fxxx
                    const / np.sqrt(3),  #fxxy
                    const / np.sqrt(3),  #fxxz
                    const / np.sqrt(3),  #fxyy
                    const,               #fxyz
                    const / np.sqrt(3),  #fxzz
                    const / np.sqrt(15), #fyyy
                    const / np.sqrt(3),  #fyyz
                    const / np.sqrt(3),  #fyzz
                    const / np.sqrt(15)  #fzzz
                ])
            elif ang == 4 and len(ao_indices) == 15: # cartesian g
                const = np.sqrt(945 / (4 * np.pi))
                ao_scaling.extend([
                    const / np.sqrt(105), #gxxxx
                    const / np.sqrt(15),  #gxxxy
                    const / np.sqrt(15),  #gxxxz
                    const / np.sqrt(9),  #gxxyy
                    const / np.sqrt(3),   #gxxyz
                    const / np.sqrt(9),  #gxxzz
                    const / np.sqrt(15), #gxyyy
                    const / np.sqrt(3),  #gxyyz
                    const / np.sqrt(3),  #gxyzz
                    const / np.sqrt(15), #gxzzz
                    const / np.sqrt(105), #gyyyy
                    const / np.sqrt(15),  #gyyyz
                    const / np.sqrt(9),  #gyyzz
                    const / np.sqrt(15), #gyzzz
                    const / np.sqrt(105), #gzzzz
                ])
            else:
                ao_scaling.extend([1.0 for _ in ao_indices])

    # reorderign MO
    C_new = np.zeros_like(C)
    for i, (ao, aoc) in enumerate(zip(ao_ordering, ao_scaling)):
        C_new[i] = C[ao] * aoc

    fchk_info = {
        'file_type' : 'fchk',
        'atoms' : atoms,
        'atomnos' : atomnos,
        'coords' : coords,
        'coords_unit' : coords_unit,
        'mo_coeff' : C_new,
        'basis' : pyscf_basis,
        'cart' : cart
    }

    return fchk_info


def load_cube(inp_file):

    with open(inp_file) as f:
        lines = f.readlines()

    natoms = abs(int(lines[2].split()[0]))
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    nx, vx = (int(lines[3].split()[0]), np.array([float(x) for x in lines[3].split()[1:4]]))
    ny, vy = (int(lines[4].split()[0]), np.array([float(x) for x in lines[4].split()[1:4]]))
    nz, vz = (int(lines[5].split()[0]), np.array([float(x) for x in lines[5].split()[1:4]]))
    grid_vecs = np.vstack([vx, vy, vz])

    atomnos = []
    coords = []
    for i in range(natoms):
        parts = lines[6 + i].split()
        atomnos.append(int((parts[0])))
        coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
    atomnos = np.array(atomnos, dtype=int)
    coords = np.array(coords, dtype=float) # bohr

    if len(lines[6 + natoms].strip().split(" ")) < 6:
        raw_values = " ".join(lines[7 + natoms :]).replace("D", "E").replace("d", "E")
    else:
        raw_values = " ".join(lines[6 + natoms :]).replace("D", "E").replace("d", "E")
    data = np.fromstring(raw_values, sep=" ")
    expected_size = nx * ny * nz
    if data.size != expected_size:
        raise ValueError(
            f"cube grid size mismatch: expected {expected_size} values for shape ({nx}, {ny}, {nz}), got {data.size}"
        )
    mo_cube = np.array(data).reshape((nx, ny, nz), order="C")
    
    cube_info = {
        'file_type' : 'cube',
        'atomnos' : atomnos,
        'coords' : coords,
        'coords_unit' : 'Bohr',
        'origin' : origin,
        'grid_vecs' : grid_vecs,
        'mo_cube' : mo_cube
    }

    return cube_info

def load_pyscfchk(inp_file):

    mo_coeff = lib.chkfile.load(inp_file, 'scf/mo_coeff')
    mol = lib.chkfile.load_mol(inp_file)

    pyscfchk_info = {
        'file_type' : 'pyscfchk',
        'atoms' : mol._atom,
        'atomnos' : mol.atom_charges(),
        'coords' : mol.atom_coords(),
        'coords_unit' : 'Bohr',
        'mo_coeff' : mo_coeff,
        'basis' : mol.basis,
        'cart' : mol.cart
    }
    
    return pyscfchk_info
