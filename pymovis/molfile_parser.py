import numpy as np
import periodictable
from pyscf import lib

###############################################################################
# fchk parser
###############################################################################

def read_fchk_array(lines : list[str], key : str, dtype=float) -> np.ndarray:
    """
    Fchk array reader.
    Search for the line containing the key, read the number of elements, 
    and then read the values until the expected number of elements is reached.

    Parameters:
        lines : list[str]
        key : str
        dtype : type

    Returns: 
        ndarray
    """

    for i, line in enumerate(lines):
        if key in line:
            n = int(line.split()[-1]) # number of elements in the array
            vals = []
            j = i + 1
            while len(vals) < n:
                vals.extend(lines[j].split())
                j += 1
            return np.array(vals[:n], dtype=dtype)

    raise ValueError(f"{key} not found")

shell_ao_map = {
     0 : [(0, [0])],
    -1 : [(0, [0]), (1, [0, 1, 2])], # SP shell, only S part is considered here, P part is handled separately
     1 : [(1, [0, 1, 2])],
    -2 : [(-2, [4, 2, 0, 1, 3])], #spherical d
     2 : [(2, [0, 3, 4, 1, 5, 2])], # cartecian d
    -3 : [(-3, [6, 4, 2, 0, 1, 3, 5])], # spherical f
     3 : [(3, [0, 4, 5, 3, 9, 6, 1, 8, 7, 2])], # cartesian f
    -4 : [(-4, [8, 6, 4, 2, 0, 1, 3, 5, 7])], # spherical g
     4 : [(4, [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0])] # cartesian g
}

d_const = np.sqrt(15 / (4 * np.pi))
f_const = np.sqrt(105 / (4 * np.pi))
g_const = np.sqrt(945 / (4 * np.pi))
ao_scaling_factors = {
    0 : [1.0], # s
    1 : [1.0, 1.0, 1.0], # p
    -2 : [1.0, 1.0, 1.0, 1.0, 1.0], # spherical d
    2 : [
        d_const / np.sqrt(3), #dxx
        d_const,              #dxy
        d_const,              #dxz
        d_const / np.sqrt(3), #dyy
        d_const,              #dyz
        d_const / np.sqrt(3)  #dzz
    ], # cartesian d
    -3 : [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], # spherical f
    3 : [
        f_const / np.sqrt(15), #fxxx
        f_const / np.sqrt(3),  #fxxy
        f_const / np.sqrt(3),  #fxxz
        f_const / np.sqrt(3),  #fxyy
        f_const,               #fxyz
        f_const / np.sqrt(3),  #fxzz
        f_const / np.sqrt(15), #fyyy
        f_const / np.sqrt(3),  #fyyz
        f_const / np.sqrt(3),  #fyzz
        f_const / np.sqrt(15)  #fzzz
    ], # cartesian f
    -4 : [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], # spherical g
    4 : [
        g_const / np.sqrt(105), #gxxxx
        g_const / np.sqrt(15),  #gxxxy
        g_const / np.sqrt(15),  #gxxxz
        g_const / np.sqrt(9),  #gxxyy
        g_const / np.sqrt(3),   #gxxyz
        g_const / np.sqrt(9),  #gxxzz
        g_const / np.sqrt(15), #gxyyy
        g_const / np.sqrt(3),  #gxyyz
        g_const / np.sqrt(3),  #gxyzz
        g_const / np.sqrt(15), #gxzzz
        g_const / np.sqrt(105), #gyyyy
        g_const / np.sqrt(15),  #gyyyz
        g_const / np.sqrt(9),  #gyyzz
        g_const / np.sqrt(15), #gyzzz
        g_const / np.sqrt(105)  #gzzzz
    ], # cartesian g
}

def load_fchk(inp_file : str) -> dict:
    """
    Load fchk file and extract necessary information for visualization.
    
    Parameters:
        inp_file : str
    Returns:
        dict with keys:
            'file_type' : 'fchk'
            'atoms' : list of [symbol, (x, y, z)]
            'atomnos' : array of atomic numbers
            'coords' : array of atomic coordinates
            'coords_unit' : 'Bohr' or 'Angstrom'
            'mo_coeff' : array of MO coefficients in AO basis, reordered and scaled for visualization
            'basis' : dict of basis information for each atom, compatible with PySCF format
            'cart' : bool, whether the basis functions are Cartesian or spherical
    """

    with open(inp_file) as f:
        lines = f.readlines()
    
    # ----------------------------------
    # parse geometries
    # ----------------------------------

    atomnos = read_fchk_array(lines, "Atomic numbers", int) # array of atomic numbers starts from 1
    coords = read_fchk_array(lines, "Current cartesian coordinates").reshape(-1, 3) # coordinates in Bohr
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
    
    # check consistency of basis information
    if not (len(shell_types) == len(nprim_per_shell) == len(shell_to_atom)):
        raise ValueError("Inconsistent shell information in fchk file")
    if not (len(primitive_exponents) == len(contraction_coeffs) == sum(nprim_per_shell)):
        raise ValueError("Inconsistent primitive information in fchk file")
    if any(stype == -1 for stype in shell_types) and p_con_coeffs is None:
        raise ValueError("SP shells detected but P(S=P) contraction coefficients are missing in fchk file")

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
    primitive_data = {}
    ao_data = {}    
    for ish, stype in enumerate(shell_types):
        nprim = nprim_per_shell[ish]
        atom_idx = shell_to_atom[ish] - 1 # 0 start index for python
        atom_symbol = atoms[atom_idx][0]
        exps = primitive_exponents[prim_ptr : prim_ptr + nprim]
        coeffs = contraction_coeffs[prim_ptr : prim_ptr + nprim]
        exps_coeffs = list(zip(exps, coeffs))

        if p_con_coeffs is not None:
            pcoeffs = p_con_coeffs[prim_ptr : prim_ptr + nprim]
            exps_pcoeffs = list(zip(exps, pcoeffs))
        else:
            pcoeffs = None
        
        prim_ptr += nprim

        # store primitive data for pyscf basis construction
        primitive_data.setdefault(atom_symbol, [])
        if stype == -1:  # SP shell
            # S shell
            s_block : list = [0]
            s_block.extend(exps_coeffs)
            primitive_data[atom_symbol].append(s_block)

            # P shell
            p_block : list = [1]
            p_block.extend(exps_pcoeffs)
            primitive_data[atom_symbol].append(p_block)

        else: # normal shells
            block = [abs(stype)]
            block.extend(exps_coeffs)
            primitive_data[atom_symbol].append(block)

        # store AO index for MO coefficient reordering and scaling
        ao_data.setdefault(atom_symbol, [])
        for shell_type, ao_map in shell_ao_map[stype]:
            ao_indices = [current_ao + i for i in ao_map]
            ao_data[atom_symbol].append((shell_type, ao_indices))
            current_ao += len(ao_map)            

    for atom_symbol in primitive_data:
        primitive_data[atom_symbol].sort(key=lambda x: x[0])

    ao_ordering = []
    ao_scaling = []
    for atom_symbol in ao_data:
        ao_data[atom_symbol].sort(key=lambda x: abs(x[0])) # sort shells by absolute angular momentum (keep S, P before D, etc.)
        for shell_type, ao_indices in ao_data[atom_symbol]:
            ao_ordering.extend(ao_indices)
            ao_scaling.extend(ao_scaling_factors[shell_type])
    
    # ----------------------------------
    # parse mo coefficients
    # ----------------------------------
    
    nao = current_ao
    coeffs = read_fchk_array(lines, "Alpha MO coefficients")
    nmo = coeffs.size // nao
    C = coeffs.reshape(nmo, nao).T # [ao, mo]

    # reorderign and scaling of MO coefficients for visualization
    C_new = np.zeros_like(C)
    for i, (ao, aoc) in enumerate(zip(ao_ordering, ao_scaling)):
        C_new[i] = C[ao] * aoc


    # return all parsed information in a dict for later use
    fchk_info = {
        'file_type' : 'fchk',
        'atoms' : atoms,
        'atomnos' : atomnos,
        'coords' : coords,
        'coords_unit' : coords_unit,
        'mo_coeff' : C_new,
        'basis' : primitive_data,
        'cart' : cart
    }

    return fchk_info


###############################################################################
# cube parser
###############################################################################

def load_cube(inp_file: str) -> dict:

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


###############################################################################
# pyscf chkfile parser
###############################################################################

def load_pyscfchk(inp_file: str) -> dict:

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
