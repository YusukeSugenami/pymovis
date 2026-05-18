from pyscf import lib

def load_pyscf(inp_file : str) -> dict:

    mo_coeff = lib.chkfile.load(inp_file, 'scf/mo_coeff')
    mol = lib.chkfile.load_mol(inp_file)

    pyscf_info = {
        'file_type' : 'pyscfchk',
        'atoms' : mol._atom,
        'atomnos' : mol.atom_charges(),
        'coords' : mol.atom_coords(),
        'coords_unit' : 'Bohr',
        'mo_coeff' : mo_coeff,
        'basis' : mol.basis,
        'cart' : mol.cart
    }
    
    return pyscf_info