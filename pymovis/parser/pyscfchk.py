"""PySCF checkpoint parser."""

from typing import Any

from pyscf import lib

###############################################################################
# pyscf chkfile parser
###############################################################################

def load_pyscfchk(inp_file: str) -> dict[str, Any]:
    """Read geometry and alpha/beta SCF orbitals from a PySCF checkpoint."""

    mo_coeff = lib.chkfile.load(inp_file, 'scf/mo_coeff')
    mo_occ = lib.chkfile.load(inp_file, 'scf/mo_occ')
    mo_energy = lib.chkfile.load(inp_file, 'scf/mo_energy')
    mol = lib.chkfile.load_mol(inp_file)
    if isinstance(mo_coeff, (tuple, list)):
        if not isinstance(mo_occ, (tuple, list)) or not isinstance(mo_energy, (tuple, list)):
            raise ValueError("spin-separated MO data is incomplete in PySCF checkpoint")
        mo_coeff_alpha, mo_coeff_beta = mo_coeff
        mo_occ_alpha, mo_occ_beta = mo_occ
        mo_energy_alpha, mo_energy_beta = mo_energy
    else:
        mo_coeff_alpha, mo_coeff_beta = mo_coeff, None
        mo_occ_alpha, mo_occ_beta = mo_occ, None
        mo_energy_alpha, mo_energy_beta = mo_energy, None

    pyscfchk_info = {
        'file_type' : 'pyscfchk',
        'atom_symbols' : [atom[0] for atom in mol._atom],
        'atomnos' : mol.atom_charges(),
        'coords' : mol.atom_coords(),
        'coords_unit' : 'Bohr',
        'mo_coeff' : mo_coeff_alpha,
        'mo_coeff_beta' : mo_coeff_beta,
        'mo_occ' : mo_occ_alpha,
        'mo_occ_beta' : mo_occ_beta,
        'mo_energy' : mo_energy_alpha,
        'mo_energy_beta' : mo_energy_beta,
        'charge' : mol.charge,
        'multiplicity' : (mol.spin or 0) + 1,
        'basis' : mol.basis,
        'shell_type' : 'cart' if mol.cart else 'sph'
    }

    return pyscfchk_info
    
