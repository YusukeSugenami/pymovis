import itertools
import numpy as np
from pyscf import gto, scf
import time

def get_shell_indices(mol, target_l):
    """
    指定 angular momentum shell の AO index を取得

    Parameters
    ----------
    target_l : int
        2=d
        3=f

    Returns
    -------
    shells : list of ndarray
    """
    ao_loc = mol.ao_loc_nr()
    shells = []

    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        if l != target_l:
            continue
        start = ao_loc[ib]

        # spherical only
        if l == 2:
            nfunc = 6
        elif l == 3:
            nfunc = 10
        else:
            continue

        end = start + nfunc
        shells.append(np.arange(start, end))

    return shells


def make_global_permutation(nao, d_shells, f_shells, d_perm, f_perm):
    """
    d/f permutation を全原子に適用
    """

    perm = np.arange(nao)

    # d shell
    for idx in d_shells:
        idx2 = idx[list(d_perm)]
        perm[idx] = idx2

    # f shell
    for idx in f_shells:
        idx2 = idx[list(f_perm)]
        perm[idx] = idx2

    return perm


def permutation_matrix(perm):
    """
    permutation vector -> matrix
    """

    n = len(perm)
    P = np.zeros((n, n))
    for i, j in enumerate(perm):
        P[i, j] = 1.0

    return P


def orthogonality_error(C, S):
    """
    || C^T S C - I ||
    """

    X = C.T @ S @ C
    I = np.eye(X.shape[0])

    return np.linalg.norm(X - I)


def search_permutation(mol, C):
    """
    d/f permutation を探索
    """

    print('search valid mo ordering')
    print(mol.basis)
    print(mol.ao_labels())
    print(C.shape)

    start = time.perf_counter()
    S = mol.intor("int1e_ovlp")
    nao = mol.nao_nr()
    d_shells = get_shell_indices(mol, 2)
    f_shells = get_shell_indices(mol, 3)

    print(d_shells)
    print(f_shells)

    # spherical 前提
    d_perms = list(itertools.permutations(range(6)))
    f_perms = list(itertools.permutations(range(10)))

    best_error = 1e100

    best_d_perm = None
    best_f_perm = None
    best_C = None

    ##################################################
    # alternating optimization
    ##################################################

    current_f_perm = tuple(range(10))

    for iteration in range(10):

        # optimize d
        print('d ordering')
        best_d_error = 1e100
        for d_perm in d_perms:
            perm = make_global_permutation(
                nao,
                d_shells,
                f_shells,
                d_perm,
                current_f_perm
            )
            P = permutation_matrix(perm)
            C_test = P.T @ C
            err = orthogonality_error(C_test, S)
            if err < best_d_error:
                best_d_error = err
                best_d_perm = d_perm

        # optimize f
        if f_shells:
            #print('f ordering')
            best_f_error = 1e100
            for f_perm in f_perms:
                perm = make_global_permutation(
                    nao,
                    d_shells,
                    f_shells,
                    best_d_perm,
                    f_perm
                )
                P = permutation_matrix(perm)
                C_test = P.T @ C
                err = orthogonality_error(C_test, S)
                if err < best_f_error:
    
                    best_f_error = err
                    best_f_perm = f_perm

            current_f_perm = best_f_perm
            perm = make_global_permutation(
                nao,
                d_shells,
                f_shells,
                best_d_perm,
                best_f_perm
            )
    
            P = permutation_matrix(perm)
            C_best = P.T @ C
            total_error = orthogonality_error(C_best, S)
    
            print(total_error)
            if abs(best_error - total_error) < 1e-10:
                break
    
            best_error = total_error
            best_C = C_best

        else:
            best_f_perm = []
            perm = make_global_permutation(
                nao,
                d_shells,
                f_shells,
                best_d_perm,
                best_f_perm
            )
            P = permutation_matrix(perm)
            C_best = P.T @ C
            best_error = best_d_error
            best_C = C_best
            break
    
    print(best_d_perm)

    end = time.perf_counter()
    print(f'reordering time: {end - start}')
    print(f'{best_C.T @ S @ best_C}')
    return {
        "d_perm": best_d_perm,
        "f_perm": best_f_perm,
        "error": best_error,
        "C_corrected": best_C
    }


def test():
    mol = gto.M(
        atom="""
        Fe 0 0 0
        """,
        basis="cc-pvtz",
        cart=False
    )
    mf = scf.RHF(mol).run()
    C_true = mf.mo_coeff.copy()
    
    ############################################################
    # fake Gaussian ordering
    ############################################################
    
    nao = mol.nao_nr()
    d_shells = get_shell_indices(mol, 2)
    f_shells = get_shell_indices(mol, 3)
    fake_d_perm = [2, 0, 1, 4, 3]
    fake_f_perm = [3, 0, 1, 6, 2, 5, 4]
    
    perm_fake = make_global_permutation(
        nao,
        d_shells,
        f_shells,
        fake_d_perm,
        fake_f_perm
    )
    
    P_fake = permutation_matrix(perm_fake)
    C_gaussian = P_fake @ C_true
    
    result = search_permutation(mol, C_gaussian)
    
    print("Recovered d permutation:")
    print(result["d_perm"])
    print("Recovered f permutation:")
    print(result["f_perm"])
    print("Final error:")
    print(result["error"])

