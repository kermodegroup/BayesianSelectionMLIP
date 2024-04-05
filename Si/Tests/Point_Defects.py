import numpy as np
from ase.optimize.precon import PreconLBFGS
from ase.optimize import BFGSLineSearch, BFGS
from ase.constraints import ExpCellFilter
import os
from ase.build import bulk
from si_models import *
from ase.lattice.cubic import Diamond
from ase.optimize.precon.precon import Exp
from ase.atom import Atom

a0 = 5.46
si_bulk = Diamond(symbol='Si', latticeconstant=a0)
tol = 3e-3
N = 3

def tetrahedral_interstitial_energy(bulk, verbose=True):
    bulk_energy = bulk.get_potential_energy()
    Nat = len(bulk)
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())

    # add an atom to introduce an interstitial
    int_struct.append(Atom('Si', (0.001, 0.002, 5.44/2.0+0.003)))

    opt = PreconLBFGS(int_struct, logfile=None)
    opt.run(tol)

    # compute formation energy as difference of bulk and int energies
    e_form = int_struct.get_potential_energy() - bulk_energy*((Nat+1.0)/Nat)

    if verbose:
        print('bulk cell energy', bulk_energy)
        print('interstitial cell energy', int_struct.get_potential_energy())
        print('interstitial formation energy', e_form)
    return e_form

def dumbbell_interstitial_energy(bulk, verbose=True):
    bulk_energy = bulk.get_potential_energy()
    Nat = len(bulk)
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())

    # add an atom to introduce an interstitial
    int_struct.append(Atom('Si', (-0.5, 0.5, 5.44/2.0+1.0)))
    p = int_struct.get_positions()
    p[149,0] -= 1.0
    p[149,1] += 1.0
    p[149,2] -= 0.5
    int_struct.set_positions(p)

    opt = PreconLBFGS(int_struct, logfile=None)
    opt.run(tol)

    # compute formation energy as difference of bulk and int energies
    e_form = int_struct.get_potential_energy() - bulk_energy*((Nat+1.0)/Nat)

    if verbose:
        print('bulk cell energy', bulk_energy)
        print('interstitial cell energy', int_struct.get_potential_energy())
        print('interstitial formation energy', e_form)
    return e_form

def hexagonal_interstitial_energy(bulk, verbose=True):
    bulk_energy = bulk.get_potential_energy()
    Nat = len(bulk)
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())
    # add an atom to introduce an interstitial
    pos = int_struct.get_positions()
    int_pos = (pos[83,:] + pos[101,:] + pos[105,:] + pos[106,:] + pos[108,:] + pos[176,:])/6.0

    int_struct.append(Atom('Si', int_pos))

    rot_mat = np.zeros( (3,3) )
    rot_mat[0,:] = (1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0), 0.0)
    rot_mat[1,:] = (1.0/np.sqrt(6.0), 1.0/np.sqrt(6.0), -2.0/np.sqrt(6.0))
    rot_mat[2,:] = (1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0))
    cell = np.dot(int_struct.get_cell(), rot_mat.T)
    int_struct.set_cell(cell, scale_atoms=True)

    int_struct.arrays['move_mask_3'] = np.ones( (len(int_struct), 3), dtype=int )
    int_struct.arrays['move_mask_3'][len(int_struct)-1,2] = 0

    opt = PreconLBFGS(int_struct, logfile=None)
    opt.run(tol)

    # compute formation energy as difference of bulk and int energies
    e_form = int_struct.get_potential_energy() - bulk_energy*((Nat+1.0)/Nat)

    if verbose:
        print('bulk cell energy', bulk_energy)
        print('interstitial cell energy', int_struct.get_potential_energy())
        print('interstitial formation energy', e_form)
    return e_form

def vacancy_energy(bulk, verbose=True):
    bulk_energy = bulk.get_potential_energy()
    Nat = len(bulk)
    vac_struct = bulk.copy()
    vac_struct.set_calculator(bulk.get_calculator())

    # delete an atom to make a vacancy
    del vac_struct[-1]

    opt = PreconLBFGS(vac_struct, logfile=None)
    opt.run(tol)

    # compute formation energy as difference of bulk and int energies
    e_form = vac_struct.get_potential_energy() - bulk_energy*((Nat-1.0)/Nat)

    if verbose:
        print('bulk cell energy', bulk_energy)
        print('vacancy cell energy', vac_struct.get_potential_energy())
        print('vacancy formation energy', e_form)
    return e_form

def test_calc(calc, verbose=False):
    at = si_bulk.copy()
    at.calc = calc

    filter = ExpCellFilter(at, mask=[True]*3 + [False]*3)

    precon = Exp(3.0)
    opt = PreconLBFGS(filter, precon=precon)

    opt.run(fmax=1e-4, smax=1e-4)

    base_bulk = at.copy() * (N, N, N)
    del base_bulk.constraints
    base_bulk.calc = calc

    E_tetra = tetrahedral_interstitial_energy(base_bulk, verbose=verbose)

    E_dumbbell = dumbbell_interstitial_energy(base_bulk, verbose=verbose)

    E_hex = hexagonal_interstitial_energy(base_bulk, verbose=verbose)

    E_vac = vacancy_energy(base_bulk, verbose=verbose)

    return E_tetra, E_dumbbell, E_hex, E_vac

calc = original_gap()

dft_vals = [
    3.91, 3.66, 3.72, 3.67
]


Et, Ed, Eh, Ev = test_calc(calc)

print("E_tetra ", Et, ": ", 100 * (Et - dft_vals[0]) / dft_vals[0], "%")
print("E_dumbbell ", Ed, ": ", 100 * (Ed - dft_vals[1]) / dft_vals[1], "%")
print("E_hex ", Eh, ": ", 100 * (Eh - dft_vals[2]) / dft_vals[2], "%")
print("E_vac ", Ev, ": ", 100 * (Ev - dft_vals[3]) / dft_vals[3], "%")