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
import json
from ase.io import read

a0 = 5.46
si_bulk = Diamond(symbol='Si', latticeconstant=a0)
tol = 3e-3
N_reps = 3


def relax_struct(ats):
    opt = PreconLBFGS(ats, logfile=None)
    opt.run(tol, steps=500)

def tetrahedral_interstitial_energy(bulk, verbose=True):
    bulk_energy = bulk.get_potential_energy()
    Nat = len(bulk)
    int_struct = bulk.copy()
    int_struct.set_calculator(bulk.get_calculator())

    # add an atom to introduce an interstitial
    int_struct.append(Atom('Si', (0.001, 0.002, 5.44/2.0+0.003)))
    
    # int_struct = read("PD_Structs/tetra.xyz")
    # int_struct.calc = bulk.calc

    relax_struct(int_struct)

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

    # int_struct = read("PD_Structs/dumbbell.xyz")
    # int_struct.calc = bulk.calc

    relax_struct(int_struct)

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

    # int_struct = read("PD_Structs/hex.xyz")
    # int_struct.calc = bulk.calc

    relax_struct(int_struct)

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
    
    # vac_struct = read("PD_Structs/dumbbell.xyz")
    # vac_struct.calc = bulk.calc

    relax_struct(vac_struct)

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

    base_bulk = at.copy() * (N_reps, N_reps, N_reps)
    del base_bulk.constraints
    base_bulk.calc = calc

    E_tetra = tetrahedral_interstitial_energy(base_bulk, verbose=verbose)

    E_dumbbell = dumbbell_interstitial_energy(base_bulk, verbose=verbose)

    E_hex = hexagonal_interstitial_energy(base_bulk, verbose=verbose)

    E_vac = vacancy_energy(base_bulk, verbose=verbose)

    return E_tetra, E_dumbbell, E_hex, E_vac


calc_type = "ACE"

method = "ACEAVGCUR"
Nc = 20

if calc_type == "GAP":
    calc_fn = test_gap
elif calc_type == "ACE":
    calc_fn = test_ace

# Try to load existing datafile
if os.path.exists(f"../Test_Results/{method}/{method}_Point_Defects_{calc_type}.json"):
    with open(f"../Test_Results/{method}/{method}_Point_Defects_{calc_type}.json", "r") as f:
        data = json.load(f)
else:
    data = {}

print(method)
for N in [11, 22, 44, 88, 121, 242]:#, 484, 726]:
    print(N)
    
    data[str(N)] = {}

    # Tetra, Dumbbell, Hexagonal Interstitials + Vacancy energies
    Ets = np.zeros((Nc))
    Eds = np.zeros_like(Ets)
    Ehs = np.zeros_like(Ets)
    Evs = np.zeros_like(Ets)

    for i in range(Nc):
        calc = calc_fn(method, N, [i])[0]
        try:
            Et, Ed, Eh, Ev = test_calc(calc)
            Ets[i] = Et
            Eds[i] = Ed
            Ehs[i] = Eh
            Evs[i] = Ev
        except RuntimeError:
            continue

    data[str(N)]["Et_mean"] = np.mean(Ets[Ets != 0])
    data[str(N)]["Ed_mean"] = np.mean(Eds[Ets != 0])
    data[str(N)]["Eh_mean"] = np.mean(Ehs[Ets != 0])
    data[str(N)]["Ev_mean"] = np.mean(Evs[Ets != 0])

    data[str(N)]["Et_std"] = np.std(Ets[Ets != 0])
    data[str(N)]["Ed_std"] = np.std(Eds[Ets != 0])
    data[str(N)]["Eh_std"] = np.std(Ehs[Ets != 0])
    data[str(N)]["Ev_std"] = np.std(Evs[Ets != 0])

    data[str(N)]["Et_raw_vals"] = list(Ets[Ets != 0])
    data[str(N)]["Ed_raw_vals"] = list(Eds[Ets != 0])
    data[str(N)]["Eh_raw_vals"] = list(Ehs[Ets != 0])
    data[str(N)]["Ev_raw_vals"] = list(Evs[Ets != 0])

os.makedirs(f"../Test_Results/{method}", exist_ok=True)

# Ensure keys are sorted as integers
Ns = sorted([int(k) for k in data.keys()])

data = {N: data[str(N)] for N in Ns}

with open(f"../Test_Results/{method}/{method}_Point_Defects_{calc_type}.json", "w") as f:
    json.dump(data, f, indent=4)