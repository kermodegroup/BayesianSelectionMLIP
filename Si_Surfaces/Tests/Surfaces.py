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
from ase.units import _e
from ase.io import write
import json

a0 = 5.46
tol = 3e-3

bulks = [
    Diamond(symbol='Si', latticeconstant=a0, directions=[[1,1,0],[1,-1,0],[0,0,1]])*(1,1,10), # 100
    Diamond(symbol="Si", latticeconstant=a0, directions=[[1,-1,0],[0,0,1],[1,1,0]])*(1,1,10), # 110
    Diamond(symbol="Si", latticeconstant=a0, directions=[[1,-1,0],[1,0,-1],[1,1,1]])*(1,1,10) # 111
]

cell_disps = [
    -10,
    -10,
    +10
]

offsets = [
    0.0, 
    1.0, 
    2.0
]

names = [
    "100", "110", "111"
]

def test_calc(calc):

    e_forms = []

    for i in range(3):
        bulk = bulks[i].copy()
        disp = cell_disps[i]
        offset = offsets[i]

        bulk.calc = calc
        Nat = len(bulk)

        bulk.positions[:,2] += offset
        bulk.wrap()

        #write("Surface_" + names[i] + "_bulk.xyz", bulk)


        # compute surface formation energy as half the difference of bulk and expanded cell
        ebulk = bulk.get_potential_energy()

        bulk.cell[2,:] += [0.0,0.0,disp]

        if i == 1:
            c = bulk.get_cell()
            t_v = c[0,:].copy()
            c[0,:] = c[1,:]
            c[1,:] = t_v
            bulk.set_cell(c)
            np.random.seed(75)

            #bulk.positions += (np.random.rand((Nat*3))*0.1).reshape([Nat,3])


        # relax expanded cell
        opt = PreconLBFGS(bulk, logfile=None)
        opt.run(tol)

        eexp  = bulk.get_potential_energy()
        
        e_forms.append(0.5*(eexp - ebulk) / np.linalg.norm(np.cross(bulk.cell[0,:],bulk.cell[1,:])))
    return np.array(e_forms) * _e * 1e20


calc_type = "ACE"

method = "ACECURMEAN"
Nc = 20

if calc_type == "GAP":
    calc_fn = test_gap
elif calc_type == "ACE":
    calc_fn = test_ace

if os.path.exists(f"../Test_Results/{method}/{method}_Surfaces_{calc_type}.json"):
    with open(f"../Test_Results/{method}/{method}_Surfaces_{calc_type}.json", "r") as f:
        data = json.load(f)
else:
    data = {}

for N in [5, 10, 20, 50, 100]:
    print(N)
    
    data[str(N)] = {}

    # dft_vals = [
    #     153.3, 56.3, 72.2
    # ]

    # alat, elastic constants
    E100s = np.zeros((Nc))
    E110s = np.zeros_like(E100s)
    E111s = np.zeros_like(E100s)

    for i in range(Nc):
        calc = calc_fn(method, N, [i])[0]
        E100, E110, E111 = test_calc(calc)
        
        
        E100s[i] = E100
        E110s[i] = E110
        E111s[i] = E111

    data[str(N)]["E100_mean"] = np.mean(E100s)
    data[str(N)]["E110_mean"] = np.mean(E110s)
    data[str(N)]["E111_mean"] = np.mean(E111s)

    data[str(N)]["E100_std"] = np.std(E100s)
    data[str(N)]["E110_std"] = np.std(E110s)
    data[str(N)]["E111_std"] = np.std(E111s)

    data[str(N)]["E100_raw_vals"] = list(E100s)
    data[str(N)]["E110_raw_vals"] = list(E110s)
    data[str(N)]["E111_raw_vals"] = list(E111s)

os.makedirs(f"../Test_Results/{method}", exist_ok=True)

# Ensure keys are sorted as integers
Ns = sorted([int(k) for k in data.keys()])

data = {N: data[str(N)] for N in Ns}
with open(f"../Test_Results/{method}/{method}_Surfaces_{calc_type}.json", "w") as f:
    json.dump(data, f, indent=4)