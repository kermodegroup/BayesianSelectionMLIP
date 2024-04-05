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
        print(names[i], cell_disps[i])
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

            bulk.positions += (np.random.rand((Nat*3))*0.1).reshape([Nat,3])

        #write("Surface_" + names[i] + "_unrelaxed.xyz", bulk)

        # relax expanded cell
        opt = PreconLBFGS(bulk, logfile=None)
        opt.run(tol)

        
        #write("Surface_" + names[i] + "_relaxed.xyz", bulk)

        eexp  = bulk.get_potential_energy()
        
        e_forms.append(0.5*(eexp - ebulk) / np.linalg.norm(np.cross(bulk.cell[0,:],bulk.cell[1,:])))
    return np.array(e_forms) * _e * 1e20

calc = original_gap()

dft_vals = [
    2.17,
    1.52,
    1.57
]


E100, E110, E111 = test_calc(calc)



print("E_100 ", E100, ": ", 100 * (E100 - dft_vals[0]) / dft_vals[0], "%")
print("E_110 ", E110, ": ", 100 * (E110 - dft_vals[1]) / dft_vals[1], "%")
print("E_111 ", E111, ": ", 100 * (E111 - dft_vals[2]) / dft_vals[2], "%")