import numpy as np
from matscipy.elasticity import fit_elastic_constants, elastic_moduli
from ase.optimize.precon import PreconLBFGS
from ase.optimize import BFGSLineSearch, BFGS
from ase.constraints import ExpCellFilter
import matplotlib.pyplot as plt
import os
from ase.units import GPa
from ase.build import bulk
from si_models import *
from ase.optimize.precon.precon import Exp


ats = bulk("Si", "diamond", cubic=True)

def test_calc(calc):
    Cs = np.zeros((6, 6))
    C_errs = np.zeros_like(Cs)

    at = ats.copy()
    at.calc = calc

    filter = ExpCellFilter(at, mask=[True]*3 + [False]*3)

    precon = Exp(3.0)
    opt = PreconLBFGS(filter, precon=precon)

    opt.run(fmax=1e-4, smax=1e-4)

    alat = np.average(np.diag(at.cell[:, :]))
    
    precon = Exp(3.0)
    
    opt = lambda atoms, **kwargs: PreconLBFGS(atoms, precon=precon, **kwargs)

    Cs, C_errs = fit_elastic_constants(
        at, optimizer=opt, fmax=1e-3, symmetry="cubic", verbose=False)#, N_steps=10, delta=1e-3)

    c = np.array([Cs[0, 0], Cs[0, 1], Cs[-1, -1]]) / GPa
    c_err = np.array([C_errs[0, 0], C_errs[0, 1], C_errs[3, 3]]) / GPa

    dat = np.array([c, c_err])

    # Fit elastic moduli
    E, nu, Gm, B, K = elastic_moduli(Cs)

    return alat, c, c_err, np.average(E)/GPa, nu, Gm, B/GPa, K


calc = original_gap()

alat, c, c_err, _, __, ___, B, ____ = test_calc(calc)

dft_vals = [
    153.3, 56.3, 72.2
]
C11, C12, C44 = c

print("Alat ", alat)
print("C11 ", C11, ": ", 100 * (C11 - dft_vals[0]) / dft_vals[0], "%")
print("C12 ", C12, ": ", 100 * (C12 - dft_vals[1]) / dft_vals[1], "%")
print("C44 ", C44, ": ", 100 * (C44 - dft_vals[2]) / dft_vals[2], "%")
print("B", B)