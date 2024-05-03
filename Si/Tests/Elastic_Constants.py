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
import json

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
        at, optimizer=opt, fmax=1e-3, symmetry="cubic", verbose=False)

    c = np.array([Cs[0, 0], Cs[0, 1], Cs[-1, -1]]) / GPa
    c_err = np.array([C_errs[0, 0], C_errs[0, 1], C_errs[3, 3]]) / GPa

    dat = np.array([c, c_err])

    # Fit elastic moduli
    E, nu, Gm, B, K = elastic_moduli(Cs)

    return alat, c, c_err, np.average(E)/GPa, nu, Gm, B/GPa, K



calc_type = "ACE"

method = "ACEAVGCUR"
Nc = 20

if calc_type == "GAP":
    calc_fn = test_gap
elif calc_type == "ACE":
    calc_fn = test_ace

if os.path.exists(f"../Test_Results/{method}/{method}_Elastic_Constants_{calc_type}.json"):
    with open(f"../Test_Results/{method}/{method}_Elastic_Constants_{calc_type}.json", "r") as f:
        data = json.load(f)
else:
    data = {}

print(method)
for N in [11, 22, 44, 88, 121, 242]:#, 484]:
    print(N)
    
    data[str(N)] = {}

    # alat, elastic constants
    alats = np.zeros((Nc))
    C11s = np.zeros_like(alats)
    C12s = np.zeros_like(alats)
    C44s = np.zeros_like(alats)

    for i in range(Nc):
        calc = calc_fn(method, N, [i])[0]
        alat, c = test_calc(calc)[:2]
        C11, C12, C44 = c
        
        
        alats[i] = alat
        C11s[i] = C11
        C12s[i] = C12
        C44s[i] = C44

    data[str(N)]["alat_mean"] = np.mean(alats)
    data[str(N)]["C11_mean"] = np.mean(C11s)
    data[str(N)]["C12_mean"] = np.mean(C12s)
    data[str(N)]["C44_mean"] = np.mean(C44s)

    data[str(N)]["alat_std"] = np.std(alats)
    data[str(N)]["C11_std"] = np.std(C11s)
    data[str(N)]["C12_std"] = np.std(C12s)
    data[str(N)]["C44_std"] = np.std(C44s)

    data[str(N)]["alat_raw_vals"] = list(alats)
    data[str(N)]["C11_raw_vals"] = list(C11s)
    data[str(N)]["C12_raw_vals"] = list(C12s)
    data[str(N)]["C44_raw_vals"] = list(C44s)

os.makedirs(f"../Test_Results/{method}", exist_ok=True)

# Ensure keys are sorted as integers
Ns = sorted([int(k) for k in data.keys()])

data = {N: data[str(N)] for N in Ns}

with open(f"../Test_Results/{method}/{method}_Elastic_Constants_{calc_type}.json", "w") as f:
    json.dump(data, f, indent=4)