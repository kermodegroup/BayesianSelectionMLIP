from ase.io import read, write
from lips_models import test_ace
import numpy as np
from ase.optimize.precon import PreconLBFGS
from ase.optimize.precon.precon import Exp
from ase.constraints import ExpCellFilter
import os
import json

# Stoichiometries from
# https://pubs.acs.org/doi/epdf/10.1021/jacs.2c01913
structs = [
    "Li3P",
    "Li7PS2",
    "Li5PS",
    "Li8P2S",
    "Li11P3S",
    "Li2S"
]

structs = {
    name : read(f"Relaxed_Structs/{name}/structs.xyz", index="-1") for name in structs
}

method = "MONTECARLO"
nc = 20

if os.path.exists(f"../Test_Results/{method}/{method}_Hull_Energies.json"):
    with open(f"../Test_Results/{method}/{method}_Hull_Energies.json", "r") as f:
        data = json.load(f)
else:
    data = {}

for N in [20, 50, 100, 200]:
    data[str(N)] = {}

    for name, struct in structs.items():
            
        data[str(N)][name] = {}

        spec = np.array(struct.get_chemical_symbols())
        data[str(N)][name]["S_frac"] = np.sum(spec == "S") / len(spec)
        
        data[str(N)][name]["Es"] = []
        data[str(N)][name]["Vols"] = []

        for i in range(nc):
            calc = test_ace(method, N, [i])[0]
            ats = struct.copy()
            ats.calc = calc
            try:
                filter = ExpCellFilter(ats, mask=[True]*6)

                precon = Exp(3.0)
                opt = PreconLBFGS(filter, precon=precon)

                opt.run(fmax=1e-3, smax=1e-3, steps=100)

                data[str(N)][name]["Es"].append(ats.get_potential_energy() / len(ats))
                data[str(N)][name]["Vols"].append(ats.get_volume() / len(ats))
            except RuntimeError:
                pass


os.makedirs(f"../Test_Results/{method}", exist_ok=True)

# Ensure keys are sorted as integers
Ns = sorted([int(k) for k in data.keys()])

data = {N: data[str(N)] for N in Ns}

with open(f"../Test_Results/{method}/{method}_Hull_Energies.json", "w") as f:
    json.dump(data, f, indent=4)