import os
import json
import numpy as np
import matplotlib.pyplot as plt

dft_vals = {
        "Li3P": {
            "E": 0.0,
            "Vol": 14.677966179334197
        },
        "Li7PS2": {
            "E": 41.0,
            "Vol": 16.39860899466704
        },
        "Li5PS": {
            "E": 39.8,
            "Vol": 14.719078222668012
        },
        "Li8P2S": {
            "E":39.1,
            "Vol": 19.061889481376785
        },
        "Li11P3S": {
            "E": 31.4,
            "Vol": 13.95505733358334
        },
        "Li2S": {
            "E": 0.0,
            "Vol": 15.509496826688414
        }
    }

props = ["Li7PS2", "Li5PS", "Li8P2S", "Li11P3S"]

models = [
    "MONTECARLO",
    "ACECURMEAN",
    "ACEAVGCUR",
    #"HALFORCE",
    "ACEAVGKMED"
]

fig, ax = plt.subplots(2, 2, figsize=(15, 15))

ax = ax.flatten()

for i, prop in enumerate(props):
    ax[i].set_title(prop)
    ax[i].set_xlabel("N_struct")
    ax[i].set_ylabel(f"Error from DFT (%)")
    #ax[i].set_yscale("log")

for m, mod in enumerate(models):
    if not os.path.exists(f"../Test_Results/{mod}/{mod}_Hull_Energies.json"):
        # Skip as data not available
        continue

    with open(f"../Test_Results/{mod}/{mod}_Hull_Energies.json", "r") as f:
        data = json.load(f)

    Ns = [int(k) for k in data.keys()]

    for i, prop in enumerate(props):

        errs = np.zeros((len(Ns)))

        for j, N in enumerate(Ns):
            ndata = data[str(N)][prop]
            
            x = ndata["S_frac"] / data[str(N)]["Li2S"]["S_frac"]
            
            Es = np.zeros((len(ndata["Es"])))
            for k in range(len(ndata["Es"])):
                Es[k] = ndata["Es"][k] - x * data[str(N)]["Li2S"]["Es"][k] - (1-x) * data[str(N)]["Li3P"]["Es"][k]

            errs[j] = np.sqrt(np.mean((Es[np.abs(Es) < 10] - dft_vals[prop]["E"])**2))

        ax[i].plot(Ns, 100 * errs / dft_vals[prop]["E"], label=mod, marker="x")

ax[0].legend()
plt.savefig("../Plots/Hull_Energies.png")

plt.clf()

################

fig, ax = plt.subplots(2, 2, figsize=(15, 15))

ax = ax.flatten()

for i, prop in enumerate(props):
    ax[i].set_title(prop)
    ax[i].set_xlabel("N_struct")
    ax[i].set_ylabel(f"Error from DFT (%)")
    #ax[i].set_yscale("log")

for m, mod in enumerate(models):
    if not os.path.exists(f"../Test_Results/{mod}/{mod}_Hull_Energies.json"):
        # Skip as data not available
        continue

    with open(f"../Test_Results/{mod}/{mod}_Hull_Energies.json", "r") as f:
        data = json.load(f)

    Ns = [int(k) for k in data.keys()]

    for i, prop in enumerate(props):

        errs = np.zeros((len(Ns)))

        for j, N in enumerate(Ns):
            ndata = data[str(N)][prop]

            vols = np.array(ndata["Vols"])

            errs[j] = np.sqrt(np.mean((vols[np.abs(vols) > 5] - dft_vals[prop]["Vol"])**2))

        ax[i].plot(Ns, 100 * errs / dft_vals[prop]["Vol"], label=mod, marker="x")

ax[0].legend()
plt.savefig("../Plots/Volumes.png")
