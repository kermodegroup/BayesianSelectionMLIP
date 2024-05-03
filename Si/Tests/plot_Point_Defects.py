import os
import json
import matplotlib.pyplot as plt
import numpy as np

# DFT vals for Et, Ed, Eh (see https://github.com/libAtoms/silicon-testing-framework/blob/master/test-results/model-CASTEP_ASE-test-bulk_diamond-properties.json)
# DFT vals for Ev (see https://github.com/libAtoms/silicon-testing-framework/blob/master/test-results/model-CASTEP_ASE-test-vacancy-energy-properties.json)
dft_values = [
    3.909861789215938, 3.661610819210182, 3.7243628192154574, 3.6732338007932412
]

props = ["Et", "Ed", "Eh", "Ev"]
full_names = ["Tetrahedral Interstitial", "Dumbbell Intersitital", "Hexagonal Intersitital", "Vacancy"]

models = [
    "MONTECARLO_ACE",
    "ACECURMEAN_ACE",
    "ACEAVGCUR_ACE",
    #"HALFORCE_ACE",
    "ACEAVGKMED_ACE"
]

fig, ax = plt.subplots(2, 2, figsize=(15, 15))

ax = ax.flatten()

for i, prop in enumerate(full_names):
    ax[i].set_title(prop)
    ax[i].set_xlabel("N_struct")
    ax[i].set_ylabel(f"Error from DFT (eV)")
    ax[i].set_yscale("log")

for m, model in enumerate(models):
    mod, pot = model.split("_")
    if not os.path.exists(f"../Test_Results/{mod}/{mod}_Point_Defects_{pot}.json"):
        # Skip as data not available
        continue

    with open(f"../Test_Results/{mod}/{mod}_Point_Defects_{pot}.json", "r") as f:
        data = json.load(f)

    Ns = [int(k) for k in data.keys()]

    for i, prop in enumerate(props):
        errs = [np.sqrt(np.mean((np.array(ndata[prop + "_raw_vals"]) - dft_values[i])**2)) if ndata[prop + "_raw_vals"] is not None else np.inf for ndata in data.values()]

        #max_err = [np.max(np.abs(np.array(ndata[prop + "_raw_vals"]) - dft_values[i]))for ndata in data.values()]
        #min_err = [np.min(np.abs(np.array(ndata[prop + "_raw_vals"]) - dft_values[i]))for ndata in data.values()]

        ax[i].plot(Ns, errs, label=model, color=f"C{m}", marker="o")
        #ax[i].plot(Ns, max_err, color=f"C{m}", marker="o", linestyle="dashed")
        #ax[i].plot(Ns, min_err, color=f"C{m}", marker="o", linestyle="dashed")

ax[0].legend()
plt.tight_layout()
plt.savefig("../Plots/Point_Defects.png", dpi=200)