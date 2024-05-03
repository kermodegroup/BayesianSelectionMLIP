import os
import json
import matplotlib.pyplot as plt
import numpy as np

# DFT vals for E100, E110, E111 (see )
dft_values = [
    2.17, 1.52, 1.57
]

props = ["E100", "E110", "E111"]

models = [
    #"MONTECARLO_GAP",
    #"SOAPCURMEAN_GAP",
    "MONTECARLO_ACE",
    "ACECURMEAN_ACE",
    "ACEAVGCUR_ACE",
    #"HALFORCE_ACE",
    "ACEAVGKMED_ACE"
]

fig, ax = plt.subplots(2, 2, figsize=(15, 15))

ax = ax.flatten()

for i, prop in enumerate(props):
    ax[i].set_title(prop)
    ax[i].set_xlabel("N_struct")
    ax[i].set_ylabel(f"Error from DFT (J/m^2)")
    ax[i].set_yscale("log")

ax[-1].axis("off")

for m, model in enumerate(models):
    mod, pot = model.split("_")
    if not os.path.exists(f"../Test_Results/{mod}/{mod}_Surfaces_{pot}.json"):
        # Skip as data not available
        continue

    with open(f"../Test_Results/{mod}/{mod}_Surfaces_{pot}.json", "r") as f:
        data = json.load(f)

    Ns = [int(k) for k in data.keys()]

    for i, prop in enumerate(props):
        errs = [np.sqrt(np.mean((np.array(ndata[prop + "_raw_vals"]) - dft_values[i])**2))for ndata in data.values()]

        max_err = [np.max(np.abs(np.array(ndata[prop + "_raw_vals"]) - dft_values[i]))for ndata in data.values()]
        min_err = [np.min(np.abs(np.array(ndata[prop + "_raw_vals"]) - dft_values[i]))for ndata in data.values()]

        ax[i].plot(Ns, errs, label=model, color=f"C{m}", marker="o")
        #ax[i].plot(Ns, max_err, color=f"C{m}", marker="o", linestyle="dashed")
        #ax[i].plot(Ns, min_err, color=f"C{m}", marker="o", linestyle="dashed")

ax[0].legend()
plt.tight_layout()
plt.savefig("../Plots/Surfaces.png", dpi=200)