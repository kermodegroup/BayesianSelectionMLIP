import json
import os
import matplotlib.pyplot as plt
import numpy as np

models = [
    #"MONTECARLO_GAP",
    #"SOAPCURMEAN_GAP",
    "MONTECARLO_ACE",
    "ACECURMEAN_ACE",
    "ACEAVGCUR_ACE",
    #"HALFORCE_ACE",
    "ACEAVGKMED_ACE"
]

data = {}

for m, model in enumerate(models):
    mod, pot = model.split("_")
    if not os.path.exists(f"../Test_Results/{mod}/{mod}_Dataset_Errs_{pot}.json"):
        # Skip as data not available
        continue

    with open(f"../Test_Results/{mod}/{mod}_Dataset_Errs_{pot}.json", "r") as f:
        model_data = json.load(f)

    for key in model_data.keys():
        if key not in data.keys():
            data[key] = {}

        data[key][model] = model_data[key]

Ns = list(data.keys())
N_ints = [int(n) for n in Ns]
models = list(data[Ns[0]].keys())
config_types = list(data[Ns[0]][models[0]].keys())

for config_type in config_types:
    plt.clf()
    fig, ax = plt.subplots(nrows=2, figsize=(10, 15), sharex=True)

    for m, model in enumerate(models):
        E_rmses = np.array([data[n][model][config_type]["E_rmse_mean"] for n in Ns if (model in data[n].keys())]) * 1000
        F_rmses = np.array([data[n][model][config_type]["F_rmse_mean"] for n in Ns if (model in data[n].keys())]) * 1000

        N_ints = [n for n in Ns if (model in data[n].keys())]
        # E_rmse_maxs = np.array([np.max(data[n][model][config_type]["E_rmse_raw_vals"]) for n in Ns if (model in data[n].keys())]) * 1000
        # F_rmse_maxs = np.array([np.max(data[n][model][config_type]["F_rmse_raw_vals"]) for n in Ns if (model in data[n].keys())]) * 1000

        # E_rmse_mins = np.array([np.min(data[n][model][config_type]["E_rmse_raw_vals"]) for n in Ns if (model in data[n].keys())]) * 1000
        # F_rmse_mins = np.array([np.min(data[n][model][config_type]["F_rmse_raw_vals"]) for n in Ns if (model in data[n].keys())]) * 1000

        ax[0].plot(N_ints, E_rmses, label=model, color=f"C{m}", marker="o")
        #ax[0].plot(N_ints, E_rmse_maxs, color=f"C{m}", marker="o", linestyle="dashed")
        #ax[0].plot(N_ints, E_rmse_mins, color=f"C{m}", marker="o", linestyle="dashed")

        ax[1].plot(N_ints, F_rmses, label=model, color=f"C{m}", marker="o")
        #ax[1].plot(N_ints, F_rmse_maxs, color=f"C{m}", marker="o", linestyle="dashed")
        #ax[1].plot(N_ints, F_rmse_mins, color=f"C{m}", marker="o", linestyle="dashed")

    ax[1].set_xlabel("N_struct")
    ax[0].set_ylabel(f"Energy RMSEs (meV/Atom)")
    ax[1].set_ylabel(f"Force RMSEs (meV/Ang)")

    ax[0].set_title("Energy Errors")
    ax[1].set_title("Force Errors")

    ax[0].set_yscale("log")
    ax[1].set_yscale("log")
    ax[0].legend()
    plt.tight_layout()

    plt.savefig(f"../Plots/Dataset_Errs/{config_type}_Dataset_Errs.png")