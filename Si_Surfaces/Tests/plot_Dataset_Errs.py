import json
import os
import matplotlib.pyplot as plt
import numpy as np
from plot_config import *
import matplotlib as mpl


###
# Plots Dataset RMSE figures (Figs. 1 & 3)


# File extension for plots, fed to plt.savefig()
plot_extension = ".eps"


plt.rcParams["axes.labelsize"] = 20
plt.rcParams["axes.titlesize"] = 24
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18
plt.rcParams["xtick.major.size"] = 12
plt.rcParams["ytick.major.size"] = 12
plt.rcParams["xtick.minor.size"] = 8
plt.rcParams["ytick.minor.size"] = 8


os.makedirs(f"../Plots/", exist_ok=True)
os.makedirs(f"../Plots/Dataset_Errs", exist_ok=True)

# with open(f"../Test_Results/2018_GAP.json", "r") as f:
#     ref_data = json.load(f)

for plot, mth in [["method", "KMED"], ["desc", "KMED"], ["desc", "FPS"], ["desc", "BLR"]]:

    data = {}

    if plot == "desc":
        models = descriptor_comparison_plots
        kwargs = desc_comp_kwargs
    else:
        models = method_comparison_plots
        kwargs = method_comp_kwargs

    for m, mod in enumerate(models):
        print(mod)

        for pot in ["ACE"]:

            if not os.path.exists(f"../Test_Results/{mod}/{mod}_Dataset_Errs_{pot}.json"):
                # Skip as data not available
                print(mod, pot, "not found")
                continue

            with open(f"../Test_Results/{mod}/{mod}_Dataset_Errs_{pot}.json", "r") as f:
                model_data = json.load(f)

            for key in model_data.keys():
                if key not in data.keys():
                    data[key] = {}

                data[key][mod + "_" + pot] = model_data[key]

    if models == descriptor_comparison_plots:
        desc_comp = True
    else:
        desc_comp = False

    Ns = list(data.keys())
    N_ints = [int(n) for n in Ns]
    config_types = list(data[Ns[0]][list(data[Ns[0]].keys())[0]].keys())

    for config_type in config_types:
        if config_type == "IsolatedAtom":
            continue
        plt.clf()

        n_plots = 2
        fig, ax = plt.subplots(ncols=n_plots, figsize=(8 * n_plots, 8), sharex=False)

        if len(ax) > 2:
            ax[-1].axis("off")

        for m, mod in enumerate(models):
            if mth not in mod and "MONTECARLO" not in mod and plot == "desc":
                continue
            if not os.path.exists(f"../Test_Results/{mod}/{mod}_Dataset_Errs_{pot}.json"):
                # Skip as data not available
                continue
            for pot in ["ACE"]:
                model = mod + "_" + pot

                E_rmses = np.array([data[n][model][config_type]["E_rmse_mean"] for n in Ns if (model in data[n].keys())]) * 1000
                F_rmses = np.array([data[n][model][config_type]["F_rmse_mean"] for n in Ns if (model in data[n].keys())]) * 1000

                #N_ints = [n for n in Ns if (model in data[n].keys())]
                
                ax[0].plot(N_ints, E_rmses, label=desc_names(mod), **kwargs[mod], marker=".", linewidth=2, markersize=10.0)

                ax[1].plot(N_ints, F_rmses, label=desc_names(mod), **kwargs[mod], marker=".", linewidth=2, markersize=10.0)

                if "MONTECARLO" in mod or "CUR" in mod:
                    E_stds = np.array([np.std(np.array(data[n][model][config_type]["E_rmse_raw_vals"])) for n in Ns])
                    E_stds *= 1000 / np.sqrt(len(data[Ns[0]][model][config_type]["E_rmse_raw_vals"]))

                    F_stds =  np.array([np.std(np.array(data[n][model][config_type]["F_rmse_raw_vals"])) for n in Ns])
                    F_stds *= 1000 / np.sqrt(len(data[Ns[0]][model][config_type]["F_rmse_raw_vals"]))

                    ax[0].errorbar(N_ints, E_rmses, yerr = E_stds, color = kwargs[mod]["color"], alpha=0.8, capsize=8.0, elinewidth=2.5, capthick=2.5)
                    ax[1].errorbar(N_ints, F_rmses, yerr = F_stds, color = kwargs[mod]["color"], alpha=0.8, capsize=8.0, elinewidth=2.5, capthick=2.5)

        # E_rmses = ref_data[config_type]["E_rmse"] * 1000
        # F_rmses = ref_data[config_type]["F_rmse"] * 1000

        # ax[0].axhline(E_rmses, label="2018 GAP Reference", color="k", linestyle="dashed", linewidth=2, markersize=10.0)

        # ax[1].axhline(F_rmses, label="2018 GAP Reference", color="k", linestyle="dashed", linewidth=2, markersize=10.0)

        ax[0].set_xlabel("Number of Surface Structures")
        ax[1].set_xlabel("Number of Surface Structures")
        ax[0].set_ylabel(f"Energy RMSE (meV/Atom)")
        ax[1].set_ylabel(f"Force RMSE (meV/Ang)")

        ax[0].set_title("Energy Errors")
        ax[1].set_title("Force Errors")

        #ax[0].set_yscale("log")
        #ax[1].set_yscale("log")
        handles, labs = ax[0].get_legend_handles_labels()
        idxs = np.argsort(labs)

        ax[-1].legend([handles[idx] for idx in idxs], [labs[idx] for idx in idxs], ncol=1, fontsize=17)

        if "Total" in config_type:
            if ax[0].get_ylim()[1] > 17.5:
                ax[0].set_ylim(None, 17.5)
            #ax[1].set_ylim(125, 210)
            pass
        plt.tight_layout()
        if desc_comp:
            plt.savefig(f"../Plots/Dataset_Errs/{config_type}_Dataset_Errs_Desc_{mth}{plot_extension}")
        else:
            plt.savefig(f"../Plots/Dataset_Errs/{config_type}_Dataset_Errs_Method{plot_extension}")

        plt.close()