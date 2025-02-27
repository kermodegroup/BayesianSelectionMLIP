import os
import json
import matplotlib.pyplot as plt
import numpy as np
from plot_config import *


###
# Plots Surface Formation Energy figures (Figs. 2 & 4)


plt.rcParams["axes.labelsize"] = 20
plt.rcParams["axes.titlesize"] = 24
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18
plt.rcParams["xtick.major.size"] = 12
plt.rcParams["ytick.major.size"] = 12
plt.rcParams["xtick.minor.size"] = 8
plt.rcParams["ytick.minor.size"] = 8


plot = "method"
plot = "desc"

mth = "KMED"
mth = "FPS"
mth = "BLR"

for plot, mth in [["method", "KMED"], ["desc", "KMED"], ["desc", "FPS"], ["desc", "BLR"]]:


    if plot == "desc":
        models = descriptor_comparison_plots
        kwargs = desc_comp_kwargs
    else:
        models = method_comparison_plots
        kwargs = method_comp_kwargs

    # DFT vals for E100, E110, E111 (see )
    dft_values = [
        2.17, 1.52, 1.57
    ]

    gap_errors = np.array([2, 1, 2]) * np.array(dft_values) / 1000
    print(gap_errors)

    props = ["E100", "E110", "E111"]

    plot_props = [0, 2]

    fig, ax = plt.subplots(1, len(plot_props), figsize=(8 * len(plot_props), 8))

    ax = ax.flatten()

    for i in range(len(plot_props)):
        prop = props[plot_props[i]]
        ax[i].set_title("(" + prop[1:] + ") Surface Formation Energy")
        ax[i].set_xlabel("Number of Surface Structures")
        ax[i].set_ylabel(f"Error from DFT (J/m$^2$)")
        ax[i].set_yscale("log")
        #ax[i].set_ylim(0, 0.4)


    #ax[-1].axis("off")

    if models == descriptor_comparison_plots:
        desc_comp = True
    else:
        desc_comp = False


    for m, mod in enumerate(models):
        if mth not in mod and "MONTECARLO" not in mod and plot == "desc":
            continue
        print(mod)
        #for pot in ["GAP", "ACE"]:
        for pot in ["ACE"]:
            if not os.path.exists(f"../Test_Results/{mod}/{mod}_Surfaces_{pot}.json"):
                # Skip as data not available
                print("Not Found")
                continue

            with open(f"../Test_Results/{mod}/{mod}_Surfaces_{pot}.json", "r") as f:
                data = json.load(f)

            Ns = [int(k) for k in data.keys()]

            for i in range(len(plot_props)):
                prop = props[plot_props[i]]
                errs = np.array([np.sqrt(np.mean((np.array(ndata[prop + "_raw_vals"]) - dft_values[i])**2)) for ndata in data.values()])

                max_err = [np.max(np.abs(np.array(ndata[prop + "_raw_vals"]) - dft_values[i]))for ndata in data.values()]
                min_err = [np.min(np.abs(np.array(ndata[prop + "_raw_vals"]) - dft_values[i]))for ndata in data.values()]
                ax[i].plot(Ns, errs, label=desc_names(mod), **kwargs[mod], marker=".", linewidth=2, markersize=10.0)

                if "MONTECARLO" in mod or "CUR" in mod:
                    stds =  np.array([np.std(np.array(ndata[prop + "_raw_vals"]))for ndata in data.values()])
                    stds /= np.sqrt(len(list(data.values())[0][prop + "_raw_vals"]))
                    #ax[i].plot(Ns, errs + 2*stds, color = kwargs[mod]["color"], linestyle = "dotted")
                    #ax[i].plot(Ns, errs - 2*stds, color = kwargs[mod]["color"], linestyle = "dotted")
                    if "MONTECARLO" in mod:
                        label = f"Error across {len(list(data.values())[0][prop + '_raw_vals'])} samples"
                    else:
                        label=None
                    label=None
                    ax[i].errorbar(Ns, errs,yerr = stds, color = kwargs[mod]["color"], alpha=0.8, capsize=8.0, elinewidth=2.5, capthick=2.5, label=label)


    for i in range(len(plot_props)):
        errs = gap_errors[plot_props[i]]
        #ax[i].axhline(errs, color="k", linestyle="dashed", label="2018 Si GAP Benchmark")
        
    handles, labels = ax[0].get_legend_handles_labels()
    idxs = np.argsort(labels)

    ax[0].legend([handles[idx] for idx in idxs], [labels[idx] for idx in idxs], fontsize=17)
    plt.tight_layout()

    if desc_comp:
        plt.savefig(f"../Plots/Surfaces_Desc_{mth}.eps")
    else:
        plt.savefig(f"../Plots/Surfaces_Method.eps")