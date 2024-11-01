import numpy as np
import os
from ase.io import read
from  plot_config import *
import matplotlib.pyplot as plt

data = {}


plot = "method"


if plot == "desc":
    models = descriptor_comparison_plots
    kwargs = desc_comp_kwargs
    labels = desc_comp_labels
else:
    models = method_comparison_plots
    kwargs = method_comp_kwargs
    labels = method_comp_labels

base_dataset = read("../Si_Core_Data.xyz", index=":")
N_base = len(base_dataset)
base_nats = sum([len(ats) for ats in base_dataset])

add_dataset = read("../Si_Surface_Structs.xyz", index=":")
struct_lens = [len(ats) for ats in add_dataset]

for mod in models:
    print(mod)
    desc, method = split_mod(mod)
    if not os.path.exists(f"../AL_Datasets/{mod}"):
        # Skip as data not available
        print(mod, "not found")
        continue
    
    data[mod] = {}

    xyzs = [file for file in os.listdir(f"../AL_Datasets/{mod}") if ".xyz" in file]

    for xyz in xyzs:
        dataset = read(f"../AL_Datasets/{mod}/" + xyz, index=":")

        i_samp = xyz.split("_")[-1].split(".")[0]
        N_structs = len(dataset) - N_base
        nats = sum([len(ats) for ats in dataset]) - base_nats

        print(xyz, len(dataset), sum([len(ats) for ats in dataset]))

        avg_nats = nats / N_structs
        
        if N_structs not in data[mod].keys():
            data[mod][N_structs] = {}

        if not len(data[mod][N_structs].keys()):
            # Empty subdict, populate fields
            data[mod][N_structs] = {
                "avg_nats" : avg_nats,
                "N" : 1
            }
        else:
            # Sum properties so they can be averaged later
            data[mod][N_structs]["avg_nats"] += avg_nats
            data[mod][N_structs]["N"] += 1

    sorted_keys = sorted(list(data[mod].keys()))
    plt.plot(sorted_keys, [data[mod][key]["avg_nats"] / data[mod][key]["N"] for key in sorted_keys], 
                label=labels[mod], **kwargs[mod])

plt.xlabel("N_struct")
plt.ylabel("Average number of atoms per structure in dataset")

lin1, lab1 = plt.gca().get_legend_handles_labels()
idxs = np.argsort(lab1)
lin1 = [lin1[idx] for idx in idxs]
lab1 = [lab1[idx] for idx in idxs]

ax2 = plt.gca().twiny()

bins = [49, 55, 93, 99, 105, 111, 141, 147]
ax2.hist(struct_lens, density=False, orientation="horizontal", color="C4", alpha=0.2, bins=bins, label="Dataset distribution")

ax2.set_xticks([], [])

lin2, lab2 = ax2.get_legend_handles_labels()

lines = lin1 + lin2
labels = lab1 + lab2

plt.legend(lines, labels, ncol=2, title="<Desc>_<Method>", loc="lower right")



ticks=[52, 96, 108, 144]
labels=["-(111) surface 3x3 das", "-(111) surface", "-(110) surface", "-(001) surface"]

for i in range(len(ticks)):
    ax2.text(0, ticks[i], labels[i], ha="left", va="center", color="C4", alpha=0.8)

plt.tight_layout()
plt.savefig("../Plots/DatasetLens.png", dpi=200)