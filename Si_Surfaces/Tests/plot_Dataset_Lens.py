import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import json
import os

colours = list(mcolors.TABLEAU_COLORS.values())
markers = [".", "x", "o", "s", "*", "+", "D", "^", "<", ">", "v", "H", "d"]

plot_colours = colours

scatter_colours_markers = []

for marker in markers:
  for colour in colours:
    scatter_colours_markers.append((colour, marker))

Ns = [5, 10, 20, 50, 100]

results_dir = "../Test_Results"

methods = [method for method in os.listdir(results_dir) 
           if os.path.exists(results_dir + os.sep + method + os.sep + "Dataset_Lens.json")]


for N in Ns:
    for i, method in enumerate(methods):
        colour, marker = scatter_colours_markers[i]
        with open(results_dir + os.sep + method + os.sep + "Dataset_Lens.json", "r") as f:
            data = json.load(f)[str(N)]
        
        keys = sorted(data.keys())
        
        plt.scatter(np.arange(len(keys)), [data[key] for key in keys], color=colour, marker=marker, label=method)

    plt.xticks(np.arange(len(keys)), keys, rotation=30)
    plt.ylabel("Number of Structures")
    plt.legend()
    plt.savefig(f"../Plots/Dataset_Lens/{N}.pdf")
