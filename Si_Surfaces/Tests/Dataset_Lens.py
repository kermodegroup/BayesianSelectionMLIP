import numpy as np
import os
from ase.io import read
from  plot_config import *
import json

Ns = [5, 10, 20, 50, 100]

surf_configs = ["surface_001", "surface_110", "surface_111", "surface_111_pandey", "surface_111_3x3_das", "decohesion"]

data_dir = "../AL_Datasets"

methods = os.listdir(data_dir)


for method in methods:
    print(method)
    ct_counts = {}
    for N in Ns:
        dataset = read(data_dir + os.sep + method + os.sep + f"Si_{method}_N_{N}_Sample_0.xyz", index=":")

        ct_counts[N] = {ct : 0 for ct in surf_configs}

        for image in dataset:
            ct = image.info["config_type"]

            if ct in surf_configs:
                ct_counts[N][ct] += 1
    
    with open(f"../Test_Results/{method}/Dataset_Lens.json", "w") as f:
        json.dump(ct_counts, f, indent=4)

    
