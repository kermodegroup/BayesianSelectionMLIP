import os
import shutil
import numpy as np

model = "ACE"

print(f"{model}:")

from Tests.plot_config import *

methods = np.unique(method_comparison_plots + descriptor_comparison_plots)

for method in methods:
    os.makedirs(f"Models/{model}s/{method}", exist_ok=True)

    os.chdir(f"Models/{model}s/{method}")

    data_src = f"../../../AL_Datasets/{method}"

    if not os.path.exists(data_src):
        os.chdir(f"../../../")
        print(method, "0 / 0 trained")
        print()
        continue

    xyz_files = [file for file in os.listdir(data_src) if ".xyz" in file]

    n_trained = 0

    for i, file in enumerate(xyz_files):
        name = file[:-4]
        os.makedirs(name, exist_ok=True)
        os.chdir(name)

        is_trained = name + ".xml" in os.listdir() or name + ".json" in os.listdir() 

        n_trained += is_trained

        #print(i+1, "/", len(xyz_files), name, is_trained)

        os.chdir("..")
    os.chdir(f"../../../")

    print(method, n_trained, "/", len(xyz_files), "trained")
    print()