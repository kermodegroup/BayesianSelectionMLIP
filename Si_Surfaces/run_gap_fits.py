import mpi4py
import os
import shutil
import numpy as np

from Tests.plot_config import *

methods = probabalistic_models

np.random.shuffle(methods)

for method in methods:
    os.makedirs(f"Models/GAPs/{method}", exist_ok=True)

    os.chdir(f"Models/GAPs/{method}")

    data_src = f"../../../AL_Datasets/{method}"

    if not os.path.exists(data_src):
        continue

    xyz_files = [file for file in os.listdir(data_src) if ".xyz" in file]

    np.random.shuffle(xyz_files)

    for i, file in enumerate(xyz_files):
        name = file[:-4]
        print(i+1, "/", len(xyz_files), name)
        os.makedirs(name, exist_ok=True)
        os.chdir(name)

        if name + ".xml" not in os.listdir():

            # Copy dataset
            shutil.copyfile("../" + data_src + os.sep + file, file)

            # Copy gap config
            shutil.copyfile("../" + data_src + os.sep + name + ".gapconfig", name + ".gapconfig")

            # Copy IP Glue core potential
            shutil.copyfile("../../../OriginalGAP/repulsive_2b.xml", "repulsive_2b.xml")

            # Run gap_fit to fit potential
            os.system(
                f"gap_fit config_file={name}.gapconfig"
            )

        os.chdir("..")
    os.chdir(f"../../../")
