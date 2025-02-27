from ase.io import read, write
import json
import numpy as np

###
# Generate the core and sampling datasets used in the paper

full_dataset = read("Si_2018_Dataset.xyz", index=":")


# Config types to take from the original dataset

core_config_types = ["isolated_atom", "dia"]
surf_config_types = ["surface_001", "surface_110", "surface_111", "surface_111_pandey", "surface_111_3x3_das", "decohesion"]

surf = []
core = []

for struct in full_dataset:
    if struct.info["config_type"] in surf_config_types:
        surf.append(struct)
    elif struct.info["config_type"] in core_config_types:
        core.append(struct)

all_structs = [struct for struct in core]
all_structs.extend(surf)


config_types = {}

# Tabulate the number of atoms and counts for each config type
for structure in all_structs:
    ct = structure.info["config_type"]
    N = len(structure)

    if ct not in config_types.keys():
        config_types[ct] = [len(structure)]
    else:
        config_types[ct].append(len(structure))

for ct in config_types.keys():
    print(ct, np.unique(config_types[ct], return_counts=True))


write("Si_Total_Dataset.xyz", all_structs)
write("Si_Surface_Structs.xyz", surf)
write("Si_Core_Data.xyz", core)