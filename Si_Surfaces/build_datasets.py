from ase.io import read, write
import json
import numpy as np

full_dataset = read("Si_2018_Dataset.xyz", index=":")

pd_config_types = [""]

excl_config_types = ["interstitial", "vacancy"]

pds = []
non_pds = []

for struct in full_dataset:
    if "surface" in struct.info["config_type"]:
        pds.append(struct)
    elif struct.info["config_type"] not in excl_config_types:
        non_pds.append(struct)

decohesion = read("decohesion.xyz", ":")

pds.extend(decohesion)
print(len(non_pds))
print(len(pds))

all_structs = [struct for struct in non_pds]
all_structs.extend(pds)


config_types = {}

for structure in all_structs:
    ct = structure.info["config_type"]
    N = len(structure)

    if ct not in config_types.keys():
        config_types[ct] = [len(structure)]
    else:
        config_types[ct].append(len(structure))

for ct in config_types.keys():
    print(ct, np.unique(config_types[ct], return_counts=True))
# write("Total_Dataset.xyz", all_structs)


# exit()
# write("Si_Surface_Structs.xyz", pds)
# write("Si_Core_Data.xyz", non_pds)