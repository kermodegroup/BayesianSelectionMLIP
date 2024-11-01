from ase.io import read, write
import json

full_dataset = read("../Si/Si_Dataset.xyz", index=":")

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
        config_types[ct] = {"Nats" : N, "Nstruct" : 1}
    else:
        config_types[ct]["Nats"] += N
        config_types[ct]["Nstruct"] += 1

for ct in config_types.keys():
    config_types[ct]["Nats"] /= config_types[ct]["Nstruct"]

with open("Total_Dataset_info.json", "w") as f:
    json.dump(config_types, f, indent=4) 

# write("Total_Dataset.xyz", all_structs)


# exit()
# write("Si_Surface_Structs.xyz", pds)
# write("Si_Core_Data.xyz", non_pds)