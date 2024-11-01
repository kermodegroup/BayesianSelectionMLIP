from ase.io import read
import numpy as np

ds = read("../Si_Surface_Structs.xyz", index=":")

data = {}


for struct in ds:
    if struct.info["config_type"] in data.keys():
        data[struct.info["config_type"]].append(len(struct))
    else:
        data[struct.info["config_type"]] = [len(struct)]


for key in data.keys():
    data[key] = np.unique(data[key])


print(data)