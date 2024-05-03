from ase.io import read, write

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

print(len(pds))

write("Si_Surface_Structs.xyz", pds)
write("Si_Core_Data.xyz", non_pds)