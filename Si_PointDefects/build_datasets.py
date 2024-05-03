from ase.io import read, write

full_dataset = read("../Si/Si_Dataset.xyz", index=":")

pd_config_types = ["interstitial"]

pds = []
non_pds = []

for struct in full_dataset:
    if struct.info["config_type"] in pd_config_types:
        pds.append(struct)
    else:
        non_pds.append(struct)

print(len(pds))

write("Si_PD_Structs.xyz", pds)
write("Si_Core_Data.xyz", non_pds)