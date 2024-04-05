from ase.io import read, write

dataset = read("gp_iter6_sparse9k.xml.xyz", index=":")

config_types = {}

ds_partitions = {}

for image in dataset:
    ct = image.info["config_type"]

    if ct not in config_types.keys():
        config_types[ct] = 1
        ds_partitions[ct] = [image]
    else:
        config_types[ct] += 1
        ds_partitions[ct].append(image)


print("ORIGINAL DATASET")
print(len(dataset))

for key, val in config_types.items():
    print(key, val)

    # write("Original_Dataset/" + key + ".xyz", ds_partitions[key])

print("----------------")
print("REDUCED DATASET")
ds_cts = ["isolated_atom", "interstitial", "vacancy", "surface_001", "surface_110", "surface_111", "surface_111_pandey", "surface_111_3x3_das", "dia"]

ds = []

for ct in ds_cts:
    ds.extend(ds_partitions[ct])

print(len(ds))
write("Si_Dataset.xyz", ds)