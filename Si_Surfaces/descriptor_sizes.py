from si_descriptors import si_ace_descriptor, si_ace_small_descriptor, si_ace_large_descriptor, si_mace_descriptor, si_total_mace_descriptor, mace_mp_descriptor
from ase.io import read

ats = read("Si_Core_Data.xyz", index="0")

descs = [si_ace_descriptor, si_ace_small_descriptor, si_ace_large_descriptor, si_mace_descriptor, si_total_mace_descriptor, mace_mp_descriptor]
labels = ["ACE (M)", "ACE (S)", "ACE (L)", "MACE (C)", "MACE (T)", "MP0"]

for i, desc in enumerate(descs): 
    d = desc(ats)

    print(labels[i], f"Len = {d.shape[-1]}")