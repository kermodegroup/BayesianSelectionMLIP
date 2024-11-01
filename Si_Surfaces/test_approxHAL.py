from ActiveLearnMLIPTests import async_hal
from ase.io import read
from si_descriptors import si_ace_descriptor


ats = read("Si_Core_Data.xyz", ":")

core_ds = ats[:10]
samp_ds = ats[11:20]

print("Starting")

samps = async_hal(samp_ds, si_ace_descriptor, [2, 5], core_ds, iso_atom_config_type="isolated_atom")
