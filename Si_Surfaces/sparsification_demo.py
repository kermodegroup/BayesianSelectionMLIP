from ActiveLearnMLIPTests import ace_descriptor, draw_avg_kmedoid_samples, mace_mp_descriptor, bayesian_selection
from ase.io import read
dataset = read("Si_Core_Data.xyz", ":")

# "Core" dataset, treated as common for all models
core_ds = dataset[:10]

# "Sampling" dataset, to sparsify
samp_ds = dataset[11:20]

# Number of structures to select
N_select = 4

# Define an ACE descriptor for silicon
ace_desc = ace_descriptor(
    species = ["Si"],
    order = 3,
    totaldegree = 16,
    rcut = 5.0)


# Sample using k-medoids. No information about the core dataset can be passed to this method
kmed_samps = draw_avg_kmedoid_samples(samp_ds, ace_desc, [N_select], iso_atom_config_type="isolated_atom" )


# Sample using Bayesian Selection, using the core dataset as a starting point
bayes_samps = bayesian_selection(samp_ds, ace_desc, [N_select], core_ds, iso_atom_config_type="isolated_atom")