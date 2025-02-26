from ActiveLearnMLIPTests import soap_descriptor, ace_descriptor, mace_mp_descriptor, mace_descriptor
from pathlib import Path


si_ace_descriptor = ace_descriptor(
    species = ["Si"],
    order = 3,
    totaldegree = 16,
    rcut = 5.0)

si_ace_small_descriptor = ace_descriptor(
    species = ["Si"],
    order = 3,
    totaldegree = 12,
    rcut = 5.0)

si_ace_large_descriptor = ace_descriptor(
    species = ["Si"],
    order = 3,
    totaldegree = 19,
    rcut = 5.0)

si_soap_descriptor = soap_descriptor(" \
    soap \
    l_max=12 \
    n_max=10 \
    atom_sigma=0.5 \
    zeta=4 \
    cutoff=5.0 \
    cutoff_transition_width=1.0 \
    central_weight=1.0 \
    delta=3.0 \
    f0=0.0 \
    covariance_type=dot_product \
")
    
mace_mpa_descriptor = mace_mp_descriptor(model="medium-mpa-0")
mace_mp_descriptor = mace_mp_descriptor(model="medium")

mace_fpath = str(Path(__file__).absolute().parent) + "/Models/BaseMACE/Si_Core_MACE_swa_compiled.model"
si_mace_descriptor = mace_descriptor(mace_fpath)

mace_fpath = str(Path(__file__).absolute().parent) + "/Models/TotalMACE/Si_Total_MACE_swa_compiled.model"
si_total_mace_descriptor = mace_descriptor(mace_fpath)