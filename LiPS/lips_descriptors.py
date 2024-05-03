from ActiveLearnMLIPTests import soap_descriptor, ace_descriptor

lips_soap_descriptor = soap_descriptor(" \
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

lips_ace_descriptor = ace_descriptor(
    species = ["Li", "P", "S"],
    order = 3,
    totaldegree = 8,
    rcut = 5.0)