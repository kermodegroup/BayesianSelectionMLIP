import ActiveLearnMLIPTests as al_strats
from ase.io import read, write
import numpy as np
import os
from lips_descriptors import lips_soap_descriptor, lips_ace_descriptor

global_params = {
    "N_samples" : 20,
    "N_structs" : [20, 50, 100, 200],
    "iso_atom_config_type" : "isolated_atom"
}

method_params = {
    "MONTECARLO": { # Pure random selection
        "method" : al_strats.draw_monte_carlo_samples,
        "write_gap_config" : True
    },

    "SOAPCURMAX" : { # CUR of SOAP Descriptor, using max score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : lips_soap_descriptor,
            "aggregate" : np.max
        },
        "write_gap_config" : True
    },

    "SOAPCURMEAN" : { # CUR of SOAP Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : lips_soap_descriptor,
            "aggregate" : np.mean
        },
        "write_gap_config" : True
    },

    "ACECURMEAN": {# CUR of ACE Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : lips_ace_descriptor,
            "aggregate" : np.mean
        },
        "write_gap_config" : False
    },
    "ACEAVGCUR" : { # CUR of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_cur_samples,
        "method_args" : {
            "desc" : lips_ace_descriptor
        },
        "write_gap_config" : False
    },

    "ACEAVGKMED" : { # KMedoids of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : lips_ace_descriptor
        },
        "write_gap_config" : False
    }
}

random_seed = 42

methods_to_gen = [
    "MONTECARLO",
    "ACEAVGKMED",
    "ACEAVGCUR",
    "ACECURMEAN"
]

dataset = read("LiPS_Dataset.xyz", index=":")

for method in methods_to_gen:
    print(method)
    os.makedirs(f"AL_Datasets/{method}", exist_ok=True)
    # Re-seed randomness to ensure reproducibility
    np.random.seed(random_seed)

    params = method_params[method]

    func = params["method"]

    func_args = global_params.copy()

    if "method_args" in params.keys():
        func_args.update(params["method_args"])

    samples = func(dataset, **func_args)

    print(method, " generated samples")
    print("Saving...")
    for i, N in enumerate(global_params["N_structs"]):
        smp = samples[i]
        for j in range(len(smp)):
            ds_name = f"LiPS_{method}_N_{N}_Sample_{j}.xyz"

            write(f"AL_Datasets/{method}/{ds_name}", smp[j])