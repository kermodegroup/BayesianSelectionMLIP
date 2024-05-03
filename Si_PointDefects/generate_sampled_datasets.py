import ActiveLearnMLIPTests as al_strats
from ase.io import read, write
import numpy as np
import os
from si_descriptors import si_soap_descriptor, si_ace_descriptor

global_params = {
    "N_samples" : 20,
    "N_structs" : [5, 10, 20, 50, 100]
}

method_params = {
    "MONTECARLO": { # Pure random selection
        "method" : al_strats.draw_monte_carlo_samples,
        "write_gap_config" : True
    },

    "SOAPCURMAX" : { # CUR of SOAP Descriptor, using max score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "aggregate" : np.max
        },
        "write_gap_config" : True
    },

    "SOAPCURMEAN" : { # CUR of SOAP Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "aggregate" : np.mean
        },
        "write_gap_config" : True
    },

    "ACECURMEAN": {# CUR of ACE Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "aggregate" : np.mean
        },
        "write_gap_config" : False
    },

    "ACEAVGCUR" : { # CUR of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_cur_samples,
        "method_args" : {
            "desc" : si_ace_descriptor
        },
        "write_gap_config" : False
    },

    "ACEAVGKMED" : { # KMedoids of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_ace_descriptor
        }
    }
}

random_seed = 42

methods_to_gen = [
    "ACEAVGCUR",
    "ACEAVGKMED",
    "ACECURMEAN",
    "MONTECARLO"
]

dataset = read("Si_PD_Structs.xyz", index=":")

core_dataset = read("Si_Core_Data.xyz", index=":")

for method in methods_to_gen:
    os.makedirs(f"AL_Datasets/{method}", exist_ok=True)
    # Re-seed randomness to ensure reproducibility
    np.random.seed(random_seed)

    params = method_params[method]

    func = params["method"]

    func_args = global_params.copy()

    if "method_args" in params.keys():
        func_args.update(params["method_args"])

    samples = func(dataset, **func_args)

    for i, N in enumerate(global_params["N_structs"]):
        smp = samples[i]
        for j in range(len(smp)):
            ds_name = f"Si_{method}_N_{N}_Sample_{j}.xyz"
            gap_name = f"Si_{method}_N_{N}_Sample_{j}.xml" 
            gap_config_name = f"Si_{method}_N_{N}_Sample_{j}.gapconfig" 

            write(f"AL_Datasets/{method}/{ds_name}", smp[j] + core_dataset)
