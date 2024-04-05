import ActiveLearnMLIPTests as al_strats
from ase.io import read, write
import numpy as np
import os
from si_desciptors import si_soap_descriptor

global_params = {
    "N_samples" : 10,
    "N_structs" : [88, 121, 242, 484, 726],
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
    }
}

random_seed = 42

methods_to_gen = [
 #"MONTECARLO", 
 #"SOAPCURMAX", 
 "SOAPCURMEAN"
]

base_gap_config_file = open("Models/GAPs/base_gap_config", "r")



dataset = read("Si_Dataset.xyz", index=":")

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
        for j in range(global_params["N_samples"]):
            ds_name = f"Si_{method}_N_{N}_Sample_{j}.xyz"
            gap_name = f"Si_{method}_N_{N}_Sample_{j}.xml" 
            gap_config_name = f"Si_{method}_N_{N}_Sample_{j}.gapconfig" 

            write(f"AL_Datasets/{method}/{ds_name}", samples[i][j])

            if params["write_gap_config"]:
                with open("Models/GAPs/base_gap_config", "r") as base_gap_config_file:
                    with open(f"AL_Datasets/{method}/{gap_config_name}", "w") as new_config:
                        new_config.writelines(base_gap_config_file) # create a copy of base config
                        new_config.writelines(
                            [
                                f"gp_file={gap_name}\n",
                                f"at_file={ds_name}"
                            ]
                        )