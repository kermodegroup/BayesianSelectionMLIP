import ActiveLearnMLIPTests as al_strats
from ase.io import read, write
import numpy as np
import os
from si_descriptors import si_soap_descriptor, si_ace_descriptor, mace_mp_descriptor, si_mace_descriptor, si_total_mace_descriptor

dataset = read("Si_Surface_Structs.xyz", index=":")

core_dataset = read("Si_Core_Data.xyz", index=":")

global_params = {
    "N_samples" : 1,
    "N_structs" : [5, 10, 20, 50, 100]
}

method_params = {
    "MONTECARLO": { # Pure random selection
        "method" : al_strats.draw_monte_carlo_samples,
        "write_gap_config" : True,
        "method_args" : {"N_samples" : 20}
    },

    "SOAPCURMAX" : { # CUR of SOAP Descriptor, using max score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "aggregate" : np.max,
            "N_samples" : 20
        },
        "write_gap_config" : True
    },

    "SOAPCURMEAN" : { # CUR of SOAP Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "aggregate" : np.mean,
            "N_samples" : 20
        },
        "write_gap_config" : True
    },
    "SOAPAVGKMED" : { # KMedoids of average SOAP Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_soap_descriptor
        }
    },
    "ACECURMEAN": {# CUR of ACE Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "aggregate" : np.mean,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "ACEAVGCUR" : { # CUR of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_cur_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "ACEAVGKMED" : { # KMedoids of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_ace_descriptor
        }
    },
    "ACEAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MP0CURMEAN": {# CUR of MACE-MP0 Descriptor, using mean score per struct
        "method" : al_strats.draw_cur_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "aggregate" : np.mean,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },
    "MP0AVGCUR" : { # CUR of average MACE-MP0 Descriptor per structure
        "method" : al_strats.draw_avg_cur_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "MP0AVGKMED" : { # KMedoids of average MACE-MP0 Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor
        },
        "write_gap_config" : False
    },

    "MP0AVGFPS" : { # FPS of average MACE-MP0 Descriptor per structure
        "method" : al_strats.draw_avg_fps_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MACEAVGCUR" : { # CUR of average MACE Descriptor per structure
        "method" : al_strats.draw_avg_cur_samples,
        "method_args" : {
            "desc" : si_mace_descriptor,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "MACEAVGKMED" : { # KMedoids of average MACE Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_mace_descriptor
        },
        "write_gap_config" : False
    },
    "MACETAVGKMED" : { # KMedoids of average MACE Descriptor per structure
        "method" : al_strats.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_total_mace_descriptor
        },
        "write_gap_config" : False
    },
    "MACEAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_mace_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MACETAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_total_mace_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "SOAPAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MP0AVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : al_strats.draw_avg_fps_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "ACEBLR" : {
        "method" : al_strats.async_hal,
        "method_args" : {
            "descriptor" : si_ace_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "MACEBLR" : {
        "method" : al_strats.async_hal,
        "method_args" : {
            "descriptor" : si_mace_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "MACETBLR" : {
        "method" : al_strats.async_hal,
        "method_args" : {
            "descriptor" : si_total_mace_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "SOAPBLR" : {
        "method" : al_strats.async_hal,
        "method_args" : {
            "descriptor" : si_soap_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "MP0BLR" : {
        "method" : al_strats.async_hal,
        "method_args" : {
            "descriptor" : mace_mp_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    }
} 

random_seed = 42


from Tests.plot_config import *
methods_to_gen = method_comparison_plots
methods_to_gen = [
    "SOAPBLR",
    "MACETBLR",
]

print(len(methods_to_gen))
# np.random.shuffle(methods_to_gen)
methods_to_gen = [methods_to_gen[0]]

for method in methods_to_gen:
    print(method)

    if os.path.exists(f"AL_Datasets/{method}"):
        pass#continue
    os.makedirs(f"AL_Datasets/{method}", exist_ok=True)
    # Re-seed randomness to ensure reproducibility
    np.random.seed(random_seed)

    params = method_params[method]

    func = params["method"]

    func_args = global_params.copy()

    if "method_args" in params.keys():
        func_args.update(params["method_args"])

    samples = func(dataset, **func_args)

    for i, N in enumerate(func_args["N_structs"]):
        smp = samples[i]
        for j in range(func_args["N_samples"]):
            ds_name = f"Si_{method}_N_{N}_Sample_{j}.xyz"
            gap_name = f"Si_{method}_N_{N}_Sample_{j}.xml" 
            gap_config_name = f"Si_{method}_N_{N}_Sample_{j}.gapconfig" 

            write(f"AL_Datasets/{method}/{ds_name}", smp[j] + core_dataset)

            with open("Models/GAPs/base_gap_config", "r") as base_gap_config_file:
                    with open(f"AL_Datasets/{method}/{gap_config_name}", "w") as new_config:
                        new_config.writelines(base_gap_config_file) # create a copy of base config
                        new_config.writelines(
                            [
                                f"gp_file={gap_name}\n",
                                f"at_file={ds_name}"
                            ]
                        )
