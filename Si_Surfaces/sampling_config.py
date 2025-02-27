import BayesianSelection as select
from ase.io import read, write
import numpy as np
import os
from si_descriptors import *

###
# Configuration of the sampling methods used to generate datasets


# Dataset to sample from
dataset = read("Si_Surface_Structs.xyz", index=":")

# Core dataset given to each sampled dataset
core_dataset = read("Si_Core_Data.xyz", index=":")

global_params = {
    "N_samples" : 1,
    "N_structs" : [5, 10, 20, 50, 100]
}

method_params = {
    "MONTECARLO": { # Pure random selection
        "method" : select.draw_monte_carlo_samples,
        "write_gap_config" : True,
        "method_args" : {"N_samples" : 20}
    },

    "SOAPCURMAX" : { # CUR of SOAP Descriptor, using max score per struct
        "method" : select.draw_cur_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "aggregate" : np.max,
            "N_samples" : 20
        },
        "write_gap_config" : True
    },

    "SOAPCURMEAN" : { # CUR of SOAP Descriptor, using mean score per struct
        "method" : select.draw_cur_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "aggregate" : np.mean,
            "N_samples" : 20
        },
        "write_gap_config" : True
    },
    "SOAPAVGKMED" : { # KMedoids of average SOAP Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_soap_descriptor
        }
    },
    "ACECURMEAN": {# CUR of ACE Descriptor, using mean score per struct
        "method" : select.draw_cur_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "aggregate" : np.mean,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "ACEAVGCUR" : { # CUR of average ACE Descriptor per structure
        "method" : select.draw_avg_cur_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "ACEAVGKMED" : { # KMedoids of average ACE Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_ace_descriptor
        }
    },
    "ACEAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_ace_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MP0CURMEAN": {# CUR of MACE-MP0 Descriptor, using mean score per struct
        "method" : select.draw_cur_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "aggregate" : np.mean,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },
    "MP0AVGCUR" : { # CUR of average MACE-MP0 Descriptor per structure
        "method" : select.draw_avg_cur_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "MP0AVGKMED" : { # KMedoids of average MACE-MP0 Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor
        },
        "write_gap_config" : False
    },

    "MP0AVGFPS" : { # FPS of average MACE-MP0 Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MACEAVGCUR" : { # CUR of average MACE Descriptor per structure
        "method" : select.draw_avg_cur_samples,
        "method_args" : {
            "desc" : si_mace_descriptor,
            "N_samples" : 20
        },
        "write_gap_config" : False
    },

    "MACEAVGKMED" : { # KMedoids of average MACE Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_mace_descriptor
        },
        "write_gap_config" : False
    },
    "MACETAVGKMED" : { # KMedoids of average MACE Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_total_mace_descriptor
        },
        "write_gap_config" : False
    },
    "MACEAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_mace_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MACETAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_total_mace_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "SOAPAVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : si_soap_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MP0AVGFPS" : { # FPS of average ACE Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : mace_mp_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "ACEBLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : si_ace_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "MACEBLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : si_mace_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "MACETBLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : si_total_mace_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "SOAPBLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : si_soap_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "MP0BLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : mace_mp_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "ACESBLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : si_ace_small_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "ACELBLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : si_ace_large_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
    "ACESAVGKMED" : { # KMedoids of average MACE Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_ace_small_descriptor
        },
        "write_gap_config" : False
    },
    "ACELAVGKMED" : { # KMedoids of average MACE Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : si_ace_large_descriptor
        },
        "write_gap_config" : False
    },
    "MPAAVGFPS" : { # FPS of average MACE MPA Descriptor per structure
        "method" : select.draw_avg_fps_samples,
        "method_args" : {
            "desc" : mace_mpa_descriptor,
            "core_dataset" : core_dataset
        },
        "write_gap_config" : False
    },
    "MPAAVGKMED" : { # KMedoids of average MACE MPA Descriptor per structure
        "method" : select.draw_avg_kmedoid_samples,
        "method_args" : {
            "desc" : mace_mpa_descriptor
        },
        "write_gap_config" : False
    },
    "MPABLR" : {
        "method" : select.bayesian_selection,
        "method_args" : {
            "descriptor" : mace_mpa_descriptor,
            "core_ds" : core_dataset,
            "iso_atom_config_type" : "isolated_atom"
        },
    },
} 