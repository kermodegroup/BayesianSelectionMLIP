import numpy as np
from .utils import extract_isoats, get_dataset_descriptors, get_dataset_avg_descriptors
from wfl.select.by_descriptor import greedy_fps_conf_global
from wfl.configset import OutputSpec, ConfigSet
from numpy.random import default_rng

def draw_avg_fps_samples(dataset, desc, N_structs, iso_atom_config_type=None, core_dataset=None, **kwargs):
    '''
    Draw a sample of N structures from dataset, based on the farthest point sampling of the average descriptor per structure.

    if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
    sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
    dataset: list of ase Atoms objects
        Dataset for structure selection
    desc: callable
        Descriptor function with signature descriptor(atoms) -> np array of shape (len(atoms), descriptor_length)
    N_structs: int or iterable of ints
        Number of structures to select per sample. Use a list or array to draw multiple different sub-dataset sizes
    iso_atom_config_type: string
        config type (see atoms.info["config_type"]) for isolated atom configs
    core_dataset: list of ase Atoms objects
        Core/previously selected structures. Used to initialise the FPS sampling

    Returns:
    samples: list
        Samples of sub datasets, for each sample and N_struct value
        Generated via [[sample(N) for i in range(N_samples)] for N in range N_structs]
    
    '''

    if np.issubdtype(type(N_structs), np.integer):
        # Convert to len 1 list
        N_structs = [N_structs]

    ds, isos = extract_isoats(dataset, iso_atom_config_type)

    avg_vecs = get_dataset_avg_descriptors(ds, desc)

    if core_dataset is not None:
        core_vecs = get_dataset_avg_descriptors(core_dataset, desc)
    else:
        core_vecs = None

    samples = []

    avg_vecs /= np.max(np.abs(avg_vecs))
    
    for N in N_structs:
        rng = default_rng(seed=42)
        samples.append([
            [item for item in greedy_fps_conf_global(ConfigSet([ats.copy() for ats in ds]), OutputSpec(), N, avg_vecs.copy(), rng=rng, O_N_sq=True, prev_selected_descs=core_vecs.copy())] + isos # Append iso ats to each sample
        ])

    return samples
