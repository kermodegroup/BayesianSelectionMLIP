import numpy as np
from .utils import extract_isoats

def draw_monte_carlo_samples(dataset, N_structs, N_samples=1, iso_atom_config_type=None):
    '''
    Draw a sample of N structures from dataset, from a uniform distribution.

    if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
    sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
    dataset: list of ase Atoms objects
        Dataset for structure selection
    N_structs: int or iterable of ints
        Number of structures to select per sample. Use a list or array to draw multiple different sub-dataset sizes
    N_samples: int
        Number of samples to draw for each N_struct value
    iso_atom_config_type: string
        config type (see atoms.info["config_type"]) for isolated atom configs

    Returns:
    samples: list
        Samples of sub datasets, for each sample and N_struct value
        Generated via [[sample(N) for i in range(N_samples)] for N in range N_structs]
    
    '''

    if np.issubdtype(type(N_structs), np.integer):
        # Convert to len 1 list
        N_structs = [N_structs]

    ds, isos = extract_isoats(dataset, iso_atom_config_type)

    idxs = np.arange(len(ds))

    samples = []
    
    for N in N_structs:
        all_sample_idxs = [np.random.choice(idxs, size=N, replace=False) for i in range(N_samples)]

        samples.append([
            [
                ds[i] for i in sample_idxs
            ] + isos # Append iso ats to each sample
            for sample_idxs in all_sample_idxs
        ])

    return samples