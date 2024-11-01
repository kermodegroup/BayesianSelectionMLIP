from sklearn_extra.cluster import KMedoids
import numpy as np
from .utils import extract_isoats, get_dataset_descriptors, get_dataset_avg_descriptors


def draw_avg_kmedoid_samples(dataset, desc, N_structs, iso_atom_config_type=None, **kwargs):
    '''
    Draw a sample of N structures from dataset, based on the kmedoid points of the average descriptor per structure.

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

    samples = []
    
    for N in N_structs:

        kmed = KMedoids(n_clusters=N, method="pam", init="k-medoids++").fit(avg_vecs)
        sample_idxs = kmed.medoid_indices_

        samples.append([
            [
                ds[i] for i in sample_idxs
            ] + isos # Append iso ats to each sample
        ])

    return samples
