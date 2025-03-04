import numpy as np
from wfl.select.by_descriptor import do_svd
from .utils import extract_isoats, get_dataset_descriptors
from scipy.linalg import svd


def draw_cur_samples(dataset, desc, N_structs, N_samples=1, aggregate=np.average, iso_atom_config_type=None):
    '''
    Draw a sample of N structures from dataset, based on the CUR scores of each individual descriptor vector.

    if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
    sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
    dataset: list of ase Atoms objects
        Dataset for structure selection
    desc: callable
        Descriptor function with signature descriptor(atoms) -> np array of shape (len(atoms), descriptor_length)
    N_structs: int or iterable of ints
        Number of structures to select per sample. Use a list or array to draw multiple different sub-dataset sizes
    N_samples: int
        Number of samples to draw for each N_struct value
    aggregate: callable
        Aggregation function for all scores corresponding to a single structure
        aggregate(struct_scores) -> final_score
        Good examples are np.average and np.max
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

    scores = cur_scores(ds, desc, aggregate)

    samples = []
    
    for N in N_structs:
        all_sample_idxs = [np.random.choice(idxs, N, p=(scores / np.sum(scores))) for i in range(N_samples)]

        samples.append([
            [
                ds[i] for i in sample_idxs
            ] + isos # Append iso ats to each sample
            for sample_idxs in all_sample_idxs
        ])

    return samples

def draw_avg_cur_samples(dataset, desc, N_structs, N_samples=1, iso_atom_config_type=None):
    '''
    Draw a sample of N structures from dataset, based on the CUR scores of the average descriptor per structure.

    if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
    sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
    dataset: list of ase Atoms objects
        Dataset for structure selection
    desc: callable
        Descriptor function with signature descriptor(atoms) -> np array of shape (len(atoms), descriptor_length)
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

    vecs, idxs = get_dataset_descriptors(ds, desc)

    N = len(ds)

    avg_vecs = np.array([np.average(vecs[idxs==i, :], axis=0) for i in range(N)])
    scores = get_cur_scores(avg_vecs)

    scores /= np.sum(scores)

    idxs = np.arange(N)
    samples = []
    
    for N in N_structs:
        all_sample_idxs = [np.random.choice(idxs, size=N, replace=False, p=scores) for i in range(N_samples)]

        samples.append([
            [
                ds[i] for i in sample_idxs
            ] + isos # Append iso ats to each sample
            for sample_idxs in all_sample_idxs
        ])

    return samples

def cur_scores(structures, desc, N, aggregate=np.average):
    '''
    Compute CUR scores, based on the CUR decomposition of the descriptor vectors.

    
    structures: list of ase Atoms objects
        Structures to compute scores of
    desc: callable
        Descriptor function with signature descriptor(atoms) -> np array of shape (len(atoms), descriptor_length)
    aggregate: callable
        Aggregation function for all scores corresponding to a single structure
        aggregate(struct_scores) -> final_score
        Good examples are np.average and np.max

    Returns:
    scores: np.array
        Scores for each structure
    
    '''
    N = len(structures)

    vecs, idxs = get_dataset_descriptors(structures, desc)

    raw_scores = get_cur_scores(vecs)

    final_scores = np.zeros((N))

    for i in range(N):
        if raw_scores[idxs==i].shape[0]:
            final_scores[i] = aggregate(raw_scores[idxs == i])

    # Normalise resulting scores
    final_scores /= np.sum(final_scores)

    return final_scores


def get_cur_scores(mat):
    '''
    Code inspired by WFL from https://github.com/libAtoms/workflow/blob/main/wfl/select/by_descriptor.py#L47
    Allows finer control over how scores are used, whilst following the same protocol as WFL
    '''

    (u, s, vt) = svd(mat.T)

    c_scores = np.sum(vt ** 2, axis=0) / vt.shape[0]

    return c_scores