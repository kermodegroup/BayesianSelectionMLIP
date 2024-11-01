from .utils import extract_isoats
import numpy as np
from tqdm import trange
import scipy

def bayesian_selection(dataset, descriptor, N_structs, core_ds=None, N_samples=1, iso_atom_config_type=None):
    '''
    Draw a sample of N structures from dataset, based on the bayesian selection scores.

    if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
    sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
    dataset: list of ase Atoms objects
        Dataset for structure selection
    N_structs: int or iterable of ints
        Number of structures to select per sample. Use a list or array to draw multiple different sub-dataset sizes
    N_samples: int
        Number of samples to draw for each N_struct value
    property: string
        Property to perform committee error selection on ("energy", "force)
    iso_atom_config_type: string
        config type (see atoms.info["config_type"]) for isolated atom configs

    Returns:
    samples: list
        Samples of sub datasets, for each sample and N_struct value
        Generated via [[sample(N) for i in range(N_samples)] for N in range N_structs]
    
    '''

    def enforce_2d(d):
        if len(d.shape) == 1:
            return d[np.newaxis, :]
        else:
            return d

    def energy_var(Psi_E, PLU):
        # Variance of energy per atom via posterior covariance

        P, L, U = PLU
        try:
            left = scipy.linalg.solve_triangular(U.T, Psi_E.T, lower=True).T
            right = scipy.linalg.solve_triangular(L, P.T @ Psi_E.T, lower=True)

            E_cov = left @ right
            C = np.ones(Psi_E.shape[0])

            return C @ E_cov @ C.T / Psi_E.shape[0]**2
        except np.linalg.LinAlgError:
            # Linalg error when L or U is singular, "skip" structure
            return -np.inf

    def max_force_comp_var(Psi_F, PLU):
        # Maximum variance of a force component
        P, L, U = PLU

        left = [scipy.linalg.solve_triangular(U.T, Psi_F[i, :[]].T, lower=True).T for i in range(Psi_F.shape[0])]
        right = [scipy.linalg.solve_triangular(L, P.T @ Psi_F[i, :].T, lower=True) for i in range(Psi_F.shape[0])]

        F_vars = [left[i] @ right[i] for i in range(Psi_F.shape[0])]

        return np.max(F_vars)

        
    if np.issubdtype(type(N_structs), np.integer):
        # Convert to len 1 list
        N_structs = [N_structs]

    samp_ds, isos = extract_isoats(dataset, iso_atom_config_type)
    
    score_func = energy_var
    prop = "E"
    
    if len(isos) == 0 and core_ds is not None:
        # No iso_ats found in sample dataset
        core_ds, isos = extract_isoats(core_ds, iso_atom_config_type)

    if core_ds is None:
        core_ds = []
    core_ds = core_ds + isos


    desc_len = descriptor(core_ds[0]).shape[-1]

    post_cov = np.zeros((desc_len, desc_len))

    for i in trange(len(core_ds)):
        structure = core_ds[i]

        d = enforce_2d(descriptor(structure))

        post_cov += d.T @ d

    post_cov /= np.max(post_cov)

    Psi_prop = []

    for i in trange(len(samp_ds)):
        structure = samp_ds[i]
        Psi_prop.append(enforce_2d(descriptor(structure)))

    ds_samples = []

    # 1-shot sampling, updating post_cov each time
    for i in trange(np.max(N_structs)):
        PLU = scipy.linalg.lu(post_cov)
        scores = [score_func(Psi, PLU) for Psi in Psi_prop]
        selected_idx = np.argmax(scores)

        Psi = Psi_prop.pop(selected_idx)
        ds_samples.append(samp_ds.pop(selected_idx))

        post_cov += Psi.T @ Psi
        post_cov /= np.max(post_cov)

    # Readout the selected datasets
    samples = [
                [ds_samples[:N] + isos] for N in N_structs
        ]

    return samples