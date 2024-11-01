# from ACEHAL.fit import fit
# from ACEHAL.basis import define_basis
from .utils import extract_isoats
import numpy as np
from sklearn.linear_model import BayesianRidge
from ase.io import write
import os
from tqdm import trange
import scipy

from ase.constraints import full_3x3_to_voigt_6_stress

# def energy_variance(committee, structures, **kwargs):
#     '''
#     Score structures based on the std of average energy per atom
#     '''

#     E_errs = np.zeros(len(structures))

#     for i, structure in enumerate(structures):
#         structure.calc = committee
        
#         E_errs[i] = committee.results_extra["err_energy"] / len(structure)

#     return E_errs

# def hal_force(committee, structures, eps=0.2, **kwargs):
#     '''
#     HAL force error metric
#     '''
#     F_errs = np.zeros(len(structures))

#     for i, structure in enumerate(structures):
#         structure.calc = committee

#         F_errs[i] = np.max(committee.results_extra["err_forces"] / 
#                         (np.linalg.norm(committee.results_extra["unbiased_forces"], axis=1) + eps))

#     return F_errs


# def select_hal_samples(dataset, basis_args, N_structs, core_ds=None, N_samples=1, score_func=hal_force, iso_atom_config_type=None, 
#                             solver=None, pot_file_root=None, pot_file_name=None, **fit_kwargs):
#     '''
#     Draw a sample of N structures from dataset, based on the HAL committee scores.

#     if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
#     sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
#     dataset: list of ase Atoms objects
#         Dataset for structure selection
#     committee: list of ase Calculator objects
#         Committee of calculators to use to determine standard deviation scores
#     N_structs: int or iterable of ints
#         Number of structures to select per sample. Use a list or array to draw multiple different sub-dataset sizes
#     N_samples: int
#         Number of samples to draw for each N_struct value
#     property: string
#         Property to perform committee std calculation on (see committee_scores fo more info)
#     iso_atom_config_type: string
#         config type (see atoms.info["config_type"]) for isolated atom configs

#     Returns:
#     samples: list
#         Samples of sub datasets, for each sample and N_struct value
#         Generated via [[sample(N) for i in range(N_samples)] for N in range N_structs]
    
#     '''


#     if np.issubdtype(type(N_structs), np.integer):
#         # Convert to len 1 list
#         N_structs = [N_structs]

#     samp_ds, isos = extract_isoats(dataset, iso_atom_config_type)

#     if len(isos) == 0 and core_ds is not None:
#         # No iso_ats found in sample dataset

#         core_ds, isos = extract_isoats(core_ds, iso_atom_config_type)

#     if core_ds is None:
#         core_ds = []
#     core_ds = core_ds + isos

#     E0s = {ats.get_chemical_symbols()[0] : ats.get_potential_energy() for ats in isos}

#     if solver is None:
#         solver = BayesianRidge

#     fit_kwargs["E0s"] = E0s
#     fit_kwargs["B_len_norm"] = define_basis(basis_args)

#     ds_samples = []

#     for isamp in trange(N_samples):
#         solver = solver()
#         samp_structs = [ats.copy() for ats in samp_ds]

#         ds = [ats.copy() for ats in core_ds]

#         curr_sample = []

#         committee = fit(ds, solver, pot_file=None, **fit_kwargs)

#         for istruct in trange(np.max(N_struct), leave=False):
#             if len(ds) - len(core_ds) + 1 in N_struct and pot_file_root is not None and pot_file_name is not None:
#                 fname = pot_file_root + f"_{len(ds)+1}_Sample_{isamp}" + os.sep + pot_file_name + f"_{len(ds)+1}_Sample_{isamp}"
#             else:
#                 fname = None

#             committee, samp_structs, ds, solver = hal_step(committee, samp_structs, dataset, solver, score_func=hal_force, model_fname=fname + ".json", score_kwargs={}, solver_kwargs=fit_kwargs)

#             if fname is not None:
#                 write(fname + ".xyz", ds)
#                 curr_sample.append([ats.copy() for ats in ds])

#         ds_samples.append(curr_sample)

#     # Rejig list ordering to get correct format, as HAL is faster to sample trajectories
#     reordered_samples = [
#         [
#             ds_samples[i][j]
#         for i in range(len(N_struct))]
#     for j in range(N_samples)]

#     return reordered_samples
            

# def hal_step(committee, samp_structures, dataset, solver, score_func=hal_force, model_fname=None, score_kwargs={}, solver_kwargs={}):
#     '''
#     Use committee and score_func to generate scores of each structure in structures.
#     Select a structure randomly weighted by the scores
#     Refit the committee using the solver
#     Optionally save the fitted model to file
#     '''

#     # Sco9re and select
#     scores = score_func(committee, structures, score_kwargs)

#     idxs = np.arange(len(structures))
#     selection = np.random.choice(idxs, scores/np.sum(scores))

#     dataset.append(samp_structures.pop(selection))

#     # Refit committee
#     committee = fit(dataset, solver, pot_file=model_fname, **solver_kwargs)

#     return committee, samp_structures, dataset, solver


def async_hal(dataset, descriptor, N_structs, core_ds=None, N_samples=1, iso_atom_config_type=None):
    '''
    Draw a sample of N structures from dataset, based on the asynchronous HAL committee scores.

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