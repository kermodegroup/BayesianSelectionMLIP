import numpy as np
from .utils import extract_isoats


def draw_committee_samples(dataset, committee, N_structs, N_samples=1, aggregate=np.average, iso_atom_config_type=None):
    '''
    Draw a sample of N structures from dataset, based on the committee scores.

    if isolated atom structures are found with iso_atom_config_type, they are separated from the rest of the dataset (and thus the 
    sampling), and are appended to each sample. This means each sample dataset will have the correct isolated atom structures.

    
    dataset: list of ase Atoms objects
        Dataset for structure selection
    committee: list of ase Calculator objects
        Committee of calculators to use to determine standard deviation scores
    N_structs: int or iterable of ints
        Number of structures to select per sample. Use a list or array to draw multiple different sub-dataset sizes
    N_samples: int
        Number of samples to draw for each N_struct value
    property: string
        Property to perform committee std calculation on (see committee_scores fo more info)
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

    scores = committee_scores(ds, committee, property)

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

def committee_scores(structures, committee, property="Energy"):
    '''
    Compute "committee scores", based on the normalised distribution of committee errors.

    s_i = std(p_j for j in committee) / sum(s_k)

    structures: list of ase Atoms objects
        Structures to compute scores of
    committee: list of ase Calculator objects
        list of committors
    property: string
        Property to use as error metric.
        "Energy" -> image.get_potential_energy()
        "Forces" -> image.get_forces()


    Returns:
    scores: np.array
        Scores for each structure
    '''
    
    def energy(image):
        return [image.get_potential_energy()]

    def force(image):
        return image.get_forces().flatten()

    if property.lower() in ["energy", "energies"]:
        prop = energy
    elif property.lower() in ["force", "forces"]:
        prop = force
    else:
        prop = energy

    nstruct = len(structures)
    
    scores = np.zeros((nstruct))

    for i, image in enumerate(structures):
        Ps = []
        for j, comm in enumerate(committee):
            image.calc = comm
            Ps[j] = prop(image)
        
        scores[i] = np.std(np.array(Ps))

    # Normalise scores
    scores /= np.sum(scores)

    return scores