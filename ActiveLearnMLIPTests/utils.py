import numpy as np

def extract_isoats(dataset, iso_atom_config_type=None):
    '''
    From a dataset of structures with some isolated atom configs, separate out into remainder + iso_ats.

    dataset: list of ase Atoms objects
        dataset of structures to parse
    iso_atom_config_type: string
        config type (see atoms.info["config_type"]) for isolated atom configs

    Returns:
    remainder: list of ase Atoms objects
        non-isolated atom structures in dataset
    iso_ats: list of ase Atoms objects
        isolated atom structures 
    
    '''
    if iso_atom_config_type is None:
        return dataset, []
    
    iso_ats = []
    remainder= []

    for image in dataset:
        is_iso = False
        if "config_type" in image.info.keys():
            if image.info["config_type"] == iso_atom_config_type:
                is_iso = True
        if is_iso:
            iso_ats.append(image)
        else:
            remainder.append(image)

    return remainder, iso_ats

def get_dataset_descriptors(dataset, descriptor):
    '''
    Gets all descriptor vectors for a full dataset

    dataset: list of ase Atoms objects
        Dataset to evaluate descriptor vectors
    descriptor: callable
        Descriptor function with signature descriptor(atoms) -> np array of shape (len(atoms), descriptor_length)

    Returns:
    vecs: np array of shape (np.sum([len(ats) for ats in dataset]), descriptor_length)
        Descriptor vecs for full dataset
    idxs: np array of shape (np.sum([len(ats) for ats in dataset]))
        Structure indexes for each descriptor vector. (idxs == i) creates a mask for all descriptor vectors for the ith structure in dataset
    '''
    nats = sum([len(struct) for struct in dataset])

    descriptor_length = descriptor(dataset[0]).shape[-1]

    vecs = np.zeros((nats, descriptor_length))
    
    struct_idxs = np.zeros((nats), dtype=int)

    nats = 0

    for i, struct in enumerate(dataset):
        res = descriptor(struct)

        vecs[nats:nats+len(struct), :] = res

        nats += len(struct)
        
        struct_idxs[nats:nats+len(struct)] = i
    
    struct_idxs = np.array(struct_idxs)

    return vecs, struct_idxs