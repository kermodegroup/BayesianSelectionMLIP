import os

def test_ace(method, N, sample_nums=range(20)):
    '''
    Get test ACE potentials corresponding to N_structs = N
    
    method: string
        Name of method to generate models for
    N: int
        N_structs value for models
    sample_nums: iterable
        Which samples to load models for. Default is all
    '''
    import pyjulip
    calcs = []

    for samp in sample_nums:
        calc_name = f"../Models/ACEs/{method}/LiPS_{method}_N_{N}_Sample_{samp}/LiPS_{method}_N_{N}_Sample_{samp}.json"
        if os.path.exists(calc_name):
            calcs.append(pyjulip.ACE1(calc_name))
    return calcs