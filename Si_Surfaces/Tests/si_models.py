from quippy.potential import Potential
import os

def original_gap():
    return Potential(param_filename="../Models/OriginalGAP/gp_iter6_sparse9k.xml")


def test_gap(method, N, sample_nums=range(10)):
    '''
    Get test GAP potentials corresponding to N_structs = N
    
    method: string
        Name of method to generate models for
    N: int
        N_structs value for models
    sample_nums: iterable
        Which samples to load models for. Default is all
    '''
    calcs = []

    for samp in sample_nums:
        xml_name = f"../Models/GAPs/{method}/Si_{method}_N_{N}_Sample_{samp}/Si_{method}_N_{N}_Sample_{samp}.xml"

        if os.path.exists(xml_name):

            calcs.append(
                Potential(param_filename=xml_name)
            )
    return calcs



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
        calc_name = f"../Models/ACEs/{method}/Si_{method}_N_{N}_Sample_{samp}/Si_{method}_N_{N}_Sample_{samp}.json"
        if os.path.exists(calc_name):
            calcs.append(pyjulip.ACE1(calc_name))
    return calcs