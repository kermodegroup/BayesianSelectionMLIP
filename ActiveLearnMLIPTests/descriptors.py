import numpy as np
from quippy.descriptors import Descriptor


def soap_descriptor(desc_string):
    '''
    Generate a callable interface for a soap descriptor defined from the given string
    
    '''
    desc_object = Descriptor(desc_string)

    def inner(atoms):
        res = desc_object.calc(atoms)
        return res["data"]
    
    return inner


def ace_descriptor(*args, **kwargs):
    '''
    Generate a callable interface for an ACE descriptor using Julia, PyCall, & ACEPotentials.jl
    
    '''

    desc_object = ...

    def inner(atoms):

        return desc_object(atoms)
    return inner