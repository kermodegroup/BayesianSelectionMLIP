from setuptools import setup

setup(
    name='BayesianSelection',
    version='1.0.0',    
    description='MLIP Dataset Bayesian selection',
    author='Thomas Rocke',
    author_email='thomas.rocke@warwick.ac.uk',
    packages=['BayesianSelection'],
    install_requires=['numpy<2.0',
                      'scipy',
                      'quippy-ase',
                      'scikit-learn-extra',
                      'ase',
                      'mpi4py',
                      'mace-torch',
                      'tqdm',
                      'wfl',
                      'julia',
                      'pyjulip @ git+https://github.com/casv2/pyjulip.git@master'
                      ]
)
