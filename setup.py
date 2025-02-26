from setuptools import setup

setup(
    name='ActiveLearnMLIPTests',
    version='0.1.0',    
    description='MLIP Active Learning Strategy Testing',
    author='Thomas Rocke',
    author_email='thomas.rocke@warwick.ac.uk',
    packages=['ActiveLearnMLIPTests'],
    install_requires=['numpy',
                      'scipy',
                      'quippy-ase',
                      'scikit-learn-extra',
                      'ase',
                      'mpi4py',
                      'mace-torch',
                      'tqdm',
                      'wfl',
                      'julia',
                      'https://github.com/casv2/pyjulip.git@master'
                      ]
)
