from setuptools import setup

setup(
    name='ActiveLearnMLIPTests',
    version='0.1.0',    
    description='MLIP Active Learning Strategy Testing',
    author='Thomas Rocke',
    author_email='thomas.rocke@warwick.ac.uk',
    packages=['ActiveLearnMLIPTests'],
    install_requires=['matscipy',
                      'numpy',
                      'quippy',
                      'scikit-learn-extra',
                      'ase',
                      'mpi4py',
                      'mace-torch',
                      "tqdm"
                      ]
)