# BayesianSelection
 
## Installation

### Prerequisites
- Python 3.9 - 3.11
- Julia 1.10
- An MPI implementation compatible with `mpi4py` (e.g. OpenMPI)

### Julia Installation
To use ACE models and descriptor features, the code requires Julia, with the ACEpotentials (v0.6), ASE, and JuLIP packages installed.

A suitable Julia environment is provided via the `Project.toml` and `Manifest.toml` files, which contain a minimal number of required packages. Switching Julia to this project can be achieved by setting the environment variable `export JULIA_PROJECT=.`, which is required so that importing the `julia` package from Python uses the correct project.

The packages can then be installed from a julia repl using:
```
using Pkg
Pkg.instantiate()
```

### Python installation
The Python requirements can be installed via:
```
pip install .
python -c "import julia; julia.install()"
```

## Reproducing the publication
The tools required to reproduce the publication are hosted in the `Si_Surfaces` directory. See the specific documentation [here](Si_Surfaces/README.md)
