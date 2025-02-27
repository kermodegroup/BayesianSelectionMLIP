# BayesianSelection
 
## Installation

### Julia Installation
To use ACE models and descriptor features, the code requires Julia, with the ACEpotentials (v0.6), ASE, and JuLIP packages installed.

A suitable Julia environment is provided via the `Project.toml` and `Manifest.toml` files, which contain a minimal number of required packages. Switching Julia to this project can be achieved by setting the environment variable `export JULIA_PROJECT=.`.


To install from scratch, follow the installation instructions below. The code is modified from the [ACEpotentials v0.6 docs](https://acesuit.github.io/ACEpotentials.jl/v0.6/gettingstarted/installation/), and is to be ran from a `julia` repl:

```
using Pkg
Pkg.activate(".")
Pkg.Registry.add("General")  # only needed when installing Julia for the first time
Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
Pkg.add(Pkg.PackageSpec(; name="ACEpotentials", version="0.6"))
Pkg.add("ASE")
Pkg.add("JuLIP")
Pkg.add("ACE1")
```

In order for the Python calls to julia to be using the correct Julia project, we require the environment variable `JULIA_PROJECT` to be set.


### Python installation
The Python requirements can be installed via:
```
pip install .
python -c "import julia; julia.install()"
```

## Reproducing the publication
The tools required to reproduce the publication are hosted in the `Si_Surfaces` directory. See the specific documentation [here](Si_Surfaces/README.md)
