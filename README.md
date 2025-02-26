# ActiveLearnMLIPTests
 
## Installation

### Julia Installation
To use ACE models and descriptor features, the code requires Julia, with the ACEpotentials (v0.6), ASE, and JuLIP packages installed.

The installation instructions below are modified from the [ACEpotentials v0.6 docs](https://acesuit.github.io/ACEpotentials.jl/v0.6/gettingstarted/installation/)

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

In order for the Python calls to julia to be using the correct Julia project, we set the environment variable `export JULIA_PROJECT=.`.

### Python installation
The Python requirements can be installed via:
```
pip install .
python -c "import julia; julia.install()"
```

## Reproducing the publication
The tools required to reproduce the publication are hosted in the `Si_Surfaces` directory. See the specific documentation [here](Si_Surfaces/README.md)
