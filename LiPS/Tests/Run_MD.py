from wfl.generate import md
from lips_models import test_ace
from wfl.configset import ConfigSet, OutputSpec

temps = []
dt = 2
nsteps = 10_000
thermo = "langevin"

method = "ACEAVGKMED"
Nc = 1

for T in temps:
    for N in [20]:
        for i in range(Nc):
            calc_spec = (test_calc, [], {'method': method, "N": N, "sample_nums" : [i]})

            outfile = OutputSpec(f"../Test_Results/{method}/MD_Traj/{method}_{N}_{i}_{T}.xyz")
