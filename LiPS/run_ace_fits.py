import os
import shutil
from julia import Main

method = "ACEAVGCUR"

os.makedirs(f"Models/ACEs/{method}", exist_ok=True)

os.chdir(f"Models/ACEs/{method}")

data_src = f"../../../AL_Datasets/{method}"

xyz_files = [file for file in os.listdir(data_src) if ".xyz" in file]

Main.eval("using ACEpotentials")

fit_ace = Main.eval('''
    function fit_ace(name)
        data_file = "./" * name * ".xyz"
        train = read_extxyz(data_file)

        model = acemodel(elements = [:Li, :P, :S],
                order = 3,
                totaldegree = 8,
                rcut = 5.0,
                Eref = [:Li => -0.02846111, :P => -0.0984573, :S => -0.09305436])

        @show length(model.basis);

        weights = Dict("default" => Dict("E" => 30, "F" => 1, "V" => 30))

        P = smoothness_prior(model; p = 4)

        solver = ACEfit.BLR(factorization=:svd)

        acefit!(model, train; solver=solver, prior=P, energy_key="energy", force_key="forces", 
               weights=weights, smoothness=5,
                export_json="./" * name * ".json")

    end
''')


for i, file in enumerate(xyz_files):
    name = file[:-4]
    print(f"{i+1}/{len(xyz_files)}", name)
    os.makedirs(name, exist_ok=True)
    os.chdir(name)

    if name + ".json" not in os.listdir():

        # Copy dataset
        shutil.copyfile("../" + data_src + os.sep + file, file)

        # Run ACE fit
        fit_ace(name)

    os.chdir("..")