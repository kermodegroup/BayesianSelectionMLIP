from ase.io import read
from si_models import *
import numpy as np
import os
import json

dataset = read("../../Si/Si_Dataset.xyz", index=":")

def test_calc(calc):

    data = {"Total" : {"E_rmse": 0.0, "N_E" : 0, "F_rmse": 0.0, "N_F" : 0}}

    for image in dataset:

        config_type = image.info["config_type"]

        if config_type not in data.keys():
            data[config_type] = {"E_rmse": 0.0, "N_E" : 0, "F_rmse": 0.0, "N_F" : 0}

        image.calc = calc

        # Energy per atom
        E_err = (image.get_potential_energy() - image.info["dft_energy"]) / len(image)

        F_err = image.get_forces() - image.arrays["dft_force"]

        data[config_type]["E_rmse"] += E_err ** 2
        data[config_type]["F_rmse"] += np.sum(F_err.flatten() ** 2)
        data[config_type]["N_E"] += 1
        data[config_type]["N_F"] += 3 * len(image)

        data["Total"]["E_rmse"] += E_err ** 2
        data["Total"]["F_rmse"] += np.sum(F_err.flatten() ** 2)
        data["Total"]["N_E"] += 1
        data["Total"]["N_F"] += 3 * len(image)


    for key in data.keys():
        # Convert sum squared errs into RMSEs
        data[key]["E_rmse"] = np.sqrt(data[key]["E_rmse"] / data[key]["N_E"])
        data[key]["F_rmse"] = np.sqrt(data[key]["F_rmse"] / data[key]["N_F"])

    return data

calc_type = "ACE"

method = "MONTECARLO"
Nc = 20

if calc_type == "GAP":
    calc_fn = test_gap
elif calc_type == "ACE":
    calc_fn = test_ace

if os.path.exists(f"../Test_Results/{method}/{method}_Dataset_Errs_{calc_type}.json"):
    with open(f"../Test_Results/{method}/{method}_Dataset_Errs_{calc_type}.json", "r") as f:
        data = json.load(f)
else:
    data = {}

for N in [5, 10, 20, 50, 100]:
    print(N)
    
    data[str(N)] = {}

    all_calc_data = {}

    for i in range(Nc):
        print(N, i)
        calc = calc_fn(method, N, [i])[0]
        all_calc_data[i] = test_calc(calc)

    for key in all_calc_data[0].keys():
        data[str(N)][key] = {
            "E_rmse_mean" : np.mean([all_calc_data[i][key]["E_rmse"] for i in range(Nc)]),
            "E_rmse_std" : np.std([all_calc_data[i][key]["E_rmse"] for i in range(Nc)]),
            "F_rmse_mean" : np.mean([all_calc_data[i][key]["F_rmse"] for i in range(Nc)]),
            "F_rmse_std" : np.std([all_calc_data[i][key]["F_rmse"] for i in range(Nc)]),

            "E_rmse_raw_vals" : list([all_calc_data[i][key]["E_rmse"] for i in range(Nc)]),
            "F_rmse_raw_vals" : list([all_calc_data[i][key]["F_rmse"] for i in range(Nc)]),

            "N_E" : all_calc_data[0][key]["N_E"],
            "N_F" : all_calc_data[0][key]["N_F"]
        }

os.makedirs(f"../Test_Results/{method}", exist_ok=True)

# Ensure keys are sorted as integers
Ns = sorted([int(k) for k in data.keys()])

data = {N: data[str(N)] for N in Ns}

with open(f"../Test_Results/{method}/{method}_Dataset_Errs_{calc_type}.json", "w") as f:
    json.dump(data, f, indent=4)