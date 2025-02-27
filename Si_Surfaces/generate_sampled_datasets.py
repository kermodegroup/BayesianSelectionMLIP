from sampling_config import *
from Tests.plot_config import *

###
# Generate sampled datasets from Si_Surface_Structs.xyz and save in ./AL_Datasets/
# methods_to_gen is a list of method names, which are parameterised by sampling_config.method_params[method_name]

random_seed = 42

# methods_to_gen = method_comparison_plots produces datasets for the models in Figs. 1 & 2
# methods_to_gen = descriptor_comparison_plots produces datasets for the models in Figs. 3 & 4
methods_to_gen = method_comparison_plots

for method in methods_to_gen:
    print(method)

    os.makedirs(f"AL_Datasets/{method}", exist_ok=True)
    # Re-seed randomness to ensure reproducibility
    np.random.seed(random_seed)

    params = method_params[method]

    func = params["method"]

    func_args = global_params.copy()

    if "method_args" in params.keys():
        func_args.update(params["method_args"])

    samples = func(dataset, **func_args)

    for i, N in enumerate(func_args["N_structs"]):
        smp = samples[i]
        for j in range(func_args["N_samples"]):
            ds_name = f"Si_{method}_N_{N}_Sample_{j}.xyz"
            gap_name = f"Si_{method}_N_{N}_Sample_{j}.xml" 
            gap_config_name = f"Si_{method}_N_{N}_Sample_{j}.gapconfig" 

            write(f"AL_Datasets/{method}/{ds_name}", smp[j] + core_dataset)

            with open("Models/GAPs/base_gap_config", "r") as base_gap_config_file:
                    with open(f"AL_Datasets/{method}/{gap_config_name}", "w") as new_config:
                        new_config.writelines(base_gap_config_file) # create a copy of base config
                        new_config.writelines(
                            [
                                f"gp_file={gap_name}\n",
                                f"at_file={ds_name}"
                            ]
                        )
print("Done")