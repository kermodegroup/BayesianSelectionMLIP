import os
import shutil

method = "SOAPCURMEAN"
#method = "MONTECARLO"

os.makedirs(f"Models/GAPs/{method}", exist_ok=True)

os.chdir(f"Models/GAPs/{method}")

data_src = f"../../../AL_Datasets/{method}"

xyz_files = [file for file in os.listdir(data_src) if ".xyz" in file]


for file in xyz_files:
    name = file[:-4]
    print(name)
    os.makedirs(name, exist_ok=True)
    os.chdir(name)

    if name + ".xyz" not in os.listdir():

        # Copy dataset
        shutil.copyfile("../" + data_src + os.sep + file, file)

        # Copy gap config
        shutil.copyfile("../" + data_src + os.sep + name + ".gapconfig", name + ".gapconfig")

        # Copy IP Glue core potential
        shutil.copyfile("../../../OriginalGAP/repulsive_2b.xml", "repulsive_2b.xml")

        # Run gap_fit to fit potential
        os.system(
            f"gap_fit config_file={name}.gapconfig"
        )

    os.chdir("..")
