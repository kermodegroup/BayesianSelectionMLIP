###
# Configuration for which models/methods go with which plot type



descriptor_comparison_plots = [ # Descriptor comparison, do kmedoids for each descriptor plus MONTECARLO as a baseline 
    "MONTECARLO",
    #"SOAPAVGKMED",
    "ACEAVGKMED",
    "ACESAVGKMED",
    "ACELAVGKMED",
    "MP0AVGKMED",
    #"MPAAVGKMED",
    "MACEAVGKMED",
    "MACETAVGKMED",
    #"SOAPAVGFPS",
    "ACEAVGFPS",
    "MP0AVGFPS",
    #"MPAAVGFPS",
    "MACEAVGFPS",
    "MACETAVGFPS",
    #"SOAPBLR",
    "ACEBLR",
    "ACESBLR",
    "ACELBLR",
    "MP0BLR",
    #"MPABLR",
    "MACEBLR",
    "MACETBLR",
]

method_comparison_plots = [ # Sparsification method comparison, do all methods using ACE descriptor plus MONTECARLO as a baseline 
    "MONTECARLO",
    "ACEAVGKMED",
    "ACEAVGFPS",
    "ACECURMEAN",
    "ACEAVGCUR",
    "ACEBLR",
]

desc_comp_kwargs = {
    "MONTECARLO" : {"color": "k", "linestyle": "solid"},
    "SOAPAVGKMED" : {"color": "C0", "linestyle": "solid"},
    "ACEAVGKMED" : {"color": "C1", "linestyle": "solid"},
    "MP0AVGKMED" : {"color": "C2", "linestyle": "solid"},
    "MACEAVGKMED" : {"color": "C3", "linestyle": "solid"},
    "MACETAVGKMED" : {"color": "C4", "linestyle": "solid"},
    "SOAPAVGFPS" : {"color": "C0", "linestyle": "solid"},
    "ACEAVGFPS" : {"color": "C1", "linestyle": "solid"},
    "MP0AVGFPS" : {"color": "C2", "linestyle": "solid"},
    "MACEAVGFPS" : {"color": "C3", "linestyle": "solid"},
    "MACETAVGFPS" : {"color": "C4", "linestyle": "solid"},
    "SOAPBLR" : {"color": "C0", "linestyle" : "solid"},
    "ACEBLR" : {"color": "C1", "linestyle" : "solid"},
    "MP0BLR" : {"color": "C2", "linestyle" : "solid"},
    "MACEBLR" : {"color": "C3", "linestyle" : "solid"},
    "MACETBLR" : {"color": "C4", "linestyle": "solid"},
    "ACESBLR" : {"color": "C5", "linestyle" : "solid"},
    "ACELBLR" : {"color": "C6", "linestyle" : "solid"},
    "ACESAVGKMED" : {"color": "C5", "linestyle" : "solid"},
    "ACELAVGKMED" : {"color": "C6", "linestyle" : "solid"},
    "MPABLR" : {"color": "C7", "linestyle" : "solid"},
    "MPAAVGFPS" : {"color": "C7", "linestyle" : "solid"},
    "MPAAVGKMED" : {"color": "C7", "linestyle" : "solid"},
}

method_names = {
    "MONTECARLO" : "Monte Carlo (SRS)",
    "AVGKMED" : " + k-medoids",
    "AVGFPS" : " + FPS",
    #"BLR" : r" + Bayesian Selection ($\mathrm{Var}(E_T)$)",
    "BLR" : r" + Bayesian Selection",
    "AVGCUR" : " + CUR (Structure Features)",
    "CURMEAN" : " + CUR (Avg of atomic scores)"
}



def desc_names(name):
    if "SOAP" in name:
        return "SOAP" + method_names[name[4:]]
    elif "MACET" in name:
        return "MACE (Total)" + method_names[name[5:]]
    elif "MACE" in name:
        return "MACE (Core)" + method_names[name[4:]]
    elif "ACE" in name:
        if name[3] == "S":
            # ACE S
            return "ACE (S)" + method_names[name[4:]]
        elif name[3] == "L":
            # ACE L
            return "ACE (L)" + method_names[name[4:]]
        else:
            return "ACE (M)" + method_names[name[3:]]
    elif "MP0" in name:
        return "MP0 Foundation" + method_names[name[3:]]
    elif "MPA" in name:
        return "MPA Foundation" + method_names[name[3:]]
    else:
        return method_names[name]

method_comp_kwargs = {
    "MONTECARLO" : {"color" : "k", "linestyle" : "solid"},
    "ACEAVGKMED" : {"color" : "C0", "linestyle" : "solid"},
    "ACEAVGFPS" : {"color" : "C1", "linestyle" : "solid"},
    "ACECURMEAN" : {"color" : "C2", "linestyle" : "solid"},
    "ACEAVGCUR" : {"color" : "C3", "linestyle" : "solid"},
    "ACEBLR" : {"color": "C4", "linestyle" : "solid"},
    "MACEBLR" : {"color": "C4", "linestyle" : "dashed"},
    "MP0BLR" : {"color": "C5", "linestyle" : "dotted"},
}


def split_mod(mod):
    if "MACET" in mod:
        desc = "MACE(T)"
        method = mod[5:]
    elif "MACE" in mod:
        desc = "MACE(C)"
        method = mod[4:]
    elif "ACE" in mod:
        print(mod, mod[3])
        if mod[3] == "S":
            # ACE S
            desc = "ACE (S)"
            method = mod[4:]
        elif mod[3] == "L":
            # ACE S
            desc = "ACE (L)"
            method = mod[4:]
        else:
            desc = "ACE (M)"
            method = mod[3:]
    elif "SOAP" in mod:
        desc = "SOAP"
        method = mod[4:]
    elif "MP0" in mod:
        desc = "MP0"
        method = mod[3:]
    else:
        desc = ""
        method = mod

    return desc, method


def get_plot_props(desc, method, pot, desc_comp=False):
    if not desc_comp:
        return {
            "color":colours[method], "marker":markers[desc], "linestyle":linestyles[pot]
        }
    else: 
        return {
            "color":desc_comp_colours[desc], "marker":"x", "linestyle":linestyles[pot]
        }
