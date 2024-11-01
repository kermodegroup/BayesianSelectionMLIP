
descriptor_comparison_plots = [ # Descriptor comparison, do kmedoids for each descriptor plus MONTECARLO as a baseline 
    "MONTECARLO",
    "SOAPAVGKMED",
    "ACEAVGKMED",
    "MP0AVGKMED",
    "MACEAVGKMED",
    "MACETAVGKMED",
    "SOAPAVGFPS",
    "ACEAVGFPS",
    "MP0AVGFPS",
    "MACEAVGFPS",
    "MACETAVGFPS",
    "SOAPBLR",
    "ACEBLR",
    "MP0BLR",
    "MACEBLR",
    "MACETBLR"
]

method_comparison_plots = [ # Sparsification method comparison, do all methods using ACE descriptor plus MONTECARLO as a baseline 
    "MONTECARLO",
    "ACEAVGKMED",
    "ACEAVGFPS",
    "ACECURMEAN",
    "ACEAVGCUR",
    "ACEBLR"
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
    "MACETBLR" : {"color": "C4", "linestyle": "solid"}
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
        return "ACE" + method_names[name[3:]]
    elif "MP0" in name:
        return "MP0 Foundation" + method_names[name[3:]]
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
        desc = "ACE"
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
