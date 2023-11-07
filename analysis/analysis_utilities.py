import numpy as np
import ROOT

TYPE_DICT = {
    "D": np.float64,
    "I": np.int32
}

def create_tree_from_dict(vars_dict,tree_name="ana_tree",tree_title="eN Analysis Tree",max_array=999):
    myvars = {}
    output_tree = ROOT.TTree(tree_name,tree_title)
    for vname, pars in vars_dict.items():

        if pars[1] not in TYPE_DICT:
            print(f"Unsupported type {pars[1]}, not registering var {vname}")
            continue
            
        if isinstance(pars[0],int):
            myvars[vname] = np.zeros(pars[0],dtype=TYPE_DICT[pars[1]])
            output_tree.Branch(vname,myvars[vname],f"{vname}/{pars[1]}")
        elif isinstance(pars[0],str):
            myvars[vname] = np.zeros(max_array,dtype=TYPE_DICT[pars[1]])
            output_tree.Branch(vname,myvars[vname],f"{vname}[{pars[0]}]/{pars[1]}")
        else:
            print(f"Cannot handle type {type(pars[0])} for var {vname}")
            continue

    return output_tree, myvars

def pt(particle):
    return np.sqrt(particle.getMomentum()[0]*particle.getMomentum()[0]+
                       particle.getMomentum()[1]*particle.getMomentum()[1])
def p(particle):
    return np.sqrt(particle.getMomentum()[0]*particle.getMomentum()[0]+
                       particle.getMomentum()[1]*particle.getMomentum()[1]+
                       particle.getMomentum()[2]*particle.getMomentum()[2])

def thetaz(particle, mom_threshold=1e-6):
    if(p(particle)<mom_threshold):
        return 0
    return np.arccos(particle.getMomentum()[2]/p(particle))
