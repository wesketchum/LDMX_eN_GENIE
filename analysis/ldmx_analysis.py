from LDMX.Framework import EventTree
import numpy as np
import ROOT
#from array import array

ELECTRON_PT_CUT = 0.
HCAL_PE_CUT = 0.
ECAL_PE_CUT = 0.

SIM_PARTICLE_P_CUT = 1.00 #MeV, for status != 1

FILENAME = 'test_Ti.root'

EVENTS_TO_PROCESS = []
#EVENTS_TO_PROCESS = [84,66]

VERBOSE = False

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

def thetaz(particle):
    if(p(particle)<1e-6):
        return 0
    return np.arccos(particle.getMomentum()[2]/p(particle))

input_tree = EventTree.EventTree(FILENAME)

#var_dict used to create the tree
# keys: names of the branches we will create
# val: tuple with size (either int array size,
# or if string the entry that will give it its length) and ROOT type
# right now, ROOT type is either D or I

var_dict = {
    "run": (1,"I"),
    "event": (1,"I"),

    #outgoing electron params
    "elec_px": (1,"D"),
    "elec_py": (1,"D"),
    "elec_pz": (1,"D"),
    "elec_pt": (1,"D"),
    "elec_p": (1,"D"),
    "elec_e": (1,"D"),
    "elec_thetaz": (1,"D"),

    "n_sim_p": (1,"I"),
    "sim_p_id": ("n_sim_p","I"),
    "sim_p_status": ("n_sim_p","I"),
    "sim_p_pdg": ("n_sim_p","I"),
    "sim_p_e": ("n_sim_p","D"),
    "sim_p_px": ("n_sim_p","D"),
    "sim_p_py": ("n_sim_p","D"),
    "sim_p_pz": ("n_sim_p","D"),
    "sim_p_pt": ("n_sim_p","D"),
    "sim_p_p": ("n_sim_p","D"),
    "sim_p_thetaz": ("n_sim_p","D"),
    "sim_p_parent": ("n_sim_p","I"),
    
    #hcal
    "hcal_e": (1,"D"),
    "hcal_et": (1,"D"),

    #ecal
    "ecal_e": (1,"D"),
    "ecal_et": (1,"D"),

    #total cal
    "cal_e": (1,"D"),
    "cal_et": (1,"D")
    }

output_file = ROOT.TFile("output_file.root","RECREATE")

output_tree, variables = create_tree_from_dict(var_dict)
    
for ie, event in enumerate(input_tree):

    if(len(EVENTS_TO_PROCESS)>0 and (event.EventHeader.getEventNumber() not in EVENTS_TO_PROCESS)): continue
    
    if(VERBOSE): event.EventHeader.Print()

    variables["run"][0] = event.EventHeader.getRun()
    variables["event"][0] = event.EventHeader.getEventNumber()
        
    n_sim_p = 0
    for i_p, particle in event.SimParticles_genie:

        mom = p(particle)

        #get electron params
        if(particle.getPdgID()==11 and particle.getGenStatus()==1):

            variables["elec_px"][0] = particle.getMomentum()[0]
            variables["elec_py"][0] = particle.getMomentum()[1]
            variables["elec_pz"][0] = particle.getMomentum()[2]
            variables["elec_pt"][0] = pt(particle)
            variables["elec_p"][0] = p(particle)
            variables["elec_e"][0] = particle.getEnergy()
            variables["elec_thetaz"][0] = thetaz(particle)
            
            if(VERBOSE): print(f'\tElectron pt is {variables["elec_pt"][0]}')

        if ( particle.getGenStatus()==1 or mom>SIM_PARTICLE_P_CUT ):
            variables["sim_p_id"][n_sim_p] = i_p
            variables["sim_p_status"][n_sim_p] = particle.getGenStatus()
            variables["sim_p_pdg"][n_sim_p] = particle.getPdgID()
            variables["sim_p_px"][n_sim_p] = particle.getMomentum()[0]
            variables["sim_p_py"][n_sim_p] = particle.getMomentum()[1]
            variables["sim_p_pz"][n_sim_p] = particle.getMomentum()[2]
            variables["sim_p_pt"][n_sim_p] =  pt(particle)
            variables["sim_p_p"][n_sim_p] = p(particle)
            variables["sim_p_e"][n_sim_p] = particle.getEnergy()
            variables["sim_p_thetaz"][n_sim_p] = thetaz(particle)
            variables["sim_p_parent"][n_sim_p] = particle.getParents()[0]            
            n_sim_p = n_sim_p+1
            
    variables["n_sim_p"][0] = n_sim_p
    
    if variables["elec_pt"][0] < ELECTRON_PT_CUT:
        continue

    hcal_e_dir = np.array([0.,0.,0.])
    total_hcal_e = 0
    for i_h, hcalhit in enumerate(event.HcalRecHits_genie):
        hcal_e_dir[0] += hcalhit.getEnergy()*hcalhit.getXPos()
        hcal_e_dir[1] += hcalhit.getEnergy()*hcalhit.getYPos()
        hcal_e_dir[2] += hcalhit.getEnergy()*hcalhit.getZPos()
        total_hcal_e += hcalhit.getEnergy()
    if(total_hcal_e>0):
        hcal_e_dir = hcal_e_dir / total_hcal_e
        hcal_e_dir = hcal_e_dir / np.sqrt(hcal_e_dir[0]*hcal_e_dir[0]+
                                              hcal_e_dir[1]*hcal_e_dir[1]+
                                              hcal_e_dir[2]*hcal_e_dir[2])
    hcal_et = total_hcal_e*np.sqrt(hcal_e_dir[0]*hcal_e_dir[0]+
                                       hcal_e_dir[1]*hcal_e_dir[1])
    if(VERBOSE): print(f'\tTotal HCAL E={total_hcal_e}, Et = {hcal_et}')

    ecal_e_dir = np.array([0.,0.,0.])
    total_ecal_e = 0
    for i_h, ecalhit in enumerate(event.EcalRecHits_genie):
        ecal_e_dir[0] += ecalhit.getEnergy()*ecalhit.getXPos()
        ecal_e_dir[1] += ecalhit.getEnergy()*ecalhit.getYPos()
        ecal_e_dir[2] += ecalhit.getEnergy()*ecalhit.getZPos()
        total_ecal_e += ecalhit.getEnergy()
    if(total_ecal_e>0):
        ecal_e_dir = ecal_e_dir / total_ecal_e
        ecal_e_dir = ecal_e_dir / np.sqrt(ecal_e_dir[0]*ecal_e_dir[0]+
                                              ecal_e_dir[1]*ecal_e_dir[1]+
                                              ecal_e_dir[2]*ecal_e_dir[2])
    ecal_et = total_ecal_e*np.sqrt(ecal_e_dir[0]*ecal_e_dir[0]+
                                       ecal_e_dir[1]*ecal_e_dir[1])
    if(VERBOSE): print(f'\tTotal ECAL E={total_ecal_e}, Et = {ecal_et}')

    total_cal_e = total_hcal_e+total_ecal_e
    cal_e_vec = np.array([total_ecal_e*ecal_e_dir[0]+total_hcal_e*hcal_e_dir[0],
                              total_ecal_e*ecal_e_dir[1]+total_hcal_e*hcal_e_dir[1],
                              total_ecal_e*ecal_e_dir[2]+total_hcal_e*hcal_e_dir[2]])
    cal_et = np.sqrt(cal_e_vec[0]*cal_e_vec[0]+cal_e_vec[1]*cal_e_vec[1])
    if(VERBOSE): print(f'\tTotal CAL E={total_cal_e}, Et = {cal_et}')

    variables["hcal_e"][0] = total_hcal_e
    variables["hcal_et"][0] = hcal_et
    variables["ecal_e"][0] = total_ecal_e
    variables["ecal_et"][0] = ecal_et
    variables["cal_e"][0] = total_cal_e
    variables["cal_et"][0] = cal_et    

    if(VERBOSE):
        print("List of particles:")
        for i_p, particle in event.SimParticles_genie:
            
            print(f'ID {i_p}: PDG {particle.getPdgID()}, Status {particle.getGenStatus()}, Energy {particle.getEnergy()}, Momentum {particle.getMomentum()}, thetaz {np.degrees(thetaz(particle))}, Parents {particle.getParents()}')

    output_tree.Fill()

output_file.Write()
output_file.Close()

