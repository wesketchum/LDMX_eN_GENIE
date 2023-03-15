from LDMX.Framework import EventTree
import numpy as np
import ROOT
from utils import *
#from array import array

#common processing function for sim calhits
#returns e_cal vector, total energy, dictionary of energy contributions by given list of IDs, and unmatched energy
def process_sim_edeps(calhits,p_ids):

    sim_cal_e_vec = np.array([0.,0.,0.])
    total_sim_cal_e = 0.0
    sim_cal_e_unmatched = 0.0

    cntrb_edep_dict = dict.fromkeys(p_ids,0.0)
    
    for i_h, calhit in enumerate(calhits):
        hit_pos = calhit.getPosition()
        sim_cal_e_vec[0] += calhit.getEdep()*hit_pos[0]
        sim_cal_e_vec[1] += calhit.getEdep()*hit_pos[1]
        sim_cal_e_vec[2] += calhit.getEdep()*hit_pos[2]
        total_sim_cal_e += calhit.getEdep()
            
        for i_c in range(calhit.getNumberOfContribs()):
            cntrb = calhit.getContrib(i_c)
            if(cntrb.incidentID in p_ids):
                cntrb_edep_dict[cntrb.incidentID] += cntrb.edep
            else:
                sim_cal_e_unmatched += cntrb.edep

    return sim_cal_e_vec, total_sim_cal_e, cntrb_edep_dict, sim_cal_e_unmatched


ELECTRON_PT_CUT = 0.
HCAL_PE_CUT = 0.
ECAL_PE_CUT = 0.

SIM_PARTICLE_P_CUT = 1.00 #MeV, for status != 1

#FILENAME = 'ldmx_genie_G18_02a_00_000_Ti_101.root'

FILES = [ 'ldmx_genie_G18_02a_00_000_Ti_101.root',
          'ldmx_genie_G18_02a_00_000_Ti_102.root',
          'ldmx_genie_G18_02a_00_000_Ti_103.root',
          'ldmx_genie_G18_02a_00_000_Ti_104.root',
          'ldmx_genie_G18_02a_00_000_Ti_105.root',
          'ldmx_genie_G18_02a_00_000_Ti_106.root',
          'ldmx_genie_G18_02a_00_000_Ti_107.root',
          'ldmx_genie_G18_02a_00_000_Ti_108.root',
          'ldmx_genie_G18_02a_00_000_Ti_109.root',
          'ldmx_genie_G18_02a_00_000_Ti_110.root']             

EVENTS_TO_PROCESS = []
#EVENTS_TO_PROCESS = [84,66]

VERBOSE = False

OUTPUT_FILE = ROOT.TFile("output_file.root","RECREATE")

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
    "sim_p_m": ("n_sim_p","D"),
    "sim_p_q": ("n_sim_p","D"),
    "sim_p_thetaz": ("n_sim_p","D"),

    #parent info
    "sim_p_parent_id": ("n_sim_p","I"),
    "sim_p_parent_status": ("n_sim_p","I"),
    "sim_p_parent_pdg": ("n_sim_p","I"),
    "sim_p_parent_e": ("n_sim_p","D"),
    "sim_p_parent_px": ("n_sim_p","D"),
    "sim_p_parent_py": ("n_sim_p","D"),
    "sim_p_parent_pz": ("n_sim_p","D"),
    "sim_p_parent_pt": ("n_sim_p","D"),
    "sim_p_parent_p": ("n_sim_p","D"),
    "sim_p_parent_m": ("n_sim_p","D"),
    "sim_p_parent_q": ("n_sim_p","D"),
    "sim_p_parent_thetaz": ("n_sim_p","D"),

    
    #edep by particle in ecal/hcal
    "sim_p_ecal_e":("n_sim_p","D"),
    "sim_p_hcal_e":("n_sim_p","D"),

    #edeps for unmatched particles
    "sim_ecal_e_um": (1,"D"),
    "sim_hcal_e_um": (1,"D"),    
    
    "sim_hcal_e": (1,"D"),
    "sim_hcal_et": (1,"D"),
    "sim_ecal_e": (1,"D"),
    "sim_ecal_et": (1,"D"),
    "sim_cal_e": (1,"D"),
    "sim_cal_et": (1,"D"),
    
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

output_tree, variables = create_tree_from_dict(var_dict)

for f in FILES:

    print(f'Processing file {f}')
    input_tree = EventTree.EventTree(f)

    for ie, event in enumerate(input_tree):

        if(ie!=0 and ie%1000==0):
            print(f'\tProcessing event {ie}')
        
        if(len(EVENTS_TO_PROCESS)>0 and (event.EventHeader.getEventNumber() not in EVENTS_TO_PROCESS)): continue
    
        if(VERBOSE): event.EventHeader.Print()

        variables["run"][0] = event.EventHeader.getRun()
        variables["event"][0] = event.EventHeader.getEventNumber()
        
        n_sim_p = 0
        sim_p_id_dict = { }
    
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

            #particle.getGenStatus()==1 or 
            if ( mom>SIM_PARTICLE_P_CUT ):
                variables["sim_p_id"][n_sim_p] = i_p
                variables["sim_p_status"][n_sim_p] = particle.getGenStatus()
                variables["sim_p_pdg"][n_sim_p] = particle.getPdgID()
                variables["sim_p_px"][n_sim_p] = particle.getMomentum()[0]
                variables["sim_p_py"][n_sim_p] = particle.getMomentum()[1]
                variables["sim_p_pz"][n_sim_p] = particle.getMomentum()[2]
                variables["sim_p_pt"][n_sim_p] =  pt(particle)
                variables["sim_p_p"][n_sim_p] = p(particle)
                variables["sim_p_e"][n_sim_p] = particle.getEnergy()
                variables["sim_p_m"][n_sim_p] = particle.getMass()
                variables["sim_p_q"][n_sim_p] = particle.getCharge()
                variables["sim_p_thetaz"][n_sim_p] = thetaz(particle)

                #initialize the edep per particle tracking ...
                sim_p_id_dict[i_p] = n_sim_p
                variables["sim_p_hcal_e"][n_sim_p] = 0.0
                variables["sim_p_ecal_e"][n_sim_p] = 0.0

                #fill parent info
                p_id = particle.getParents()[0]
                variables["sim_p_parent_id"][n_sim_p] = p_id
                if(p_id in sim_p_id_dict.keys()):
                    p_my_id = sim_p_id_dict[p_id]
                    variables["sim_p_parent_status"][n_sim_p] = variables["sim_p_status"][p_my_id]
                    variables["sim_p_parent_pdg"][n_sim_p] = variables["sim_p_pdg"][p_my_id]
                    variables["sim_p_parent_px"][n_sim_p] = variables["sim_p_px"][p_my_id]
                    variables["sim_p_parent_py"][n_sim_p] = variables["sim_p_py"][p_my_id]
                    variables["sim_p_parent_pz"][n_sim_p] = variables["sim_p_pz"][p_my_id]
                    variables["sim_p_parent_pt"][n_sim_p] = variables["sim_p_pt"][p_my_id]
                    variables["sim_p_parent_p"][n_sim_p] = variables["sim_p_p"][p_my_id]
                    variables["sim_p_parent_e"][n_sim_p] = variables["sim_p_e"][p_my_id]
                    variables["sim_p_parent_m"][n_sim_p] = variables["sim_p_m"][p_my_id]
                    variables["sim_p_parent_q"][n_sim_p] = variables["sim_p_q"][p_my_id]
                    variables["sim_p_parent_thetaz"][n_sim_p] = variables["sim_p_thetaz"][p_my_id]                    
                
                n_sim_p = n_sim_p+1
            
        variables["n_sim_p"][0] = n_sim_p

        sim_hcal_e_vec, total_sim_hcal_e, sim_hcal_cntrb_edep_dict, sim_hcal_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys())
        sim_ecal_e_vec, total_sim_ecal_e, sim_ecal_cntrb_edep_dict, sim_ecal_e_um = process_sim_edeps(event.EcalSimHits_genie,sim_p_id_dict.keys())

        for pid, edep_sum in sim_hcal_cntrb_edep_dict.items():
            variables["sim_p_hcal_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_ecal_cntrb_edep_dict.items():
            variables["sim_p_ecal_e"][sim_p_id_dict[pid]] += edep_sum
    
        variables["sim_hcal_e"][0] = total_sim_hcal_e
        variables["sim_hcal_et"][0] = np.sqrt(sim_hcal_e_vec[0]*sim_hcal_e_vec[0]+sim_hcal_e_vec[1]*sim_hcal_e_vec[1])
        variables["sim_hcal_e_um"][0] = sim_hcal_e_um
        variables["sim_ecal_e"][0] = total_sim_ecal_e
        variables["sim_ecal_et"][0] = np.sqrt(sim_ecal_e_vec[0]*sim_ecal_e_vec[0]+sim_ecal_e_vec[1]*sim_ecal_e_vec[1])
        variables["sim_ecal_e_um"][0] = sim_ecal_e_um
        variables["sim_cal_e"][0] = total_sim_hcal_e+total_sim_ecal_e
        sim_cal_e_vec = sim_hcal_e_vec+sim_ecal_e_vec
        variables["sim_cal_et"][0] = np.sqrt(sim_cal_e_vec[0]*sim_cal_e_vec[0]+sim_cal_e_vec[1]*sim_cal_e_vec[1])    

        if variables["elec_pt"][0] < ELECTRON_PT_CUT:
            output_tree.Fill()
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


OUTPUT_FILE.Write()
OUTPUT_FILE.Close()

