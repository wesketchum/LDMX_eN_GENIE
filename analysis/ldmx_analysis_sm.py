from LDMX.Framework import EventTree
import numpy as np
import ROOT
from analysis_utilities import *
#from utils import *
#from array import array

#in ldmx-env
import libDetDescr

#common processing function for sim calhits
#returns e_cal vector, total energy, dictionary of energy contributions by given list of IDs, and unmatched energy
def process_sim_edeps(calhits,p_ids,section=-1):

    #print(f"Processing section {section}")

    sim_cal_e_vec = np.array([0.,0.,0.])
    total_sim_cal_e = 0.0
    sim_cal_e_unmatched = 0.0

    cntrb_edep_dict = dict.fromkeys(p_ids,0.0)
    
    for i_h, calhit in enumerate(calhits):
        this_section = (calhit.getID() >> 18) & 0x7
        #print(f"\thit {i_h}: e={calhit.getEdep()}, section={this_section}")
        if(section!=-1 and this_section!=section):
            continue
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

    #print(f"Total energy={total_sim_cal_e}, Total unmatched={sim_cal_e_unmatched}")
    return sim_cal_e_vec, total_sim_cal_e, cntrb_edep_dict, sim_cal_e_unmatched

def process_sim_edeps_lw(ecalhits,p_ids):
    
    total_sim_ecal_e = 0.0
    sim_ecal_e_unmatched = 0.0

    cntrb_edep_dict = dict.fromkeys(p_ids,0.0)

    for i_h, calhit in enumerate(ecalhits):
        #bit shifting
        layer_index = ( (calhit.getID() >> 17) & 0x3f)
        #in ldmx-env
        #layer_index = libDetDescr.EcalId(p_ids).layer()
        layer_weight = layerWeights[layer_index]
        total_sim_ecal_e += (1.+ layer_weight/mip_si_energy)*secondOrderEnergyCorrection*calhit.getEdep()
#        print(f'calhit is {calhit.getEdep()}')
#        print(f'prefactor is {(1.+ layer_weight/mip_si_energy)*secondOrderEnergyCorrection}')

        for i_c in range(calhit.getNumberOfContribs()):
            cntrb = calhit.getContrib(i_c)
            if(cntrb.incidentID in p_ids):
                cntrb_edep_dict[cntrb.incidentID] += (1.+layer_weight/mip_si_energy)*secondOrderEnergyCorrection*cntrb.edep
            else:
                sim_ecal_e_unmatched += cntrb.edep

    return total_sim_ecal_e, cntrb_edep_dict, sim_ecal_e_unmatched

def recoil_nhits(recoiltracks,p_ids):

    id_hit_dict = dict.fromkeys(p_ids,0)

    for pids in p_ids:
        id_hit_dict[pids] = 0
        for track in recoiltracks:
            if(track.getPdgID()!= 11 and track.getTrackID() == pids):
            #if(track.getTrackID() == pids):
                id_hit_dict[pids] += track.getNhits()

    return id_hit_dict

# pi0 photon energy correction funcitions
ECAL_CF_LW = 1.60
#HCAL_CF = [9.25,10.27,9.81,9.83,9.74,9.57,8.99,9.91]
HCAL_CF = [22., 11.1, 8.9, 8.4, 8.]

def ecal_corrected_energy(e):
#def ecal_corrected_energy(e,ph):
    #return ECAL_CF_LW*e[ph[0]]
    return ECAL_CF_LW*e

#def hcal_corrected_energy(e,ph):
#def hcal_corrected_energy(e):
#    if(e>=0 and e<100):
#        return HCAL_CF[0]*e
#    if(e>=100 and e<200):
#        return HCAL_CF[1]*e
#    if(e>=200 and e<300):
#        return HCAL_CF[2]*e
#    if(e>=300 and e<400):
#        return HCAL_CF[3]*e
#    if(e>=400 and e<500):
#        return HCAL_CF[4]*e
#    if(e>=500 and e<600):
#        return HCAL_CF[5]*e
#    if(e>=600 and e<800):
#        return HCAL_CF[6]*e
#    if(e>=800):
#        return HCAL_CF[7]*e

# def hcal energy as a function of deposited energy
def hcal_corrected_energy(e):
    if(e<1.5):
        return 25.*e
    if(e>=1.5 and e<10):
        return HCAL_CF[0]*e
    if(e>=10 and e<20):
        return HCAL_CF[1]*e
    if(e>=20 and e<30):
        return HCAL_CF[2]*e
    if(e>=30 and e<40):
        return HCAL_CF[3]*e
    if(e>=40):
        return HCAL_CF[4]*e
    
ELECTRON_PT_CUT = 0.
HCAL_PE_CUT = 0.
ECAL_PE_CUT = 0.

SIM_PARTICLE_P_CUT = 1.00 #MeV, for status != 1

# pi0 photon constants
ECAL_ENERGY_FRAC = 0.9
ECAL_ENERGY_MIN = 0.0
HCAL_ENERGY_FRAC = 0.9
HCAL_ENERGY_MIN = 0.500

def pi0_photon_sorter(ecal_e_lw,hcal_e):
    total_corrected_energy = ecal_corrected_energy(ecal_e_lw)+hcal_corrected_energy(hcal_e)
    if(ecal_e_lw > ECAL_ENERGY_MIN and ecal_corrected_energy(ecal_e_lw)/total_corrected_energy > ECAL_ENERGY_FRAC):
        return 1
    elif(hcal_e > HCAL_ENERGY_MIN and hcal_corrected_energy(hcal_e)/total_corrected_energy > HCAL_ENERGY_FRAC):
        return 2
    else:
        return 0

FILES = ['/Users/wketchum/Data/LDMX/production_02Aug2024/ldmx_genie_output_run_110000/ldmx_genie_G18_02a_02_11b_Ti_8GeV_110000_reco.root']
#	'ldmx_genie_output_run_101011/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101011_reco.root',
#	'ldmx_genie_output_run_101012/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101012_reco.root',
#	'ldmx_genie_output_run_101013/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101013_reco.root',
#	'ldmx_genie_output_run_101014/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101014_reco.root',
#	'ldmx_genie_output_run_101015/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101015_reco.root',
#	'ldmx_genie_output_run_101016/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101016_reco.root',
#	'ldmx_genie_output_run_101017/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101017_reco.root',
#	'ldmx_genie_output_run_101018/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101018_reco.root',
#	'ldmx_genie_output_run_101019/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101019_reco.root',
#	'ldmx_genie_output_run_101020/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101020_reco.root',
#	'ldmx_genie_output_run_101021/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101021_reco.root',
#	'ldmx_genie_output_run_101022/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101022_reco.root',
#	'ldmx_genie_output_run_101023/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101023_reco.root',
#	'ldmx_genie_output_run_101024/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101024_reco.root',
#	'ldmx_genie_output_run_101025/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101025_reco.root',
#	'ldmx_genie_output_run_101026/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101026_reco.root',
#	'ldmx_genie_output_run_101027/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101027_reco.root',
#	'ldmx_genie_output_run_101028/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101028_reco.root',
#	'ldmx_genie_output_run_101029/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101029_reco.root',
#	'ldmx_genie_output_run_101030/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101030_reco.root',
#	'ldmx_genie_output_run_101031/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101031_reco.root',
#	'ldmx_genie_output_run_101032/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101032_reco.root',
#	'ldmx_genie_output_run_101034/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101034_reco.root',
#	'ldmx_genie_output_run_101035/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101035_reco.root',
#	'ldmx_genie_output_run_101036/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101036_reco.root',
#	'ldmx_genie_output_run_101037/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101037_reco.root',
#	'ldmx_genie_output_run_101038/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101038_reco.root',
#	'ldmx_genie_output_run_101039/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101039_reco.root',
#	'ldmx_genie_output_run_101040/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101040_reco.root',
#	'ldmx_genie_output_run_101042/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101042_reco.root',
#	'ldmx_genie_output_run_101043/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101043_reco.root',
#	'ldmx_genie_output_run_101044/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101044_reco.root',
#	'ldmx_genie_output_run_101045/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101045_reco.root',
#	'ldmx_genie_output_run_101046/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101046_reco.root',
#	'ldmx_genie_output_run_101047/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101047_reco.root',
#	'ldmx_genie_output_run_101048/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101048_reco.root',
#	'ldmx_genie_output_run_101049/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101049_reco.root',
#	'ldmx_genie_output_run_101050/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101050_reco.root',
#	'ldmx_genie_output_run_101051/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101051_reco.root',
#	'ldmx_genie_output_run_101052/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101052_reco.root',
#	'ldmx_genie_output_run_101053/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101053_reco.root',
#	'ldmx_genie_output_run_101054/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101054_reco.root',
#	'ldmx_genie_output_run_101055/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101055_reco.root',
#	'ldmx_genie_output_run_101056/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101056_reco.root',
#	'ldmx_genie_output_run_101057/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101057_reco.root',
#	'ldmx_genie_output_run_101058/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101058_reco.root',
#	'ldmx_genie_output_run_101059/ldmx_genie_G18_02a_02_11b_Ti_8GeV_101059_reco.root']

EVENTS_TO_PROCESS = []

VERBOSE = False

OUTPUT_FILE = ROOT.TFile("output_file_31July_test.root","RECREATE")

#layer weight imported
layerWeights = [2.312, 4.312, 6.522, 7.490, 8.595, 10.253, 10.915, 10.915, 10.915, 10.915, 10.915,
                10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915,
                10.915, 10.915, 14.783, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539,
                18.539, 18.539, 9.938]
mip_si_energy = 0.13 #MeV
secondOrderEnergyCorrection = 4000./3940.5

#hadron_dict = {
#    "proton": [2212,-2212],
#    "neutron": [2112,-2112],
#    "piplus": [211],
#    "piminus": [-211],
#    "pi0": [111],
#    "Kplus": [321],
#    "Kminus": [-321],
#    "K0": [311,-311] }

hadron_list = {2212,-2212,2112,-2112,211,-211,111,321,-321,311,-311}

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

    #incoming electron params (to be updated later)
    "elec_pxi": (1,"D"),
    "elec_pyi": (1,"D"),
    "elec_pzi": (1,"D"),

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

    #endpoint info (mm)
    "sim_p_end_x": ("n_sim_p","D"),
    "sim_p_end_y": ("n_sim_p","D"),
    "sim_p_end_z": ("n_sim_p","D"),

    #vertex info (mm)
    "sim_p_vertex_x": ("n_sim_p","D"),
    "sim_p_vertex_y": ("n_sim_p","D"),
    "sim_p_vertex_z": ("n_sim_p","D"),

    #ecal hit positions (mm)
#    "ecal_x": ("n_sim_p","D"),
#    "ecal_y": ("n_sim_p","D"),
#    "ecal_z": ("n_sim_p","D"),

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

    #layer weight (lw) edep by particle in ecal
    "sim_p_ecal_e_lw":("n_sim_p","D"),

#    #hcal back
    "hcal_back_e": ("n_sim_p","D"),
#
#    #hcal side IDs
#    #Top = 1, Bottom = 2, Left = 4, Right = 3
    "hcal_top_e": ("n_sim_p","D"),
    "hcal_bottom_e": ("n_sim_p","D"),
    "hcal_left_e": ("n_sim_p","D"),
    "hcal_right_e": ("n_sim_p","D"),
#
#    #hcal section
#    "hcal_section": ("n_sim_p","D"),

    #splitting pi0 photon kinematics
    "n_sim_pi0": (1,"I"),
    "pi0_idx": ("n_sim_pi0","I"),
    "pi0_photon1_idx": ("n_sim_pi0","I"),
    "pi0_photon2_idx": ("n_sim_pi0","I"),
    "pi0_photon1_det": ("n_sim_pi0","I"), # 1 for mostly Ecal, 2 for mostly Hcal, 0 for neither
    "pi0_photon2_det": ("n_sim_pi0","I"), # 1 for mostly Ecal, 2 for mostly Hcal, 0 for neither

    #distance between ph t ns at face of ecal
    "sim_p_ecalx": ("n_sim_p","D"),
    "sim_p_ecaly": ("n_sim_p","D"),
    "ecald_xy": ("n_sim_p","D"),
    "t_param": ("n_sim_p","D"),

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
    "cal_et": (1,"D"),

    #delta pT (one per event?)
    "delta_pt": (1,"D"),
    "hsum_px": (1,"D"),
    "hsum_py": (1,"D"),
    "hsum_pz": (1,"D"),
    "delta_px": (1,"D"),
    "delta_py": (1,"D"),
    "delta_pz": (1,"D"),

    #recoil hits in tracker for charged particles
    "recoil_nhits":("n_sim_p","I"),

    #particle multiplicity
    "n_sim_prot": (1,"I"),
    "n_sim_neut": (1,"I"),
    "n_sim_pip": (1,"I"),
    "n_sim_pim": (1,"I")
    
#    #findable tracks
##    "reco_id":("n_rec_p","I"),
#    "n_rec_p":(1,"I"),
#    "reco_pdg":("n_rec_p","I")
}

output_tree, variables = create_tree_from_dict(var_dict)

for f in FILES:

    print(f'Processing file {f}')
    input_tree = EventTree.EventTree(f)

    odd_num_ph = []

    for ie, event in enumerate(input_tree):#        if(event.EventHeader.getEventNumber()>=15):
        #if(event.EventHeader.getEventNumber()>=30):
        #    break
        if(ie!=0 and ie%1000==0):
            print(f'\tProcessing event {ie}')
        
        if(len(EVENTS_TO_PROCESS)>0 and (event.EventHeader.getEventNumber() not in EVENTS_TO_PROCESS)): continue
    
        if(VERBOSE): event.EventHeader.Print()

        variables["run"][0] = event.EventHeader.getRun()
        variables["event"][0] = event.EventHeader.getEventNumber()
        
        n_sim_p = 0
        sim_p_id_dict = { }
        n_pi0 = 0
        pi0_dict = {}
        pi0_idx = []
        n_prot = 0
        n_neut = 0
        n_pip = 0
        n_pim = 0

        #initial electron variable initialization (to be updated later)
        variables["elec_pxi"][0] = 0.0
        variables["elec_pyi"][0] = 0.0
        variables["elec_pzi"][0] = 4000.0 #units are MeV (essentially beam energy )
        #set hsum to 0 for every event
        variables["hsum_px"][0] = 0.0
        variables["hsum_py"][0] = 0.0
        variables["hsum_pz"][0] = 0.0

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
        
#                #initialize nhits to zero (to be looped over later)
#                variables["recoil_nhits"][n_sim_p] = 0

                #endpoint info
                variables["sim_p_end_x"][n_sim_p] = particle.getEndPoint()[0]
                variables["sim_p_end_y"][n_sim_p] = particle.getEndPoint()[1]
                variables["sim_p_end_z"][n_sim_p] = particle.getEndPoint()[2]

                variables["sim_p_vertex_x"][n_sim_p] = particle.getVertex()[0]
                variables["sim_p_vertex_y"][n_sim_p] = particle.getVertex()[1]
                variables["sim_p_vertex_z"][n_sim_p] = particle.getVertex()[2]

                #distance at ecal face variables
                variables["t_param"][n_sim_p] = (240 - particle.getVertex()[2])/particle.getEndPoint()[2]
                variables["sim_p_ecalx"][n_sim_p] = particle.getVertex()[0] + variables["t_param"][n_sim_p]*particle.getEndPoint()[0]

                variables["sim_p_ecaly"][n_sim_p] = particle.getVertex()[1] + variables["t_param"][n_sim_p]*particle.getEndPoint()[1]

                #delta pT
                if(particle.getPdgID() in hadron_list):
                    variables["hsum_px"][0]+=particle.getMomentum()[0]
                    variables["hsum_py"][0]+=particle.getMomentum()[1]
                    variables["hsum_pz"][0]+=particle.getMomentum()[2]

                #hcalhit section for detector ID (not correct but bitshifting is)
                #for i_h, calhits in enumerate(event.HcalSimHits_genie):
                    #variables["hcal_section"][n_sim_p] = ((calhits.getID() >> 18) & 0x7)

                #initialize the edep per particle tracking ...
                sim_p_id_dict[i_p] = n_sim_p
                variables["sim_p_hcal_e"][n_sim_p] = 0.0
                variables["sim_p_ecal_e"][n_sim_p] = 0.0
                variables["sim_p_ecal_e_lw"][n_sim_p] = 0.0

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

                if(particle.getPdgID()==111):
                    #print(f'event is {event.EventHeader.getEventNumber()}')
                    variables["pi0_idx"][n_pi0] = n_sim_p#index of pi0 in our list
                    pi0_dict[i_p] = []
                    n_pi0 +=1
                    pi0_idx.append(i_p)
                    #print(f'pi0 index is: {pi0_idx}')

                variables["n_sim_pi0"][0] = n_pi0

#		#particle multiplicity
#                if(particle.getPdgID()==2212):
#                    n_prot +=1
#                variables["n_sim_prot"][0] = n_prot
#
#                if(particle.getPdgID()==2112):
#                    n_neut +=1
#                variables["n_sim_neut"][0] = n_neut
#
#                if(particle.getPdgID()==211):
#                    n_pip +=1
#                variables["n_sim_pip"][0] = n_pip
#
#                if(particle.getPdgID()==-211):
#                    n_pim +=1
#                variables["n_sim_pim"][0] = n_pim
#
                n_sim_p = n_sim_p+1
                
        if(n_pi0>0):
            for x in range(n_sim_p):
                if(variables["sim_p_pdg"][x]==22 and variables["sim_p_parent_pdg"][x]==111):# and variables["sim_p_parent_id"][n_sim_p]==x):
                    if(variables["sim_p_parent_id"][x] in pi0_dict.keys()):
                        pi0_dict[variables["sim_p_parent_id"][x]].append(x)

        if(n_pi0>0):
            for my_pi0_id in range(len(pi0_idx)):
                pi0idx = pi0_idx[my_pi0_id]
                if((len(pi0_dict[pi0idx]))%2==0): # some weird length issue to resolve later
                    if(variables["sim_p_e"][pi0_dict[pi0idx][0]] > variables["sim_p_e"][pi0_dict[pi0idx][1]]):
                        variables["pi0_photon1_idx"][my_pi0_id] = pi0_dict[pi0idx][0]
                        variables["pi0_photon2_idx"][my_pi0_id] = pi0_dict[pi0idx][1]
                    elif(variables["sim_p_e"][pi0_dict[pi0idx][0]] < variables["sim_p_e"][pi0_dict[pi0idx][1]]):
                        variables["pi0_photon1_idx"][my_pi0_id] = pi0_dict[pi0idx][1]
                        variables["pi0_photon2_idx"][my_pi0_id] = pi0_dict[pi0idx][0]

                elif((len(pi0_dict[pi0idx]))%2==1):
                    odd_num_ph.append(event.EventHeader.getEventNumber())
#                    #print(odd_num_ph)


        #event level
        #delta pT continuation
        variables["delta_px"][0] = variables["elec_px"][0] - variables["elec_pxi"][0] + variables["hsum_px"][0]
        variables["delta_py"][0] = variables["elec_py"][0] - variables["elec_pyi"][0] + variables["hsum_py"][0]
        variables["delta_pz"][0] = variables["elec_pz"][0] - variables["elec_pzi"][0] + variables["hsum_pz"][0]
#        print(variables["delta_px"])
        variables["delta_pt"][0] = np.sqrt(variables["delta_px"][0]*variables["delta_px"][0]+variables["delta_py"][0]*variables["delta_py"][0])
#        print(f'\tdelta pt for event {event.EventHeader.getEventNumber()} is {variables["delta_pt"][0]}')
#        print(f'\tdelta px for the event is {variables["delta_px"][0]}')
#        print(f'\thsum px for the event is {variables["hsum_px"][0]}')


        variables["n_sim_p"][0] = n_sim_p

        sim_hcal_e_vec, total_sim_hcal_e, sim_hcal_cntrb_edep_dict, sim_hcal_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys())
        sim_ecal_e_vec, total_sim_ecal_e, sim_ecal_cntrb_edep_dict, sim_ecal_e_um = process_sim_edeps(event.EcalSimHits_genie,sim_p_id_dict.keys())

        sim_hcal_back_e_vec, total_sim_hcal_back_e, sim_hcal_back_cntrb_edep_dict, sim_hcal_back_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),0)
        sim_hcal_top_e_vec, total_sim_hcal_top_e, sim_hcal_top_cntrb_edep_dict, sim_hcal_top_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),1)
        sim_hcal_bottom_e_vec, total_sim_hcal_bottom_e, sim_hcal_bottom_cntrb_edep_dict, sim_hcal_bottom_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),2)
        sim_hcal_right_e_vec, total_sim_hcal_right_e, sim_hcal_right_cntrb_edep_dict, sim_hcal_right_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),3)
        sim_hcal_left_e_vec, total_sim_hcal_left_e, sim_hcal_left_cntrb_edep_dict, sim_hcal_left_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),4)

        #ecal layer weights
        total_sim_ecal_e_lw, sim_ecal_cntrb_edep_dict_lw, sim_ecal_e_um_lw = process_sim_edeps_lw(event.EcalSimHits_genie,sim_p_id_dict.keys())

        for pid, edep_sum in sim_hcal_cntrb_edep_dict.items():
            variables["sim_p_hcal_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_ecal_cntrb_edep_dict.items():
            variables["sim_p_ecal_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_ecal_cntrb_edep_dict_lw.items():
            variables["sim_p_ecal_e_lw"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_hcal_back_cntrb_edep_dict.items():
            variables["hcal_back_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_hcal_top_cntrb_edep_dict.items():
            variables["hcal_top_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_hcal_bottom_cntrb_edep_dict.items():
            variables["hcal_bottom_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_hcal_right_cntrb_edep_dict.items():
            variables["hcal_right_e"][sim_p_id_dict[pid]] += edep_sum

        for pid, edep_sum in sim_hcal_left_cntrb_edep_dict.items():
            variables["hcal_left_e"][sim_p_id_dict[pid]] += edep_sum
            #print(sim_hcal_left_cntrb_edep_dict)

        #print(sim_hcal_top_cntrb_edep_dict)
        #print(sim_hcal_bottom_cntrb_edep_dict)
        #print(sim_hcal_right_cntrb_edep_dict)
        #print(sim_hcal_left_cntrb_edep_dict)
        recoil_hits_dict = recoil_nhits(event.RecoilTracks_genie,sim_p_id_dict.keys())

        for pid, nhits in recoil_hits_dict.items():
            variables["recoil_nhits"][sim_p_id_dict[pid]] = nhits
#            if(variables["recoil_nhits"][sim_p_id_dict[pid]]!=0):
#                print(f'event is {event.EventHeader.getEventNumber()}')
#                print(f'nhits is {variables["recoil_nhits"][sim_p_id_dict[pid]]}')
#                print(f'particle id is {variables["sim_p_id"][sim_p_id_dict[pid]]}')
#
#        n_rec_p = 0
#        
#        for rechit in event.RecoilTracks_genie_reco:
##            variables["reco_id"][n_rec_p] = i_p
#            variables["reco_pdg"][n_rec_p] = rechit.getPdgID()
#
#            n_rec_p = n_rec_p + 1
#
#        variables["n_rec_p"][0] = n_rec_p

# photon sorter
        if(n_pi0>0):
            for my_pi0_id in range(len(pi0_idx)):
               if((len(pi0_dict[pi0idx]))%2==0):
                    variables["pi0_photon1_det"][my_pi0_id] = pi0_photon_sorter(variables["sim_p_ecal_e_lw"][variables["pi0_photon1_idx"][my_pi0_id]],variables["sim_p_hcal_e"][variables["pi0_photon1_idx"][my_pi0_id]])
                    variables["pi0_photon2_det"][my_pi0_id] = pi0_photon_sorter(variables["sim_p_ecal_e_lw"][variables["pi0_photon2_idx"][my_pi0_id]],variables["sim_p_hcal_e"][variables["pi0_photon2_idx"][my_pi0_id]])
#                    print(variables["pi0_photon1_det"][my_pi0_id])

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

    if(len(odd_num_ph)>0):
        file = open("odd_num_ph_eventnumber.txt","w+")
        events = str(odd_num_ph)
        file.write(events)
        file.close()
#        np.savetxt("odd_num_ph_eventnumber.txt",odd_num_ph)

OUTPUT_FILE.Write()
OUTPUT_FILE.Close()

