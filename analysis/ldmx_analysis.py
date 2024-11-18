from LDMX.Framework import EventTree
import numpy as np
import ROOT
from analysis_utilities import *
import glob
#from array import array

import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True)

arg = parser.parse_args()

#import libDetDescr
LAYER_WEIGHTS = [2.312, 4.312, 6.522, 7.490, 8.595, 10.253, 10.915, 10.915, 10.915, 10.915, 10.915,
                 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915, 10.915,
                 10.915, 10.915, 14.783, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539, 18.539,
                 18.539, 18.539, 9.938]
MIP_SI_ENERGY = 0.13 #MeV
ENERGY_CORRECTION = 4000./3940.5

#common processing function for sim calhits
#returns e_cal vector, total energy, dictionary of energy contributions by given list of IDs, and unmatched energy
def process_sim_edeps(calhits,p_ids,section=-1,use_layer_weights=False):

    sim_cal_e_vec = np.array([0.,0.,0.])
    total_sim_cal_e = 0.0
    sim_cal_e_unmatched = 0.0

    cntrb_edep_dict = dict.fromkeys(p_ids,0.0)
    
    for i_h, calhit in enumerate(calhits):
        this_section = (calhit.getID() >> 18) & 0x7

        if(section!=-1 and this_section!=section):
            continue

        hit_pos = calhit.getPosition()
        this_edep_e = calhit.getEdep()
        this_edep_e_corr = 1.0;

        if use_layer_weights:
            this_layer_index = (calhit.getID() >> 17) & 0x3f
            this_layer_weight = LAYER_WEIGHTS[this_layer_index]
            this_edep_e_corr = (1+this_layer_weight/MIP_SI_ENERGY)*ENERGY_CORRECTION

        sim_cal_e_vec[0] += this_edep_e*this_edep_e_corr*hit_pos[0]
        sim_cal_e_vec[1] += this_edep_e*this_edep_e_corr*hit_pos[1]
        sim_cal_e_vec[2] += this_edep_e*this_edep_e_corr*hit_pos[2]
        total_sim_cal_e += this_edep_e*this_edep_e_corr
            
        for i_c in range(calhit.getNumberOfContribs()):
            cntrb = calhit.getContrib(i_c)
            if(cntrb.incidentID in p_ids):
                cntrb_edep_dict[cntrb.incidentID] += cntrb.edep*this_edep_e_corr
            else:
                sim_cal_e_unmatched += cntrb.edep*this_edep_e_corr

    return sim_cal_e_vec, total_sim_cal_e, cntrb_edep_dict, sim_cal_e_unmatched

# pi0 photon energy correction funcitions

def ecal_corrected_energy(e,ECAL_CF_LW=1.60):
    return ECAL_CF_LW*e

# def hcal energy as a function of deposited energy
def hcal_corrected_energy(e,HCAL_CF=[50,21.85,10.92,8.95,8.25,8.05]):
    if(e<1.5):
        return HCAL_CF[0]*e
    if(e>=1.5 and e<10):
        return HCAL_CF[1]*e
    if(e>=10 and e<20):
        return HCAL_CF[2]*e
    if(e>=20 and e<30):
        return HCAL_CF[3]*e
    if(e>=30 and e<40):
        return HCAL_CF[4]*e
    if(e>=40):
        return HCAL_CF[5]*e

def pi0_photon_sorter(ecal_e_lw,hcal_e,ECAL_ENERGY_FRAC=0.9,ECAL_ENERGY_MIN=0.0,HCAL_ENERGY_FRAC=0.9,HCAL_ENERGY_MIN=0.500):
    total_corrected_energy = ecal_corrected_energy(ecal_e_lw)+hcal_corrected_energy(hcal_e)
    if(ecal_e_lw > ECAL_ENERGY_MIN and ecal_corrected_energy(ecal_e_lw)/total_corrected_energy > ECAL_ENERGY_FRAC):
        return 1
    elif(hcal_e > HCAL_ENERGY_MIN and hcal_corrected_energy(hcal_e)/total_corrected_energy > HCAL_ENERGY_FRAC):
        return 2
    else:
        return 0

ELECTRON_PT_CUT = 0.
HCAL_PE_CUT = 0.
ECAL_PE_CUT = 0.

SIM_PARTICLE_P_CUT = 1.00 #MeV, for status != 1

#EVENTS_TO_PROCESS = []

VERBOSE = False

def string_to_4vector(istring):
    return np.array(istring.split(" "),dtype='f')

def unpack_genie_interaction_attribute(hepmc_event,attr_name):
    try:
        if attr_name=="FSLepton":
            pdg = 11
        else:
            pdg = int(hepmc_event.attribute_as_string(f"GENIE.Interaction.{attr_name}PDG"))
        p4 = string_to_4vector(hepmc_event.attribute_as_string(f"GENIE.Interaction.{attr_name}P4"))
        p4 = 1000*p4
        return pdg, p4[0], p4[1], p4[2], p4[3]
    except:
        return 0, 0.0, 0.0, 0.0, 0.0

GENIE_KVAR_SELX = 11
GENIE_KVAR_SELY = 12
GENIE_KVAR_SELQ2 = 13
GENIE_KVAR_SELW = 15

def process_file(ifile):

    ofile = ifile[:-5]+"_ana.root"
    OUTPUT_FILE = ROOT.TFile(ofile,"RECREATE")
    output_tree, variables = create_tree_from_dict(var_dict)

    print(f'Processing file {ifile}. Output: {ofile}')
    input_tree = EventTree.EventTree(ifile)

    for ie, event in enumerate(input_tree):

        #if(ie!=0 and ie%1000==0):
        #    print(f'\tProcessing event {ie}')
        
        #if(len(EVENTS_TO_PROCESS)>0 and (event.EventHeader.getEventNumber() not in EVENTS_TO_PROCESS)): continue
    
        if(VERBOSE): event.EventHeader.Print()

        variables["run"][0] = event.EventHeader.getRun()
        variables["event"][0] = event.EventHeader.getEventNumber()

        hepmc_event = event.SimHepMC3Events_genie[0].getHepMCGenEvent()

        (variables["genie_lep_i_pdg"][0],
         variables["genie_lep_i_px"][0],
         variables["genie_lep_i_py"][0],
         variables["genie_lep_i_pz"][0],
         variables["genie_lep_i_e"][0]) = unpack_genie_interaction_attribute(hepmc_event,"Probe")

        (variables["genie_lep_f_pdg"][0],
         variables["genie_lep_f_px"][0],
         variables["genie_lep_f_py"][0],
         variables["genie_lep_f_pz"][0],
         variables["genie_lep_f_e"][0]) = unpack_genie_interaction_attribute(hepmc_event,"FSLepton")

        (variables["genie_tgt_pdg"][0],
         variables["genie_tgt_px"][0],
         variables["genie_tgt_py"][0],
         variables["genie_tgt_pz"][0],
         variables["genie_tgt_e"][0]) = unpack_genie_interaction_attribute(hepmc_event,"Target")

        (variables["genie_hnuc_pdg"][0],
         variables["genie_hnuc_px"][0],
         variables["genie_hnuc_py"][0],
         variables["genie_hnuc_pz"][0],
         variables["genie_hnuc_e"][0]) = unpack_genie_interaction_attribute(hepmc_event,"HitNucleon")

        variables["genie_scatter"][0] = int(hepmc_event.attribute_as_string("GENIE.Interaction.ScatteringType"))

        kvar_labels = hepmc_event.attribute_as_string("GENIE.Interaction.KineVarLabels").split(" ")
        kvar_vals = hepmc_event.attribute_as_string("GENIE.Interaction.KineVarValues").split(" ")

        kvar_dict = {}
        for i_k in range(len(kvar_labels)):
            kvar_dict[int(kvar_labels[i_k])] = float(kvar_vals[i_k])

        if (GENIE_KVAR_SELX in kvar_dict.keys()): variables["genie_x_bj"][0] = kvar_dict[GENIE_KVAR_SELX]
        if (GENIE_KVAR_SELY in kvar_dict.keys()): variables["genie_y_inel"][0] = kvar_dict[GENIE_KVAR_SELY]
        if (GENIE_KVAR_SELQ2 in kvar_dict.keys()): variables["genie_Q2"][0] = kvar_dict[GENIE_KVAR_SELQ2]*1000
        if (GENIE_KVAR_SELW in kvar_dict.keys()): variables["genie_W"][0] = kvar_dict[GENIE_KVAR_SELW]*1000


        variables["genie_qel"][0] = (variables["genie_scatter"][0]==1)
        variables["genie_dis"][0] = (variables["genie_scatter"][0]==3)
        variables["genie_res"][0] = (variables["genie_scatter"][0]==4)
        variables["genie_mec"][0] = (variables["genie_scatter"][0]==10)

        if hasattr(event,"genieEventWeights_genie"):
            weights = event.genieEventWeights_genie.getWeights()
            variables["n_weights"][0] = len(weights)
            for i_w, w in enumerate(weights):
                variables["weights"][i_w] = w
        else:
            variables["n_weights"][0] = 0
            variables["weights"][0] = 1.0

        n_sim_p = 0
        sim_p_id_dict = { }
        n_pi0 = 0
        pi0_dict = {}
        pi0_idx = []
        n_prot = 0
        n_neut = 0
        n_pip = 0
        n_pim = 0

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

                #initialize the edep per particle tracking ...
                sim_p_id_dict[i_p] = n_sim_p
                variables["sim_p_hcal_e"][n_sim_p] = 0.0
                variables["sim_p_ecal_e"][n_sim_p] = 0.0
                variables["sim_p_ecal_e_lw"][n_sim_p] = 0.0
                variables["hcal_back_e"][n_sim_p] = 0.0
                variables["hcal_right_e"][n_sim_p] = 0.0
                variables["hcal_top_e"][n_sim_p] = 0.0
                variables["hcal_bottom_e"][n_sim_p] = 0.0
                variables["hcal_left_e"][n_sim_p] = 0.0

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
                    variables["pi0_idx"][n_pi0] = n_sim_p#index of pi0 in our list
                    pi0_dict[i_p] = []
                    n_pi0 +=1
                    pi0_idx.append(i_p)
                elif(particle.getPdgID()==211):
                    n_pip +=1
                elif(particle.getPdgID()==-211):
                    n_pim +=1
                elif(particle.getPdgID()==2112):
                    n_prot +=1
                elif(particle.getPdgID()==2212):
                    n_neut +=1

                n_sim_p = n_sim_p+1
            
        variables["n_sim_p"][0] = n_sim_p
        variables["n_sim_pi0"][0] = n_pi0
        variables["n_sim_pip"][0] = n_pip
        variables["n_sim_pim"][0] = n_pim
        variables["n_sim_prot"][0] = n_prot
        variables["n_sim_neut"][0] = n_neut

        if(n_pi0>0):
            for x in range(n_sim_p):
                if(variables["sim_p_pdg"][x]==22 and variables["sim_p_parent_pdg"][x]==111):# and variables["sim_p_parent_id"][n_sim_p]==x):
                    if(variables["sim_p_parent_id"][x] in pi0_dict.keys()):
                        pi0_dict[variables["sim_p_parent_id"][x]].append(x)

            for my_pi0_id in range(len(pi0_idx)):
                pi0idx = pi0_idx[my_pi0_id]
                if((len(pi0_dict[pi0idx]))%2==0): # some weird length issue to resolve later
                    if(variables["sim_p_e"][pi0_dict[pi0idx][0]] > variables["sim_p_e"][pi0_dict[pi0idx][1]]):
                        variables["pi0_photon1_idx"][my_pi0_id] = pi0_dict[pi0idx][0]
                        variables["pi0_photon2_idx"][my_pi0_id] = pi0_dict[pi0idx][1]
                    elif(variables["sim_p_e"][pi0_dict[pi0idx][0]] < variables["sim_p_e"][pi0_dict[pi0idx][1]]):
                        variables["pi0_photon1_idx"][my_pi0_id] = pi0_dict[pi0idx][1]
                        variables["pi0_photon2_idx"][my_pi0_id] = pi0_dict[pi0idx][0]

                #elif((len(pi0_dict[pi0idx]))%2==1):
                #    odd_num_ph.append(event.EventHeader.getEventNumber())

        sim_hcal_e_vec, total_sim_hcal_e, sim_hcal_cntrb_edep_dict, sim_hcal_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys())
        sim_ecal_e_vec, total_sim_ecal_e, sim_ecal_cntrb_edep_dict, sim_ecal_e_um = process_sim_edeps(event.EcalSimHits_genie,sim_p_id_dict.keys())

        sim_hcal_back_e_vec, total_sim_hcal_back_e, sim_hcal_back_cntrb_edep_dict, sim_hcal_back_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),section=0)
        sim_hcal_top_e_vec, total_sim_hcal_top_e, sim_hcal_top_cntrb_edep_dict, sim_hcal_top_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),section=1)
        sim_hcal_bottom_e_vec, total_sim_hcal_bottom_e, sim_hcal_bottom_cntrb_edep_dict, sim_hcal_bottom_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),section=2)
        sim_hcal_right_e_vec, total_sim_hcal_right_e, sim_hcal_right_cntrb_edep_dict, sim_hcal_right_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),section=3)
        sim_hcal_left_e_vec, total_sim_hcal_left_e, sim_hcal_left_cntrb_edep_dict, sim_hcal_left_e_um = process_sim_edeps(event.HcalSimHits_genie,sim_p_id_dict.keys(),section=4)

        #ecal layer weights
        sim_ecal_e_vec_lw, total_sim_ecal_e_lw, sim_ecal_cntrb_edep_dict_lw, sim_ecal_e_um_lw = process_sim_edeps(event.EcalSimHits_genie,sim_p_id_dict.keys(),use_layer_weights=True)

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

        # photon sorter
        if(n_pi0>0):
            for my_pi0_id in range(len(pi0_idx)):
                if((len(pi0_dict[pi0idx]))%2==0):
                    variables["pi0_photon1_det"][my_pi0_id] = pi0_photon_sorter(variables["sim_p_ecal_e_lw"][variables["pi0_photon1_idx"][my_pi0_id]],variables["sim_p_hcal_e"][variables["pi0_photon1_idx"][my_pi0_id]])
                    variables["pi0_photon2_det"][my_pi0_id] = pi0_photon_sorter(variables["sim_p_ecal_e_lw"][variables["pi0_photon2_idx"][my_pi0_id]],variables["sim_p_hcal_e"][variables["pi0_photon2_idx"][my_pi0_id]])

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

        max_layer = 0
        for i_s, sum in enumerate(event.ecalTrigSums_genie):
            if not sum.energy()>0.0: continue
            if sum.layer()>=max_layer:
                max_layer = sum.layer()+1
            for i_layer in range(sum.layer()):
                variables["trig_Ecal_e_afterLayer"][i_layer] += sum.energy()
        variables["trig_Ecal_nLayer"][0] = max_layer

        max_layer = 0
        for i_s, sum in enumerate(event.hcalTrigQuadsBackLayerSums_genie):
            if not sum.hwEnergy()>0.0: continue
            if sum.layer()>=max_layer:
                max_layer = sum.layer()+1
            for i_layer in range(sum.layer()):
                variables["trig_Hcal_e_afterLayer"][i_layer] += sum.hwEnergy()
        variables["trig_Hcal_nLayer"][0] = max_layer

        variables["trig_n_Electron"][0] = len(event.trigElectrons_genie)
        variables["trig_maxE"][0] = -1
        variables["trig_maxPt"][0] = -1
        max_e_val = 0.0
        max_pt_val = 0.0
        for i_tE, tEle in enumerate(event.trigElectrons_genie):
            if tEle.energy()>max_e_val:
                max_e_val = tEle.energy()
                variables["trig_maxE"][0] = i_tE
            if tEle.pt()>max_pt_val:
                max_pt_val = tEle.pt()
                variables["trig_maxPt"][0] = i_tE
            variables["trig_Electron_e"][i_tE] = tEle.energy()
            variables["trig_Electron_eClus"][i_tE] = tEle.getClusEnergy()
            variables["trig_Electron_zClus"][i_tE] = tEle.endz()
            variables["trig_Electron_px"][i_tE] = tEle.px()
            variables["trig_Electron_py"][i_tE] = tEle.py()
            variables["trig_Electron_pz"][i_tE] = tEle.pz()
            variables["trig_Electron_dx"][i_tE] = tEle.endx() - tEle.vx()
            variables["trig_Electron_dy"][i_tE] = tEle.endy() - tEle.vy()
            variables["trig_Electron_x"][i_tE] = tEle.vx()
            variables["trig_Electron_y"][i_tE] = tEle.vy()
            variables["trig_Electron_tp"][i_tE] = tEle.getClusTP()
            variables["trig_Electron_depth"][i_tE] = tEle.getClusDepth()

        if(VERBOSE):
            print("List of particles:")
            for i_p, particle in event.SimParticles_genie:
            
                print(f'ID {i_p}: PDG {particle.getPdgID()}, Status {particle.getGenStatus()}, Energy {particle.getEnergy()}, Momentum {particle.getMomentum()}, thetaz {np.degrees(thetaz(particle))}, Parents {particle.getParents()}')

        output_tree.Fill()


    OUTPUT_FILE.Write()
    OUTPUT_FILE.Close()




#OUTPUT_FILE = ROOT.TFile("/Users/wketchum/Data/LDMX/ldmx_genie_G18_02a_02_11b_Ti_3_ana.root","RECREATE")

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

    "genie_lep_i_pdg": (1,"I"),
    "genie_lep_i_px": (1,"D"),
    "genie_lep_i_py": (1,"D"),
    "genie_lep_i_pz": (1,"D"),
    "genie_lep_i_e": (1,"D"),

    "genie_lep_f_pdg": (1,"I"),
    "genie_lep_f_px": (1,"D"),
    "genie_lep_f_py": (1,"D"),
    "genie_lep_f_pz": (1,"D"),
    "genie_lep_f_e": (1,"D"),

    "genie_tgt_pdg": (1,"I"),
    "genie_tgt_px": (1,"D"),
    "genie_tgt_py": (1,"D"),
    "genie_tgt_pz": (1,"D"),
    "genie_tgt_e": (1,"D"),

    "genie_hnuc_pdg": (1,"I"),
    "genie_hnuc_px": (1,"D"),
    "genie_hnuc_py": (1,"D"),
    "genie_hnuc_pz": (1,"D"),
    "genie_hnuc_e": (1,"D"),

    "genie_scatter": (1,"I"),
    "genie_qel": (1,"I"),
    "genie_dis": (1,"I"),
    "genie_res": (1,"I"),
    "genie_mec": (1,"I"),

    "genie_x_bj": (1,"D"),
    "genie_y_inel": (1,"D"),
    "genie_Q2": (1,"D"),
    "genie_W": (1,"D"),

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

    #particle multiplicity
    "n_sim_prot": (1,"I"),
    "n_sim_neut": (1,"I"),
    "n_sim_pip": (1,"I"),
    "n_sim_pim": (1,"I"),

    #trigger info
    "trig_Ecal_nLayer": (1,"I"),
    "trig_Ecal_e_afterLayer": ("trig_Ecal_nLayer","D"),

    "trig_Hcal_nLayer": (1,"I"),
    "trig_Hcal_e_afterLayer": ("trig_Ecal_nLayer","D"),

    "trig_n_Electron": (1,"I"),
    "trig_maxE": (1,"I"),
    "trig_maxPt": (1,"I"),
    "trig_Electron_e": ("trig_n_Electron","D"),
    "trig_Electron_eClus": ("trig_n_Electron","D"),
    "trig_Electron_zClus": ("trig_n_Electron","D"),
    "trig_Electron_px": ("trig_n_Electron","D"),
    "trig_Electron_py": ("trig_n_Electron","D"),
    "trig_Electron_pz": ("trig_n_Electron","D"),
    "trig_Electron_dx": ("trig_n_Electron","D"),
    "trig_Electron_dy": ("trig_n_Electron","D"),
    "trig_Electron_x": ("trig_n_Electron","D"),
    "trig_Electron_y": ("trig_n_Electron","D"),
    "trig_Electron_tp": ("trig_n_Electron","D"),
    "trig_Electron_depth": ("trig_n_Electron","D"),

    #reweighting
    "n_weights": (1,"I"),
    "weights": ("n_weights","D")

}

FILES = glob.glob(arg.input)
#print(FILES)

#import concurrent.futures
#MAX_WORKERS=10

#with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
#        future_to_ifile = {executor.submit(process_file,ifile): ifile for ifile in FILES[:12]}
#        for future in concurrent.futures.as_completed(future_to_ifile):
#                ifile = future_to_ifile[future]
#                res = future.result()
for f in FILES:
    process_file(f)
                    
