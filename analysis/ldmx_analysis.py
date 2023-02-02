from LDMX.Framework import EventTree
import numpy as np
import ROOT
#from array import array

ELECTRON_PT_CUT = 0.
HCAL_PE_CUT = 0.
ECAL_PE_CUT = 0.
FILENAME = 'test_Ti.root'

EVENTS_TO_PROCESS = []
#EVENTS_TO_PROCESS = [84,66]

VERBOSE = False


def pt(particle):
    return np.sqrt(particle.getMomentum()[0]*particle.getMomentum()[0]+
                       particle.getMomentum()[1]*particle.getMomentum()[1])
def p(particle):
    return np.sqrt(particle.getMomentum()[0]*particle.getMomentum()[0]+
                       particle.getMomentum()[1]*particle.getMomentum()[1]+
                       particle.getMomentum()[2]*particle.getMomentum()[2])

def thetaz(particle):
    return np.arccos(particle.getMomentum()[2]/p(particle))

input_tree = EventTree.EventTree(FILENAME)

variables = { "elec_px" : np.zeros(1,dtype=float),
              "elec_py" : np.zeros(1,dtype=float),
              "elec_pz" : np.zeros(1,dtype=float),
              "elec_pt" : np.zeros(1,dtype=float),
              "elec_p" : np.zeros(1,dtype=float),
              "hcal_e": np.zeros(1,dtype=float),
              "hcal_et": np.zeros(1,dtype=float),
              "ecal_e": np.zeros(1,dtype=float),
              "ecal_et": np.zeros(1,dtype=float),
              "cal_e": np.zeros(1,dtype=float),
              "cal_et": np.zeros(1,dtype=float)
              
    }

output_file = ROOT.TFile("output_file.root","RECREATE")

output_tree = ROOT.TTree("ana_tree","eN Analysis Tree")
output_tree.Branch("elec_px",variables["elec_px"],'elec_px/D')
output_tree.Branch("elec_py",variables["elec_py"],'elec_py/D')
output_tree.Branch("elec_pz",variables["elec_pz"],'elec_pz/D')
output_tree.Branch("elec_pt",variables["elec_pt"],'elec_pt/D')
output_tree.Branch("elec_p",variables["elec_p"],'elec_p/D')

output_tree.Branch("hcal_e",variables["hcal_e"],'hcal_e/D')
output_tree.Branch("hcal_et",variables["hcal_et"],'hcal_et/D')
output_tree.Branch("ecal_e",variables["ecal_e"],'ecal_e/D')
output_tree.Branch("ecal_et",variables["ecal_et"],'ecal_et/D')
output_tree.Branch("cal_e",variables["cal_e"],'cal_e/D')
output_tree.Branch("cal_et",variables["cal_et"],'cal_et/D')
    
for ie, event in enumerate(input_tree):

    if(len(EVENTS_TO_PROCESS)>0 and (event.EventHeader.getEventNumber() not in EVENTS_TO_PROCESS)): continue
    
    if(VERBOSE): event.EventHeader.Print()

    e_pt = 0
    for i_p, particle in event.SimParticles_genie:
        if(particle.getPdgID()==11 and particle.getGenStatus()==1):
            e_pt = pt(particle)
            if(VERBOSE): print(f'\tElectron pt is {e_pt}')

            variables["elec_px"][0] = particle.getMomentum()[0]
            variables["elec_py"][0] = particle.getMomentum()[1]
            variables["elec_pz"][0] = particle.getMomentum()[2]
            variables["elec_pt"][0] = pt(particle)
            variables["elec_p"][0] = p(particle)

    if e_pt < ELECTRON_PT_CUT:
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

    variables["hcal_e"][0] = total_hcal_e
    variables["hcal_et"][0] = hcal_et

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

        
    if(e_pt > ecal_et and VERBOSE):
        print("List of particles:")
        for i_p, particle in event.SimParticles_genie:
            
            print(f'ID {i_p}: PDG {particle.getPdgID()}, Status {particle.getGenStatus()}, Energy {particle.getEnergy()}, Momentum {particle.getMomentum()}, thetaz {np.degrees(thetaz(particle))}, Parents {particle.getParents()}')

    output_tree.Fill()

output_file.Write()
output_file.Close()

