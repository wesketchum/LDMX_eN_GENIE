from ROOT import *
import numpy as np
import random

f_4GeV = TFile("output_file_29Aug_4GeV_400k.root")
ana_tree = f_4GeV.Get("ana_tree")

#calibration factors
cf_lw = 1.60
cf_byslices = [50,21.85,10.92,8.95,8.25,8.05]

# pT selection
pT_4GeV = 400.
pT_8GeV = 600.

# event selection
ecal_min = 134.97 - 30.
ecal_max = 134.97 + 30.
hcal_min = 134.97 - 50.
hcal_max = 134.97 + 50.

# Eventrate Table
er_table_se_4GeV = np.zeros((9,3),dtype=double)

# confusion matrix
Npi0_X = np.zeros((4,4),dtype=double)
Np_X = np.zeros((4,4),dtype=double)
Npip_X = np.zeros((4,4),dtype=double)
Npim_X = np.zeros((4,4),dtype=double)
Npipm_X = np.zeros((4,4),dtype=double)

#print(er_table)

# function definitions
def cos_theta(px, py, pz, ph1, ph2,npi0):
    return (px[ph1[npi0]]*px[ph2[npi0]]+py[ph1[npi0]]*py[ph2[npi0]]+pz[ph1[npi0]]*pz[ph2[npi0]])/(np.sqrt(px[ph1[npi0]]*px[ph1[npi0]]+py[ph1[npi0]]*py[ph1[npi0]]+pz[ph1[npi0]]*pz[ph1[npi0]])*np.sqrt(px[ph2[npi0]]*px[ph2[npi0]]+py[ph2[npi0]]*py[ph2[npi0]]+pz[ph2[npi0]]*pz[ph2[npi0]]))

def ecal_energy(e,ph,npi0):
    return cf_lw*e[ph[npi0]]

def hcal_energy(e,ph,npi0):
    if(e[ph[npi0]]<1.5):
        return cf_byslices[0]*e[ph[npi0]]
    if(e[ph[npi0]]>=1.5 and e[ph[npi0]]<10):
        return cf_byslices[1]*e[ph[npi0]]
    if(e[ph[npi0]]>=10 and e[ph[npi0]]<20):
        return cf_byslices[2]*e[ph[npi0]]
    if(e[ph[npi0]]>=20 and e[ph[npi0]]<30):
        return cf_byslices[3]*e[ph[npi0]]
    if(e[ph[npi0]]>=30 and e[ph[npi0]]<40):
        return cf_byslices[4]*e[ph[npi0]]
    if(e[ph[npi0]]>=40):
        return cf_byslices[5]*e[ph[npi0]]

def invarmass(e1,e2,px,py,pz,ph1,ph2,npi0):
    return np.sqrt(2.*e1*e2*(1-cos_theta(px,py,pz,ph1,ph2,npi0)))
    
npi0_1_counter = 0
npi0_2_counter = 0
npi0_3_counter = 0
np_1_counter = 0
np_2_counter = 0
np_3_counter = 0
npip_1_counter = 0
npip_2_counter = 0
npip_3_counter = 0
npim_1_counter = 0
npim_2_counter = 0
npim_3_counter = 0

for entryNum in range(0,ana_tree.GetEntries()):
    ana_tree.GetEntry(entryNum)
    # electron
    e_pt = getattr(ana_tree,"elec_pt")
   
    n_sim_p = getattr(ana_tree,"n_sim_p")
    p_pdg = getattr(ana_tree,"sim_p_pdg")
    p_e_true = getattr(ana_tree,"sim_p_e")
    p_p = getattr(ana_tree,"sim_p_p")
    p_e_lw = getattr(ana_tree,"sim_p_ecal_e_lw")
    p_h = getattr(ana_tree,"sim_p_hcal_e")
    p_px = getattr(ana_tree,"sim_p_px")
    p_py = getattr(ana_tree,"sim_p_py")
    p_pz = getattr(ana_tree,"sim_p_pz")
    p_thetaz = getattr(ana_tree,"sim_p_thetaz")
    n_sim_pi0 = getattr(ana_tree,"n_sim_pi0")
    pi0_idx = getattr(ana_tree,"pi0_idx")
    ph1_idx = getattr(ana_tree,"pi0_photon1_idx")
    ph2_idx = getattr(ana_tree,"pi0_photon2_idx")
    ph1_det = getattr(ana_tree,"pi0_photon1_det")
    ph2_det = getattr(ana_tree,"pi0_photon2_det")
    hcal_top_e = getattr(ana_tree,"hcal_top_e")
    hcal_bottom_e = getattr(ana_tree,"hcal_bottom_e")
    hcal_left_e = getattr(ana_tree,"hcal_left_e")
    hcal_right_e = getattr(ana_tree,"hcal_right_e")
    hcal_back_e = getattr(ana_tree,"hcal_back_e")

    n_event = 0
    n_p = 0
    n_n = 0
    n_pi0 = 0
    n_pip = 0
    n_pim = 0
    n_p_thresh = 0
    n_n_thresh = 0
    n_pi0_thresh = 0
    n_pi0_thresh_ee = 0
    n_pi0_thresh_eh = 0
    n_pi0_thresh_he = 0
    n_pi0_thresh_hh = 0
    n_pip_thresh = 0
    n_pim_thresh = 0
    n_Kp = 0
    n_Kp_thresh = 0
    n_Km = 0
    n_Km_thresh = 0
    n_K0s = 0
    n_K0s_thresh = 0
    n_K0l = 0
    n_K0l_thresh = 0

    # Particle multiplicity with thresholds
    # Trigger Events requires e_pT >400MeV @ 4GeV/ >500MeV @ 8GeV
    if(e_pt>pT_4GeV):
        #print(e_pt)
        n_event += 1
        for particle in range(n_sim_p):
            if(p_pdg[particle] == 2212):
                n_p += 1
                if(p_p[particle]>=800 and p_thetaz[particle]<=80*3.1415/180):
                    n_p_thresh += 1
            if(p_pdg[particle] == 2112):
                n_n += 1
                if(p_e_true[particle]>=1000 and p_thetaz[particle]<=80*3.1415/180):
                    n_n_thresh += 1
            if(p_pdg[particle] == 111):
                n_pi0 += 1
            if(p_pdg[particle] == 211):
                n_pip += 1
                if(p_p[particle]>=100 and p_thetaz[particle]<=80*3.1415/180):
                    n_pip_thresh += 1
            if(p_pdg[particle] == -211):
                n_pim += 1
                if(p_p[particle]>=100 and p_thetaz[particle]<=80*3.1415/180):
                    n_pim_thresh += 1
            if(p_pdg[particle] == 321):
                n_Kp += 1
                if(p_p[particle] >= 800 and p_thetaz[particle]<=80*3.1415/180):
                    n_Kp_thresh += 1
            if(p_pdg[particle] == -321):
                n_Km += 1
                if(p_p[particle] >= 800 and p_thetaz[particle]<=80*3.1415/180):
                    n_Km_thresh += 1
            if(p_pdg[particle] == 310):
                n_K0s += 1
                if(p_p[particle] >= 800 and p_thetaz[particle]<=80*3.1415/180):
                    n_K0s_thresh += 1
            if(p_pdg[particle] == 130):
                n_K0l += 1
                if(p_p[particle] >= 800 and p_thetaz[particle]<=80*3.1415/180):
                    n_K0l_thresh += 1
    
        # pi0 threshold
        if(n_pi0>=1):
            for x in range(n_sim_pi0):
                # ee
                if(ph1_det[x] == 1 and ph2_det[x]== 1 and len(p_thetaz)>ph1_idx[x] and len(p_thetaz)>ph2_idx[x] and invarmass(ecal_energy(p_e_lw,ph1_idx,x),ecal_energy(p_e_lw,ph2_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) >= ecal_min and invarmass(ecal_energy(p_e_lw,ph1_idx,x),ecal_energy(p_e_lw,ph2_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) <= ecal_max):
                    n_pi0_thresh_ee += 1
                # eh
                if(ph1_det[x] == 1 and ph2_det[x] == 2 and len(p_thetaz)>ph1_idx[x] and len(p_thetaz)>ph2_idx[x] and invarmass(ecal_energy(p_e_lw,ph1_idx,x),hcal_energy(p_h,ph2_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) >= hcal_min and invarmass(ecal_energy(p_e_lw,ph1_idx,x),hcal_energy(p_h,ph2_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) <= hcal_max):
                    n_pi0_thresh_eh += 1
                # he
                if(ph1_det[x] == 2 and ph2_det[x] == 1 and len(p_thetaz)>ph1_idx[x] and len(p_thetaz)>ph2_idx[x] and invarmass(ecal_energy(p_e_lw,ph2_idx,x),hcal_energy(p_h,ph1_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) >= hcal_min and invarmass(ecal_energy(p_e_lw,ph2_idx,x),hcal_energy(p_h,ph1_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) <= hcal_max):
                    n_pi0_thresh_he += 1
                # hh
                if(ph1_det[x] == 2 and ph2_det[x] == 2 and len(p_thetaz)>ph1_idx[x] and len(p_thetaz)>ph2_idx[x] and invarmass(hcal_energy(p_h,ph1_idx,x),hcal_energy(p_h,ph2_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) >= hcal_min and invarmass(hcal_energy(p_h,ph1_idx,x),hcal_energy(p_h,ph2_idx,x),p_px,p_py,p_pz,ph1_idx,ph2_idx,x) <= hcal_max):
                    n_pi0_thresh_hh += 1

        n_pi0_thresh = n_pi0_thresh_ee + n_pi0_thresh_eh + n_pi0_thresh_he + n_pi0_thresh_hh

        # Eventrate Table
        # Trigger Events |  Reconstructable Events | Reconstructable Acceptance
        # semi-exclusive topologies
        
        # Npi0 + X
        if(n_pi0>=1):
            er_table_se_4GeV[0][0] += 1
            if(n_pi0_thresh>=1):
                er_table_se_4GeV[0][1] += 1
            er_table_se_4GeV[0][2] = round(100*er_table_se_4GeV[0][1]/er_table_se_4GeV[0][0],1)
        # Npi p/m + X
        if(n_pip>=1 or n_pim>=1):
            er_table_se_4GeV[1][0] += 1
            if(n_pip_thresh>=1 or n_pim_thresh>=1):
                er_table_se_4GeV[1][1] += 1
            er_table_se_4GeV[1][2] = round(100*er_table_se_4GeV[1][1]/er_table_se_4GeV[1][0],1)
        # Npi0 + Np + X
        if(n_pi0>=1 and n_p>=1):
            er_table_se_4GeV[2][0] += 1
            if(n_pi0_thresh>=1 and n_p_thresh>=1):
                er_table_se_4GeV[2][1] += 1
            er_table_se_4GeV[2][2] = round(100*er_table_se_4GeV[2][1]/er_table_se_4GeV[2][0],1)
        # Npi p/m + Np + X
        if((n_pip>=1 or n_pim>=1) and n_p>=1):
            er_table_se_4GeV[3][0] += 1
            if((n_pip_thresh>=1 or n_pim_thresh>=1) and n_p>=1):
                er_table_se_4GeV[3][1] += 1
            er_table_se_4GeV[3][2] = round(100*er_table_se_4GeV[3][1]/er_table_se_4GeV[3][0],1)
        # Npi0 + Nn + X
        if(n_pi0>=1 and n_n>=1):
            er_table_se_4GeV[4][0] += 1
            if(n_pi0_thresh>=1 and n_n_thresh>=1):
                er_table_se_4GeV[4][1] += 1
            er_table_se_4GeV[4][2] = round(100*er_table_se_4GeV[4][1]/er_table_se_4GeV[4][0],1)
        # Npi p/m + Nn + X
        if((n_pip>=1 or n_pim>=1) and n_n>=1):
            er_table_se_4GeV[5][0] += 1
            if((n_pip_thresh>=1 or n_pim_thresh>=1) and n_n>=1):
                er_table_se_4GeV[5][1] += 1
            er_table_se_4GeV[5][2] = round(100*er_table_se_4GeV[5][1]/er_table_se_4GeV[5][0],1)
        # NK p/m + X
        if(n_Kp>=1 or n_Km>=1):
            er_table_se_4GeV[6][0] += 1
            if(n_Kp_thresh>=1 or n_Km_thresh>=1):
                er_table_se_4GeV[6][1] += 1
            er_table_se_4GeV[6][2] = round(100*er_table_se_4GeV[6][1]/er_table_se_4GeV[6][0],1)
        # NK0s + X
        if(n_K0s>=1):
            er_table_se_4GeV[7][0] += 1
            if(n_K0s>=1):
                er_table_se_4GeV[7][1] += 1
            er_table_se_4GeV[7][2] = round(100*er_table_se_4GeV[7][1]/er_table_se_4GeV[7][0],1)
        # NK0l + X
        if(n_K0l>=1):
            er_table_se_4GeV[8][0] += 1
            if(n_K0l>=1):
                er_table_se_4GeV[8][1] += 1
            er_table_se_4GeV[8][2] = round(100*er_table_se_4GeV[8][1]/er_table_se_4GeV[8][0],1)

        # confusion matrix Npi0 + X
        if(n_pi0==0):
            Npi0_X[3][0] += 1
        if(n_pi0==1):
            npi0_1_counter += 1
            if(n_pi0_thresh==1):
                Npi0_X[2][1] += 1
            Npi0_X[2][0] = npi0_1_counter - Npi0_X[2][1]
        if(n_pi0==2):
            npi0_2_counter += 1
            if(n_pi0_thresh==2):
                Npi0_X[1][2] += 1
            if(n_pi0_thresh==1):
                Npi0_X[1][1] += 1
            Npi0_X[1][0] = npi0_2_counter - Npi0_X[1][1] - Npi0_X[1][2]
        if(n_pi0>2):
            npi0_3_counter += 1
            if(n_pi0_thresh>2):
                Npi0_X[0][3] += 1
            if(n_pi0_thresh==2):
                Npi0_X[0][2] += 1
            if(n_pi0_thresh==1):
                Npi0_X[0][1] += 1
            Npi0_X[0][0] = npi0_3_counter - Npi0_X[0][1] - Npi0_X[0][2] - Npi0_X[0][3]
    
        # confusion matrix Np + X
        if(n_p==0):
            Np_X[3][0] += 1
        if(n_p==1):
            np_1_counter += 1
            if(n_p_thresh==1):
                Np_X[2][1] += 1
            Np_X[2][0] = np_1_counter - Np_X[2][1]
        if(n_p==2):
            np_2_counter += 1
            if(n_p_thresh==2):
                Np_X[1][2] += 1
            if(n_p_thresh==1):
                Np_X[1][1] += 1
            Np_X[1][0] = np_2_counter - Np_X[1][1] - Np_X[1][2]
        if(n_p>2):
            np_3_counter += 1
            if(n_p_thresh>2):
                Np_X[0][3] += 1
            if(n_p_thresh==2):
                Np_X[0][2] += 1
            if(n_p_thresh==1):
                Np_X[0][1] += 1
            Np_X[0][0] = np_3_counter - Np_X[0][1] - Np_X[0][2] - Np_X[0][3]
    
        # confusion matrix Npip + X
        # to be added with Npim + X
        if(n_pip==0):
            Npip_X[3][0] += 1
        if(n_pip==1):
            npip_1_counter += 1
            if(n_pip_thresh==1):
                Npip_X[2][1] += 1
            Npip_X[2][0] = npip_1_counter - Npip_X[2][1]
        if(n_pip==2):
            npip_2_counter += 1
            if(n_pip_thresh==2):
                Npip_X[1][2] += 1
            if(n_pip_thresh==1):
                Npip_X[1][1] += 1
            Npip_X[1][0] = npip_2_counter - Npip_X[1][1] - Npip_X[1][2]
        if(n_pip>2):
            npip_3_counter += 1
            if(n_pip_thresh>2):
                Npip_X[0][3] += 1
            if(n_pip_thresh==2):
                Npip_X[0][2] += 1
            if(n_pip_thresh==1):
                Npip_X[0][1] += 1
            Npip_X[0][0] = npip_3_counter - Npip_X[0][1] - Npip_X[0][2] - Npip_X[0][3]
    
        # confusion matrix Npim + X
        if(n_pim==0):
            Npim_X[3][0] += 1
        if(n_pim==1):
            npim_1_counter += 1
            if(n_pim_thresh==1):
                Npim_X[2][1] += 1
            Npim_X[2][0] = npim_1_counter - Npim_X[2][1]
        if(n_pim==2):
            npim_2_counter += 1
            if(n_pim_thresh==2):
                Npim_X[1][2] += 1
            if(n_pim_thresh==1):
                Npim_X[1][1] += 1
            Npim_X[1][0] = npim_2_counter - Npim_X[1][1] - Npim_X[1][2]
        if(n_pim>2):
            npim_3_counter += 1
            if(n_pim_thresh>2):
                Npim_X[0][3] += 1
            if(n_pim_thresh==2):
                Npim_X[0][2] += 1
            if(n_pim_thresh==1):
                Npim_X[0][1] += 1
            Npim_X[0][0] = npim_3_counter - Npim_X[0][1] - Npim_X[0][2] - Npim_X[0][3]

        Npipm_X = np.add(Npip_X,Npim_X)

#print(er_table_se_4GeV)

# LaTeX printout
# ensure \usepackage{multirow}

latex_table_4GeV = "\\begin{center}\n"
latex_table_4GeV += "\\begin{tabular}{c | c c c}\n"
latex_table_4GeV += " & Trigger Events & Reconstructable Events & Reconstructable Efficiency\\\ \n"
latex_table_4GeV += "\\hline\n"
latex_table_4GeV += "$N\pi^{0} + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[0][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$N\pi^{\pm} + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[1][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$N\pi^{0} + Np + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[2][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$N\pi^{\pm} + Np + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[3][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$N\pi^{0} + Nn + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[4][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$N\pi^{\pm} + Nn + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[5][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$NK^{\pm} + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[6][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$NK^{0}_{s} + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[7][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$NK^{0}_{l} + X$"
for num in range(0,3):
    latex_table_4GeV += " & " + str(er_table_se_4GeV[8][num])
    if(num==2):
        latex_table_4GeV += "\%"
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "\\end{tabular}\n"
latex_table_4GeV += "\\end{center}\n"

latex_table_4GeV += "\\begin{center}\n"
latex_table_4GeV += "\\begin{tabular}{c | r r r r}\n"
latex_table_4GeV += "\multicolumn{1}{}{} & \multicolumn{4}{c}{Reco}\\\ \n"
latex_table_4GeV += " & $0\pi^{0}+X$ & $1\pi^{0}+X$ & $2\pi^{0}+X$ & $\geq3\pi^{0}+X$\\\ \n"
latex_table_4GeV += "\\hline\n"
latex_table_4GeV += "\\rotatebox{90}{\multirow{4}{*}{\\vspace{35pt}\hspace{-30pt}Truth}}"
latex_table_4GeV += "$\geq3\pi^{0}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npi0_X[0][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$2\pi^{0}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npi0_X[1][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$1\pi^{0}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npi0_X[2][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$0\pi^{0}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npi0_X[3][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "\\end{tabular}\n"
latex_table_4GeV += "\\end{center}\n"

latex_table_4GeV += "\\begin{center}\n"
latex_table_4GeV += "\\begin{tabular}{c | r r r r}\n"
latex_table_4GeV += "\multicolumn{1}{}{} & \multicolumn{4}{c}{Reco}\\\ \n"
latex_table_4GeV += " & $0p+X$ & $1p+X$ & $2p+X$ & $\geq3p+X$\\\ \n"
latex_table_4GeV += "\\hline\n"
latex_table_4GeV += "\\rotatebox{90}{\multirow{4}{*}{\\vspace{35pt}\hspace{-30pt}Truth}}"
latex_table_4GeV += "$\geq3p+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Np_X[0][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$2p+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Np_X[1][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$1p+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Np_X[2][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$0p+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Np_X[3][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "\\end{tabular}\n"
latex_table_4GeV += "\\end{center}\n"

latex_table_4GeV += "\\begin{center}\n"
latex_table_4GeV += "\\begin{tabular}{c | r r r r}\n"
latex_table_4GeV += "\multicolumn{1}{}{} & \multicolumn{4}{c}{Reco}\\\ \n"
latex_table_4GeV += " & $0\pi^{\pm}+X$ & $1\pi^{\pm}+X$ & $2\pi^{\pm}+X$ & $\geq3\pi^{\pm}+X$\\\ \n"
latex_table_4GeV += "\\hline\n"
latex_table_4GeV += "\\rotatebox{90}{\multirow{4}{*}{\\vspace{35pt}\hspace{-30pt}Truth}}"
latex_table_4GeV += "$\geq3\pi^{\pm}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npipm_X[0][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$2\pi^{\pm}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npipm_X[1][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$1\pi^{\pm}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npipm_X[2][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "$0\pi^{\pm}+X$"
for reco in range(0,4):
    latex_table_4GeV += " & " + str(Npipm_X[3][reco])
latex_table_4GeV += "\\\ \n"
latex_table_4GeV += "\\end{tabular}\n"
latex_table_4GeV += "\\end{center}\n"

#print(latex_table_4GeV)
with open('er_table_4GeV_se_with_cm.txt','w') as f:
    f.write(latex_table_4GeV)
