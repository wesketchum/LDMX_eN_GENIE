#!/usr/bin/env python3

import os,sys
sys.path.append("./")
from analysis_common import *
from root_common import *

LEPTON_PT_TRIGGER = 0.4 #GeV

RTRACKER_ACCEPT_P = 0.05 #GeV, From Snowmass White Paper
RTRACKER_ACCEPT_ANGLE = radians(40)

#49cm is roughly the height (width is 51 CM), 24 cm from target to front face, from gdml v13
ECAL_MAXDET_ANGLE = ROOT.TMath.ATan2(49./2.,24.0) 

#3.1x3.1 m height/width at same position to target from front face, from gdml v13
HCAL_MAXDET_ANGLE = ROOT.TMath.ATan2(310./2.,24.0)

ECAL_ACCEPT_ANGLE = radians(40)
#HCAL_ACCEPT_ANGLE = radians(65)
HCAL_ACCEPT_ANGLE = radians(40)

ELECTRON_ACCEPT_KE = 0.06
CHPION_ACCEPT_KE = 0.06
PROTON_ACCEPT_KE = 0.06
NEUTRON_ACCEPT_KE = 0.06
PHOTON_ACCEPT_KE = 0.06

def create_gst_chain(files,verbose=False):
    gst_chain = ROOT.TChain("gst")
    for f in files:
        gst_chain.Add(f)
    if(verbose):
        print(f'Created gst chain from {len(files)} files with {gst_chain.GetEntries()} total events.')
    return gst_chain


leading_nucleon_ke_code='''
using namespace ROOT::VecOps;
using namespace ROOT::Math;
double leading_nucleon_ke(const RVec<double> &kefa_proton,const RVec<double> &kefa_neutron)
{

  double max_proton_ke = -999;
  double max_neutron_ke = -999;

  if(kefa_proton.size()!=0){
    max_proton_ke = Max(kefa_proton);
  }
  if(kefa_neutron.size()!=0){
    max_neutron_ke = Max(kefa_neutron);
  }

  if(max_proton_ke > max_neutron_ke) return max_proton_ke;
  else return max_neutron_ke;

  return -999;

}
'''
ROOT.gInterpreter.Declare(leading_nucleon_ke_code)

subleading_nucleon_ke_code='''
using namespace ROOT::VecOps;
using namespace ROOT::Math;
double subleading_nucleon_ke(const RVec<double> &kefa_proton,const RVec<double> &kefa_neutron)
{

  double max_proton_ke = -999;
  double max_neutron_ke = -999;
  double sub_max_proton_ke = -999;
  double sub_max_neutron_ke = -999;

  if(kefa_proton.size()>0){
    max_proton_ke = Max(kefa_proton);
    if(kefa_proton.size()>1){
       sub_max_proton_ke = Reverse(Sort(kefa_proton)).at(1);
    }
  }

  if(kefa_neutron.size()>0){
    max_neutron_ke = Max(kefa_neutron);
    if(kefa_neutron.size()>1){
       sub_max_neutron_ke = Reverse(Sort(kefa_neutron)).at(1);
    }
  }

  if(max_proton_ke > max_neutron_ke){
    if(sub_max_proton_ke > max_neutron_ke) return sub_max_proton_ke;
    else return max_neutron_ke;
  }
  if(max_neutron_ke > max_proton_ke){
    if(sub_max_neutron_ke > max_proton_ke) return sub_max_neutron_ke;
    else return max_proton_ke;
  }

  return -999;

}
'''
ROOT.gInterpreter.Declare(subleading_nucleon_ke_code)


getP4vec_code='''
using namespace ROOT::VecOps;
RVec<TLorentzVector> getP4vec(const RVec<double> &vpx,const RVec<double> &vpy,const RVec<double> &vpz,const RVec<double> &ve)
{
   auto p_4vec = [](const double &px,const double &py,const double &pz,const double &e)
   { return TLorentzVector(px,py,pz,e); };
   return Map(vpx,vpy,vpz,ve,p_4vec);
}
'''

ROOT.gInterpreter.Declare(getP4vec_code)

decayPi0_photon1_code='''
const double PI0_DECAY_E_PHOTON = 0.1395704*0.5;
using namespace ROOT::VecOps;
using namespace ROOT::Math;
RVec<TLorentzVector> decayPi0_photon1(const RVec<double> &vpx,const RVec<double> &vpy,const RVec<double> &vpz,const RVec<double> &ve)
{

   auto photon1_4vec = [](const double& px, const double& py, const double& pz, const double& e)
   {
      TRandom PI0_DECAY_RAND(ULong_t(e*1e14)+ULong_t(px*1e13));
   
      TLorentzVector pi0_4vec = TLorentzVector(px,py,pz,e);
      //photon1 cos_theta and phi
      double cos_theta = 2.*PI0_DECAY_RAND.Uniform(0,1)-1.;
      double sin_theta = sqrt(1.-cos_theta*cos_theta);
      double phi = 2*Pi()*PI0_DECAY_RAND.Uniform(0,1);
   
      auto photon1_4vec = TLorentzVector(PI0_DECAY_E_PHOTON*sin_theta*cos(phi),
                                         PI0_DECAY_E_PHOTON*sin_theta*sin(phi),
                                         PI0_DECAY_E_PHOTON*cos_theta,
                                         PI0_DECAY_E_PHOTON);
      photon1_4vec.Boost(pi0_4vec.BoostVector());
      return photon1_4vec;
   };
   return Map(vpx,vpy,vpz,ve,photon1_4vec);
}
'''
ROOT.gInterpreter.Declare(decayPi0_photon1_code)

decayPi0_photon2_code='''
using namespace ROOT::VecOps;
RVec<TLorentzVector> decayPi0_photon2(const RVec<double> &vpx,const RVec<double> &vpy,const RVec<double> &vpz,const RVec<double> &ve, const RVec<double> &ph1px, const RVec<double> &ph1py, const RVec<double> &ph1pz, const RVec<double> &ph1e)
{
   auto photon2_4vec = [](const double& px, const double& py, const double& pz, const double& e, const double& px1, const double& py1, const double& pz1, const double& e1)
   { return TLorentzVector(px,py,pz,e)-TLorentzVector(px1,py1,pz1,e1); };
   return Map(vpx,vpy,vpz,ve,ph1px,ph1py,ph1pz,ph1e,photon2_4vec);
}
'''
ROOT.gInterpreter.Declare(decayPi0_photon2_code)

getPx_code='''
using namespace ROOT::VecOps;
RVec<double> getPx(const RVec<TLorentzVector> &p4vec)
{
   auto px = [](const TLorentzVector &p4vec)
   { return p4vec.Px(); };
   return Map(p4vec,px);
}
'''
getPy_code='''
using namespace ROOT::VecOps;
RVec<double> getPy(const RVec<TLorentzVector> &p4vec)
{
   auto py = [](const TLorentzVector &p4vec)
   { return p4vec.Py(); };
   return Map(p4vec,py);
}
'''
getPz_code='''
using namespace ROOT::VecOps;
RVec<double> getPz(const RVec<TLorentzVector> &p4vec)
{
   auto pz = [](const TLorentzVector &p4vec)
   { return p4vec.Pz(); };
   return Map(p4vec,pz);
}
'''
getE_code='''
using namespace ROOT::VecOps;
RVec<double> getE(const RVec<TLorentzVector> &p4vec)
{
   auto e = [](const TLorentzVector &p4vec)
   { return p4vec.E(); };
   return Map(p4vec,e);
}
'''
ROOT.gInterpreter.Declare(getPx_code)
ROOT.gInterpreter.Declare(getPy_code)
ROOT.gInterpreter.Declare(getPz_code)
ROOT.gInterpreter.Declare(getE_code)



#define a function that makes new lepton variables
def define_df_gst_lep_vars(df_gst):
    df_gst = df_gst.Define("ptl","sqrt(pxl*pxl+pyl*pyl)") #pt lepton
    df_gst = df_gst.Define("thetazl","atan2(ptl,pzl)") #theta_z of lepton
    df_gst = df_gst.Define("energy_transfer","Ev-El")
    df_gst = df_gst.Define("pv","sqrt(pxv*pxv+pyv*pyv+pzv*pzv)")

    for var in ["px","py","pz"]:
        df_gst = df_gst.Define(f"{var}dl",f"{var}l-{var}v")
    df_gst = df_gst.Define(f"pdl",f"sqrt(pxdl*pxdl+pydl*pydl+pzdl*pzdl)")

    return df_gst

#define a function to make a bunch of hadron related variables
#sfx is the suffix on the vars (initial, final)
def define_df_gst_hadron_vars(df_gst,sfx=["i","f"]):
    
    for s in sfx:
        
        df_gst = df_gst.Define(f"thetaxz{s}",f"atan2(px{s},pz{s})")
        df_gst = df_gst.Define(f"thetayz{s}",f"atan2(py{s},pz{s})")
        df_gst = df_gst.Define(f"pt{s}",f"sqrt(px{s}*px{s}+py{s}*py{s})")
        df_gst = df_gst.Define(f"thetaz{s}",f"atan2(pt{s},pz{s})")
        
        #note, only for 'i' is total momentum not already in the gst tree
        if s=="i":
            df_gst = df_gst.Define(f"p{s}",f"sqrt(px{s}*px{s}+py{s}+py{s}+pz{s}*pz{s})")
            
        df_gst = df_gst.Define(f"mass{s}",f"sqrt(abs(E{s}*E{s}-p{s}*p{s}))")
        df_gst = df_gst.Define(f"ke{s}",f"E{s}-mass{s}")
        
        #for the hadrons, get the indices sorted by KE
        #df_gst = df_gst.Define(f"idx_ke{s}",f"Reverse(Argsort(ke{s}))")
    
    return df_gst


def define_df_gst_hadrons_by_pdg(df_gst,
                                 hvars=["E","p","px","py","pz",
                                        "pt","mass","ke",
                                        "thetaxz","thetayz","thetaz"],
                                 sfx=["f","i"]):
    
    particle_dict = { "proton":[2212,-2212],
                      "neutron":[2112,-2112],
                      "piplus":[211],
                      "piminus":[-211],
                      "pi0":[111],
                      "K0":[311,-311],
                      "Kplus":[321],
                      "Kminus":[-321] }
    
    for s in sfx:
        hvars_sfx = [ f'{hv}{s}' for hv in hvars ]
        
        for pname,pdgcodes in particle_dict.items():    
            pdgstr='('
            for pdgcode in pdgcodes:
                if len(pdgstr)==1:
                    pdgstr+=f'pdg{s}=={pdgcode}'
                else:
                    pdgstr+=f'||pdg{s}=={pdgcode}'
            pdgstr+=")"
            
            for hv in hvars_sfx:
                df_gst = df_gst.Define(f'{hv}_{pname}',f'{hv}[{pdgstr}]')        
                
    return df_gst

def define_df_gst_hadron_acceptance(df_gst,
                                    hvars=["E","p","px","py","pz",
                                        "pt","mass","ke",
                                        "thetaxz","thetayz","thetaz"],
                                    particles=["proton","neutron","piplus","piminus"],
                                    ke_min=[PROTON_ACCEPT_KE,NEUTRON_ACCEPT_KE,
                                            CHPION_ACCEPT_KE,CHPION_ACCEPT_KE],
                                    thetaz_max=[RTRACKER_ACCEPT_ANGLE,HCAL_ACCEPT_ANGLE,
                                                RTRACKER_ACCEPT_ANGLE,RTRACKER_ACCEPT_ANGLE],
                                    sfx=["f"]):

    for i in range(len(particles)):
        pname=particles[i]
        ke_min_val=ke_min[i]
        thetaz_max_val=thetaz_max[i]

        for s in sfx:
            df_gst = df_gst.Define(f"Accept{s}_{pname}",f"ke{s}_{pname}>{ke_min_val}&&thetaz{s}_{pname}<{thetaz_max_val}")
            for hv in hvars:
                df_gst = df_gst.Define(f"{hv}{s}a_{pname}",f"{hv}{s}_{pname}[Accept{s}_{pname}==1]")
                df_gst = df_gst.Define(f"{hv}{s}m_{pname}",f"{hv}{s}_{pname}[Accept{s}_{pname}!=1]")
        
        df_gst = df_gst.Define(f"nfa_{pname}",f"Sum(Acceptf_{pname})")
        df_gst = df_gst.Define(f"nfm_{pname}",f"(int)Acceptf_{pname}.size()-nfa_{pname}")
        
        

    return df_gst
    
def define_df_gst_photon_acceptance(df_gst,
                                    phvars=["E","px","py","pz","pt","thetaz"],
                                    ke_min=PHOTON_ACCEPT_KE,thetaz_max=ECAL_ACCEPT_ANGLE):

    for phname in ["ph1","ph2"]:
        df_gst = df_gst.Define(f"Acceptf_pi0_{phname}",
                               f"Ef_pi0_{phname}>{ke_min}&&thetazf_pi0_{phname}<{thetaz_max}")
        for var in phvars:
            df_gst = df_gst.Define(f"{var}fa_pi0_{phname}",f"{var}f_pi0_{phname}[Acceptf_pi0_{phname}==1]")
            df_gst = df_gst.Define(f"{var}fm_pi0_{phname}",f"{var}f_pi0_{phname}[Acceptf_pi0_{phname}!=1]")

    df_gst = df_gst.Define("Acceptf_pi0_phAll","Acceptf_pi0_ph1&&Acceptf_pi0_ph2")        
    df_gst = df_gst.Define("Acceptf_pi0_phNone","!Acceptf_pi0_ph1&&!Acceptf_pi0_ph2")

    df_gst = df_gst.Define(f"nfa_pi0_ph",f"Sum(Acceptf_pi0_ph1)+Sum(Acceptf_pi0_ph2)")
    df_gst = df_gst.Define(f"nfm_pi0_ph",f"(int)Acceptf_pi0_ph1.size()+(int)Acceptf_pi0_ph2.size()-nfa_pi0_ph")

    df_gst = df_gst.Define(f"nfa_pi0",f"Sum(Acceptf_pi0_phAll)")
    df_gst = df_gst.Define(f"nfm_pi0",f"(int)Acceptf_pi0_phAll.size()-nfa_pi0")
    
#    for hv in ["pxf","pyf","pzf","Ef"]:
#        df_gst = df_gst.Define(f"sum_{hv}_pi0_ph",f"{hv}_pi0_ph1+{hv}_pi0_ph2")

    
    return df_gst


def define_df_gst_pi0decay(df_gst):
    df_gst = df_gst.Define("p4vec_pi0_ph1","decayPi0_photon1(pxf_pi0,pyf_pi0,pzf_pi0,Ef_pi0)")
    df_gst = df_gst.Define("pxf_pi0_ph1","getPx(p4vec_pi0_ph1)")
    df_gst = df_gst.Define("pyf_pi0_ph1","getPy(p4vec_pi0_ph1)")
    df_gst = df_gst.Define("pzf_pi0_ph1","getPz(p4vec_pi0_ph1)")
    df_gst = df_gst.Define("Ef_pi0_ph1","getE(p4vec_pi0_ph1)")
    df_gst = df_gst.Define("p4vec_pi0_ph2",
                           "decayPi0_photon2(pxf_pi0,pyf_pi0,pzf_pi0,Ef_pi0,pxf_pi0_ph1,pyf_pi0_ph1,pzf_pi0_ph1,Ef_pi0_ph1)")
    df_gst = df_gst.Define("pxf_pi0_ph2","getPx(p4vec_pi0_ph2)")
    df_gst = df_gst.Define("pyf_pi0_ph2","getPy(p4vec_pi0_ph2)")
    df_gst = df_gst.Define("pzf_pi0_ph2","getPz(p4vec_pi0_ph2)")
    df_gst = df_gst.Define("Ef_pi0_ph2","getE(p4vec_pi0_ph2)")
    
    df_gst = df_gst.Define("ptf_pi0_ph1","sqrt(pxf_pi0_ph1*pxf_pi0_ph1+pyf_pi0_ph1*pyf_pi0_ph1)")
    df_gst = df_gst.Define("ptf_pi0_ph2","sqrt(pxf_pi0_ph2*pxf_pi0_ph2+pyf_pi0_ph2*pyf_pi0_ph2)")
    
    df_gst = df_gst.Define("thetazf_pi0_ph1","atan2(ptf_pi0_ph1,pzf_pi0_ph1)")
    df_gst = df_gst.Define("thetazf_pi0_ph2","atan2(ptf_pi0_ph2,pzf_pi0_ph2)")
        
    return df_gst


def define_df_gst_hadron_sums(df_gst,
                              particles=["proton","neutron","piplus","piminus"],
                              sfx=["i","f"],do_pi0=True):
    
    for pname in particles:
        for s in sfx:
            for hv in ["px","py","pz","E","ke"]:
                df_gst = df_gst.Define(f"sum_{hv}{s}_{pname}",f"Sum({hv}{s}_{pname})")
            
            df_gst = df_gst.Define(f"sum_pt{s}_{pname}",
                                   f"sqrt(sum_px{s}_{pname}*sum_px{s}_{pname}+sum_py{s}_{pname}*sum_py{s}_{pname})")
            df_gst = df_gst.Define(f"sum_p{s}_{pname}",
                                   f"sqrt(sum_pt{s}_{pname}*sum_pt{s}_{pname}+sum_pz{s}_{pname}*sum_pz{s}_{pname})")
            
    if do_pi0:
        for s in sfx:
            for hv in ["px","py","pz","E"]:
                if s=="i":
                    df_gst = df_gst.Define(f"sum_{hv}{s}_pi0",f"Sum({hv}{s}_{pname})")
                else:
                    df_gst = df_gst.Define(f"sum_{hv}{s}_pi0_ph",f"Sum({hv}{s}_pi0_ph1)+Sum({hv}{s}_pi0_ph2)")
                    
            if s=="i":
                df_gst = df_gst.Define(f"sum_pt{s}_pi0",
                                       f"sqrt(sum_px{s}_pi0*sum_px{s}_pi0+sum_py{s}_pi0*sum_py{s}_pi0)")
                df_gst = df_gst.Define(f"sum_p{s}_pi0",
                                       f"sqrt(sum_pt{s}_pi0*sum_pt{s}_pi0+sum_pz{s}_pi0*sum_pz{s}_pi0)")
            else:
                df_gst = df_gst.Define(f"sum_pt{s}_pi0_ph",
                                       f"sqrt(sum_px{s}_pi0_ph*sum_px{s}_pi0_ph+sum_py{s}_pi0_ph*sum_py{s}_pi0_ph)")
                df_gst = df_gst.Define(f"sum_p{s}_pi0_ph",
                                       f"sqrt(sum_pt{s}_pi0_ph*sum_pt{s}_pi0_ph+sum_pz{s}_pi0_ph*sum_pz{s}_pi0_ph)")
                

                    
    for s in sfx:
        if s!="i" and s!="f": continue
        
        for hv in ["px","py","pz","E","ke"]:
            df_gst = df_gst.Define(f"sum_{hv}{s}",f"Sum({hv}{s})")
            
        df_gst = df_gst.Define(f"sum_pt{s}",f"sqrt(sum_px{s}*sum_px{s}+sum_py{s}*sum_py{s})")
        df_gst = df_gst.Define(f"sum_p{s}",f"sqrt(sum_pt{s}*sum_pt{s}+sum_pz{s}*sum_pz{s})")
        
    
    return df_gst


def define_df_gst_momentum_imbalance(df_gst,suffix_list,cname):
    
    for var in ["px","py","pz"]:
        df_gst = df_gst.Define(f"hsum_{var}{cname}","+".join(f"sum_{var}{s}" for s in suffix_list))
        df_gst = df_gst.Define(f"delta_{var}{cname}",f"{var}dl+hsum_{var}{cname}")
    
    df_gst = df_gst.Define(f"hsum_pt{cname}",f"sqrt(hsum_px{cname}*hsum_px{cname}+hsum_py{cname}*hsum_py{cname})")
    df_gst = df_gst.Define(f"hsum_p{cname}",f"sqrt(hsum_pt{cname}*hsum_pt{cname}+hsum_pz{cname}*hsum_pz{cname})")
        
    df_gst = df_gst.Define(f"delta_pt{cname}",f"sqrt(delta_px{cname}*delta_px{cname}+delta_py{cname}*delta_py{cname})")
    df_gst = df_gst.Define(f"delta_cosalphat{cname}",f"-1*(pxl*delta_px{cname}+pyl*delta_py{cname})/(ptl*delta_pt{cname})")
    df_gst = df_gst.Define(f"delta_cosphit{cname}",f"-1*(pxl*hsum_px{cname}+pyl*hsum_py{cname})/(ptl*hsum_pt{cname})")
    df_gst = df_gst.Define(f"delta_alphat{cname}",f"acos(delta_cosalphat{cname})")
    df_gst = df_gst.Define(f"delta_phit{cname}",f"acos(delta_cosphit{cname})")
    
    df_gst = df_gst.Define(f"delta_p{cname}",
                        f"sqrt(delta_px{cname}*delta_px{cname}+delta_py{cname}*delta_py{cname}+delta_pz{cname}*delta_pz{cname})")
    df_gst = df_gst.Define(f"delta_cosalpha{cname}",
                        f"-1*(pxdl*delta_px{cname}+pydl*delta_py{cname}+pzdl*delta_pz{cname})/(pdl*delta_p{cname})")
    df_gst = df_gst.Define(f"delta_cosphi{cname}",
                        f"-1*(pxdl*hsum_px{cname}+pydl*hsum_py{cname}+pzdl*hsum_pz{cname})/(pdl*hsum_p{cname})")
    df_gst = df_gst.Define(f"delta_alpha{cname}",f"acos(delta_cosalpha{cname})")
    df_gst = df_gst.Define(f"delta_phi{cname}",f"acos(delta_cosphi{cname})")
    
    return df_gst


#@ROOT.Numba.Declare(["float"], "float")
#def pt_res(pt):
#    return pt*0.1
