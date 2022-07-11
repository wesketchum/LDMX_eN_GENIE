# LDMX_eN_GENIE
eN studies for LDMX

### Examples

See https://github.com/wesketchum/LDMX_eN_GENIE/blob/main/notebooks/examples/Analyze_GST_RDataFrame.ipynb for an example notebook that uses the basic functions.

### Useful documentation
- [GENIE Event Generator](http://www.genie-mc.org/)
  - [Model tunes list](https://hep.ph.liv.ac.uk/~costasa/genie/tunes.html)
  - [User Manual](https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=2)

### Common helper functions

`create_gst_chain(files,verbose=False)`  
_Returns a TChain from a list of GENIE gst files._

`define_df_gst_lep_vars(df_gst)`  
_Returns an RDataframe with leptonic kinematic variables added/calculated._

`define_df_gst_hadron_vars(df_gst,sfx=["i","f"])`  
_Returns an RDataFrame with hadronic kinematic variables added/calculated. Does this for the given suffixes in gst files ('i' for initial state hadrons, 'f' for final state hadrons)._

`define_df_gst_hadrons_by_pdg(df_gst,hvars=["E","p","px","py","pz","pt","mass","ke","thetaxz","thetayz","thetaz"],sfx=["f","i"])`  
_Returns an RDataFrame with hadronic variables broken down by particle type. `hvars` is the list of hadronic variables to do, and `sfx` is the list of suffixes to consider._

`define_df_gst_pi0decay(df_gst)`  
_Returns an RDataFrame that performs a simulation of_ $\pi^{0}\rightarrow\gamma\gamma$ _decays._

`define_df_gst_hadron_acceptance(df_gst,  
                                 hvars=["E","p","px","py","pz",  
                                        "pt","mass","ke",  
                                        "thetaxz","thetayz","thetaz"],  
                                 particles=["proton","neutron","piplus","piminus"],  
                                 ke_min=[PROTON_ACCEPT_KE,NEUTRON_ACCEPT_KE,  
                                         CHPION_ACCEPT_KE,CHPION_ACCEPT_KE],  
                                 thetaz_max=[RTRACKER_ACCEPT_ANGLE,HCAL_ACCEPT_ANGLE,  
                                             RTRACKER_ACCEPT_ANGLE,RTRACKER_ACCEPT_ANGLE],  
                                 sfx=["f"]):`  
_Returns an RDataFrame with hadronic variables selected for detector acceptance. `hvars` is the list of hadronic variables to create, `particles` is the list of particles to consider, `ke_min` is a list of kinetic energy thresholds to apply for each listed particle (in GeV), `thetaz_max` is a list of angular thresholds to apply for each listed particle (radians), and `sfx` is the list of suffixes to consider (generally, just final state particles or `f`). New variables will have suffix `fa` for those accepted, and `fm` for those missed._

`define_df_gst_photon_acceptance(df_gst,phvars=["E","px","py","pz","pt","thetaz"],ke_min=PHOTON_ACCEPT_KE,thetaz_max=ECAL_ACCEPT_ANGLE)`  
                                    
_Returns an RDataFrame with pi0 decay photon variables selected for detector acceptance. `phvars` is the list of photon variables to create, `ke_min` is the kinetic energy threshold (GeV), and `thataz_max` is the angular threshold (radians)._

`define_df_gst_hadron_sums(df_gst,particles=["proton","neutron","piplus","piminus"],sfx=["i","f"],do_pi0=True)`  
_Returns an RDataFrame that calculates hadronic system sums used for momentum imbalance calculations. `particles` is the list of particles to consider, `sfx` is the list of suffixes to calculate for, and `do_pi0` is a boolean on whether to consider neutral pions or not._

`define_df_gst_momentum_imbalance(df_gst,suffix_list,cname)`  
_Returns an RDataFrame with momentum imbalance variables. `suffix_list` is the list of suffixes to consider, and `cname` is a custom name to attach to the end of momentum imblance variable column names._

### List of variables in the common DataFrame

GENIE `gst` variables can be found in the [GENIE user manual](https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=2) (section 9.5.2.1).

Description of other generated variables below.

#### Event-level kinematic variables
- `ptl`: Transverse momentum of final state lepton (GeV/c)
- `energy_transfer`: $E_{\nu} - E_{l}$ (GeV)

#### Hadron/other final state particle variables
Generally, variables try to follow a convention like `{variable_name}{suffix}_{particle_type}`. Variables are generally all `RVec` vector types (except in cases like `pi0_ph1` and `pi0_ph2`).

##### Common suffixes
Suffixes are generally:
- `i`: initial state particles (pre-FSI)
- `f`: final state particles (post-FSI, but before any detector acceptance/resolution effects)
- `fa`: final state accepted particles (with some detector acceptance applied)
- `fm`: final state not accepted ("missed") particles (those that fail detector acceptance)

##### Common particle types
Particle types are:
- `proton`: both proton (`pdg=2212`) and antiprotons (`pdg=-2212`)
- `neutron`: both neutrons (`pdg=2112`) and antineutrons (`pdg=-2112`)
- `piplus`: $\pi^{+}$ (`pdg=211`)
- `piminus`: $\pi^{-}$ (`pdg=-211`)
- `pi0`: $\pi^{0}$ (`pdg=111`)
- `K0`: $K^{0}$ (`pdg=311`) and $\bar{K^{0}}$ (`pdg=-311`)
- `Kplus`: $K^{+}$ (`pdg=321`)
- `Kminus`: $K^{-}$ (`pdg=-321`)
- `pi0_ph1`,`pi0_ph2`: photons from decay of $\pi^{0}$

##### Kinematic variables
- `E`: total energy (GeV)
- `ke`: kinetic energy (GeV)
- `mass`: particle mass (GeV/c/c)
- `pdg`: PDG code
- `px`,`py`,`pz`: components of three-momentum (GeV/c)
- `pt`: transverse momentum (GeV/c)
- `p`: total momentum (GeV/c)
- `sum_{var}`: sum of variable (generally for a given suffix and particle type)
- `thetaxz`: angle from `z=0` in the _xz_-plane (radians)
- `thetayz`: angle from `z=0` in the _yz_-plane (radians)
- `thetaz`: _zenith_ angle from `z=0` (radians)

##### Counting/acceptance variables
- `n`: number of particles in hadronic system
- `Acceptf`: whether or not (final state) particle is within detector acceptance cuts

##### Momentum imbalance variables
Most momentum imbalance variables may also have a final `_{optional_descriptor}`, which can used to speficy which final state particles are included in the calculation, like `all` or `chhad` for charged hadrons only.

- `hsum_{var}{suffix}`: sum of hadronic elements for momentum imbalance (see kinematic variable list above)

- `delta_alpha{suffix}`: 3D-momentum imbalance angle $\delta\alpha$
- `delta_alphat{suffix}`: transverse (2D) momentum imbalance angle $\delta\alpha_{T}$
- `delta_cosalpha{suffix}`: $cos(\delta\alpha)$
- `delta_cosalphat{suffix}`: $cos(\delta\alpha_{T})$

- `delta_phi{suffix}`: 3D-momentum imbalance angle $\delta\phi$
- `delta_phit{suffix}`: transverse (2D) momentum imbalance angle $\delta\phi_{T}$
- `delta_cosphi{suffix}`: $cos(\delta\phi)$
- `delta_cosphi{suffix}`: $cos(\delta\phi_{T})$

- `delta_p{suffix}`: 3D-momentum imbalance magnitude $\delta p$
- `delta_pt{suffix}`: transverse (2D) momentum imbalance magnitude $\delta p_{T}$
