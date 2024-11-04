import argparse
import sys

parser = argparse.ArgumentParser(f'ldmx fire {sys.argv[0]}')
parser.add_argument('-n','--n_events',default=1000,type=int)
parser.add_argument('-t','--target',default='Ti')
parser.add_argument('-r','--run',default=100,type=int)
parser.add_argument('--tune',default='G18_02a_02_11b')
parser.add_argument('-v','--verbosity',default=0,type=int)
parser.add_argument('-e', '--energy', default=4,type=int)
parser.add_argument('-o','--output',default=None)
parser.add_argument('--genie_splines',default='./')
parser.add_argument('--genie_messenger_xml',default='./Messenger_ErrorOnly.xml')
parser.add_argument('--prop_map',default='./propagationMap.root')
parser.add_argument('--sim_only',default=False)
parser.add_argument('--reco_only',default=False)

arg = parser.parse_args()

if arg.sim_only and arg.reco_only:
    print(f'Cannot use both sim_only and reco_only flags!')
    sys.exit()

#basic params for running
N_EVENTS = arg.n_events
TARGET=arg.target
RUN=arg.run
TUNE=arg.tune
VERBOSITY=arg.verbosity
OUTPUT_FILE_NAME=arg.output
ENERGY=arg.energy

if OUTPUT_FILE_NAME is None:
    OUTPUT_FILE_NAME= f'ldmx_genie_{TUNE}_{TARGET}_{ENERGY}GeV_{RUN}.root'

if OUTPUT_FILE_NAME[-5:]!=".root":
    OUTPUT_FILE_NAME=OUTPUT_FILE_NAME+".root"

output_file_name_local=OUTPUT_FILE_NAME.split("/")[-1]
HIST_OUTPUT_FILE_NAME = OUTPUT_FILE_NAME[:len(OUTPUT_FILE_NAME)-len(output_file_name_local)]+output_file_name_local[:-5]+"_hist.root"

#locations of things we need...
PATH_TO_GENIE_SPLINES=arg.genie_splines
GENIE_MESSENGER_XML_FILE=arg.genie_messenger_xml

PROPAGATION_MAP=arg.prop_map

DO_SIM = not arg.reco_only;
DO_RECO = not arg.sim_only;

#function for filling targets and abundances
def get_targets(target_name):

    targets = []
    abundances = []

    if(target_name=='Ti'):
        #targets = [ 1000220460, 1000220470, 1000220480, 1000220490, 1000220500 ]
        #abundances = [ 0.0825, 0.0744, 0.7372, 0.0541, 0.0581 ]
        targets = [ 1000220460, 1000220470, 1000220480, 1000220490 ]
        abundances = [ 0.0825, 0.0744, 0.7372, 0.0541 ]
    elif(target_name=='W'):
        targets = [ 1000741820, 1000741830, 1000741840, 1000741860 ]
        abundances = [ 0.2650, 0.1431, 0.3064, 0.2843 ]
    else:
        print("Unsupported target name: ",TARGET)
        print("Default to Ti...")
        return get_targets('Ti')

    return targets, abundances

if ENERGY==4:
    DET_NAME="ldmx-det-v14"
elif ENERGY==8:
    DET_NAME="ldmx-det-v14-8gev"
else:
    print("Energy must be 4 or 8 GeV!")
    sys.exit()

from LDMX.Framework import ldmxcfg

p=ldmxcfg.Process("genie")

if DO_SIM:
    from LDMX.SimCore import generators
    from LDMX.SimCore import simulator
    from LDMX.SimCore import genie_reweight
    import LDMX.Ecal.EcalGeometry
    import LDMX.Hcal.HcalGeometry
    from LDMX.DQM import dqm

    sim = simulator.simulator('sim')
    sim.setDetector(det_name=DET_NAME,include_scoring_planes=True)

    targets, abundances = get_targets(TARGET)

    genie = generators.genie(name=f'genie_{TUNE}',
                             energy = float(ENERGY),
                             targets = targets,
                             target_thickness = 0.3504,
                             abundances = abundances,
                             time = 0.0,
                             position = [0.,0.,0.],
                             beam_size = [ 20., 80. ],
                             direction = [0.,0.,1.],
                             tune=TUNE,
                             spline_file=f'{PATH_TO_GENIE_SPLINES}/gxspl_emode_GENIE_v3_04_00.xml',
                             message_threshold_file=GENIE_MESSENGER_XML_FILE,
                             verbosity=VERBOSITY)

    genie_rw = genie_reweight.GenieReweightProducer(name='genie_reweight')
    genie_rw.hepmc3CollName="SimHepMC3Events"
    genie_rw.hepmc3PassName=""
    genie_rw.var_types = ["GENIE_INukeTwkDial_MFP_pi","GENIE_INukeTwkDial_MFP_N"]
    genie_rw.verbosity = VERBOSITY

    sim.generators = [ genie ]

    p.sequence.append(sim)
    p.sequence.append(genie_rw)
    #p.sequence.append( dqm.GenieTruthDQM(coll_name="SimHepMC3Events") )

if DO_RECO:
    ###reco parts
    import LDMX.Ecal.EcalGeometry
    import LDMX.Ecal.ecal_hardcoded_conditions
    import LDMX.Hcal.HcalGeometry
    import LDMX.Hcal.hcal_hardcoded_conditions
    import LDMX.Ecal.digi as ecal_digi
    import LDMX.Ecal.vetos as ecal_vetos
    import LDMX.Hcal.digi as hcal_digi
    import LDMX.Hcal.vetos as hcal_vetos

    from LDMX.Hcal import hcal_trig_digi
    from LDMX.Ecal import ecal_trig_digi
    from LDMX.Trigger import trigger_energy_sums
    from LDMX.Recon import pfReco



    p.sequence.extend([
        ecal_digi.EcalDigiProducer(),
        ecal_trig_digi.EcalTrigPrimDigiProducer(),
        ecal_digi.EcalRecProducer(),
        hcal_digi.HcalDigiProducer(),
        hcal_trig_digi.HcalTrigPrimDigiProducer(),
        hcal_digi.HcalRecProducer(),
        trigger_energy_sums.EcalTPSelector(),
        trigger_energy_sums.TrigEcalEnergySum(),
        trigger_energy_sums.TrigHcalEnergySum(),
        trigger_energy_sums.TrigEcalClusterProducer(),
        trigger_energy_sums.TrigElectronProducer(propMapName=PROPAGATION_MAP),
        pfReco.pfTruthProducer()
    ])

    ###tracking parts

    from LDMX.Tracking import tracking
    from LDMX.Tracking import geo

    # Truth seeder
    truth_tracking           = tracking.TruthSeedProcessor()
    truth_tracking.debug             = True
    truth_tracking.trk_coll_name     = "RecoilTruthSeeds"
    #truth_tracking.pdgIDs            = [11]
    truth_tracking.scoring_hits      = "TargetScoringPlaneHits"
    truth_tracking.z_min             = 0.
    truth_tracking.track_id          = -1
    truth_tracking.p_cut             = 0.05 # In GeV
    truth_tracking.pz_cut            = 0.03
    truth_tracking.p_cutEcal         = 0. # In GeV

    #smearings
    uSmearing = 0.006       #mm #could bump up to 10 micron if we want
    vSmearing = 0.000001    #mm #~unused

    # Smearing Processor - Tagger
    # Runs G4 hit smearing producing measurements in the Tagger tracker.
    # Hits that belong to the same sensor with the same trackID are merged together to reduce combinatorics
    digiTagger = tracking.DigitizationProcessor("DigitizationProcessor")
    digiTagger.hit_collection = "TaggerSimHits"
    digiTagger.out_collection = "DigiTaggerSimHits"
    digiTagger.merge_hits = True
    digiTagger.sigma_u = uSmearing
    digiTagger.sigma_v = vSmearing

    # Smearing Processor - Recoil
    digiRecoil = tracking.DigitizationProcessor("DigitizationProcessorRecoil")
    digiRecoil.hit_collection = "RecoilSimHits"
    digiRecoil.out_collection = "DigiRecoilSimHits"
    digiRecoil.merge_hits = True
    digiRecoil.sigma_u = uSmearing
    digiRecoil.sigma_v = vSmearing

    #Seed finder processor - Recoil
    seederRecoil = tracking.SeedFinderProcessor("SeedRecoil")
    seederRecoil.perigee_location = [0.,0.,0.]
    seederRecoil.input_hits_collection =  digiRecoil.out_collection
    seederRecoil.out_seed_collection = "RecoilRecoSeeds"
    seederRecoil.bfield = 1.5
    seederRecoil.pmin  = 0.1
    seederRecoil.pmax  = 8.
    seederRecoil.d0min = -0.5
    seederRecoil.d0max = 0.5
    seederRecoil.z0max = 10.

    #Recoil tracking
    #CKF Options
    tracking_recoil  = tracking.CKFProcessor("Recoil_TrackFinder")
    tracking_recoil.dumpobj = False
    tracking_recoil.debug = False
    tracking_recoil.propagator_step_size = 1000.  #mm
    #tracking_recoil.bfield = -1.5  #in T #From looking at the BField map #not used if below is False
    tracking_recoil.const_b_field = False

    #Target location for the CKF extrapolation
    #tracking_recoil.seed_coll_name = seederRecoil.out_seed_collection
    tracking_recoil.seed_coll_name = "RecoilTruthSeeds"
    tracking_recoil.out_trk_collection = "RecoilTracks"

    #smear the hits used for finding/fitting
    tracking_recoil.trackID = -1 #1
    tracking_recoil.pdgID = -9999 #11
    tracking_recoil.measurement_collection = digiRecoil.out_collection
    tracking_recoil.min_hits = 6  #truth runs with at least 7 hits, so require the same here too

    from LDMX.Tracking import dqm
    digi_dqm = dqm.TrackerDigiDQM()
    tracking_dqm = dqm.TrackingRecoDQM()

    #seed_recoil_dqm = dqm.TrackingRecoDQM("SeedRecoilDQM")
    #seed_recoil_dqm.track_collection = seederRecoil.out_seed_collection
    #seed_recoil_dqm.truth_collection = "RecoilTruthTracks"
    #seed_recoil_dqm.title = ""

    recoil_dqm = dqm.TrackingRecoDQM("RecoilDQM")
    recoil_dqm.track_collection = tracking_recoil.out_trk_collection
    recoil_dqm.truth_collection = "RecoilTruthTracks"
    recoil_dqm.title = ""

    p.sequence.extend([ digiRecoil,
                        truth_tracking,
                        #seederTagger, seederRecoil,
                        # #tracking_tagger,
                        tracking_recoil ]),
#                      recoil_dqm ])#, seed_recoil_dqm]


p.outputFiles=[OUTPUT_FILE_NAME]
p.maxEvents = N_EVENTS
p.run = RUN
p.logFrequency = 1
p.histogramFile = HIST_OUTPUT_FILE_NAME

#uses the run number as a seed by default. But if you want to set differently...
#p.randomNumberSeedService.external(RUN)
