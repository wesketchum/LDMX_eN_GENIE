import argparse
import sys

parser = argparse.ArgumentParser(f'ldmx fire {sys.argv[0]}')
parser.add_argument('-i','--input',required=True)
parser.add_argument('-o','--output',default=None)

arg = parser.parse_args()
OUTPUT_FILE_NAME=arg.output

if OUTPUT_FILE_NAME is None:
    OUTPUT_FILE_NAME= f'{arg.input[:-5]}_reco.root'

if OUTPUT_FILE_NAME[-5:]!=".root":
    OUTPUT_FILE_NAME=OUTPUT_FILE_NAME+".root"

output_file_name_local=OUTPUT_FILE_NAME.split("/")[-1]
HIST_OUTPUT_FILE_NAME = OUTPUT_FILE_NAME[:len(OUTPUT_FILE_NAME)-len(output_file_name_local)]+output_file_name_local[:-5]+"_hist.root"

from LDMX.Framework import ldmxcfg

p=ldmxcfg.Process("genie_reco")

import LDMX.Ecal.EcalGeometry
import LDMX.Hcal.HcalGeometry

#p.inputFiles=["/Users/wketchum/Data/LDMX/ldmx_genie_G18_02a_00_000_Ti_501.root"]
p.inputFiles=[arg.input]
p.outputFiles=[OUTPUT_FILE_NAME]
p.logFrequency = 1
p.histogramFile = HIST_OUTPUT_FILE_NAME

#uses the run number as a seed by default. But if you want to set differently...
#p.randomNumberSeedService.external(RUN) 



###reco parts
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions
import LDMX.Ecal.digi as ecal_digi
import LDMX.Hcal.digi as hcal_digi

p.sequence = [
    ecal_digi.EcalDigiProducer(),
    ecal_digi.EcalRecProducer(),
    hcal_digi.HcalDigiProducer(),
    hcal_digi.HcalRecProducer()
]

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
                      #tracking_tagger,
                      tracking_recoil])
#                      recoil_dqm ])#, seed_recoil_dqm]

