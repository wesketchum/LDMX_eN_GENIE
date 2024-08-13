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

arg = parser.parse_args()

#basic params for running
N_EVENTS = arg.n_events
TARGET=arg.target
RUN=arg.run
TUNE=arg.tune
VERBOSITY=arg.verbosity
OUTPUT_FILE_NAME=arg.output
ENERGY=arg.energy

if OUTPUT_FILE_NAME is None:
    OUTPUT_FILE_NAME= f'ldmx_genie_{TUNE}_{TARGET}_{ENERGY}GeV_{RUN}_100MeVphotongun.root'

if OUTPUT_FILE_NAME[-5:]!=".root":
    OUTPUT_FILE_NAME=OUTPUT_FILE_NAME+".root"

output_file_name_local=OUTPUT_FILE_NAME.split("/")[-1]
HIST_OUTPUT_FILE_NAME = OUTPUT_FILE_NAME[:len(OUTPUT_FILE_NAME)-len(output_file_name_local)]+output_file_name_local[:-5]+"_hist.root"

#locations of things we need...
PATH_TO_GENIE_SPLINES=arg.genie_splines
GENIE_MESSENGER_XML_FILE=arg.genie_messenger_xml

#function for filling targets and abundances
def get_targets(target_name):

    targets = []
    abundances = []

    if(target_name=='Ti'):
        targets = [ 1000220460, 1000220470, 1000220480, 1000220490, 1000220500 ]
        abundances = [ 0.0825, 0.0744, 0.7372, 0.0541, 0.0581 ]
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
from LDMX.SimCore import generators
from LDMX.SimCore import simulator

p=ldmxcfg.Process("singlephoton")

import LDMX.Ecal.EcalGeometry
import LDMX.Hcal.HcalGeometry

sim = simulator.simulator('sim')
sim.setDetector(det_name=DET_NAME,include_scoring_planes=True)

myGun = generators.gun( 'myGun' )
myGun.particle = 'gamma'
myGun.energy = 0.1
myGun.direction = [ 0., 0., 1. ]
myGun.position = [ 0., 0., 0. ]

sim.generators = [ myGun ]
p.sequence.append(sim)
p.outputFiles=[OUTPUT_FILE_NAME]
p.maxEvents = N_EVENTS
p.run = RUN
p.logFrequency = 1
p.histogramFile = HIST_OUTPUT_FILE_NAME

###reco parts
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions
import LDMX.Ecal.digi as ecal_digi
import LDMX.Hcal.digi as hcal_digi

p.sequence.extend([
    ecal_digi.EcalDigiProducer(),
    ecal_digi.EcalRecProducer(),
    hcal_digi.HcalDigiProducer(),
    hcal_digi.HcalRecProducer()
])

