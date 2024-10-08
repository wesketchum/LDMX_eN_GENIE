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
    OUTPUT_FILE_NAME= f'ldmx_genie_{TUNE}_{TARGET}_{ENERGY}GeV_{RUN}.root'

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
from LDMX.SimCore import genie_reweight

p=ldmxcfg.Process("genie")

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
#genie_rw.hepmc3RunInfoCollName="SimHepMC3RunInfo"
genie_rw.hepmc3PassName=""
genie_rw.var_types = ["GENIE_INukeTwkDial_MFP_pi","GENIE_INukeTwkDial_MFP_N"]
genie_rw.verbosity = 5

sim.generators = [ genie ]
p.sequence.append(sim)
p.sequence.append(genie_rw)
p.sequence.append( dqm.GenieTruthDQM(coll_name="SimHepMC3Events") )
p.outputFiles=[OUTPUT_FILE_NAME]
p.maxEvents = N_EVENTS
p.run = RUN
p.logFrequency = 1
p.histogramFile = HIST_OUTPUT_FILE_NAME

