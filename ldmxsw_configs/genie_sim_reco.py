#basic params for running
N_EVENTS = 100
TARGET='Ti'
OUTPUT_FILE_NAME = f'test_{TARGET}.root'
SEED=10
VERBOSITY=0

#locations of things we need...
PATH_TO_GENIE_SPLINES='/Users/wketchum/ldmx-genie-splines'
GENIE_MESSENGER_XML_FILE='./ldmxsw_configs/Messenger.xml'

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

from LDMX.Framework import ldmxcfg
from LDMX.SimCore import generators
from LDMX.SimCore import simulator

p=ldmxcfg.Process("genie")

import LDMX.Ecal.EcalGeometry
import LDMX.Hcal.HcalGeometry

sim = simulator.simulator('sim')
sim.setDetector('ldmx-det-v14')

targets, abundances = get_targets(TARGET)

genie_G18_02a_00_000 = generators.genie(name='genie_G18_02a_00_000',
                                        energy = 4.0,
                                        targets = targets,
                                        abundances = abundances,
                                        time = 0.0,
                                        position = [0.,0.,0.],
                                        direction = [0.,0.,1.],
                                        tune='G18_02a_00_000',
                                        spline_file=f'{PATH_TO_GENIE_SPLINES}/G18_02a_00/gxspl_emode_G18_02a_00.xml',
                                        message_threshold_file=GENIE_MESSENGER_XML_FILE,
                                        verbosity=VERBOSITY,
                                        seed=SEED)

sim.generators = [ genie_G18_02a_00_000 ]
p.sequence.append(sim)
p.outputFiles=[OUTPUT_FILE_NAME]
p.maxEvents = N_EVENTS
p.logFrequency = 1


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
