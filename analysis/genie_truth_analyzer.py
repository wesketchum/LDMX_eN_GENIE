from LDMX.Framework import EventTree
import numpy as np
import ROOT
from analysis_utilities import *
import glob
#from array import array

import argparse
import sys

import pyHepMC3
import tempfile

def string_to_4vector(istring):
    return np.array(istring.split(" "),dtype='f')


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True)

arg = parser.parse_args()

HepMC3_Print = pyHepMC3.HepMC3.Print()

tmp = tempfile.NamedTemporaryFile()
HepMC3_Reader = pyHepMC3.HepMC3.ReaderAscii(tmp.name)

FILES = glob.glob(arg.input)
for ifile in FILES:

    input_tree = EventTree.EventTree(ifile)

    for ie, event in enumerate(input_tree):

        sim_events = event.SimHepMC3Events_genie
        print(len(sim_events))
        for sim_event in sim_events:

            #this grabs the hepmc3 object out of the event, but it's not easy to use with pyHepMC3 as is
            #can use the various C++ functions though
            hepmc_event = sim_event.getHepMCGenEvent()
            e_initial_p4 = string_to_4vector(hepmc_event.attribute_as_string("GENIE.Interaction.ProbeP4"))
            N_initial_p4 = string_to_4vector(hepmc_event.attribute_as_string("GENIE.Interaction.TargetP4"))
            e_final_p4 = string_to_4vector(hepmc_event.attribute_as_string("GENIE.Interaction.FSLeptonP4"))

            #use print function from ldmx-sw definition
            sim_event.Print()

            #as a hack, we can conver the object to a string, write to a file, and then read back up in pyHepMC3
            ev_string = sim_event.get_as_string()
            print(ev_string)
            ev_regen = pyHepMC3.HepMC3.GenEvent()
            with open(tmp.name, 'w') as f:
                f.write(ev_string)
            pyHepMC3.HepMC3.ReaderAscii(tmp.name).read_event(ev_regen)
            HepMC3_Print.content(ev_regen)





