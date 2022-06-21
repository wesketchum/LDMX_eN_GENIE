#!/usr/bin/env python3

import ROOT

def create_gst_chain(files,verbose=False):
    gst_chain = ROOT.TChain("gst")
    for f in files:
        gst_chain.Add(f)
    if(verbose):
        print(f'Created gst chain from {len(files)} files with {gst_chain.GetEntries()} total events.')
    return gst_chain