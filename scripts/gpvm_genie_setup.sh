#!/bin/sh

source /cvmfs/larsoft.opensciencegrid.org/products/setup
setup log4cpp v1_1_3c -qe20:prof
setup libxml2 v2_9_10a
setup lhapdf v6_3_0 -qe20:p392:prof
setup root v6_22_08d -qe20:p392:prof

export GENIE=`pwd`/Generator
export GENIE_REWEIGHT=`pwd`/Reweight

export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE_REWEIGHT/lib:$LD_LIBRARY_PATH
export PATH=$GENIE_REWEIGHT/bin:$PATH
