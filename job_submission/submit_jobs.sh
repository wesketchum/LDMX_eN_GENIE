#!/usr/bin/bash

#n_events=$1
#energy=$2
#spline=$3
#tune=$4
#out_name=$5
#out_loc=$6

n_events=100000


#setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc
setup jobsub_client

#out_name="ldmx_eTi_4GeV_GEM21_11b_00"
out_name="ldmx_eTi_4GeV_G18_02a_00"


out_loc="/pnfs/dune/scratch/users/wketchum/eN_Ti_GENIE_v3_2_0/${out_name}"
mkdir -p $out_loc


#jobsub_submit -N 100 -G dune file:///dune/app/users/wketchum/GENIE_v3_2_0/job_submission/run_gevgen.sh $n_events 4.0 gxspl_emode_Ti48_GEM21_11b_00.xml GEM21_11b_00_000 $out_name $out_loc --expected-lifetime=24h --OS=SL7 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --disk=20GB --memory=2000MB

jobsub_submit -N 100 -G dune --expected-lifetime=24h --OS=SL7 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --disk=20GB --memory=2000MB file:///dune/app/users/wketchum/LDMX_eN_GENIE/job_submission/run_gevgen.sh $n_events 4.0 gxspl_emode_Ti48_default.xml G18_02a_00_000 $out_name $out_loc
