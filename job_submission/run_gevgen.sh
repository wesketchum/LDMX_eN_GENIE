#!/bin/bash

#inputs
run_number=`expr $PROCESS + 1`
seed=$run_number
n_events=$1
energy=$2
spline=$3
tune=$4
out_name=$5
out_loc=$6

#setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup log4cpp v1_1_3c -qe20:prof
setup libxml2 v2_9_10a
setup lhapdf v6_3_0 -qe20:p392:prof
setup root v6_22_08d -qe20:p392:prof
setup ifdhc

echo "INITIAL IFDH SETUP COMPLETE"

ifdh cp /pnfs/dune/resilient/users/wketchum/GENIE_v3_2_0_21April2022.tar GENIE_v3_2_0.tar
#cp /pnfs/dune/resilient/users/wketchum/GENIE_v3_2_0_21April2022.tar GENIE_v3_2_0.tar

tar -xvf GENIE_v3_2_0.tar

export GENIE=`pwd`/GENIE_v3_2_0/Generator
export GENIE_REWEIGHT=`pwd`/GENIE_v3_2_0/Reweight

export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE_REWEIGHT/lib:$LD_LIBRARY_PATH
export PATH=$GENIE_REWEIGHT/bin:$PATH

echo "GENIE SETUP COMPLETE"

out_ghep_filename="${out_name}_run${run_number}.ghep.root"
out_gevgen_logname="${out_name}_run${run_number}.gevgen.log"

gevgen -r $run_number -n $n_events -p 11 -t 1000220480 --cross-sections `pwd`/GENIE_v3_2_0/GENIE_splines/$spline --seed $seed -e $energy --event-generator-list EM --tune $tune -o $out_ghep_filename >& $out_gevgen_logname

##make the gst file
out_gst_filename="${out_name}_run${run_number}.gst.root"
gntpc -i $out_ghep_filename -o $out_gst_filename -f gst &> /dev/null

##copyback
ifdh cp ./$out_ghep_filename $out_loc/$out_ghep_filename
ifdh cp ./$out_gevgen_logname $out_loc/$out_gevgen_logname
ifdh cp ./$out_gst_filename $out_loc/$out_gst_filename

#cleanup
rm -rf GENIE_v3_2_0
rm -rf *.tar
rm -rf *.root
