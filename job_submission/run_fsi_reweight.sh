#!/bin/bash

#inputs
inputfile=$1
out_loc=$2

#setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup genie v3_00_06p -qe20:prof
setup ifdhc

echo "IFDHC SETUP COMPLETE"

ifdh cp /pnfs/dune/resilient/users/wketchum/genie_rw_fsi_10Oct2021.tar ./genie_rw_fsi_10Oct2021.tar
tar -xf genie_rw_fsi_10Oct2021.tar
export GENIE_INC=`pwd`/install/include/
export GENIE_LIB=`pwd`/install/lib
#ifdh cp /pnfs/dune/resilient/users/wketchum/Messenger.xml ./Messenger.xml

export PYTHONUSERBASE=`pwd`/install/python_libs
export PYTHONPATH=$PYTHONUSERBASE/bin:$PYTHONPATH
export PATH=$PYTHONUSERBASE/bin:$PATH

echo "GENIE SETUP COMPLETE"
echo $GENIE_INC
echo $GENIE_LIB
echo $PYTHONPATH

#copy in input file
myinput=`basename ${inputfile}`
ifdh cp $inputfile ./$myinput

echo "File ${myinput} copied."

##first make the gst file
out_filename_gst="${myinput%%.*}_gst.root"
gntpc -i $myinput -o $out_filename_gst -f gst &> /dev/null


#genie RW pars
seed=1000
ntwk=5
min_twk=-2
max_twk=2

knobs=(
    "FrAbs_N"
    "FrInel_N"
    "FrPiProd_N"
    "FrCEx_N"
    "FrAbs_pi"
    "FrInel_pi"
    "FrPiProd_pi"
    "FrCEx_pi"
)

for k in "${knobs[@]}"; do
    echo "Launch knob $k"
    out_filename="${myinput%%.*}_weights_${k}.root"
    install/bin/grwght1p -f $myinput -s $k -t $ntwk --min-tweak $min_twk --max-tweak $max_twk --seed $seed -p 11 -o $out_filename --event-record-print-level 0 --message-thresholds Messenger.xml >& /dev/null
done

#do the python part now
#python install/python/make_gst_dataframes.py

##copyback
ifdh cp ./$out_filename_gst $out_loc/$out_filename_gst

for f in ./*_weights_*.root; do
    ifdh cp $f $out_loc/$f
done

for f in ./*.hdf; do
    ifdh cp $f $out_loc/$f
done


#cleanup
rm -rf install
rm -rf *.tar
rm -rf *.root
rm -rf *.xml
rm -rf *.hdf
