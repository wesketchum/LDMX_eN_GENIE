#!/bin/sh

Help()
{
    echo "Usage:"
    echo "run_genie_sim_and_reco.sh [-h] -n number_of_events -r run -t target -T tune -e energy (GeV)"
    echo
}

N_EVENTS=100
RUN=999
TARGET=Ti
TUNE=G18_02a_02_11b
ENERGY=4

while getopts ":hn:r:t:T:e:" option; do
    case $option in
	h) #display Help
	    Help
	    exit;;
	n) N_EVENTS=$OPTARG;;
	r) RUN=$OPTARG;;
	t) TARGET=$OPTARG;;
	T) TUNE=$OPTARG;;
	e) ENERGY=$OPTARG;;
	\?) #Invalid
	    echo "Error: Invalid option"
	    Help
	    exit;;
    esac
done

TMP_OUTDIR=./ldmx_genie_output_run_${RUN}
mkdir -p $TMP_OUTDIR

#sim_file=$TMP_OUTDIR/ldmx_genie_${TUNE}_${TARGET}_${ENERGY}GeV_${RUN}.root
#echo $sim_file

#fire /LDMX_eN_GENIE/ldmxsw_configs/genie_sim.py -n $N_EVENTS -r $RUN --target $TARGET --tune $TUNE --energy $ENERGY --genie_splines /ldmx-genie-splines --output $sim_file --genie_messenger_xml /LDMX_eN_GENIE/ldmxsw_configs/Messenger_ErrorOnly.xml

reco_file=$TMP_OUTDIR/ldmx_genie_${TUNE}_${TARGET}_${ENERGY}GeV_${RUN}_reco.root
echo $reco_file
fire /LDMX_eN_GENIE/ldmxsw_configs/genie_sim_reco.py -n $N_EVENTS -r $RUN --target $TARGET --tune $TUNE --energy $ENERGY --genie_splines /ldmx-genie-splines --output $reco_file --genie_messenger_xml /LDMX_eN_GENIE/ldmxsw_configs/Messenger_ErrorOnly.xml --prop_map /LDMX_eN_GENIE/propagationMap.root

python3 /LDMX_eN_GENIE/analysis/ldmx_analysis.py -i $reco_file

#rm $sim_file

tar -czf ldmx_genie_output_run_${RUN}.tgz $TMP_OUTDIR
rm -rf $TMP_OUTDIR
