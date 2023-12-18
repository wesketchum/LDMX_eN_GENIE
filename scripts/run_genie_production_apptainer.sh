#!/bin/bash

export EOS_MGM_URL=root://eosuser.cern.ch
export APPTAINER_CACHEDIR="/tmp/$(whoami)/apptainer"

N_EVENTS=100
TARGET=Ti
TUNE=G18_02a_02_11b
ENERGY=4
CLUSTER_ID=100
PROC_ID=0
EOS_OUTDIR=./

while getopts ":hn:t:T:e:c:p:o:" option; do
    case $option in
	h) #display Help
	    Help
	    exit;;
	n) N_EVENTS=$OPTARG;;
	t) TARGET=$OPTARG;;
	T) TUNE=$OPTARG;;
	e) ENERGY=$OPTARG;;
	c) CLUSTER_ID=$OPTARG;;
	p) PROC_ID=$OPTARG;;
	o) EOS_OUTDIR=$OPTARG;;
	\?) #Invalid
	    echo "Error: Invalid option"
	    Help
	    exit;;
    esac
done

RUN=$((CLUSTER_ID+PROC_ID))
TMP_OUTDIR=./ldmx_genie_output_run_${RUN}_${CLUSTER_ID}_${PROC_ID}

mkdir -p $TMP_OUTDIR

apptainer run -B $TMP_OUTDIR:/output:rw docker://ghcr.io/wesketchum/ldmxsw_genie_prod:prod -n $N_EVENTS -r $RUN -t $TARGET -T $TUNE -e $ENERGY

mv *.root $TMP_OUTDIR

tar -czf ldmx_genie_output_run_${RUN}_${CLUSTER_ID}_${PROC_ID}.tgz $TMP_OUTDIR
eos cp ldmx_genie_output_run_${RUN}_${CLUSTER_ID}_${PROC_ID}.tgz $EOS_OUTDIR
