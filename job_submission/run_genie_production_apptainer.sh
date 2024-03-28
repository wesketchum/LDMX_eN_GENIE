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

apptainer run docker://ghcr.io/wesketchum/ldmxsw_genie_prod:latest -n $N_EVENTS -r $RUN -t $TARGET -T $TUNE -e $ENERGY

eos cp ldmx_genie_output_run_${RUN}.tgz $EOS_OUTDIR/
rm ldmx_genie_output_run_${RUN}.tgz
