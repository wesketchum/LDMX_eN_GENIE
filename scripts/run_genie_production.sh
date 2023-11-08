#!/bin/bash

export EOS_MGM_URL=root://eosuser.cern.ch

N_EVENTS_PER_JOB=$1
RUN_START=$2
TARGET=$3
TUNE=$4
ENERGY=$5
CLUSTER_ID=$6
PROC_ID=$7
EOS_OUTDIR=$8

RUN=$((RUN_START+PROC_ID))
TMP_OUTDIR=ldmx_genie_output_run_${RUN}_${CLUSTER_ID}_${PROC_ID}

mkdir $TMP_OUTDIR

docker run -e N_EVENTS_PER_JOB=100 -e RUN=9999 -e TARGET=TARGET -e TUNE=$TUNE -e ENERGY=$ENERGY -v $TMP_OUTDIR:/output ghcr.io/wesketchum/ldmxsw_genie_prod

tar -czf ldmx_genie_output_run_${RUN}_${CLUSTER_ID}_${PROC_ID}.tgz $TMP_OUTDIR
eos cp output.tgz $EOS_OUTDIR
