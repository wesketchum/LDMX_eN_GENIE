#!/bin/bash

N_EVENTS_PER_JOB=20000
N_JOBS=100
N_MAX_RUNNING=12
START_RUN_NUMBER=200
TARGET=Ti
TUNE=G18_02a_02_11b

LDMX_DOCKER_TAG=ghcr.io/wesketchum/ldmxsw_genie_prod:latest

SPLINES_DIR=/Users/wketchum/ldmx-genie-splines
EN_GENIE_DIR=/Users/wketchum/LDMX_eN_GENIE
OUTPUT_DIR=/Users/wketchum/Data/LDMX/production_27Oct2023

for (( c=0; c<$N_JOBS; c++ ))
do
    running=`docker ps | grep $LDMX_DOCKER_TAG | wc -l | awk '{print $1}'`
    while [ $running -ge $N_MAX_RUNNING ]
    do
	sleep 60
	running=`docker ps | grep $LDMX_DOCKER_TAG | wc -l | awk '{print $1}'`
    done
    RUN=$((START_RUN_NUMBER+c))
    echo "Starting run ${RUN} -- ${N_EVENTS_PER_JOB} events."
    docker run -d -e N_EVENTS_PER_JOB=$N_EVENTS_PER_JOB -e RUN=$RUN -e TARGET=$TARGET -e TUNE=$TUNE -v $HOME:$HOME -v $OUTPUT_DIR:/output -v $SPLINES_DIR:/ldmx-genie-splines -v $EN_GENIE_DIR:/LDMX_eN_GENIE $LDMX_DOCKER_TAG $HOME $EN_GENIE_DIR/scripts/run_genie_sim_and_reco.sh
done

