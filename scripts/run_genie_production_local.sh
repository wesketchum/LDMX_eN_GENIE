#!/bin/bash

N_EVENTS_PER_JOB=20000
N_JOBS=100
N_MAX_RUNNING=12
START_RUN_NUMBER=1000
TARGET=Ti
ENERGY=8
TUNE=G18_02a_02_11b

LDMX_DOCKER_TAG=ghcr.io/wesketchum/ldmxsw_genie_prod:latest

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
    docker run -v $OUTPUT_DIR:/output -d $LDMX_DOCKER_TAG -n $N_EVENTS_PER_JOB -r $RUN -t $TARGET -T $TUNE -e $ENERGY
done

