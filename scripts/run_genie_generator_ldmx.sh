#!/bin/sh

MY_LDMX_BASE=~/ldmx-sw
N_EVENTS_PER_JOB=20000
N_JOBS=100
N_MAX_RUNNING=10
START_RUN_NUMBER=300

source $MY_LDMX_BASE/scripts/ldmx-env.sh
ldmx use dev latest

#alias ldmx_detached='docker run --rm -d -e LDMX_BASE -v $LDMX_BASE:$LDMX_BASE ${LDMX_DOCKER_TAG} $(pwd)'

alias ldmx_detached='docker run -d -e LDMX_BASE -v $LDMX_BASE:$LDMX_BASE ${LDMX_DOCKER_TAG} $(pwd)'

for (( c=1; c<=$N_JOBS; c++ ))
do
    running=`docker ps | grep $LDMX_DOCKER_TAG | wc -l | awk '{print $1}'`
    while [ $running -ge $N_MAX_RUNNING ]
    do
	sleep 60
	running=`docker ps | grep $LDMX_DOCKER_TAG | wc -l | awk '{print $1}'`
    done
    run=$((START_RUN_NUMBER+c))
    echo "Starting run ${run} -- ${N_EVENTS_PER_JOB} events."
    ldmx_detached fire ~/LDMX_eN_GENIE/ldmxsw_configs/genie_sim_reco.py -n $N_EVENTS_PER_JOB -r $run --output_dir ~/Data/LDMX --genie_splines ~/ldmx-genie-splines --genie_messenger_xml ~/LDMX_eN_GENIE/ldmxsw_configs/Messenger_ErrorOnly.xml
done

