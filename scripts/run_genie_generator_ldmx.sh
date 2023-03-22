#!/bin/sh

source ~/ldmx-sw/scripts/ldmx-env.sh
ldmx use dev latest

alias ldmx_detached='docker run --rm -d -e LDMX_BASE -v $LDMX_BASE:$LDMX_BASE ${LDMX_DOCKER_TAG} $(pwd)'

ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 111 -s 11
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 112 -s 12
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 113 -s 13
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 114 -s 14
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 115 -s 15
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 116 -s 16
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 117 -s 17
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 118 -s 18
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 119 -s 19
ldmx_detached fire ldmxsw_configs/genie_sim_reco.py -n 100000 -r 120 -s 20
