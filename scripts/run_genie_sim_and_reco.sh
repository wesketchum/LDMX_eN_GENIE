#!/bin/sh

sim_file=/output/ldmx_genie_${TUNE}_${TARGET}_${ENERGY}GeV_${RUN}.root
echo $sim_file

fire /LDMX_eN_GENIE/ldmxsw_configs/genie_sim.py -n $N_EVENTS_PER_JOB -r $RUN --target $TARGET --tune $TUNE --genie_splines /ldmx-genie-splines --output_dir /output --genie_messenger_xml /LDMX_eN_GENIE/ldmxsw_configs/Messenger_ErrorOnly.xml

fire /LDMX_eN_GENIE/ldmxsw_configs/genie_reco_only.py -i $sim_file

python3 /LDMX_eN_GENIE/analysis/ldmx_analysis.py

rm $sim_file
