## Preparing a local-code Docker image

You can create a image that has locally built code, and an executable to run automatically if you wish. An example is below.

Important notes:
* you don't necessarily need to pull the two repos here, clearly. But just to show you can.
* The entrypoint of this container is changed to be one of the scripts that's in the LDMX_eN_GENIE repo. So, that means you can call the container and directly give it arguments to that script, which is very useful for simplifying running container jobs. Do something similar, or not at all.

```
###############################################################################
# This dockerfile is meant to build the production image of ldmx-sw
#   for the development image, look at the LDMX-Software/docker repo
###############################################################################

FROM ldmx/dev:4.1.0

ARG NPROC=1

# install ldmx-sw into the container at /usr/local
COPY . /code
RUN mkdir /code/build &&\
    /etc/entry.sh /code/build cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. &&\
    /etc/entry.sh /code/build make -j${NPROC} install &&\
    rm -rf code &&\
    ldconfig /usr/local/lib

ENV LDMX_SW_INSTALL=/usr/local

#need this for git
RUN sudo apt-get update&& sudo apt-get install git

WORKDIR /ldmx-genie-splines
RUN git clone https://github.com/wesketchum/ldmx-genie-splines.git

WORKDIR /LDMX_eN_GENIE
RUN git clone https://github.com/wesketchum/LDMX_eN_GENIE.git

VOLUME /output

ENTRYPOINT ["/bin/bash","/LDMX_eN_GENIE/scripts/run_genie_sim_and_reco.sh"]
```

## Job submission script
You'll want a script that can run on the worker node, pulling down the appropriate container and then doing something with it. Here's an example script [https://github.com/wesketchum/LDMX_eN_GENIE/blob/main/job_submission/run_genie_production_apptainer.sh]:
```
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
```
It sets up EOS copying, an APPTAINER cache directory, some arguments to be handled in the submission, then calls `apptainer run` (rather than `docker run`) on our image. docker, even when run in priveleged mode, gives me problems on CERN htcondor, but apptainer (/singularity) appears to work just fine. Then, move the outputs to a directory, tar them up, and copy them out.

## htcondor submission file
Finally, you need an htcondor `.sub` file. Something like:
```
universe                = vanilla

executable              = run_genie_production_apptainer.sh
arguments               = -n 20000 -p $(ProcId) -c $(ClusterId) -t Ti -T G18_02a_02_11b -e 8 -o production_7Nov2023
output                  = logs/ldmx_genie_prod.$(ClusterId).$(ProcId).out
error                   = logs/ldmx_genie_prod.$(ClusterId).$(ProcId).err
log                     = logs/ldmx_genie_prod.$(ClusterId).$(ProcId).log

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
output_destination      = root://eosuser.cern.ch//eos/user/w/wketchum/ldmx_en_samples/production_7Nov2023/
MY.XRDCP_CREATE_DIR     = True

request_cpus   = 1
#request_memory = 2000M
#request_disk   = 10240M

#+MaxRuntime   = 7200
+JobFlavour  = "workday"

should_transfer_files = yes
queue 100
```
See [https://batchdocs.web.cern.ch/local/submit.html] for a lot of the details. And I won't claim that there's not a better way to do this that handles automated transfers and such better.

But, this is likely straightforward or can be found in the docs. 
`condor_submit your_sub_file.sub` will do the submission...obviously it's good to do a test submission first.
