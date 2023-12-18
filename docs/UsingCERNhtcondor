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
