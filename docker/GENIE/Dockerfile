FROM ubuntu:22.04

ARG ROOT_VERSION=6.26.06
ARG CMAKE_CXX_STANDARD=17
ARG LHAPDF_VERSION=6.5.2

ARG GENIE_VERSION=3_02_00
ARG GENIE_REWEIGHT_VERSION=1_02_00

ENV GENIE_BASE=/usr/local/GENIE
ENV GENIE=$GENIE_BASE/Generator
ENV GENIE_REWEIGHT=$GENIE_BASE/Reweight

ENV PYTHIA6_LIBRARY=/usr/pythia6
ENV ROOTSYS=/usr/local


# Start by stealing from ROOT team
#LABEL maintainer.name="ROOT team"
#LABEL maintainer.email="root-dev@cern.ch"

ENV LANG=C.UTF-8

COPY packages_ROOT packages_ROOT
COPY packages_GENIE packages_GENIE

RUN apt-get update -qq && \
    ln -sf /usr/share/zoneinfo/UTC /etc/localtime && \
    apt-get -y install $(cat packages_ROOT) wget git && \
    apt-get -y install $(cat packages_GENIE) && \
    apt-get autoremove -y && \
    apt-get clean -y && \
    rm -rf /var/cache/apt/archives/* && \
    rm -rf /var/lib/apt/lists/*

RUN yes | unminimize

# clone build and install

RUN mkdir /tmp/source && \
    cd /tmp/source && \
    wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz -O LHAPDF.tar.gz&& \
    tar xf LHAPDF.tar.gz && \
    cd LHAPDF-$LHAPDF_VERSION && \
    ./configure --prefix=/usr && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /tmp/source


RUN mkdir /usr/pythia6 && \
    wget https://root.cern.ch/download/pythia6.tar.gz && \
    tar zxvf pythia6.tar.gz -C /usr/ && rm -rf pythia6.tar.gz && \
    wget --no-check-certificate https://pythia.org/download/pythia6/pythia6428.f && \
    mv pythia6428.f /usr/pythia6/pythia6428.f && rm -rf /usr/pythia6/pythia6416.f && \
    cd /usr/pythia6 && \
    sed -i 's/int py/extern int py/g' pythia6_common_address.c && \
    sed -i 's/extern int pyuppr/int pyuppr/g' pythia6_common_address.c && \
    sed -i 's/char py/extern char py/g' pythia6_common_address.c && \
    cat pythia6_common_address.c && \
    echo 'void MAIN__() {}' >main.c && \
    gcc -c -m64 -fPIC -shared main.c -lgfortran && \
    gcc -c -m64 -fPIC -shared pythia6_common_address.c -lgfortran && \
    gfortran -c -m64 -fPIC -shared pythia*.f && \
    gfortran -c -m64 -fPIC -shared -fno-second-underscore tpythia6_called_from_cc.F && \
    gfortran -m64 -shared -Wl,-soname,libPythia6.so -o libPythia6.so main.o  pythia*.o tpythia*.o


RUN ROOT_GIT_URL=https://github.com/root-project/root.git \
    && if [ -z "$(git ls-remote --heads $ROOT_GIT_URL $ROOT_VERSION)" ] ; then \
    export ROOT_GIT_REVISION="v$(echo $ROOT_VERSION | cut -d. -f1)-$(echo $ROOT_VERSION | cut -d. -f2)-$(echo $ROOT_VERSION | cut -d. -f3)" \
    ; else \
    export ROOT_GIT_REVISION=$ROOT_VERSION \
    ; fi \
    # Above lines will set ROOT_GIT_REVISION to ROOT_VERSION argument if it corresponds to a valid branch name (such as ROOT_VERSION=master),
    # otherwise it will assume its a semantic version string and try to convert it into the tag format (such as 6.26.00 -> v6-26-00)
    && git clone --branch $ROOT_GIT_REVISION --depth=1 $ROOT_GIT_URL /tmp/source \
    && cd /tmp/source \
    && mkdir -p /tmp/build &&  cd /tmp/build \
    && cmake /tmp/source \
    -DCMAKE_CXX_STANDARD=$CMAKE_CXX_STANDARD \
    -Dgnuinstall=OFF \
    -Dsoversion=ON \
    -Druntime_cxxmodules=OFF \
    -Dgsl_shared=ON \
    -Dpythia6=ON \
    -Dmathmore=ON \
    # For ROOT version 6.26.00 it is necessary to set `-Druntime_cxxmodules=OFF` (https://github.com/root-project/root/pull/10198)
    -DCMAKE_INSTALL_PREFIX=$ROOTSYS \
    -DPYTHIA6_LIBRARY=$PYTHIA6_LIBRARY/libPythia6.so \
    && make -j$(nproc) install \
    && rm -rf /tmp/build /tmp/source

RUN /bin/bash -c 'source $ROOTSYS/bin/thisroot.sh'

RUN mkdir -p $GENIE_BASE

RUN git clone --branch R-$GENIE_VERSION https://github.com/GENIE-MC/Generator.git $GENIE_BASE/Generator

RUN cd $GENIE && \
    ./configure --enable-lhapdf6 --disable-lhapdf5 --enable-gfortran --with-gfortran-lib=/usr/x86_64-linux-gnu/ --disable-pythia8 --with-pythia6-lib=$PYTHIA6_LIBRARY && \
    gmake -j$(nproc)

#ENV GENIE_LIB_PATH=$GENIE/lib
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENIE_LIB_PATH
ENV PATH=$PATH:$GENIE/bin


#RUN echo $LD_LIBRARY_PATH
#RUN git clone --branch R-$GENIE_REWEIGHT_VERSION https://github.com/GENIE-MC/Reweight.git $GENIE_BASE/Reweight
#RUN echo $GENIE_REWEIGHT
#RUN git clone --branch wketchum_GHepParticleToPointer https://github.com/wesketchum/Reweight.git $GENIE_BASE/Reweight

#RUN cd $GENIE_REWEIGHT && \
#    gmake

#RUN ls $GENIE_REWEIGHT

