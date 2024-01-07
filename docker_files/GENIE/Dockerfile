FROM ubuntu:22.04
MAINTAINER Wesley Ketchum <wketchum@fnal.gov>

ARG NPROC=1

ARG ROOT_VERSION=6.26.06
ARG CMAKE_CXX_STANDARD=17
ARG LHAPDF_VERSION=6.5.2

ARG GENIE_VERSION=3_02_00
ARG GENIE_REWEIGHT_VERSION=1_02_00

ENV GENIE_BASE=/usr/local/GENIE
ENV GENIE=$GENIE_BASE/Generator
ENV GENIE_REWEIGHT=$GENIE_BASE/Reweight

ENV PYTHIA6_LIBRARY=/usr/pythia6
ENV ROOTSYS=/usr/local/root
#ENV LANG=C.UTF-8

COPY packages packages

RUN apt-get update -qq && \
    ln -sf /usr/share/zoneinfo/UTC /etc/localtime && \
    apt-get -y install $(cat packages) wget git && \
    apt-get autoremove -y && \
    apt-get clean -y && \
    rm -rf /var/cache/apt/archives/* && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove --purge &&\
    apt-get clean all &&\
    python3 -m pip install --upgrade pip &&\
    python3 -m pip install --upgrade --no-cache-dir cmake
#RUN yes | unminimize

RUN mkdir /tmp/source && \
    cd /tmp/source && \
    wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-$LHAPDF_VERSION.tar.gz -O LHAPDF.tar.gz&& \
    tar xf LHAPDF.tar.gz && \
    cd LHAPDF-$LHAPDF_VERSION && \
    ./configure --prefix=/usr && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /tmp/source

#RUN mkdir /tmp/source && \
#    cd /tmp/source && \
#    wget https://sourceforge.net/projects/log4cpp/files/latest/download -O log4cpp.tar.gz && \
#    tar xf log4cpp.tar.gz &&\
#    cd log4cpp &&\
#    ./autogen.sh &&\
#    ./configure &&\
#    make check &&\
#    make install -j$NPROC &&\
#    rm -rf /tmp/source

RUN apt-get update -qq && \
    apt-get -y install liblog4cpp5-dev

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

RUN mkdir -p /tmp/source && \
    cd /tmp && \
    wget https://root.cern/download/root_v$ROOT_VERSION.source.tar.gz -O root.tar.gz && \
    tar xf root.tar.gz --strip-components=1 --directory source && \
    cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_INSTALL_PREFIX=$ROOTSYS \
      -Dgminimal=ON \
      -Dasimage=ON \
      -Dgdml=ON \
      -Dopengl=ON \
      -Dpyroot=ON \
      -Dgnuinstall=ON \
      -Dxrootd=OFF \
      -Dgsl_shared=ON \ 
      -Dmathmore=ON \   
      -Dpythia6=ON \    
      -DPYTHIA6_LIBRARY=$PYTHIA6_LIBRARY/libPythia6.so \
      -B build \
      -S source && \ 
    cmake --build build --target install -j$NPROC && \
    rm -rf build source

RUN cd /usr/local/root/lib &&\
    ln -s /usr/local/root/lib/root/* . &&\
    mkdir -p $GENIE &&\
    cd $GENIE_BASE &&\
    wget -q -O - https://github.com/GENIE-MC/Generator/archive/refs/tags/R-${GENIE_VERSION}.tar.gz |\
    tar -xz -C Generator --strip-components 1 &&\
    cd Generator &&\
    export LD_LIBRARY_PATH=/usr/local/lib:$ROOTSYS/lib/root &&\
    export PATH=$PATH:$ROOTSYS/bin &&\
    ./configure --enable-lhapdf6 --disable-lhapdf5 \
                --enable-gfortran --with-gfortran-lib=/usr/x86_64-linux-gnu/ \
                --disable-pythia8 --with-pythia6-lib=$PYTHIA6_LIBRARY \
                --enable-test && \
    make -j$NPROC && \
    make -j$NPROC install

