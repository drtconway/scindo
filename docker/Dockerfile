FROM ubuntu:groovy

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y \
        cmake \
        g++ \
        git \
        pkg-config \
        samtools \
        catch2 \
        nlohmann-json3-dev \
        libboost-log-dev \
        libdocopt-dev \
        libhts-dev \
        libz-dev

RUN git clone https://github.com/taocpp/PEGTL.git && \
    cd PEGTL && \
    git checkout 3.2.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j4 && \
    make install && \
    cd ../.. && \
    rm -rf PEGTL

RUN git clone https://github.com/simongog/sdsl-lite.git && \
    cd sdsl-lite && \
    ./install.sh /usr/local && \
    cd .. && \
    rm -rf sdsl-lite
