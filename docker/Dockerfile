FROM ubuntu:groovy AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y \
        cmake \
        g++ \
        gdb \
        git \
        pkg-config \
        samtools \
        catch2 \
        nlohmann-json3-dev \
        libboost-iostreams-dev \
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

ADD . /source
ENV CMAKE_MODULE_PATH=/source/cmake/
RUN cd /source && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make && \
    make install

FROM ubuntu:groovy
RUN apt update && \
    apt install -y \
        libboost-iostreams1.71.0 \
        libboost-log1.71.0 \
        libdocopt0 \
        libhts3 && \
    rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/bin/scindo /usr/local/bin/
COPY --from=builder /usr/local/bin/hecil /usr/local/bin/
