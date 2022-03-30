# Dockerfile for blast2gff
#
FROM debian:stretch

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    curl \
    ca-certificates \
    build-essential \
    perl \
    procps

RUN curl -L https://github.com/guigolab/blast2gff/archive/refs/heads/master.tar.gz | \
    tar xz && \
    cd blast2gff-master && \
    make && \
    cp ./bin/blast2gff /usr/local/bin/


ENTRYPOINT [ "/usr/local/bin/blast2gff" ]
