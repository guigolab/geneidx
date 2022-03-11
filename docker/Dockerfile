# Dockerfile for geneid
#
FROM guigolab/geneid:1.4.5 as geneidbin
FROM debian:stretch

COPY --from=geneidbin /usr/local/bin/geneid /usr/local/bin/

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    curl \
    ca-certificates \
    build-essential \
    exonerate \
    perl \
    procps

RUN curl -L https://github.com/guigolab/blast2gff/archive/refs/heads/master.tar.gz | \
    tar xz && \
    cd blast2gff-master && \
    make && \
    cp ./bin/blast2gff /usr/local/bin/

COPY ./sgp_getHSPSR.pl /usr/local/bin/
RUN chmod +x /usr/local/bin/sgp_getHSPSR.pl


ENTRYPOINT [ "/usr/local/bin/geneid" ]
ENTRYPOINT [ "/usr/local/bin/blast2gff" ]
ENTRYPOINT [ "/usr/local/bin/sgp_getHSPSR.pl" ]
