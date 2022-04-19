# Dockerfile for computing the Geneid parameter file matrices
#
FROM debian:stretch

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    curl \
    ca-certificates \
    build-essential \
    procps \
    gawk \
    emboss

ADD ./FastaToTbl /usr/local/bin/
RUN chmod +x /usr/local/bin/FastaToTbl
ENTRYPOINT [ "/usr/local/bin/FastaToTbl" ]

ADD ./TblToFasta /usr/local/bin/
RUN chmod +x /usr/local/bin/TblToFasta
ENTRYPOINT [ "/usr/local/bin/TblToFasta" ]

