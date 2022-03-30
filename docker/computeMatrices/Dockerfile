# Dockerfile for computing the Geneid parameter file matrices
#
FROM debian:stretch

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    curl \
    ca-certificates \
    build-essential \
    gawk

ADD MarkovMatrices.awk scripts/
ADD MarkovMatrices-noframe.awk scripts/
ADD pro2log_ini.awk scripts/
ADD pro2log_tran.awk scripts/

ADD ./FastaToTbl /usr/local/bin/
RUN chmod +x /usr/local/bin/FastaToTbl
ENTRYPOINT [ "/usr/local/bin/FastaToTbl" ]
