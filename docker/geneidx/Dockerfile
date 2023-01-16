# Dockerfile for GeneidX
#
FROM guigolab/geneid:1.4.5 as geneidbin

FROM python:3.9.14

ADD MarkovMatrices.awk scripts/
ADD MarkovMatrices-noframe.awk scripts/
ADD pro2log_ini.awk scripts/
ADD pro2log_tran.awk scripts/

ADD ./FastaToTbl /usr/local/bin/
ADD ./TblToFasta /usr/local/bin/


COPY --from=geneidbin /usr/local/bin/geneid /usr/local/bin/

COPY ./sgp_getHSPSR.pl /usr/local/bin/

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    curl \
    ca-certificates \
    build-essential \
    perl \
    gawk \
    exonerate \
    wget \
    procps

RUN python3 -m pip install --upgrade pip
RUN pip install requests
RUN pip install lxml
RUN pip install pandas
RUN pip install numpy

RUN curl -L https://github.com/guigolab/blast2gff/archive/refs/heads/master.tar.gz | \
    tar xz && \
    cd blast2gff-master && \
    make && \
    cp ./bin/blast2gff /usr/local/bin/

RUN chmod +x /usr/local/bin/FastaToTbl

RUN chmod +x /usr/local/bin/sgp_getHSPSR.pl

