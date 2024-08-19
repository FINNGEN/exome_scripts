FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8

RUN apt-get update && apt-get upgrade --yes
RUN apt-get install -qqy python3 python3-venv python3-pip && pip3 install numpy

#Latest plink2 version
RUN mkdir plink2 && cd plink2 && curl -O s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip && \
    unzip plink2_linux_x86_64_latest.zip && mv plink2 /usr/local/bin && rm * && \
    curl -O http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip && \
    unzip plink2_linux_avx2_latest.zip && mv plink2 /usr/local/bin && cd .. && rm -rf plink2

ADD Scripts ../Scripts