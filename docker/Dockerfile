FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

RUN apt-get update && apt-get upgrade --yes
RUN apt-get install -qqy python3 python3-venv python3-pip parallel pv pigz zstd && \
    pip3 install numpy pandas matplotlib tqdm

#Latest plink2 version
RUN mkdir plink2 && cd plink2 && curl -O s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip && \
    unzip plink2_linux_x86_64_latest.zip && mv plink2 /usr/local/bin && rm * && \
    curl -O http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip && \
    unzip plink2_linux_avx2_latest.zip && mv plink2 /usr/local/bin && cd .. && rm -rf plink2

#KING
RUN wget https://www.kingrelatedness.com/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz -C /usr/local/bin && rm Linux*king*.tar.gz

# Latest bcftools (built from source)
RUN apt-get install -qqy --no-install-recommends \
        gcc make autoconf zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev
	
RUN BCFTOOLS_VERSION=1.23.1 &&  \
    curl -fsSL "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" | tar -xj -C /tmp && \
    cd "/tmp/bcftools-${BCFTOOLS_VERSION}" && \
    ./configure --prefix=/usr/local && \
    make -j"$(nproc)" && \
    make install && \
    rm -rf "/tmp/bcftools-${BCFTOOLS_VERSION}"

COPY scripts/ /scripts/

