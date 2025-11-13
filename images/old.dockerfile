# Use Ubuntu as the base image
FROM ubuntu:22.04

# Set environment variables
ENV R_VERSION=4.4.0 \
    DEBIAN_FRONTEND=noninteractive \
    PATH="/usr/local/miniconda/bin:/usr/local/bin:$PATH"

# Update and install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl \
    parallel \
    fastqc \
    multiqc \
    fastp \
    kallisto \
    less \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libagg-dev \
    libproj-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libopenblas-dev \
    python3-pip \
    ca-certificates \
    build-essential \
    gfortran \
    libreadline-dev \
    xorg-dev \
    libbz2-dev \
    liblzma-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libsodium-dev \
    libxt-dev \
    libpng-dev \
    wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Download and install R
#https://stackoverflow.com/questions/78749769/install-specific-version-of-r-langauge-r4-x-x-on-a-docker-image-based-on-ubunt
RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
    && tar -xf R-${R_VERSION}.tar.gz \
    && cd R-${R_VERSION} \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

# Install Miniconda and Snakemake 8.14.0
#RUN curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh && \
#    bash miniconda.sh -b -p /usr/local/miniconda && \
#    rm miniconda.sh && \
#    /usr/local/miniconda/bin/conda install -c conda-forge -c bioconda -y snakemake=8.14.0 && \
#    /usr/local/miniconda/bin/conda clean --all -y

# Install R libraries
RUN R -e "install.packages(c('BiocManager', 'openxlsx', 'readxl', 'data.table', 'ggplot2', 'cowplot', 'fastmatch', 'scales', 'BH', 'Cairo'), repos = 'https://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('fgsea', 'edgeR', 'EnhancedVolcano', 'DESeq2', 'flextable', 'ComplexHeatmap', 'tximport', 'tidyverse', 'cowplot', 'biomaRt', 'pathview', 'DOSE', 'enrichplot', 'clusterProfiler', 'ReactomePA', 'HDO.db'))"

# Install required Python packages
RUN pip install --no-cache-dir pandas pillow python-docx openpyxl

# Set working directory
WORKDIR /opt/

#cp necessary file
ADD rna_ref /opt/rna_ref
#ADD rna_script /opt/rna_script

RUN mkdir -p rna_input rna_output /root/.local/share/GOSemSim/ && \
    cp /opt/rna_ref/HDO.sqlite.gz /root/.local/share/GOSemSim/

#https://ask.csdn.net/questions/8134929

