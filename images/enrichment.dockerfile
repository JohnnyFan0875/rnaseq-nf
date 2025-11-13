# Use Bioconductor base image with R 4.3 and Bioconductor 3.20
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install specific R packages with versions
RUN R -e "install.packages(c('optparse', 'jsonlite', 'tidyverse', 'data.table', 'ggridges'), repos='https://cloud.r-project.org/')"

RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"

RUN R -e "BiocManager::install(c('clusterProfiler', 'enrichplot', 'org.Hs.eg.db', 'ReactomePA', 'DOSE', 'pathview'), version='3.20', ask=FALSE)"

# Create working directory
WORKDIR /opt/