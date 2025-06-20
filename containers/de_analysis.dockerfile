FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install R packages
RUN R -e " \
    pkgs <- c('tximport', 'DESeq2', 'edgeR', 'openxlsx', 'tidyverse', 'ComplexHeatmap', 'RColorBrewer', 'circlize', 'jsonlite', 'optparse', 'rhdf5', 'AnnotationDbi', 'org.Hs.eg.db', 'matrixStats', 'EnhancedVolcano'); \
    BiocManager::install(pkgs, ask=FALSE, update=FALSE); \
    "

# Set working directory
WORKDIR /opt/
