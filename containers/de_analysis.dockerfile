FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install specific CRAN packages
RUN R -e "install.packages(c('optparse', 'jsonlite', 'tidyverse', 'openxlsx', 'RColorBrewer', 'matrixStats', 'circlize'), repos='https://cloud.r-project.org/')"

# Ensure BiocManager is installed
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"

# Install Bioconductor packages using BiocManager for release 3.20
RUN R -e "BiocManager::install(c('tximport', 'DESeq2', 'edgeR', 'rhdf5', 'AnnotationDbi', 'org.Hs.eg.db', 'EnhancedVolcano', 'ComplexHeatmap'), version='3.20', ask=FALSE, update=FALSE)"

# Set working directory
WORKDIR /opt/
