FROM zavolab/gffread:0.11.7-slim AS gffread_stage

FROM quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0

USER root

ENV LD_LIBRARY_PATH=/usr/local/lib

COPY --from=gffread_stage /usr/bin/gffread /usr/bin/gffread

RUN gffread --help >/dev/null 2>&1 || test $? -eq 1

WORKDIR /opt/
