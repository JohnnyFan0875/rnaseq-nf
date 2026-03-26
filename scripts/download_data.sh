TARGET="SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 SRR1039520 SRR1039521"

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP033351&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
| tail -n +2 \
| while IFS=$'\t' read -r srr paths; do
    if echo "$TARGET" | grep -qw "$srr"; then
        echo "$paths" | tr ';' '\n' | grep -E "_[12]\.fastq\.gz$" | while read -r path; do
            wget "ftp://${path}"
        done
    fi
done