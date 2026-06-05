process prepare_reference_assets {

    tag { ref_dir }

    conda "${moduleDir}/environment.yml"

    input:
    val ref_dir
    path existing_reference_files
    val gtf_path
    val genome_fasta_path
    val transcripts_fasta_path
    val kallisto_index_path
    val gtf_url
    val genome_fasta_url
    val auto_download

    output:
    path "Homo_sapiens.GRCh38.113.gtf", emit: gtf
    path "Homo_sapiens.GRCh38.dna.primary_assembly.fa", emit: genome_fasta
    path "transcripts.GRCh38.113.fa", emit: transcripts_fasta
    path "transcripts.GRCh38.113.kallisto.idx", emit: kallisto_index

    script:
    """
    set -eu

    tmp_gtf="Homo_sapiens.GRCh38.113.gtf.gz"
    tmp_fasta="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

    if [ ! -f Homo_sapiens.GRCh38.113.gtf ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing GTF and reference_auto_download=false." >&2
            echo "  GTF: ${gtf_path}" >&2
            exit 1
        fi

        wget -O "\$tmp_gtf" "${gtf_url}"
        gunzip -c "\$tmp_gtf" > Homo_sapiens.GRCh38.113.gtf
        rm -f "\$tmp_gtf"
    fi

    if [ ! -f Homo_sapiens.GRCh38.dna.primary_assembly.fa ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing genome FASTA and reference_auto_download=false." >&2
            echo "  Genome FASTA: ${genome_fasta_path}" >&2
            exit 1
        fi

        wget -O "\$tmp_fasta" "${genome_fasta_url}"
        gunzip -c "\$tmp_fasta" > Homo_sapiens.GRCh38.dna.primary_assembly.fa
        rm -f "\$tmp_fasta"
    fi

    if [ ! -f transcripts.GRCh38.113.fa ]; then
        gffread -w transcripts.GRCh38.113.fa -g Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.113.gtf
    fi

    if [ ! -f transcripts.GRCh38.113.kallisto.idx ]; then
        kallisto index -i transcripts.GRCh38.113.kallisto.idx transcripts.GRCh38.113.fa
    fi

    if [ ! -f Homo_sapiens.GRCh38.113.gtf ] || [ ! -f Homo_sapiens.GRCh38.dna.primary_assembly.fa ] || [ ! -f transcripts.GRCh38.113.fa ] || [ ! -f transcripts.GRCh38.113.kallisto.idx ]; then
        echo "ERROR: Failed to materialize one or more reference assets." >&2
        exit 1
    fi
    """

    stub:
    """
    touch Homo_sapiens.GRCh38.113.gtf
    touch Homo_sapiens.GRCh38.dna.primary_assembly.fa
    touch transcripts.GRCh38.113.fa
    touch transcripts.GRCh38.113.kallisto.idx
    """
}
