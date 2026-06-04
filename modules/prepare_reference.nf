process download_reference_files {

    tag { ref_dir }

    input:
    val ref_dir
    val gtf_path
    val genome_fasta_path
    val gtf_url
    val genome_fasta_url
    val auto_download

    output:
    path "resolved_reference/reference.gtf", emit: gtf
    path "resolved_reference/genome.fa",     emit: genome_fasta

    script:
    """
    set -eu

    mkdir -p "${ref_dir}"
    mkdir -p "\$(dirname "${gtf_path}")"
    mkdir -p "\$(dirname "${genome_fasta_path}")"

    if [ ! -f "${gtf_path}" ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing GTF file and reference_auto_download=false: ${gtf_path}" >&2
            exit 1
        fi

        tmp_gtf="${gtf_path}.gz"
        wget -O "\$tmp_gtf" "${gtf_url}"
        gunzip -c "\$tmp_gtf" > "${gtf_path}"
        rm -f "\$tmp_gtf"
    fi

    if [ ! -f "${genome_fasta_path}" ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing genome FASTA file and reference_auto_download=false: ${genome_fasta_path}" >&2
            exit 1
        fi

        tmp_fasta="${genome_fasta_path}.gz"
        wget -O "\$tmp_fasta" "${genome_fasta_url}"
        gunzip -c "\$tmp_fasta" > "${genome_fasta_path}"
        rm -f "\$tmp_fasta"
    fi

    mkdir -p resolved_reference
    ln -sf "${gtf_path}" resolved_reference/reference.gtf
    ln -sf "${genome_fasta_path}" resolved_reference/genome.fa
    """
}

process prepare_tx2gene {

    tag { tx2gene_path }

    input:
    val tx2gene_path
    val auto_download
    path script_dir

    output:
    path "resolved_reference/tx2gene.tsv", emit: tx2gene

    script:
    """
    set -eu

    mkdir -p "\$(dirname "${tx2gene_path}")"

    if [ ! -f "${tx2gene_path}" ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing tx2gene file and reference_auto_download=false: ${tx2gene_path}" >&2
            exit 1
        fi

        Rscript ${script_dir}/tx2gene.R --output "${tx2gene_path}"
    fi

    if [ ! -f "${tx2gene_path}" ]; then
        echo "ERROR: tx2gene file was not created: ${tx2gene_path}" >&2
        exit 1
    fi

    mkdir -p resolved_reference
    ln -sf "${tx2gene_path}" resolved_reference/tx2gene.tsv
    """
}

process prepare_transcript_fasta {

    tag { transcripts_fasta_path }

    input:
    path gtf
    path genome_fasta
    val transcripts_fasta_path
    val auto_download

    output:
    path "resolved_reference/transcripts.fa", emit: transcripts_fasta

    script:
    """
    set -eu

    mkdir -p "\$(dirname "${transcripts_fasta_path}")"

    if [ ! -f "${transcripts_fasta_path}" ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing transcripts FASTA and reference_auto_download=false: ${transcripts_fasta_path}" >&2
            exit 1
        fi

        gffread -w "${transcripts_fasta_path}" -g "${genome_fasta}" "${gtf}"
    fi

    if [ ! -f "${transcripts_fasta_path}" ]; then
        echo "ERROR: Transcript FASTA was not created: ${transcripts_fasta_path}" >&2
        exit 1
    fi

    mkdir -p resolved_reference
    ln -sf "${transcripts_fasta_path}" resolved_reference/transcripts.fa
    """
}

process prepare_kallisto_index {

    tag { kallisto_index_path }

    input:
    path transcripts_fasta
    val kallisto_index_path
    val auto_download

    output:
    path "resolved_reference/kallisto_index.idx", emit: kallisto_index

    script:
    """
    set -eu

    mkdir -p "\$(dirname "${kallisto_index_path}")"

    if [ ! -f "${kallisto_index_path}" ]; then
        if [ "${auto_download}" != "true" ]; then
            echo "ERROR: Missing kallisto index and reference_auto_download=false: ${kallisto_index_path}" >&2
            exit 1
        fi

        kallisto index -i "${kallisto_index_path}" "${transcripts_fasta}"
    fi

    if [ ! -f "${kallisto_index_path}" ]; then
        echo "ERROR: Kallisto index was not created: ${kallisto_index_path}" >&2
        exit 1
    fi

    mkdir -p resolved_reference
    ln -sf "${kallisto_index_path}" resolved_reference/kallisto_index.idx
    """
}

workflow prepare_reference {

    take:
    ref_dir
    gtf_path
    genome_fasta_path
    transcripts_fasta_path
    kallisto_index_path
    tx2gene_path
    gtf_url
    genome_fasta_url
    auto_download
    script_dir

    main:
    downloaded_reference = download_reference_files(
        ref_dir,
        gtf_path,
        genome_fasta_path,
        gtf_url,
        genome_fasta_url,
        auto_download
    )

    generated_tx2gene = prepare_tx2gene(
        tx2gene_path,
        auto_download,
        script_dir
    )

    generated_transcripts = prepare_transcript_fasta(
        downloaded_reference.gtf,
        downloaded_reference.genome_fasta,
        transcripts_fasta_path,
        auto_download
    )

    generated_index = prepare_kallisto_index(
        generated_transcripts.transcripts_fasta,
        kallisto_index_path,
        auto_download
    )

    emit:
    kallisto_index = generated_index.kallisto_index
    tx2gene = generated_tx2gene.tx2gene
}
