nextflow.enable.dsl = 2

/*
 * RNA-seq multi-omics integration pipeline
 *
 * Main steps:
 * 1. Download paired-end reads from SRA
 * 2. Trim reads using fastp
 * 3. Align reads using HISAT2
 * 4. Generate gene-level counts using featureCounts
 * 5. Annotate proteins using eggNOG-mapper
 * 6. Integrate expression, annotation, and optional metabolic model data
 */

process DOWNLOAD_FASTQ {

    tag "${srr_id}"
    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    val srr_id

    output:
    tuple val(srr_id), path("${srr_id}_1.fastq"), path("${srr_id}_2.fastq")

    script:
    """
    prefetch ${srr_id}
    fasterq-dump ${srr_id} --split-files --threads ${task.cpus}
    """
}

process FASTP {

    tag "${sample_id}"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_clean_1.fastq"), path("${sample_id}_clean_2.fastq")

    script:
    """
    fastp \
        -i ${read1} \
        -I ${read2} \
        -o ${sample_id}_clean_1.fastq \
        -O ${sample_id}_clean_2.fastq \
        --thread ${task.cpus}
    """
}

process HISAT2_ALIGN {

    tag "${sample_id}"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam")

    script:
    """
    hisat2 \
        -x ${params.index} \
        -1 ${read1} \
        -2 ${read2} \
        -p ${task.cpus} | \
    samtools view -@ ${task.cpus} -bS - > ${sample_id}.aligned.bam
    """
}

process FEATURECOUNTS {

    tag "${sample_id}"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.counts.txt"), path("${sample_id}.counts.txt.summary")

    script:
    """
    featureCounts \
        -p \
        -T ${task.cpus} \
        -a ${params.gtf} \
        -o ${sample_id}.counts.txt \
        ${bam}
    """
}

process EGGNOG {

    tag "eggnog_annotation"
    publishDir "${params.outdir}/eggnog", mode: 'copy'

    input:
    path protein

    output:
    path "eggnog_output.emapper.annotations"

    script:
    """
    emapper.py \
        -i ${protein} \
        --itype proteins \
        --output eggnog_output \
        --data_dir ${params.eggnog_db} \
        --cpu ${task.cpus} \
        --sensmode fast \
        -m diamond \
        --tax_scope Bacteria \
        --override
    """
}

process INTEGRATION {

    tag "multiomics_integration"
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path counts
    path eggnog

    output:
    path "final_complete_multiomics.csv"

    script:
    def metabolicArg = params.metabolic_model ? "--metabolic-model ${params.metabolic_model}" : ""

    """
    python3 ${projectDir}/bin/integrate.py \
        --counts ${counts} \
        --eggnog ${eggnog} \
        ${metabolicArg} \
        --output final_complete_multiomics.csv \
        --locus-prefix ${params.locus_prefix}
    """
}

workflow {

    protein_ch = Channel.fromPath(params.protein, checkIfExists: true)

    reads_ch  = DOWNLOAD_FASTQ(params.srr)
    clean_ch  = FASTP(reads_ch)
    bam_ch    = HISAT2_ALIGN(clean_ch)
    counts_ch = FEATURECOUNTS(bam_ch)

    eggnog_ch = EGGNOG(protein_ch)

    counts_for_integration = Channel.fromPath(params.counts, checkIfExists: true)

    INTEGRATION(
        counts_for_integration,
        eggnog_ch
    )
}
