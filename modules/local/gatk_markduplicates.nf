process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_low_multi'

    publishDir "$params.outdir/mapping", pattern: '*markdup.bam', mode: 'copy'
    publishDir "$params.outdir/mapping", pattern: '*markdup.bai', mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*markdup.bam"), emit: bam
    tuple val(meta), path("*markdup.bai"), emit: bai
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml",                emit: versions


    script:
    def avail_mem = 16
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.markdup.bam \\
        --METRICS_FILE ${prefix}.metrics \\
        --CREATE_INDEX true \\
        --TMP_DIR . \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --VALIDATION_STRINGENCY LENIENT \\
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 800 \\

    rm ${bam}

    samtools index ${prefix}.markdup.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS

    """
}