process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low_multi'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path("*.stats")     , emit: stats
    path "versions.yml" , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        stats \\
        --threads 2 \\
        ${bam} \\
        > ${prefix}_samtools.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
