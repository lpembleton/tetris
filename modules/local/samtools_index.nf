process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    if( params.skip_markdup == 'true' ) {
        publishDir "$params.outdir/mapping", pattern: '*.bai', mode: 'copy'
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai") , emit: bai
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}