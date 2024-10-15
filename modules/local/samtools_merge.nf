process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    if( params.skip_markdup == true ) {
        publishDir "$params.outdir/mapping", pattern: '*.merged.bam', mode: 'copy'
    }   

    input:
    tuple val(meta), path(input_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${prefix}.merged.bam") , emit: bam
    path  "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        merge \\
        --threads 1 \\
        $args \\
        ${prefix}.merged.bam \\
        $input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}