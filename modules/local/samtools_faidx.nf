process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    input:
    path(fasta)

    output:
    path ("*.fai")          , emit: fai
    path "versions.yml"     , emit: versions

    script:
    """
    samtools \\
        faidx \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
