process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    input:
    path(fasta)

    output:
    path(bwa) , emit: index
    path "versions.yml"        , emit: versions
    
    script:
    """
    mkdir bwa
    bwa \\
        index \\
        -p bwa/\$(basename ${fasta}) \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
