process BCFTOOLS_MPILEUP {
    tag "mpileup $meta.name"
    label 'process_medium'

    publishDir "$params.outdir/mpileup", pattern: '*.bcf.gz', mode: 'copy'
    publishDir "$params.outdir/mpileup", pattern: '*.bcf.gz.csi', mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path region

    output:
    tuple val(meta), path("*.bcf.gz") , emit: bcf
    tuple val(meta), path("*.bcf.gz.csi") , emit: csi
    path "versions.yml"        , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.name}"
    // consider adding in -B option when using a SNP list
    """
    bcftools mpileup \\
        --output-type z \\
        --fasta-ref ${fasta} \\
        --output ${prefix}.bcf.gz \\
        --max-depth 40000 \\
        --regions-file ${region} \\
        --annotate AD,DP \\
        ${bam}

    bcftools index \\
        ${prefix}.bcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
