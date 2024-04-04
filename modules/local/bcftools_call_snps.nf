process BCFTOOLS_CALL_SNPS {
    tag "calling $meta.name"
    label 'process_medium'

    publishDir "$params.outdir/variants", pattern: '*.vcf.gz', mode: 'copy'

    input:
    tuple val(meta), path(bcf)
    tuple val(meta), path(csi)
    path region
    path target_index

    output:
    tuple val(meta), path("*.vcf.gz") , emit: vcf
    path "versions.yml"        , emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.name}"
    """
    bcftools call \\
        $args \\
        --multiallelic-caller \\
        --output-type z \\
        --targets-file ${region} \\
        --constrain alleles \\
        --group-samples - \\
        --skip-variants indels \\
        --output ${prefix}.vcf.gz \\
        ${bcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
