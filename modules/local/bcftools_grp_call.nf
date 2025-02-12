process BCFTOOLS_GRP_CALL {
    tag "grouped call"
    label 'process_medium'

    publishDir "$params.outdir/variants", pattern: '*.vcf.gz', mode: 'copy'

    input:
    path(bam)
    path(bai)
    path fasta
    path fai
    tuple val(region), path(bed)

    output:
    path("*.vcf.gz") , emit: vcf
    path "versions.yml"        , emit: versions
    
    script:
      // consider adding in -B option when using a SNP list
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$params.outPrefix"
    """
    bcftools mpileup \\
        --output-type u \\
        --fasta-ref ${fasta} \\
        --max-depth 40000 \\
        --regions-file ${bed} \\
        --annotate AD,DP \\
        ${bam} | bcftools call \\
            $args \\
            --multiallelic-caller \\
            --output-type z \\
            --targets-file ${bed} \\
            --group-samples - \\
            --skip-variants indels \\
            --output ${prefix}_${region}.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
