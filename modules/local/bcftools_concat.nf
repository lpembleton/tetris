process BCFTOOLS_CONCAT {
    tag "concat vcf"
    label 'process_medium'

    publishDir "$params.outdir/variants", pattern: '*.vcf.gz', mode: 'copy'

    input:
    path(vcf)

    output:
    path("*.vcf.gz") , emit: vcf
    path "versions.yml"        , emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$params.outPrefix"
    """
	bcftools concat *vcf.gz \\
		--output-type z \\
        --output ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
