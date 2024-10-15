process BWA_MEM {
    tag "mapping $meta.id"
    label 'process_medium_multi'

    if( params.split_fastq == 0 && params.skip_markdup == true) {
        publishDir "$params.outdir/mapping", pattern: '*.bam', mode: 'copy'
    }    

    input:
    tuple val(meta), path(reads)
    path(index)
    val(sort_bam)
    
    output:
    tuple val(meta), path("*.bam")              , emit: bam
    tuple val(meta), path("*bai"), optional: true   , emit: bai
    path  "versions.yml"                        , emit: versions


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    if (params.dedup == true) {
        """
            INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

            bwa mem \\
                $args \\
                -t 4 \\
                \$INDEX \\
                $reads \\
                | samtools $samtools_command $args2 --threads 4 -o ${prefix}.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
        """
    } else{
        """
            INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

            bwa mem \\
                $args \\
                -t 4 \\
                \$INDEX \\
                $reads \\
                | samtools $samtools_command $args2 --threads 4 -o ${prefix}.bam
            
            samtools \\
                index \\
                -@ 1 \\
                ${prefix}.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
                samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            END_VERSIONS
        """
    }
}
// todo: add in $task.cpus as a variable