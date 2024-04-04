process MOSDEPTH {
    tag "$meta.id"
    label 'process_low_multi'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path('*.global.dist.txt')   , emit: global_txt
    path('*.summary.txt')       , emit: summary_txt
    path  "versions.yml"        , emit: versions    

    script:
    """
    mosdepth -n --fast-mode --by 500 $meta.id $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS        
    """
    
}