process AWK_SPLITBED {
    tag "${regions}"
    label 'process_single'

    input:
    path(regions)

    output:
    path('*.bed'), emit: subregion
    
    shell:
    """
   
    awk '{print \$0 > \$4".bed"}' ${regions}

    rm !{regions}

    """
}
