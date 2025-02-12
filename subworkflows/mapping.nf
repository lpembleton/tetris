/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    READ MAPPING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP } from '../modules/local/fastp'
include { FASTQC } from '../modules/local/fastqc'
include { SAMTOOLS_FAIDX } from '../modules/local/samtools_faidx'
include { SAMTOOLS_INDEX } from '../modules/local/samtools_index'
include { SAMTOOLS_STATS } from '../modules/local/samtools_stats'
include { SAMTOOLS_MERGE } from '../modules/local/samtools_merge'
include { BWA_INDEX } from '../modules/local/bwa_index'
include { BWA_MEM } from '../modules/local/bwa_mem'
include { MOSDEPTH } from '../modules/local/mosdepth'
include { GATK4_MARKDUPLICATES } from '../modules/local/gatk_markduplicates'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow MAPPING {

	take:
		input_sample
		reference
		skip_fastp
        split_fastq
		skip_markdup
		markdup_tool

	main:

    ch_versions = Channel.empty()
    ch_reports = Channel.empty()
    
    // ==== P R E P R O C E S S I N G   S T E P S ====

	if (skip_fastp) {
		//skip fastp
		post_fp_reads = input_sample
	} else {
		FASTP(input_sample)
		post_fp_reads = FASTP.out.reads
		// Gather used softwares versions and reports
		ch_versions = ch_versions.mix(FASTP.out.versions)
		ch_reports = ch_reports.mix(FASTP.out.json.collect{ meta, json -> json })
		ch_reports = ch_reports.mix(FASTP.out.html.collect{ meta, html -> html })
	}


    // TODO: code this to also work with single read data
    if (split_fastq > 0) {
            reads_for_alignment = post_fp_reads.map{ meta, reads ->
                read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                [ meta + [ size:read_files.size() ], read_files ]
            }.transpose()
    } else {
        reads_for_alignment = post_fp_reads
    }
    //reads_for_alignment.view()

    // QC on trimmed reads
    FASTQC(reads_for_alignment)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })

    // ==== M A P P I N G   S T E P S ====

    // Prepare reference fasta index
    //SAMTOOLS_FAIDX(reference)
    BWA_INDEX(reference)
    // Gather used softwares versions and reports
    //ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)


    // Map reads to reference fasta
    sort = true // sort after mapping
    BWA_MEM(reads_for_alignment, BWA_INDEX.out.index, sort)
    // Gather used softwares versions and reports
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    if (split_fastq > 0) { // if reads were split merge them back together before next steps
        SAMTOOLS_MERGE(BWA_MEM.out.bam.groupTuple())
        aligned_bams_ch = SAMTOOLS_MERGE.out.bam  
    } else {
        aligned_bams_ch = BWA_MEM.out.bam
    }

    if (params.skip_markdup) { // skip duplicate read marking (e.g. RE-GBS or amplicon seq)

        SAMTOOLS_INDEX(aligned_bams_ch)
        //join the bams channel with the newly created bai index channel using a matching key, meta, with the join() command
        //https://www.nextflow.io/docs/latest/operator.html#join
        indexed_bams_ch = aligned_bams_ch.join(SAMTOOLS_INDEX.out.bai)

    } else {
        if (params.markdup_tool == 'gatk') { // mark duplicate reads using gatk4 (current default tooL)
            GATK4_MARKDUPLICATES(aligned_bams_ch)
            //GATK markduplicates will do the indexexing, so no need to call samtools_index after
            indexed_bams_ch = GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai)
        } 
        //if (param.markdup_samtools) { //TODO: add in tooling for samtools duplicate read marking
        
        //}
    }

    // Get bam stats
    SAMTOOLS_STATS(indexed_bams_ch)
    // Gather used softwares versions and reports
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats)


    // WGS coverage stats
    MOSDEPTH(indexed_bams_ch)
    // Gather used softwares versions and reports
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)    
    ch_reports = ch_reports.mix(MOSDEPTH.out.global_txt)
    ch_reports = ch_reports.mix(MOSDEPTH.out.summary_txt)

	emit:
		bam = indexed_bams_ch
		versions = ch_versions
		reports = ch_reports

}
