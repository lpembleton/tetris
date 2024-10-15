/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
if (params.input) { csv_file = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
 
if (params.reference == null) error "Please specify a reference genome fasta file with --reference"
// demo reference can be found in docs/example_data

// TODO: need add more parameter checks

params.split_fastq = 0 // number of reads to split fastq files by (default: 0, none)
params.skip_fastp = false // T/F skip fastp processing [currently skipping is not supported]
params.skip_markdup = false // T/F skip duplicate read marking
params.markdup_tool = 'gatk' // use GATK markduplicates tool for duplicate read marking, samtools option to be added
params.grouped_call = false // T/F run mpileup and call on all bams together

//params.variants_only = false // see config for the interation with this flag

params.constrain_alleles = false // true is currently only supported in non grouped call
params.targets_index = null
// if constrain_alleles = true then a indexed set of target positions and alleles must be provided 
//  bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' file.vcf | bgzip -c > als.tsv.gz && tabix -s1 -b2 -e2 als.tsv.gz

// print parameters
log.info """\
    =======================================================================
    T E T R I S
    A   D N A   M A P P I N G   
    A N D   V A R I A N T   C A L L I N G   P I P E L I N E
    ======================================================================
    samplesheet: ${params.input}
    reference: ${params.reference}
    skip fastp: ${params.skip_fastp}
    split fastq (0 to skip): ${params.split_fastq}
    skip markduplicates: ${params.skip_markdup}
    markduplicate tool: ${params.markdup_tool}
    call variants only: ${params.variants_only}
    grouped call: ${params.grouped_call}
    constrain alleles: ${params.constrain_alleles}
    output directory: ${params.outdir}
    genomic regions: ${params.regions}
    output prefix: ${params.outPrefix}
    ======================================================================
    """
    .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
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
include { BCFTOOLS_MPILEUP } from '../modules/local/bcftools_mpileup'
include { BCFTOOLS_CALL } from '../modules/local/bcftools_call'
include { BCFTOOLS_CALL_SNPS } from '../modules/local/bcftools_call_snps'
include { BCFTOOLS_GRP_CALL } from '../modules/local/bcftools_grp_call'
include { GATK4_MARKDUPLICATES } from '../modules/local/gatk_markduplicates'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom/custom_dumpsoftwareversions'
include { MULTIQC } from '../modules/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TETRIS {
    // To gather all QC reports for MultiQC
    versions = Channel.empty()
    // To gather used softwares versions for MultiQC
    reports = Channel.empty()

    input_sample = extract_csv(file(csv_file))

    // ==== P R E P R O C E S S I N G   S T E P S ====

	if (params.skip_fastp) {
		//skip fastp
		post_fp_reads = input_sample
	} else {
		FASTP(input_sample)
		post_fp_reads = FASTP.out.reads
		// Gather used softwares versions and reports
		versions = versions.mix(FASTP.out.versions)
		reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
		reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })
	}


    // TODO: code this to also work with single read data
    if (params.split_fastq) {
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
    versions = versions.mix(FASTQC.out.versions)
    reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })

    // ==== M A P P I N G   S T E P S ====

    // Prepare reference fasta index
    SAMTOOLS_FAIDX(params.reference)
    BWA_INDEX(params.reference)
    // Gather used softwares versions and reports
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(BWA_INDEX.out.versions)


    // Map reads to reference fasta
    sort = true // sort after mapping
    BWA_MEM(reads_for_alignment, BWA_INDEX.out.index, sort)
    // Gather used softwares versions and reports
    versions = versions.mix(BWA_MEM.out.versions)

    if (params.split_fastq > 0) { // if reads were split merge them back together before next steps
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
    versions = versions.mix(SAMTOOLS_STATS.out.versions)
    reports = reports.mix(SAMTOOLS_STATS.out.stats)


    // WGS coverage stats
    MOSDEPTH(indexed_bams_ch)
    // Gather used softwares versions and reports
    versions = versions.mix(MOSDEPTH.out.versions)    
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.summary_txt)



    // ==== V A R I A N T   C A L L I N G ====
 
    if (params.grouped_call) { // perform variant calling on all samples together (e.g. variant discovery, mapping families)

        //BCFTOOLS_GRP_CALL(BWA_MEM.out.bam.map { tuple -> tuple[1] }.collect(), BWA_MEM.out.bai.map { tuple -> tuple[1] }.collect(), params.reference, SAMTOOLS_FAIDX.out.fai, params.regions)
        BCFTOOLS_GRP_CALL(indexed_bams_ch.map { tuple -> tuple[1] }.collect(), indexed_bams_ch.map { tuple -> tuple[2] }.collect(), params.reference, SAMTOOLS_FAIDX.out.fai, params.regions)
        // Gather used softwares versions and reports
        versions = versions.mix(BCFTOOLS_GRP_CALL.out.versions)  

    } else { // perform variant calling seperately for each sample (e.g. genotyping pipeline)

        // Groups multi bam from the same sample 'name' to run mpileup on
        //BCFTOOLS_MPILEUP(BWA_MEM.out.bam.map{ meta, bam -> [ meta.subMap('name'), bam ] }.groupTuple(), BWA_MEM.out.bai.map{ meta, bai -> [ meta.subMap('name'), bai ] }.groupTuple(), params.reference, SAMTOOLS_FAIDX.out.fai, params.regions)
        BCFTOOLS_MPILEUP(indexed_bams_ch.map{ meta, bam, bai -> [ meta.subMap('name'), bam, bai ] }.groupTuple(), params.reference, SAMTOOLS_FAIDX.out.fai, params.regions)
        // Gather used softwares versions and reports
        versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)  

        if(params.constrain_alleles) { // constrain genotype calls to those alleles provided an indexed als.tsv type file

            BCFTOOLS_CALL_SNPS(BCFTOOLS_MPILEUP.out.bcf, BCFTOOLS_MPILEUP.out.csi, params.regions, params.targets_index)
            // Gather used softwares versions and reports
            versions = versions.mix(BCFTOOLS_CALL_SNPS.out.versions) 

        } else{

            BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.bcf, BCFTOOLS_MPILEUP.out.csi, params.regions)
            // Gather used softwares versions and reports
            versions = versions.mix(BCFTOOLS_CALL.out.versions) 

        }

    }

    // ==== R E P O R T I N G ====

    version_yaml = Channel.empty()
    CUSTOM_DUMPSOFTWAREVERSIONS(versions.unique().collectFile(name: 'collated_versions.yml'))
    version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    multiqc_files = Channel.empty()
    multiqc_files = multiqc_files.mix(version_yaml)
    multiqc_files = multiqc_files.mix(reports.collect().ifEmpty([]))

    MULTIQC(multiqc_files.collect())

    multiqc_report = MULTIQC.out.report.toList()
    versions = versions.mix(MULTIQC.out.versions)
    

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {

    Channel.of(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by name and id and counting how many entries there are for this combination
        .map{ row ->
            //sample_count_all++
            if (!(row.name && row.seqid)) {
                error("Missing field in csv file header. The csv file must have fields named 'name' and 'seqid'.")
            }
            else if (row.name.contains(" ") || row.seqid.contains(" ")) {
                error("Invalid value in csv file. Values for 'name' and 'seqid' can not contain space.")
            }
            [ [ row.name.toString(), row.seqId.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map{ row, num_lanes -> // from here do the usual thing for csv parsing

        def meta = [:]

        if (row.name) meta.name = row.name.toString()
        if (row.seqId)  meta.seqid  = row.seqid.toString()

        if (row.seq_type) meta.seq_type = row.seq_type.toString()
        else meta.seq_type = 'paired'
       

        meta.id         = "${row.name}_X_${row.seqid}".toString()
        def fastq_1     = file(row.fastq_1, checkIfExists: true)
        if(row.fastq_2) fastq_2     = file(row.fastq_2, checkIfExists: true)

        def flowcell    = flowcellLaneFromFastq(fastq_1)
        // Don't use a random element for ID, it breaks resuming
        // the below read_group needs to be revised once all meta params are accessible
        def read_group  = "\"@RG\\tID:${flowcell}.${meta.seqid}\\tSM:${meta.name}\\tDS:${params.reference}\""

        //meta.num_lanes  = num_lanes.toInteger()
        meta.read_group = read_group.toString()
        meta.data_type  = 'fastq'

        if (row.fastq_2) return [ meta, [fastq_1, fastq_2] ]
        else {
            return [ meta, [fastq_1] ]
        }
    }
}


// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
