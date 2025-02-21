params.help = false // default


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


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
params.mapping_only = false // T/F whether to only perform read mapping and no variant calling
params.calling_only = false // T/F whether to only perform read mapping and no variant calling
params.regions = 'not_provided'

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
    perform only mapping: ${params.mapping_only}
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

include { MAPPING } from './subworkflows/mapping.nf'
include { VARIANT_CALLING } from './subworkflows/variant_calling.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/local/custom/custom_dumpsoftwareversions'
include { MULTIQC } from './modules/local/multiqc'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {
    // To gather all QC reports for MultiQC
    versions = Channel.empty()
    // To gather used softwares versions for MultiQC
    reports = Channel.empty()

	if (params.mapping_only) {
		input_sample = extract_csv(file(csv_file))
		MAPPING(input_sample, params.reference, params.skip_fastp, params.split_fastq, params.skip_markdup, params.markdup_tool)
		bam_ch = MAPPING.out.bam
		//versions = versions.mix(MAPPING.out.versions)
		reports = reports.mix(MAPPING.out.reports)
	}

	if (params.calling_only) {

		// Create the channel from the CSV file
		Channel
			.fromPath(csv_file)
			.splitCsv(header: true)
			.map { row -> 
				def meta = [name: row.name]
				tuple(meta, file(row.bam), file(row.bai))
			}
			.set { input_bams }

		VARIANT_CALLING(input_bams, params.reference, params.regions, params.grouped_call, params.target_index)
		//versions = versions.mix(VARIANT_CALLING.out.versions)
		reports = reports.mix(VARIANT_CALLING.out.reports)

	}
    
	if (!params.mapping_only && !params.calling_only) {
		// run both mapping and variant calling

	}
	   

    // ==== R E P O R T I N G ====

	MULTIQC(reports.collect())
    

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
