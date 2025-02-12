/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT CALLING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX } from '../modules/local/samtools_faidx'
include { SAMTOOLS_INDEX } from '../modules/local/samtools_index'
include { BCFTOOLS_MPILEUP } from '../modules/local/bcftools_mpileup'
include { BCFTOOLS_CALL } from '../modules/local/bcftools_call'
include { BCFTOOLS_CALL_SNPS } from '../modules/local/bcftools_call_snps'
include { BCFTOOLS_GRP_CALL } from '../modules/local/bcftools_grp_call'
include { BCFTOOLS_CONCAT } from '../modules/local/bcftools_concat'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom/custom_dumpsoftwareversions'
include { MULTIQC } from '../modules/local/multiqc'
include { AWK_SPLITBED } from '../modules/local/awk_splitbed'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow VARIANT_CALLING {

	take:
		indexed_bams
		reference
		regions
		grouped_call
		targets_index

	main:

		ch_versions = Channel.empty()
		ch_reports = Channel.empty()

		// Prepare reference fasta index
    	SAMTOOLS_FAIDX(reference)
    
    	// ==== V A R I A N T   C A L L I N G ====

		if (grouped_call) { // perform variant calling on all samples together (e.g. variant discovery, mapping families)

            split_beds = AWK_SPLITBED(regions)

			// Create a channel with region name and bed file
    		region_beds = split_beds.flatten()
        		.map { bed -> 
            	def region = bed.baseName
            	return tuple(region, bed)
        		}

            BCFTOOLS_GRP_CALL(indexed_bams.map { tuple -> tuple[1] }.collect(), indexed_bams.map { tuple -> tuple[2] }.collect(), reference, SAMTOOLS_FAIDX.out.fai, region_beds)

			//BCFTOOLS_GRP_CALL(indexed_bams.map { tuple -> tuple[1] }.collect(), indexed_bams.map { tuple -> tuple[2] }.collect(), reference, SAMTOOLS_FAIDX.out.fai, AWK_SPLITBED.out.subregion.flatten())
			// Gather used softwares versions and reports
			ch_versions = ch_versions.mix(BCFTOOLS_GRP_CALL.out.versions) 

            BCFTOOLS_CONCAT(BCFTOOLS_GRP_CALL.out.vcf.collect())
            // Gather used softwares versions and reports
			ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions) 


		} else { // perform variant calling seperately for each sample (e.g. genotyping pipeline)

			// Groups multi bam from the same sample 'name' to run mpileup on
			BCFTOOLS_MPILEUP(indexed_bams.map{ meta, bam, bai -> [ meta.subMap('name'), bam, bai ] }.groupTuple(), reference, SAMTOOLS_FAIDX.out.fai, regions)
			// Gather used softwares versions and reports
			ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)  

			if(params.constrain_alleles) { // constrain genotype calls to those alleles provided an indexed als.tsv type file

				BCFTOOLS_CALL_SNPS(BCFTOOLS_MPILEUP.out.bcf, BCFTOOLS_MPILEUP.out.csi, regions, targets_index)
				// Gather used softwares versions and reports
				ch_versions = ch_versions.mix(BCFTOOLS_CALL_SNPS.out.versions) 

			} else{

				BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.bcf, BCFTOOLS_MPILEUP.out.csi, regions)
				// Gather used softwares versions and reports
				ch_versions = ch_versions.mix(BCFTOOLS_CALL.out.versions) 

			}

		}

	emit:
		versions = ch_versions
		reports = ch_reports

}
