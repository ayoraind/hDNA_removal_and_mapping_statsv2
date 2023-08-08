#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()
// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)
final_params = check_params(merged_params)
// starting pipeline
pipeline_start_message(version, final_params)

// include processes
include { MINIMAP2_INDEX; MINIMAP2_SAM; SAM_SORT_AND_INDEX; EXTRACT_MICROBIAL_READS; MINIMAP2_INDEX_SECOND_ITERATION; MINIMAP2_SAM_SECOND_ITERATION; SAM_SORT_AND_INDEX_SECOND_ITERATION; EXTRACT_MICROBIAL_READS_SECOND_ITERATION; BAM_STATISTICS;  COMBINE_BAM_STATISTICS; BAM_STATISTICS_SECOND_ITERATION;  COMBINE_BAM_STATISTICS_SECOND_ITERATION } from './modules/processes.nf' addParams(final_params)

workflow  {
         fasta_ch = channel
                          .fromPath( final_params.reference_fasta, checkIfExists: true )
                          .map { file -> tuple(file.baseName, file) }
         reads_ch = channel
                          .fromPath( final_params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }


         MINIMAP2_INDEX(fasta_ch)
	 
	 ch_index  = MINIMAP2_INDEX.out.index_ch

         // combine method used instead of join method as both channels (reads_ch and ch_index) share no key
         combine_ch = reads_ch.combine(ch_index)

         MINIMAP2_SAM ( combine_ch )

         SAM_SORT_AND_INDEX(MINIMAP2_SAM.out.sam_ch)

         EXTRACT_MICROBIAL_READS(SAM_SORT_AND_INDEX.out.bam_ch)
	 
	 BAM_STATISTICS(SAM_SORT_AND_INDEX.out.bam_ch)


         collected_bam_statistics_ch = BAM_STATISTICS.out.bam_stats_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_BAM_STATISTICS(collected_bam_statistics_ch)
	 
	 
	 if (final_params.second_reference_fasta){
	 
	 second_fasta_ch = channel
                          .fromPath( final_params.second_reference_fasta, checkIfExists: true )
                          .map { file -> tuple(file.baseName, file) }
			  
	 MINIMAP2_INDEX_SECOND_ITERATION(second_fasta_ch)
	 
	 ch_index  = MINIMAP2_INDEX_SECOND_ITERATION.out.index_ch

         // combine method used instead of join method as both channels (reads_ch and ch_index) share no key
         combine_ch = EXTRACT_MICROBIAL_READS.out.microbial_reads_ch.combine(ch_index)

         MINIMAP2_SAM_SECOND_ITERATION ( combine_ch )

         SAM_SORT_AND_INDEX_SECOND_ITERATION(MINIMAP2_SAM_SECOND_ITERATION.out.sam_ch)

         EXTRACT_MICROBIAL_READS_SECOND_ITERATION(SAM_SORT_AND_INDEX_SECOND_ITERATION.out.bam_ch)
	 
	 BAM_STATISTICS_SECOND_ITERATION(SAM_SORT_AND_INDEX_SECOND_ITERATION.out.bam_ch)


         collected_bam_statistics_ch = BAM_STATISTICS_SECOND_ITERATION.out.bam_stats_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_BAM_STATISTICS_SECOND_ITERATION(collected_bam_statistics_ch)
	 
	 }

         
}

workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
