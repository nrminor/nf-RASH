#!/usr/bin/env nextflow

include { RUN_HIFIASM_HIFI_ONLY } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HIFI_ONLY_WGS {
    
    take:
        ch_hifi_only_samples

    main:
        RUN_HIFIASM_HIFI_ONLY (
            ch_hifi_only_samples
                .map { 
                    sample_id, pb_fastq, _ont_fastq -> 
                        tuple( file(pb_fastq), sample_id, "whole_genome" )
                }
                .filter {
                    pb_fastq, _sample_id, _region ->
                        file(pb_fastq).countFastq() > params.min_reads
                }
        )

        CONVERT_CONTIGS_TO_FASTA (
            RUN_HIFIASM_HIFI_ONLY.out
        )

}
