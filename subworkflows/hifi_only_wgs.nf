#!/usr/bin/env nextflow

include { RUN_HIFIASM_HIFI_ONLY } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HIFI_ONLY_WGS {
    
    take:
        ch_pb_reads

    main:
        RUN_HIFIASM_HIFI_ONLY (
            ch_pb_reads
                .map { 
                    pb_fastq, basename, platform -> 
                        tuple( file(pb_fastq), basename, "whole_genome" )
                }
                .filter {
                    pb_fastq, basename, region ->
                        file(pb_fastq).countFastq() > params.min_reads
                }
        )

        CONVERT_CONTIGS_TO_FASTA (
            RUN_HIFIASM_HIFI_ONLY.out
        )

}