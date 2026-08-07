#!/usr/bin/env nextflow

include { RUN_HIFIASM } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HYBRID_WGS {

    take:
        ch_hybrid_samples

    main:
        RUN_HIFIASM (
            ch_hybrid_samples
                .map { 
                    sample_id, pb_fastq, ont_fastq -> 
                        tuple( file(pb_fastq), file(ont_fastq), sample_id, "whole_genome" )
                }
        )

        CONVERT_CONTIGS_TO_FASTA (
            RUN_HIFIASM.out
        )

}
