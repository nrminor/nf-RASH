#!/usr/bin/env nextflow

include { HYBRID_WGS } from "../subworkflows/hybrid_wgs"
include { HIFI_ONLY_WGS } from "../subworkflows/hifi_only_wgs"

workflow WHOLE_GENOME {
    
    take:
        ch_pb_reads

    main:
        // if an ont fastq is provided, run hybrid assembly
        if ( params.ont_fastq ) {

            // raise an error if the provided ont FASTQ path doesn't exist
            assert file(params.ont_fastq).exists() : "Provided path to Nanopore FASTQ does not exist."

            // create the ont channel tuple
            ch_ont_reads = Channel
                .fromPath ( params.ont_fastq )
                .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "ont" )}

            // run hybrid assembly
            HYBRID_WGS (
                ch_pb_reads,
                ch_ont_reads
            )
        
        // otherwise, just use the provided PacBio HiFi reads
        } else {

            // run hifi-only assembly
            HIFI_ONLY_WGS (
                ch_pb_reads
            )

        }

}
