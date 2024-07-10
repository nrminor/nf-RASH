#!/usr/bin/env nextflow

include { HIFI_ONLY_REGIONAL } from "../subworkflows/hifi_only_regional"
include { HYBRID_REGIONAL } from "../subworkflows/hybrid_regional"

workflow REGIONAL_BAITING {

    take:
        ch_pb_reads,
        ch_ref,
        ch_desired_regions

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
            HYBRID_REGIONAL (
                ch_pb_reads,
                ch_ont_reads,
                ch_ref,
                ch_desired_regions
            )
        
        // otherwise, just use the provided PacBio HiFi reads
        } else {

            // run hifi-only assembly
            HIFI_ONLY_REGIONAL (
                ch_pb_reads,
                ch_ref,
                ch_desired_regions
            )

        }

}
