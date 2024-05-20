#!/usr/bin/env nextflow

include { QUICK_SPLIT_FASTQ } from '../modules/quick_split_fastqs'
include { MAP_TO_REF } from '../modules/map_to_ref'
include { EXTRACT_REGIONS } from '../modules/extract_regions'
include { MERGE_PACBIO_FASTQS } from '../modules/merge_fastqs'
include { RUN_HIFIASM_HIFI_ONLY } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HIFI_ONLY {
    
    take:
        ch_pb_reads
        ch_ref
        ch_desired_regions

    main:
        QUICK_SPLIT_FASTQ (
            ch_pb_reads
        )

        MAP_TO_REF (
            QUICK_SPLIT_FASTQ.out
                .map { fastqs, platform -> fastqs }
                .flatten ( )
                .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "pacbio" ) },
                ch_ref
        )

        EXTRACT_REGIONS (
            MAP_TO_REF.out,
            ch_desired_regions
        )

        MERGE_PACBIO_FASTQS (
            EXTRACT_REGIONS.out
                .groupTuple ( by: [ 1, 2, 3 ] )
        )

        RUN_HIFIASM_HIFI_ONLY (
            MERGE_PACBIO_FASTQS.out
                .map { 
                    pb_fastq, basename, platform, region -> 
                        tuple( file(pb_fastq), basename, region )
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