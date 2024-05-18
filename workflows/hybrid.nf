#!/usr/bin/env nextflow

include { QUICK_SPLIT_FASTQ } from '../modules/quick_split_fastqs'
include { MAP_TO_REF } from '../modules/map_to_ref'
include { EXTRACT_REGIONS } from '../modules/extract_regions'
include { MERGE_PACBIO_FASTQS } from '../modules/merge_fastqs'
include { MERGE_ONT_FASTQS } from '../modules/merge_fastqs'
include { RUN_HIFIASM } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HYBRID {

    take:
        ch_pb_reads
        ch_ont_reads
        ch_ref
        ch_desired_regions

    main:
        QUICK_SPLIT_FASTQ (
            ch_pb_reads
                .mix ( ch_ont_reads )
        )

        MAP_TO_REF (
            QUICK_SPLIT_FASTQ.out
                .filter { fastq, platform -> platform == "pacbio" }
                .map { fastqs, platform -> fastqs }
                .flatten ( )
                .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "pacbio" ) }
                .mix (

                    QUICK_SPLIT_FASTQ.out
                        .filter { fastq, platform -> platform == "ont" }
                        .map { fastqs, platform -> fastqs }
                        .flatten ( )
                        .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "ont" ) }

                ),
                ch_ref
        )

        EXTRACT_REGIONS (
            MAP_TO_REF.out,
            ch_desired_regions
        )

        MERGE_PACBIO_FASTQS (
            EXTRACT_REGIONS.out
                .filter { x -> x[2] == "pacbio" }
                .groupTuple ( by: [ 1, 2, 3 ] )
        )

        MERGE_ONT_FASTQS (
            EXTRACT_REGIONS.out
                .filter { x -> x[2] == "ont" }
                .groupTuple ( by: [ 1, 2, 3 ] )
        )

        RUN_HIFIASM (
            MERGE_PACBIO_FASTQS.out
                .join ( MERGE_ONT_FASTQS.out, by: [ 1, 3 ] )
                .map { 
                    basename, region, pb_fastq, pacbio, ont_fastq, ont -> 
                        tuple( file(pb_fastq), file(ont_fastq), basename, region )
                }
                .filter {
                    pb_fastq, ont_fastq, basename, region ->
                        file(pb_fastq).countFastq() > params.min_reads &&
                        file(ont_fastq).countFastq() > params.min_reads
                }
        )

        CONVERT_CONTIGS_TO_FASTA (
            RUN_HIFIASM.out
        )

}
