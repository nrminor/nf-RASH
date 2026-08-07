#!/usr/bin/env nextflow

include { QUICK_SPLIT_FASTQ } from '../modules/quick_split_fastqs'
include { MAP_TO_REF } from '../modules/map_to_ref'
include { EXTRACT_REGIONS } from '../modules/extract_regions'
include { MERGE_PACBIO_FASTQS } from '../modules/merge_fastqs'
include { MERGE_ONT_FASTQS } from '../modules/merge_fastqs'
include { RUN_HIFIASM } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HYBRID_REGIONAL {

    take:
        ch_hybrid_samples
        ch_ref
        ch_desired_regions

    main:
        QUICK_SPLIT_FASTQ (
            ch_hybrid_samples
                .flatMap { sample_id, pb_fastq, ont_fastq ->
                    [
                        tuple( file(pb_fastq), sample_id, "pacbio" ),
                        tuple( file(ont_fastq), sample_id, "ont" )
                    ]
                }
        )

        MAP_TO_REF (
            QUICK_SPLIT_FASTQ.out
                .flatMap { fastqs, sample_id, platform ->
                    fastqs.collect { fastq -> tuple( file(fastq), sample_id, platform ) }
                },
                ch_ref
        )

        EXTRACT_REGIONS (
            MAP_TO_REF.out
                .combine ( ch_desired_regions )
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
                    sample_id, region, pb_fastq, _pacbio, ont_fastq, _ont -> 
                        tuple( file(pb_fastq), file(ont_fastq), sample_id, region )
                }
                .filter {
                    pb_fastq, ont_fastq, _sample_id, _region ->
                        file(pb_fastq).countFastq() > params.min_reads &&
                        file(ont_fastq).countFastq() > params.min_reads
                }
        )

        CONVERT_CONTIGS_TO_FASTA (
            RUN_HIFIASM.out
        )

}
