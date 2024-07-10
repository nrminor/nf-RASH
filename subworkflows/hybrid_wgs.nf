#!/usr/bin/env nextflow

include { QUICK_SPLIT_FASTQ } from '../modules/quick_split_fastqs'
include { MAP_TO_REF } from '../modules/map_to_ref'
include { EXTRACT_REGIONS } from '../modules/extract_regions'
include { MERGE_PACBIO_FASTQS } from '../modules/merge_fastqs'
include { MERGE_ONT_FASTQS } from '../modules/merge_fastqs'
include { RUN_HIFIASM } from '../modules/hifiasm'
include { CONVERT_CONTIGS_TO_FASTA } from '../modules/convert_to_fasta'

workflow HYBRID_WGS {

    take:
        ch_pb_reads
        ch_ont_reads

    main:
        RUN_HIFIASM (
            ch_pb_reads
                .join ( ch_ont_reads, by: 2 )
                .map { 
                    basename, pb_fastq, pacbio, ont_fastq, ont -> 
                        tuple( file(pb_fastq), file(ont_fastq), basename, "whole_genome" )
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
