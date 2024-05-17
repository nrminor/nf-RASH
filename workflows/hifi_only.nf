#!/usr/bin/env nextflow

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

        EXTRACT_DESIRED_REGIONS (
            MAP_TO_REF.out,
            ch_desired_regions
        )

        MERGE_PACBIO_FASTQS (
            EXTRACT_DESIRED_REGIONS.out
                .groupTuple ( by: [ 1, 2, 3 ] )
        )

        RUN_HIFIASM_HIFI_ONLY (
            MERGE_PACBIO_FASTQS.out
                .map { 
                    basename, region, pb_fastq, pacbio, ont_fastq, ont -> 
                        tuple( file(pb_fastq), file(ont_fastq), basename, region )
                }
                .filter {
                    pb_fastq, ont_fastq, basename, region ->
                        file(pb_fastq).countFastq() > params.min_reads
                }
        )

        CONVERT_CONTIGS_TO_FASTA (
            RUN_HIFIASM.out
        )

}
// -------------------------------------------------------------------------- //