#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HYBRID } from './workflows/hybrid'
include { HIFI_ONLY } from './workflows/hifi_only'

// log out some of the information provided by the user
log.info    """
            RASH: Regional ASsembly Helper
            ===================================
            RASH is a containerized Nextflow pipeline for extracting genome regions
            of interest from PacBio HiFi and Oxford Nanopore sequencing reads and 
            running them through high-accuracy hybrid assembly using Hifiasm.
            RASH also supports HiFi-only assembly through the hifi_only workflow,
            which will be invoked when a Nanopore FASTQ isn't provided by the user.
            (version 0.1.2)
            ===================================

            Inputs and Outputs:
            ----------------------------------
            PacBio FASTQ           : ${params.pb_fastq}
            ONT FASTQ (optional)   : ${params.ont_fastq_log}
            Reference FASTA        : ${params.ref_fasta}
            Regions TSV            : ${params.desired_regions}
            results_dir            : ${params.results}

            Run settings:
            -----------------------------------
            Reads per split FASTQ : ${params.split_max}
            Min reads per region  : ${params.min_reads}
            cleanup               : ${params.cleanup}
            cpus per task         : ${params.cpus}
            """
            .stripIndent()

// define the main workflow
workflow {

    // make sure the user provided inputs exist
    assert params.pb_fastq : "Please provide a PacBio HiFi CCS FASTQ.gz file with the --pb_fastq argument."
    assert file(params.pb_fastq).exists() : "Provided path to PacBio FASTQ does not exist."
	assert params.ref_fasta : "Please provide a reference FASTA with the --ref_fasta argument."
    assert file(params.ref_fasta).exists() : "Provided path to reference FASTA does not exist."

	// input channels shared by both workflows
    ch_pb_reads = Channel
        .fromPath ( params.pb_fastq )
        .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "pacbio" )}
    
    ch_ref = Channel
        .fromPath ( params.ref_fasta )
    
    ch_desired_regions = Channel
        .fromPath ( params.desired_regions )
        .splitCsv ( header: true, sep: "\t", strip: true )
        .map { 
            row -> tuple( 
                "${row.chromosome}:${row.start}-${row.stop}", row.region, row.merge_key
            ) 
        }

    // if an ont fastq is provided, run hybrid assembly
    if ( params.ont_fastq ) {

        // raise an error if the provided ont FASTQ path doesn't exist
        assert file(params.ont_fastq).exists() : "Provided path to Nanopore FASTQ does not exist."

        // create the ont channel tuple
        ch_ont_reads = Channel
            .fromPath ( params.ont_fastq )
            .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "ont" )}

        // run hybrid assembly
        HYBRID (
            ch_pb_reads,
            ch_ont_reads,
            ch_ref,
            ch_desired_regions
        )
    
    // otherwise, just use the provided PacBio HiFi reads
    } else {

        // run hifi-only assembly
        HIFI_ONLY (
            ch_pb_reads,
            ch_ref,
            ch_desired_regions
        )

    }

}

