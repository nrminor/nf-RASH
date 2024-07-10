#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*

Welcome to RASH!

Like many Nextflow pipelines, nf-RASH is a metaphorical Russian nesting doll of
workflows that functions like a binary tree of switchboards. The core switchboard that
is handled by main.nf is:
    > Did the user request regional assembly by providing a TSV file of genome
    > coordinates in a reference sequence?

If the user does (e.g., the TSV's provided in the resources/ directory), the pipeline
switches into its regional assembly mode. If the user doesn't provide the region TSV,
the pipeline switches into whole-genome assembly mode.

Each of those two modes themselves have a second switchboard:
    > Did the user provide a path to Oxford Nanopore ultra-long reads?

If they did, regional assembly mode and whole genome mode will switch into hybrid
assembly modes. If not, the two workflow modes will work only with the provided Hifi
reads.

In sum, nf-RASH allows users to do regional assembly or whole genome assembly and do so
with Hifi reads or Hifi and Nanopore reads. Perhaps a diagram will help:

                                PacBio Hifi reads
                                      |
                                      |
                                      | Desired regions?
                                      |
                         Yes:         |          No:
                       Regional       |      Whole-genome
                       Assembly       |       Assembly
                        ————————————————————————————
                        |                          |
             ONT reads? |                          | ONT reads?
                        |                          |
               No:      |     Yes:        No:      |     Yes:
                ————————————————           ————————————————
                |              |           |              |
                |              |           |              |
               Hifi          HiFi        HiFi            HiFi
               only         and ONT      only           and ONT
                           (hybrid)                    (hybrid)


*/

include { WHOLE_GENOME } from './workflows/whole_genome'
include { REGIONAL_BAITING } from './workflows/regional_baiting'

// log out some of the information provided by the user
log.info    """
            RASH: Regional ASsembly Helper
            ===================================
            RASH is a containerized Nextflow pipeline for extracting genome regions
            of interest from PacBio HiFi and Oxford Nanopore sequencing reads and 
            running them through high-accuracy hybrid assembly using Hifiasm.
            RASH also supports HiFi-only assembly through the hifi_only workflow,
            which will be invoked when a Nanopore FASTQ isn't provided by the user.
            (version 0.2.3)
            ===================================

            Inputs and Outputs:
            ----------------------------------
            PacBio FASTQ           : ${params.pb_fastq}
            ONT FASTQ (optional)   : ${params.ont_fastq ?: ""}
            Reference FASTA        : ${params.ref_fasta ?: ""}
            Regions TSV            : ${params.desired_regions ?: ""}
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

    // make sure the minimum user provided input, PacBio HiFi reads, exists
    assert params.pb_fastq : "Please provide a PacBio HiFi CCS FASTQ.gz file with the --pb_fastq argument."
    assert file(params.pb_fastq).exists() : "Provided path to PacBio FASTQ does not exist."

	// input channels shared by both workflows
    ch_pb_reads = Channel
        .fromPath ( params.pb_fastq )
        .map { fastq -> tuple( file(fastq), file(fastq).getSimpleName(), "pacbio" )}

    // if desired regions are provided, assemble just those regions
    if ( params.desired_regions ) {

        // raise an error if the provided TSV path doesn't exist or if the reference doesn't exist
        assert file(params.desired_regions).exists() : "Provided path to desired region TSV does not exist."
        assert params.ref_fasta : "Please provide a reference FASTA with the --ref_fasta argument."
        assert file(params.ref_fasta).exists() : "Provided path to reference FASTA does not exist."
    
        // launch input channels for the reference FASTA and the desired region coordinates
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

        // run regional baiting
        REGIONAL_BAITING (
            ch_pb_reads,
            ch_ref,
            ch_desired_regions
        )
    
    // otherwise, run whole-genome assembly
    } else {

        WHOLE_GENOME (
            ch_pb_reads
        )

    }

}

