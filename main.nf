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

// define the main workflow
workflow {

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
                Samplesheet            : ${params.samplesheet}
                Reference FASTA        : ${params.ref_fasta ?: ""}
                Regions TSV            : ${params.desired_regions ?: ""}
                results_dir            : ${params.results}

                Run settings:
                -----------------------------------
                Reads per split FASTQ : ${params.split_max}
                Min reads per region  : ${params.min_reads}
                cleanup               : ${params.cleanup}
                """
                .stripIndent()

    // parse one row per sample and route it according to whether ONT reads were provided
    assert params.samplesheet : "Please provide a CSV samplesheet with the --samplesheet argument."
    assert file(params.samplesheet).exists() : "Provided path to samplesheet does not exist."

    ch_samples_by_type = channel
        .fromPath ( params.samplesheet )
        .splitCsv ( header: true, sep: ",", strip: true )
        .map { row ->
            assert row.containsKey("sample_id") : "Samplesheet is missing the sample_id column."
            assert row.containsKey("pb_fastq") : "Samplesheet is missing the pb_fastq column."
            assert row.containsKey("ont_fastq") : "Samplesheet is missing the ont_fastq column."

            def sample_id = row.sample_id?.trim()
            assert sample_id : "Every samplesheet row must include a sample_id."
            assert sample_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/ : "Invalid sample_id '${sample_id}'. Use only letters, numbers, periods, underscores, and hyphens."

            assert row.pb_fastq : "Sample '${sample_id}' is missing a PacBio FASTQ path."
            def pb_fastq = file(row.pb_fastq)
            assert pb_fastq.exists() : "PacBio FASTQ for sample '${sample_id}' does not exist: ${row.pb_fastq}"

            def ont_fastq = row.ont_fastq ? file(row.ont_fastq) : null
            assert !ont_fastq || ont_fastq.exists() : "ONT FASTQ for sample '${sample_id}' does not exist: ${row.ont_fastq}"

            tuple( sample_id, pb_fastq, ont_fastq )
        }
        .ifEmpty { error "Samplesheet must contain at least one sample." }
        .collect ( flat: false )
        .flatMap { samples ->
            def sample_ids = samples.collect { sample -> sample[0] }
            assert sample_ids.size() == sample_ids.toSet().size() : "Samplesheet sample_id values must be unique."
            samples
        }
        .branch { sample ->
            hybrid: sample[2] != null
            hifi_only: true
        }

    // if desired regions are provided, assemble just those regions
    if ( params.desired_regions ) {

        // raise an error if the provided TSV path doesn't exist or if the reference doesn't exist
        assert file(params.desired_regions).exists() : "Provided path to desired region TSV does not exist."
        assert params.ref_fasta : "Please provide a reference FASTA with the --ref_fasta argument."
        assert file(params.ref_fasta).exists() : "Provided path to reference FASTA does not exist."
    
        // launch input channels for the reference FASTA and the desired region coordinates
        ch_ref = channel
            .fromPath ( params.ref_fasta )

        ch_desired_regions = channel
            .fromPath ( params.desired_regions )
            .splitCsv ( header: true, sep: "\t", strip: true )
            .map { 
                row -> tuple( 
                    "${row.chromosome}:${row.start}-${row.stop}", row.region, row.merge_key
                ) 
            }

        // run regional baiting
        REGIONAL_BAITING (
            ch_samples_by_type.hybrid,
            ch_samples_by_type.hifi_only,
            ch_ref,
            ch_desired_regions
        )
    
    // otherwise, run whole-genome assembly
    } else {

        WHOLE_GENOME (
            ch_samples_by_type.hybrid,
            ch_samples_by_type.hifi_only
        )

    }

}
