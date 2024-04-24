#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// -------------------------------------------------------------------------- //
workflow {


    // make sure the user provided inputs that exist
    assert params.pb_fastq : "Please provide a PacBio HiFi CCS FASTQ.gz file with the --pb_fastq argument."
    assert file(params.pb_fastq).exists() : "Provided path to PacBio FASTQ does not exist."
	assert params.ont_fastq : "Please provide a Oxford Nanopore FASTQ.gz file file with the --ont_fastq argument."
    assert file(params.ont_fastq).exists() : "Provided path to Nanopore FASTQ does not exist."
	assert params.ref_fasta : "Please provide a reference FASTA with the --ref_fasta argument."
    assert file(params.ref_fasta).exists() : "Provided path to reference FASTA does not exist."


	// input channels
    ch_pb_reads = Channel
        .fromPath ( params.pb_fastq )
        .map { fastq -> tuple( file(fastq), "pacbio" )}

    ch_ont_reads = Channel
        .fromPath ( params.ont_fastq )
        .map { fastq -> tuple( file(fastq), "ont" )}
    
    ch_ref = Channel
        .fromPath ( params.ref_fasta )
    
    ch_desired_regions = Channel
        .fromPath ( params.desired_regions )
        .splitCsv( header: true, sep: "\t", strip: true )
        .map { row -> tuple( row.samtools_expression, row.file_label, row.description ) }


	// Workflow steps
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

    EXTRACT_DESIRED_REGIONS (
        MAP_TO_REF.out,
        ch_desired_regions
    )

    MERGE_PACBIO_FASTQS (
        EXTRACT_DESIRED_REGIONS.out
            .filter { x -> x[2] == "pacbio" }
            .groupTuple ( by: [ 1, 2, 3 ] )
    )

    MERGE_ONT_FASTQS (
        EXTRACT_DESIRED_REGIONS.out
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
    )

    CONVERT_CONTIGS_TO_FASTA (
        RUN_HIFIASM.out
    )

}
// -------------------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// -------------------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// where to place the PacBio and ONT FASTQ for each extracted region
params.extracted = params.results + "/01_extracted_regions"

// where to place hifiasm results, including a FASTA of contigs
params.assembly = params.results + "/02_hifiasm_assembly"

// -------------------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// -------------------------------------------------------------------------- //

process QUICK_SPLIT_FASTQ {

	/*
    Rather than using Nextflow's fastq-splitting API (which is very cool but 
    also slow), we use Seqkit to split the FASTQs into smaller pieces based on
    the `split_max` parameter in `nextflow.config`. By default, this parameter
    tells Seqkit to split the large input FASTQ's into as many FASTQs as it
    takes for each to have no more than 500,000 reads. In our experience, this
    results in FASTQ's that are 3 to 5 GB in size each.

    Each split FASTQ is mapped to a reference downstream, and then in parallel,
    regions based on that reference are extracted from each split FASTQ. Because
    the input FASTQ's are often very large, splitting makes mapping and region
    extraction considerable faster and less computationally intensive.
    */

    tag "${platform}, ${params.split_max} reads per file"
    label "seqkit"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

    input:
    tuple path(big_ol_fastq), val(platform)

    output:
    tuple path("split/*.fastq.gz"), val(platform)

    script:
    """
    seqkit split2 \
    --by-size ${params.split_max} \
    --extension ".gz" \
    --out-dir split/ \
    --threads ${task.cpus} \
    ${big_ol_fastq}
    """

}

process MAP_TO_REF {

	/*
    Users provide a reference with the parameter `ref_fasta` in
    `nextflow.config`, which is used to bait out reads that are likely to 
    contain regions of interest.
    */

	tag "${basename}, ${platform}"
    label "map_and_extract"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path(fastq), val(basename), val(platform)
    each path(ref_fasta)

	output:
    path "*.bam"

	script:
    minimap2_preset = platform == "pacbio" ? "map-hifi" : "map-ont"
	"""
    minimap2 -t ${task.cpus} -L --eqx -ax ${minimap2_preset} \
    ${ref_fasta} \
    ${fastq} \
    | samtools view -Sbt ${ref_fasta} \
    | samtools sort - -o ${basename}_${platform}.bam
	"""

}

process EXTRACT_DESIRED_REGIONS {

	/* 
    Users provide a TSV of regions of interest with the parameter 
    `desired_regions` in `nextflow.config`. Each row of that TSV is split into
    its own task in queue. Every input FASTQ to this process is thus
    multiplexed across regions called for in the TSV; if there are N rows, each
    FASTQ will be run through this process N times to extract reads mapping to
    those N regions.
    */

	tag "${basename}, ${platform}, ${file_label}"
    label "map_and_extract"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    each path(bam)
    tuple val(expression), val(file_label), val(description)

	output:
    tuple path("${basename}_${platform}_${file_label}.fastq.gz"), val(basename), val(platform), val(file_label)

	script:
    bam_components = bam.toString().replace(".bam", "").split("_")
    assert bam_components.size() == 2 : "Necessary information could not be parsed from $bam.toString()."
    basename = bam_components[0]
    platform = bam_components[1]
	"""
    samtools index ${bam}
    samtools view -b ${bam} ${expression} \
    | samtools fastq - \
    | reformat.sh qin=33 int=f in=stdin.fq \
    out=${basename}_${platform}_${file_label}.fastq.gz
	"""

}

process MERGE_PACBIO_FASTQS {

	/*
    Now that we've split, mapped, and extracted regions from our input FASTQs,
    we need to merge them back together such that we have one PacBio FASTQ and
    one Oxford Nanopore FASTQ for each region. This process uses Seqkit to
    run these merges performantly on PacBio reads.
    */

    tag "${basename}, ${platform}, ${file_label}"
    label "seqkit"
	publishDir params.extracted, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path("to_merge/???.fastq.gz"), val(basename), val(platform), val(file_label)

	output:
    tuple path("${basename}_${platform}_${file_label}.fastq.gz"), val(basename), val(platform), val(file_label)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq ./to_merge/ \
    | gzip -c > ${basename}_${platform}_${file_label}.fastq.gz
	"""

}

process MERGE_ONT_FASTQS {

	/*
    Now that we've split, mapped, and extracted regions from our input FASTQs,
    we need to merge them back together such that we have one PacBio FASTQ and
    one Oxford Nanopore FASTQ for each region. This process uses Seqkit to
    run these merges performantly on Oxford Nanopore reads.
    */

    tag "${basename}, ${platform}, ${file_label}"
    label "seqkit"
	publishDir params.extracted, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path("to_merge/???.fastq.gz"), val(basename), val(platform), val(file_label)

	output:
    tuple path("${basename}_${platform}_${file_label}.fastq.gz"), val(basename), val(platform), val(file_label)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq ./to_merge/ \
    | gzip -c > ${basename}_${platform}_${file_label}.fastq.gz
	"""

}

process RUN_HIFIASM {

	/*
    Now that we have PacBio and Oxford Nanopore FASTQs prepared for all regions
    requested in the `desired_regions` parameter TSV (see above), we're ready
    to run hybrid assembly on them with Hifiasm. This process handles assembly,
    making sure that only FASTQs with matching regions are assembled together.
    */

	tag "${basename}, ${region}"
	publishDir "${params.assembly}/${basename}_${region}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path(ont_fastq), path(pacbio_fastq), val(basename), val(region)

	output:
    tuple path("*"), val(basename), val(region)

	script:
    assert ont_fastq.toString().toLowerCase().contains("ont")
    assert pb_fastq.toString().toLowerCase().contains("pacbio")
	"""
    hifiasm -o ${basename}_${region} -t ${task.cpus} --ul ${ont_fastq} ${pb_fastq}
	"""

}

process CONVERT_CONTIGS_TO_FASTA {

	/*
    A cute trick from https://hifiasm.readthedocs.io/en/latest/faq.html#id1
    allows us to convert our contigs into a FASTA file that we can then
    inspect in Geneious (or wherever).
    */

	tag "${basename}, ${platform}, ${file_label}"
    label "map_and_extract"
	publishDir "${params.assembly}/${basename}_${region}", mode: 'copy', overwrite: true

    cpus 3

	input:
    tuple path("hifiasm_files/*"), val(basename), val(region)

	output:
    path "${basename}_${region}.p_contigs.fasta"

	shell:
	'''
    awk '/^S/{print ">"$2"n"$3}' hifiasm_files/!{basename}_!{region}.bp.p_ctg.gfa \
    | fold > !{basename}_!{region}.p_contigs.fasta
	'''

}

// -------------------------------------------------------------------------- //
