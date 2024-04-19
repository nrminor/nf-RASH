#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
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
        .splitFastq ( by: params.split_max, file: true, compress: true )
        .map { fastq -> tuple( file(fastq), file(fastq).getBaseName(), "pacbio" ) }

    ch_ont_reads = Channel
        .fromPath ( params.ont_fastq )
        .splitFastq ( by: params.split_max, file: true, compress: true )
        .map { fastq -> tuple( file(fastq), file(fastq).getBaseName(), "ont" ) }
    
    ch_ref = Channel
        .fromPath ( params.ref_fasta )
    
    ch_desired_regions = Channel
        .fromPath ( params.desired_regions )
        .splitCsv( header: true, sep: "\t", strip: true )
        .map { row -> tuple( row.samtools_expression, row.file_label, row.description ) }
	
	
	// Workflow steps
    MAP_TO_REF (
        ch_pb_reads
            .mix ( ch_ont_reads ),
        ch_ref
    )

    EXTRACT_DESIRED_REGIONS (
        MAP_TO_REF.out,
        ch_desired_regions
    )

    MERGE_PACBIO_FASTQS (
        EXTRACT_DESIRED_REGIONS.out
            .filter { x[2] == "pacbio" }
            .groupTuple( by: [ 1, 2, 3] )
    )

    MERGE_ONT_FASTQS (
        EXTRACT_DESIRED_REGIONS.out
            .filter { x[2] == "ont" }
            .groupTuple( by: [ 1, 2, 3] )
    )

    RUN_HIFIASM (
        MERGE_PACBIO_FASTQS.out,
        MERGE_ONT_FASTQS.out
    )

    CONVERT_CONTIGS_TO_FASTA (
        RUN_HIFIASM.out
    )

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

params.extracted = params.results + "/01_extracted_regions"
params.assembly = params.results + "/02_hifiasm_assembly"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process MAP_TO_REF {

	/* */

	tag "${basename}, ${platform}"
    label "map_and_extract"

    cpus 8

	input:
    tuple path(fastq), val(basename), val(platform)
    each path(ref_fasta)

	output:
    path "*.bam"

	script:
    minimap2_preset = platform == "pacbio" ? "map-hifi" : "map-ont"
	"""
    minimap2 -t 6 -L --eqx -ax ${minimap2_preset} \
    ${ref_fasta} \
    ${fastq} \
    | samtools view -Sbt ${ref_fasta} \
    | samtools sort - -o ${basename}_{platform}.bam
	"""

}

process EXTRACT_DESIRED_REGIONS {

	/* */

	tag "${basename}, ${platform}, ${file_label}"
    label "map_and_extract"

	input:
    each path(bam)
    tuple val(expression), val(file_label), val(description)

	output:
    tuple path("${basename}_${platform}_${file_label}.fastq.gz"), val(basename), val(platform), val(file_label)

	script:
    basename = file(bam).getSimpleName().split("_")[0]
    platform = file(bam).getSimpleName().split("_")[1]
	"""
    samtools index ${bam}
    samtools view -b ${bam} ${expression} \
    | samtools fastq - \
    | reformat.sh qin=33 int=f in=stdin.fq \
    out=${basename}_${platform}_${file_label}.fastq.gz
	"""

}

process MERGE_PACBIO_FASTQS {

	/* */

    tag "${basename}, ${platform}, ${file_label}"
	publishDir params.extracted, mode: 'copy', overwrite: true

    cpus 10

	input:
    tuple path("to_merge/*"), val(basename), val(platform), val(file_label)

	output:
    tuple path("${basename}_${platform}_${file_label}.fastq.gz"), val(basename), val(platform), val(file_label)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq
    to_merge/ | gzip -c > ${basename}_${platform}_${file_label}.fastq.gz
	"""

}

process MERGE_ONT_FASTQS {

	/* */

    tag "${basename}, ${platform}, ${file_label}"
	publishDir params.extracted, mode: 'copy', overwrite: true

    cpus 10

	input:
    tuple path("to_merge/*"), val(basename), val(platform), val(file_label)

	output:
    tuple path("${basename}_${platform}_${file_label}.fastq.gz"), val(basename), val(platform), val(file_label)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq
    to_merge/ | gzip -c > ${basename}_${platform}_${file_label}.fastq.gz
	"""

}

process RUN_HIFIASM {

	/* */

	tag "${basename}, ${file_label}"
	publishDir "${params.assembly}/${basename}_${file_label}", mode: 'copy', overwrite: true

    cpus 8

	input:
    tuple path(pb_fastq), val(basename), val(platform), val(file_label)
    tuple path(ont_fastq), val(basename), val(platform), val(file_label)

	output:
    tuple path("*"), val(basename), val(platform), val(file_label)

	script:
	"""
    hifiasm -o ${basename}_${file_label} -t6 --ul ${ont_fastq} ${pb_fastq}
	"""

}

process CONVERT_CONTIGS_TO_FASTA {

	/* */

	tag "${basename}, ${platform}, ${file_label}"
    label "map_and_extract"
	publishDir params.assembly, mode: 'copy', overwrite: true

    cpus 8

	input:
    tuple path("hifiasm_files/*"), val(basename), val(platform), val(file_label)

	output:
    path "${basename}_${file_label}.p_contigs.fasta"

	shell:
	'''
    awk '/^S/{print ">"$2"n"$3}' hifiasm_files/!{basename}_!{file_label}.bp.p_ctg.gfa \
    | fold > ${basename}_${file_label}.p_contigs.fasta
	'''

}

// --------------------------------------------------------------- //