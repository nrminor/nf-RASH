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
    tuple path(pb_fastq), path(ont_fastq), val(basename), val(region)

	output:
    tuple path("*"), val(basename), val(region)

	script:
    assert pb_fastq.toString().toLowerCase().contains("pacbio")
    assert ont_fastq.toString().toLowerCase().contains("ont")
	"""
    hifiasm -o ${basename}_${region} -t ${task.cpus} --ul ${ont_fastq} ${pb_fastq}
	"""

}

process RUN_HIFIASM_HIFI_ONLY {

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
    tuple path(pb_fastq), val(basename), val(region)

	output:
    tuple path("*"), val(basename), val(region)

	script:
    assert pb_fastq.toString().toLowerCase().contains("pacbio")
	"""
    hifiasm -o ${basename}_${region} -t ${task.cpus} ${pb_fastq}
	"""

}
