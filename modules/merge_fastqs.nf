process MERGE_PACBIO_FASTQS {

	/*
    Now that we've split, mapped, and extracted regions from our input FASTQs,
    we need to merge them back together such that we have one PacBio FASTQ and
    one Oxford Nanopore FASTQ for each region. This process uses Seqkit to
    run these merges performantly on PacBio reads.
    */

    tag "${basename}, ${platform}, ${region}"
    label "seqkit"
	publishDir params.extracted, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path("to_merge/???.fastq.gz"), val(basename), val(platform), val(region)

	output:
    tuple path("${basename}_${platform}_${region}.fastq.gz"), val(basename), val(platform), val(region)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq ./to_merge/ \
    | gzip -c > ${basename}_${platform}_${region}.fastq.gz
	"""

}

process MERGE_ONT_FASTQS {

	/*
    Now that we've split, mapped, and extracted regions from our input FASTQs,
    we need to merge them back together such that we have one PacBio FASTQ and
    one Oxford Nanopore FASTQ for each region. This process uses Seqkit to
    run these merges performantly on Oxford Nanopore reads.
    */

    tag "${basename}, ${platform}, ${region}"
    label "seqkit"
	publishDir params.extracted, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path("to_merge/???.fastq.gz"), val(basename), val(platform), val(region)

	output:
    tuple path("${basename}_${platform}_${region}.fastq.gz"), val(basename), val(platform), val(region)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq ./to_merge/ \
    | gzip -c > ${basename}_${platform}_${region}.fastq.gz
	"""

}
