process MERGE_PACBIO_FASTQS {

	/*
    Now that we've split, mapped, and extracted regions from our input FASTQs,
    we need to merge them back together such that we have one PacBio FASTQ and
    one Oxford Nanopore FASTQ for each region. This process uses Seqkit to
    run these merges performantly on PacBio reads.
    */

    tag "${sample_id}, ${platform}, ${region}"
    label "seqkit"
	publishDir params.extracted, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
    tuple path("to_merge/???.fastq.gz"), val(sample_id), val(platform), val(region)

	output:
    tuple path("${sample_id}_${platform}_${region}.fastq.gz"), val(sample_id), val(platform), val(region)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq ./to_merge/ \
    | gzip -c > ${sample_id}_${platform}_${region}.fastq.gz
	"""

}

process MERGE_ONT_FASTQS {

	/*
    Now that we've split, mapped, and extracted regions from our input FASTQs,
    we need to merge them back together such that we have one PacBio FASTQ and
    one Oxford Nanopore FASTQ for each region. This process uses Seqkit to
    run these merges performantly on Oxford Nanopore reads.
    */

    tag "${sample_id}, ${platform}, ${region}"
    label "seqkit"
	publishDir params.extracted, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
    tuple path("to_merge/???.fastq.gz"), val(sample_id), val(platform), val(region)

	output:
    tuple path("${sample_id}_${platform}_${region}.fastq.gz"), val(sample_id), val(platform), val(region)

	script:
	"""
    seqkit scat \
    --threads ${task.cpus} \
    --find-only \
    --out-format fastq ./to_merge/ \
    | gzip -c > ${sample_id}_${platform}_${region}.fastq.gz
	"""

}
