process EXTRACT_REGIONS {

	/* 
    Users provide a TSV of regions of interest with the parameter 
    `desired_regions` in `nextflow.config`. Each row of that TSV is split into
    its own task in queue. Every input FASTQ to this process is thus
    multiplexed across regions called for in the TSV; if there are N rows, each
    FASTQ will be run through this process N times to extract reads mapping to
    those N regions.
    */

	tag "${sample_id}, ${platform}, ${region}"
    label "map_and_extract"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    tuple path(bam), val(sample_id), val(platform), val(expression), val(region), val(merge_key)

	output:
    tuple path("*_${merge_key}.fastq.gz"), val(sample_id), val(platform), val(merge_key)

	script:
	bam_id = bam.getSimpleName()
	"""
    samtools index ${bam}
    samtools view -b ${bam} ${expression} \
    | samtools fastq - \
    | reformat.sh qin=33 int=f in=stdin.fq \
    out=${bam_id}_${merge_key}.fastq.gz
	"""

}
