process EXTRACT_REGIONS {

	/* 
    Users provide a TSV of regions of interest with the parameter 
    `desired_regions` in `nextflow.config`. Each row of that TSV is split into
    its own task in queue. Every input FASTQ to this process is thus
    multiplexed across regions called for in the TSV; if there are N rows, each
    FASTQ will be run through this process N times to extract reads mapping to
    those N regions.
    */

	tag "${basename}, ${platform}, ${region}"
    label "map_and_extract"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

	input:
    each path(bam)
    tuple val(expression), val(region), val(merge_key)

	output:
    tuple path("${basename}_${platform}_${merge_key}.fastq.gz"), val(basename), val(platform), val(merge_key)

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
    out=${basename}_${platform}_${merge_key}.fastq.gz
	"""

}
