process MAP_TO_REF {

	/*
    Users provide a reference with the parameter `ref_fasta` in
    `nextflow.config`, which is used to bait out reads that are likely to 
    contain regions of interest.
    */

	tag "${sample_id}, ${platform}"
    label "map_and_extract"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
    tuple path(fastq), val(sample_id), val(platform)
    each path(ref_fasta)

	output:
    tuple path("*.bam"), val(sample_id), val(platform)

	script:
    minimap2_preset = platform == "pacbio" ? "map-hifi" : "map-ont"
	chunk_id = fastq.getSimpleName()
	"""
    minimap2 -t ${task.cpus} -L --eqx -ax ${minimap2_preset} \
    `realpath ${ref_fasta}` \
    `realpath ${fastq}` \
    | samtools view -Sbt `realpath ${ref_fasta}` \
    | samtools sort - -o ${sample_id}_${chunk_id}_${platform}.bam
	"""

}
