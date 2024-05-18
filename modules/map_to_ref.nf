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
    `realpath ${ref_fasta}` \
    `realpath ${fastq}` \
    | samtools view -Sbt `realpath ${ref_fasta}` \
    | samtools sort - -o ${basename}_${platform}.bam
	"""

}
