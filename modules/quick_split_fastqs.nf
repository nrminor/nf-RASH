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

    storeDir "$launchDir/${sample_id}_split_fastqs"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus params.cpus

    input:
    tuple path(big_ol_fastq), val(sample_id), val(platform)

    output:
    tuple path("split/*.fastq.gz"), val(platform)

    script:
    """
    seqkit split2 \
    --by-size ${params.split_max} \
    --extension ".gz" \
    --out-dir split/ \
    --threads ${task.cpus} \
    --force \
    ${big_ol_fastq}
    """

}