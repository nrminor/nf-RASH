process CONVERT_CONTIGS_TO_FASTA {

	/*
    A cute trick from https://hifiasm.readthedocs.io/en/latest/faq.html#id1
    allows us to convert our contigs into a FASTA file that we can then
    inspect in Geneious (or wherever).
    */

	tag "${sample_id}, ${region}"
    label "map_and_extract"
	publishDir "${params.assembly}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
    tuple path("hifiasm_files/*"), val(sample_id), val(region)

	output:
    path "${sample_id}_${region}.p_contigs.fasta"

	shell:
	'''
    awk '/^S/{print ">"$2;print $3}' hifiasm_files/!{sample_id}_!{region}.bp.p_ctg.gfa \
    | fold > !{sample_id}_!{region}.p_contigs.fasta
	'''

}
