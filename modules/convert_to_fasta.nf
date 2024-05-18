process CONVERT_CONTIGS_TO_FASTA {

	/*
    A cute trick from https://hifiasm.readthedocs.io/en/latest/faq.html#id1
    allows us to convert our contigs into a FASTA file that we can then
    inspect in Geneious (or wherever).
    */

	tag "${basename}, ${region}"
    label "map_and_extract"
	publishDir "${params.assembly}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 3

	input:
    tuple path("hifiasm_files/*"), val(basename), val(region)

	output:
    path "${basename}_${region}.p_contigs.fasta"

	shell:
	'''
    awk '/^S/{print ">"$2;print $3}' hifiasm_files/!{basename}_!{region}.bp.p_ctg.gfa \
    | fold > !{basename}_!{region}.p_contigs.fasta
	'''

}
