#!/usr/bin/env nextflow

/*	This is nextflow pipeline to process Illumina paired-end data from Mycobacterium bovis isolates  
*	It will first deduplicate the reads using fastuniq, trim them using trimmomatic and then map to the reference genome
*	
*
*	written by ellisrichardj, based on pipeline developed by Javier Nunez
*
*	Version 0.1.0	31/07/18	Initial version
*/


params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"
params.ref = "/home/richard/MyScripts/references/AF2122.fna"

ref = file(params.ref)


/*	Collect pairs of fastq files and infer sample names */
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs } 

/* remove duplicates from raw data */
process Deduplicate {
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file(forward), file(reverse) from read_pairs

	output:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") into dedup_read_pairs

	"""
	gunzip -c ${forward} > ${pair_id}_R1.fastq 
	gunzip -c ${reverse} > ${pair_id}_R2.fastq
	echo '${pair_id}_R1.fastq\n${pair_id}_R2.fastq' > fqin.lst
	~/FastUniq/source/fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
	"""
}	

/* trim adapters and low quality bases from fastq data */
process Trim {
	errorStrategy 'ignore'

	maxForks 1

	input:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") from dedup_read_pairs

	output:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_read_pairs

	"""
	java -jar ~/MyScripts/trimmomatic-0.30.jar PE -threads 4 -phred33 ${pair_id}_uniq_R1.fastq ${pair_id}_uniq_R2.fastq  ${pair_id}_trim_R1.fastq ${pair_id}_fail1.fastq ${pair_id}_trim_R2.fastq ${pair_id}_fail2.fastq ILLUMINACLIP:/home/richard/ReferenceSequences/adapter.fasta:2:30:10 SLIDINGWINDOW:10:20 MINLEN:36
rm ${pair_id}_fail1.fastq
rm ${pair_id}_fail2.fastq
	"""
}

/* map to reference sequence */
process Map2Ref {
	errorStrategy 'ignore'

	maxForks 1

	input:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from trim_read_pairs

	output:
	set pair_id, file("${pair_id}.mapped.sorted.bam") into mapped_bam

	"""
	bwa mem -T10 -t4 $ref  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq |
	 samtools view -@4 -ShuF 4 - |
	 samtools sort -@4 - -o ${pair_id}.mapped.sorted.bam
	"""
}

/* Variant calling */
process VarCall {
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.mapped.sorted.bam") from mapped_bam

	output:
	set pair_id, file("${pair_id}.pileup.vcf.gz") into vcf

	"""
	samtools index ${pair_id}.mapped.sorted.bam
	samtools mpileup -q 10 -uvf $ref ${pair_id}.mapped.sorted.bam |
	 bcftools call --ploidy Y -cf GQ - -Oz -o ${pair_id}.pileup.vcf.gz
	"""
}

/* Consensus calling */
process VCF2Consensus {
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.pileup.vcf.gz") from vcf

	output:
	set pair_id, file("${pair_id}_consensus.fas") into consensus

	"""
	bcftools index ${pair_id}.pileup.vcf.gz
	bcftools consensus -f $ref -o ${pair_id}_consensus.fas ${pair_id}.pileup.vcf.gz
	"""
}

/* Mapping Statistics */


/* SNP filtering and annotation */



workflow.onComplete {
		log.info "Completed sucessfully:	$workflow.success"		
		log.info "Nextflow Version:	$workflow.nextflow.version"
		log.info "Duration:		$workflow.duration"
		log.info "Output Directory:	$params.outdir"
}
