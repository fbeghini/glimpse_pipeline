#!/usr/bin/env nextflow

/*
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/
ref = file(params.referenceGenome)

process align {
	cache "lenient"
	executor "slurm"

	input:
		tuple val(pair_id), path(reads)
		path(ref_fasta)
		path(ref_index)

	output:
		tuple val(pair_id), path("${pair_id}_aligned_reads.bam"), path("${pair_id}_aligned_reads.bam.bai")
		path("${pair_id}_fq2bam.log")
		path("${pair_id}_duplicate_metrics.txt") optional true

	publishDir "${params.out}/results/aligned_reads/", pattern: "*bam*", mode: "symlink"
	publishDir "${params.out}/results/logs/align/", pattern: "*{_fq2bam.log,_duplicate_metrics.txt}", mode: "copy"

	script:
	def readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:illumina\\tPM:novaseq\\tSM:${pair_id}"
	def markdups  = params.markDuplicates ? "--out-duplicate-metrics ${pair_id}_duplicate_metrics.txt" : "--no-markdups"
	def lowmem    = params.lowMemory ? "--low-memory" : ""
	"""
	mkdir -p pbrun_tmp

	pbrun fq2bam \\
		--ref ${ref_fasta} \\
		--in-fq ${reads[0]} ${reads[1]} \"${readGroup}\" \\
		--out-bam ${pair_id}_aligned_reads.bam \\
		--bwa-options=\"-K 100000000 -Y\" \\
		${markdups} ${lowmem} \\
		--num-gpus ${params.gpus} \\
		--tmp-dir pbrun_tmp \\
		--logfile ${pair_id}_fq2bam.log

	# fq2bam writes the BAM index itself; only index here if it somehow did not.
	if [ ! -f ${pair_id}_aligned_reads.bam.bai ]; then
		samtools index -@ ${task.cpus} ${pair_id}_aligned_reads.bam
	fi
	"""
}

process chunk {
	cache "lenient"
	executor "slurm"
	module 'BCFtools/1.16-GCCcore-10.2.0'
	input:
		tuple path(sites_vcf), path(sites_vcf_index)

	output:
		path("*.chunk.*.txt")
	
	publishDir "${params.out}/results/logs/chunk/", pattern: "*.chunks.log", mode: "move"

	"""
	n_chrom=`bcftools index -s ${sites_vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${sites_vcf} | cut -f1`

 	${params.chunk_exec} --input ${sites_vcf} --region \${chrom} --sequential --map ${params.glimpse_maps}/\${chrom}.b[0-9]*.gmap.gz --output \${chrom}.chunks.txt --threads ${task.cpus} > \${chrom}.chunks.log
	split -l 1 -d --additional-suffix=.txt \${chrom}.chunks.txt \${chrom}.chunk.
	"""	
}


process reference_by_chrom {
	// executor "local" // DEBUG
	cpus 1
	memory "5 GB" // DEBUG
	module 'BCFtools/1.16-GCCcore-10.2.0'
	input:
	tuple path(vcf), path(vcf_index)

	output:
	tuple stdout, path(vcf), path(vcf_index)

	"""
	n_chrom=`bcftools index -s ${vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
	"""
}

process split_reference {
	executor "slurm"
	memory "5 GB"
	cpus 1
		storeDir "${params.out}/results/split_reference"
	input:
		tuple val(chrom), path(chunk), path(ref_vcf), path(ref_vcf_index)

	output:
		tuple val(chrom), path("${chunk.baseName}_*.bin")

	"""
	# split -l 1 -d --additional-suffix=.txt ${chunk} ${chrom}.chunk.
	while IFS= read -r line
	do
	IRG=`cut -f3 <<< \$line`
	ORG=`cut -f4 <<< \$line`
	GLIMPSE2_split_reference --threads ${task.cpus} -R ${ref_vcf} --map ${params.glimpse_maps}/${chrom}.b[0-9]*.gmap.gz --input-region \${IRG} --output-region \${ORG} -O ${chunk.baseName}
	done < ${chunk}
	"""
}

process impute_chunks {
	errorStrategy "retry"
	maxRetries 3
	cache "lenient"
	executor "slurm"
	// executor "local" //DEBUG
	module 'parallel/20210322-GCCcore-10.2.0'
	memory { 25.GB * task.attempt }
	time "24h"

	input:
	tuple val(chrom), path(bins), val(pair_id), path(bam), path(bam_index)
	
	output:
		tuple val(pair_id), val(chrom), path("*.imputed.bcf")
		path("*.imputed.log")

	publishDir "results/logs/impute/", pattern: "*.imputed.log", mode: "move"
	publishDir "results/imputed_chunks/", pattern: "*.imputed.bcf", mode: "symlink"
	
	"""
	${params.parallel} --joblog log --resume-failed -j ${task.cpus} --trim rl ${params.phase_exec} --bam-file ${bam} --reference {} --threads 1 --output ${pair_id}.{}.imputed.bcf --log ${pair_id}.{}.imputed.log ::: ${bins}
	"""
}


process ligate_chunks {
	//errorStrategy "retry"
	//maxRetries 3
	cache	 "lenient"
	executor "slurm"
	memory "4 GB"
	time "12h"
	
	input:
		tuple val(pair_id), val(chrom), path(imputed_bcf), path(imputed_bcf_index)
	
	output:
		tuple val(pair_id), path("${pair_id}_${chrom}.imputed.bcf")
		path("${pair_id}_${chrom}.imputed.bcf.csi")

	"""
	for f in ${imputed_bcf}; do echo "\${f}"; done | sort -V > files_list.txt
	${params.ligate_exec} --input files_list.txt --output ${pair_id}_${chrom}.imputed.bcf --thread ${task.cpus}
	"""
}

process merge_chrom_sample {
	cache "lenient"
	memory "20 GB"
	module 'BCFtools/1.16-GCCcore-10.2.0'
	input:
		tuple val(pair_id), path(imputed_bcf)
	
	output:
		tuple val(pair_id), path('*.imputed.bcf*')

	publishDir "${params.out}/results/imputed", pattern: "${pair_id}.imputed.bcf*", mode: "symlink"

	"""
	for f in ${imputed_bcf}; do echo "\${f}"; done | sort -V > files_list.txt
	bcftools concat -f files_list.txt -Ob -o ${pair_id}.imputed.bcf
	bcftools index ${pair_id}.imputed.bcf
	"""
}

workflow {
	// Reference FASTA + its BWA/samtools index files are staged into the task dir so the
	// Parabricks container can see them (fq2bam needs *.amb, *.ann, *.bwt, *.pac, *.sa, *.fai).
	ref_index = Channel.fromPath("${params.referenceGenome}.{amb,ann,bwt,pac,sa,fai}", checkIfExists: true).collect()
	aligned_reads = align(Channel.fromFilePairs( params.reads, followLinks: true), Channel.value(ref), ref_index)[0]
	// aligned_reads = Channel.fromPath( params.study_bams, followLinks: true) | map { bam -> [bam.getBaseName(), bam, bam+".bai"]}
	chunks = Channel.fromPath(params.reference_sites_vcfs) | map{ vcf -> [vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".csi")] } | chunk
	reference_vcfs = Channel.fromPath(params.reference_vcfs).map{ file -> [file, file + (file.getExtension() == "bcf" ? ".csi" : ".tbi")] } | reference_by_chrom
	//chunks | transpose | flatten  | map { chr_chunk -> [chr_chunk.getSimpleName(), chr_chunk]} | combine(reference_vcfs, by: 0) | view
	bins = chunks | transpose | flatten | map { chr_chunk -> [chr_chunk.getSimpleName(), chr_chunk]} | combine(reference_vcfs, by: 0) | split_reference 
	// bins | groupTuple() | combine(aligned_reads) | view
	// aligned_reads | combine(bins) | groupTuple(by: [0,3]) | view
	// bins | groupTuple(by: 0, sort: true) | combine(aligned_reads) | view
	imputed_chunks = bins | groupTuple(by: 0, sort: true) | combine(aligned_reads) | impute_chunks
	// imputed_chunks[0] | view //transpose | map{pair_id, chr, imputed_bcf -> [pair_id, chr, imputed_bcf, imputed_bcf + ".csi"]} | groupTuple(by: [0,1], sort: true) | view
	ligated_sample = imputed_chunks[0] | transpose | map{pair_id, chr, imputed_bcf -> [pair_id, chr, imputed_bcf, imputed_bcf + ".csi"]} | groupTuple(by: [0,1]) | ligate_chunks
	merged_chroms = ligated_sample[0] | groupTuple | merge_chrom_sample
}

