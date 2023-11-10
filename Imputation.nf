#!/usr/bin/env nextflow

/*
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/
ref = file(params.referenceGenome)

process align {
	module 'BWA/0.7.17-GCCcore-10.2.0'
	module 'SAMtools/1.16-GCCcore-10.2.0'
	// executor "local"
	cpus 12

	input:
		tuple val(pair_id), path(reads)

	output:
		tuple val(pair_id), path("${pair_id}_aligned_reads.bam"), path("${pair_id}_aligned_reads.bam.bai")
	
	publishDir "${params.out}/results/aligned_reads/", pattern: "*bam*", mode: "symlink"

    script:
    readGroup = \
	"@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:illumina\\tPM:novaseq\\tSM:${pair_id}"
    """
    zcat ${reads[0]} ${reads[1]} | bwa mem \
	-K 100000000 \
	-v 3 \
	-t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$ref \
	- | samtools view -Sb - | samtools sort --threads ${task.cpus} -O BAM -o ${pair_id}_aligned_reads.bam

	samtools index -@ ${task.cpus} ${pair_id}_aligned_reads.bam
    """

}

process chunk {
	cache "lenient"
	// executor "local"
	module 'BCFtools/1.16-GCCcore-10.2.0'
	input:
	tuple path(sites_vcf), path(sites_vcf_index)

	output:
		tuple stdout, path("*.chunk.*.txt")
		path("*.chunks.log")
	
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
	printf "\${chrom}"
	"""	
}


process reference_by_chrom {
	// executor "local"
	cpus 1
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
	// executor "local"
	cpus 1
	publishDir "${params.out}/results/split_reference" , pattern: "*.bin", mode: "copy"
	input:
		tuple val(chrom), path(chunk), path(ref_vcf), path(ref_vcf_index)

	output:
		tuple val(chrom), path("*bin")

	"""
	IRG=`cut -f3 ${chunk}`
	ORG=`cut -f4 ${chunk}`
	GLIMPSE2_split_reference --threads ${task.cpus} -R ${ref_vcf} --map ${params.glimpse_maps}/${chrom}.b[0-9]*.gmap.gz --input-region \${IRG} --output-region \${ORG} -O bin
	"""
}

process impute_chunks {
	//errorStrategy "retry"
	//maxRetries 3
	cache "lenient"
	// executor "local"
	memory "12 GB"
	time "1h"
	
	input:
	tuple val(pair_id), path(bam), path(bam_index), val(chrom), path(bin)
	
	output:
	tuple val(pair_id), val(chrom), path("*.imputed.bcf")
	path("*.imputed.log")

	publishDir "results/logs/impute/", pattern: "*.imputed.log", mode: "move"
	publishDir "results/imputed_chunks/", pattern: "*.imputed.bcf", mode: "symlink"
	
	"""
	${params.phase_exec} --bam-file ${bam} --reference ${bin} --threads ${task.cpus} --output ${pair_id}.${bin.getBaseName()}.imputed.bcf --log ${pair_id}.${bin.getBaseName()}.imputed.log
	"""
}


process ligate_chunks {
	//errorStrategy "retry"
	//maxRetries 3
	cache "lenient"
	executor "local"
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
	memory "4 GB"
	module 'BCFtools/1.16-GCCcore-10.2.0'
	input:
		tuple val(pair_id), path(imputed_bcf)
	
	output:
		tuple val(pair_id), path('*.imputed.bcf*')

	publishDir "${params.out}/results/imputed", pattern: "${pair_id}.imputed.bcf*", mode: "symlink"

	"""
	for f in ${imputed_bcf}; do echo "\${f}"; done | sort -V > files_list.txt
	bcftools concat -f files.txt -Ob -o ${pair_id}.imputed.bcf
	bcftools index ${pair_id}.imputed.bcf
	"""
}

workflow {
	aligned_reads = align(Channel.fromFilePairs( params.reads, followLinks: true))
	chunks = chunk(Channel.fromPath(params.reference_sites_vcfs).map{ vcf -> [vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".csi")] })
	reference_vcfs = reference_by_chrom(Channel.fromPath(params.reference_vcfs).map{ file -> [file, file + (file.getExtension() == "bcf" ? ".csi" : ".tbi")] })
	bins = split_reference(chunks[0].take(2).transpose().combine(reference_vcfs, by: 0))

	imputed_chunks = impute_chunks(aligned_reads.combine(bins))
	ligated_sample = ligate_chunks(imputed_chunks[0].map{pair_id, chr, imputed_bcf -> [pair_id, chr, imputed_bcf, imputed_bcf + ".csi"]}.groupTuple(by: [0,1]))
	merged_chroms = merge_chrom_sample(ligated_sample[0].groupTuple())
}

