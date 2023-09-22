for sample in /gpfs/ycga/scratch60/christakis/fb343/gatk/out_hg38/aligned_reads/*sam
do
for vcf in ~/project/databases/20201028_3202_phased/*.sites.vcf.gz
do
samplename=$(basename ${sample})
samplename=${sample%%.*}
chr=${vcf%%.*}
chr=$(rev <<< $chr | cut -f1 -d'_' | rev)
echo "bcftools mpileup --threads 12 -f /gpfs/ycga/datasets/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -I -E -a 'FORMAT/DP' -T ${vcf} ${sample} -Ou | bcftools call -Aim -C alleles -T ${vcf/vcf/tsv} -Oz -o ${samplename}_${chr}.vcf.gz && bcftools index -f ${samplename}_${chr}.vcf.gz" >> dsq_mpileup
done
done