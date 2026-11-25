# Genotype imputation using low depth sequencing data

## 1. Description
Pipeline for genotype imputation from low depth sequencing data using [GLIMPSE1](https://github.com/odelaneau/GLIMPSE). This pipeline uses [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) for estimating genotype probabilities at reference sites prior to imputation. The pipeline was developed using [Nextflow](https://www.nextflow.io/) and was tested on [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler.

## 2. Prerequisites
The following software is required:
- Singularity
- Nextflow
- bcftools
- One or more NVIDIA GPUs with a recent driver (required by the `align` step, which uses NVIDIA Clara Parabricks `pbrun fq2bam`)

## 3. Installation
To run this pipline you will need:
1. Download and install [GLIMPSE1](https://github.com/odelaneau/GLIMPSE)
2. Build GATK singularity container: `singularity build gatk_VERSION.sif docker://broadinstitute/gatk:VERSION`
3. Build the Clara Parabricks singularity container (used for alignment):
   `singularity build parabricks_4.8.0-1.sif docker://nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1`
4. Index the reference genome for BWA and samtools, next to the FASTA file:
   ```
   bwa index genome.fasta      # creates genome.fasta.{amb,ann,bwt,pac,sa}
   samtools faidx genome.fasta # creates genome.fasta.fai
   ```
   `pbrun fq2bam` needs all of these; the pipeline stages them into the task directory alongside the FASTA.
5. Clone this repo: `git clone https://github.com/CERC-Genomic-Medicine/glimpse_pipeline.git`

## 4. Execution
1. Modify `nextflow.config` configuration file:
* `params.reference_vcfs` -- path to VCF/BCF files with phased reference panel genotypes. Each VCF/BCF file must have the corresponding tbi/csi index.
* `params.reference_sites_vcfs` -- path to sites-only VCF/BCF files of the reference panel. Each VCF/BCF file must have the corresponding tbi/csi index.
* `params.study_bams` -- path to BAM/CRAM files. One BAM/CRAM file per study participant. Each BAM/CRAM file must have the corresponding bai/crai index.
* `params.referenceGenome` -- name of the reference genome *.fa file (e.g. hs37d5.fa).
* `params.parabricksContainer` -- path to the Clara Parabricks singularity image file (.sif) used by the `align` step.
* `params.gpus` -- number of GPUs per alignment task. Passed to `pbrun fq2bam --num-gpus` and requested from SLURM via `process.withName:'align'.clusterOptions`.
* `params.markDuplicates` -- `true` (default) lets `fq2bam` mark duplicates and write a duplicate-metrics file; `false` adds `--no-markdups`, reproducing the previous BWA-MEM behaviour of leaving duplicates unmarked.
* `params.lowMemory` -- `true` adds `fq2bam --low-memory`, which reduces the number of streams per GPU. Use it if alignment runs out of GPU memory.
* `params.window_size` -- imputation window size in base-pairs. This is a parameter to the `GLIMPSE_chunk` executable. See [GLIMPSE1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_chunk) for more details.
* `params.buffer_size` -- imputation window buffer size in base-pairs. This is a parameter to the `GLIMPSE_chunk` executable. See [GLIMPSE1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_chunk) for more details.
* `params.chunk_exec` -- path to the `GLIMPSE_chunk` executable.
* `params.phase_exec` -- path to the `GLIMPSE_phase` executable.
* `params.ligate_exec` -- path to the `GLIMPSE_ligate` executable.
* `params.glimpse_maps` -- path to the GLIMPSE's genetic maps folder with the corresponding human genome build version.
* `process.*` and `executor.*` -- set this arguments according to your compute cluster configuration. In particular, adjust `withName: 'align'`'s `clusterOptions` (GPU partition name and `--gpus`/`--gres=gpu:` syntax) to match your scheduler: Nextflow's `accelerator` directive is **not** honoured by the SLURM executor, so GPUs have to be requested through `clusterOptions`.

2. Run pipleine. Example of interactive SLURM job:
```
salloc --time=1-0 -p day -cpus-per-task 10 --ntasks=1 --mem-per-cpu=16G
module load nextflow
module load bcftools
nextflow run Imputation.nf
```
