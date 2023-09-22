#! /bin/bash
module load miniconda
conda activate nextflow
cd /home/fb343/git/glimpse_pipeline
nextflow run Imputation.nf -resume
