#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 32                      # Number of CPUs. Equivalent to the -pe whole_nodes 1 option in SGE
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=charliew@mit.edu  # Email to which notifications will be sent.

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate nf-core_v24
module add singularity/3.10.4

nextflow run nf-core/rnaseq -r 3.19.0 -profile singularity \
--input 082526_samples.csv \
--fasta /net/ostrom/data/bcc/charliew/Genomes/mm39_ens113/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--gtf /net/ostrom/data/bcc/charliew/Genomes/mm39_ens113/Mus_musculus.GRCm39.113.gtf \
--star_index '/net/ostrom/data/bcc/charliew/Genomes/nfcore_rnaseq/mm39_ens113/genome/index/star' \
--salmon_index '/net/ostrom/data/bcc/charliew/Genomes/nfcore_rnaseq/mm39_ens113/genome/index/salmon' \
--rsem_index '/net/ostrom/data/bcc/charliew/Genomes/nfcore_rnaseq/mm39_ens113/genome/rsem' \
--multiqc_title fionaac_082526 \
--aligner 'star_salmon' \
--pseudo_aligner 'salmon' \
--skip_stringtie \
--outdir fionaac_082526 \
