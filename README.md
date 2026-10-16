# ovca-dc-tcell-response-bulk-rnaseq

Bulk RNA-seq analysis of T cells from an ovarian cancer / type 2 dendritic cell
study (Spranger lab; sequencing done at the Whitehead Institute Genome Technology Core).
Samples were sequenced in two batches about a year apart and analysed separately
(preprocessing) and together (differential expression, PCA, GSEA).

> **Status:** work in progress. The primary contrast (IL12-MSA vs delta123) is
> confounded with sequencing batch; see [Known issues](#known-issues).

## Experimental design

| Condition   | n | Batch | Notes                                      |
| ----------- | - | ----- | ------------------------------------------ |
| WT (WT_IP)  | 3 | 1     | role (baseline vs active wild-type) to confirm |
| delta123    | 2 | 1     | labelled `Delta3_IP` in some files         |
| Control     | 3 | 2     | Control_3 is a possible outlier            |
| Flt3L      | 3 | 2     |                                            |
| IL12-MSA    | 3 | 2     | high within-group variance                 |

Batch 1 also includes lymph node / tumor samples (`LN`, `Tumor`) processed in the
`WI_fionaac` run. Sample metadata carries a `Dataset` field (`IntratumorT`, `IP`).

## Repository layout

```text
072225_rnaseq/    Batch 1 (July 2025 run)
  nf-core_rnaseq_3.19.0.sh   SLURM + Nextflow launcher
  WI_fionaac.csv             nf-core samplesheet
  WI_fionaac/                nf-core output directory (output not tracked)
    tximport_mm/             per-batch Rmds (IP, LN/Tumor), figures, sample sheet
    gsea.IP/, gsea.LN.Tumor/ GSEA inputs (.rnk, .gmt) and prerank scripts
  documents/                 gene signature tables, .gmx sets, manuscript draft
  geo_submission/            GEO metadata sheets and md5 lists (no data)
082526_rnaseq/    Batch 2 (August 2026 run)
  nf-core_rnaseq_3.19.0.sh, 082526_samples.csv
  fionaac_082526/            nf-core output directory (output not tracked)
    tximport_mm/             intratumorT Rmd, DE tables, figures, sample sheet
    gsea/                    GSEA inputs and prerank scripts
  documents/                 gene signature tables
  geo_submission/            GEO metadata sheet and md5s
combined/         Cross-batch analysis
  analysis/*.Rmd             combined_processing, intratumorT_processing, WI_fionaac_IP_processing
  analysis/                  result tables (xlsx) and figures (PDF)
  opus_summary.md            design notes and open questions
Rcode/ (inside each tximport_mm/ and combined/analysis/)
                             ssGSEA scripts from the Broad Institute (own license)
```

Raw FASTQ files, nf-core output, GSEA run folders, GEO data directories and R
libraries are **not** tracked (about 300 GB); see `.gitignore`.

## Workflow

1. **Preprocess** (per batch): `nf-core/rnaseq` 3.19.0, STAR + Salmon, Singularity.
   Reference: mouse GRCm39, Ensembl release 113.

   ```bash
   cd 072225_rnaseq && sbatch nf-core_rnaseq_3.19.0.sh
   cd 082526_rnaseq && sbatch nf-core_rnaseq_3.19.0.sh
   ```

   Edit the genome paths in the scripts for your environment.
2. **Import and normalise** (R, tximport + DESeq2): `combined/analysis/*.Rmd`.
3. **Pathway analysis:** GSEA preranked (v4.4.0), fgsea on PCA loadings, g:Profiler on DE gene lists.
   Custom signatures are in the `documents/` folders.

R is run in a RStudio Singularity container (`combined/singularity_Rstudio_451.sh`, R 4.5.1).
Package versions are recorded in `combined/analysis/*_sessionInfo.txt`.

## Data availability

Raw and processed data are deposited at GEO (accession: _TBD_). Submission
metadata is in the `geo_submission/` folders.

## Known issues

- IL12-MSA (batch 2) and delta123 (batch 1) share no batch, so `~ batch + condition`
  is not full rank. Resolving this depends on whether WT_IP is an untreated
  baseline or active wild-type IL-12. Details in `combined/opus_summary.md`.
- Whether "IP" denotes administration route or the anatomical source of the T cells must be confirmed.
- Transcript-to-gene tables are named `mm10_ens113_*` while the pipeline uses GRCm39 (mm39); confirm this is intended.

## Contact

Charlie Whittaker, MIT Bioinformatics and Computing Core (<charliew@mit.edu>)

## License

See [LICENSE](LICENSE). Third-party code in `combined/analysis/Rcode/` is covered by its own license.
