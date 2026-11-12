# ovca-dc-tcell-response-bulk-rnaseq

Bulk RNA-seq analysis of T cells from an ovarian cancer / type 2 dendritic cell
study (Spranger lab; sequencing done at the Whitehead Institute Genome Technology Core).
Samples were sequenced in two batches about a year apart and analysed separately
(preprocessing) and together (differential expression, PCA, GSEA).

## Experimental design

| Condition  | n | Batch | Notes                                         |
| ---------- | - | ----- | --------------------------------------------- |
| WT (WT_IP) | 3 | 1     | role: baseline vs active wild-type; confirm   |
| delta123   | 2 | 1     | labelled `Delta3_IP` in some files            |
| Control    | 3 | 2     | Control_3 is a possible outlier               |
| Flt3L      | 3 | 2     |                                               |
| IL12-MSA   | 3 | 2     | high within-group variance                    |

Batch 1 also includes lymph node / tumor samples (`LN_1-3`, `Tumor_1-3`; 11 samples
in total) processed in the `WI_fionaac` run. The 14 samples in the table above are
the ones used in the cross-batch T-cell analysis. Sample metadata carries a
`Dataset` field (`IntratumorT`, `IP`).

## Repository layout

```text
072225_rnaseq/    Batch 1 (July 2025 run)
  nf-core_rnaseq_3.19.0.sh   SLURM + Nextflow launcher
  WI_fionaac.csv             nf-core samplesheet (11 samples: LN, Tumor, WT_IP, Delta3_IP)
  250711_WIGTC-AVITI1A_*/    sequencing core lane annotation and library fragment sizes
  FASTQ/                     md5sums.txt tracked; FASTQ files are not
  WI_fionaac/                nf-core output directory (output not tracked)
    tximport_mm/             per-dataset Rmds (IP, LN/Tumor), figures, sample sheet
    gsea.IP/, gsea.LN.Tumor/ GSEA inputs (.rnk, .gmt) and prerank scripts
  documents/                 gene signature tables and .gmx sets
  geo_submission/            GEO metadata sheets and md5 lists (no data); also holds
                             metadata for the companion single-cell submission
082526_rnaseq/    Batch 2 (August 2026 run)
  nf-core_rnaseq_3.19.0.sh, 082526_samples.csv
  260820_WIGTC-AVITI1B_*/    sequencing core lane annotation and library fragment sizes
  fionaac_082526/            nf-core output directory (output not tracked)
    tximport_mm/             intratumorT Rmd, DE tables, figures, sample sheet
    gsea/                    GSEA inputs (.rnk, .gmt) and prerank scripts
  documents/                 gene signature tables
  geo_submission/            GEO metadata sheet and md5s
combined/         Cross-batch analysis
  analysis/*.Rmd             combined_processing, intratumorT_processing, WI_fionaac_IP_processing
  analysis/                  result tables (xlsx, txt) and figures (PDF)
  analysis/gsea/             GSEA inputs and scripts (run folders not tracked)
  opus_summary.md            design notes and open questions
Rcode/            ssGSEA / GenePattern scripts from the Broad Institute (own license),
                  sourced by the Rmds as 'Rcode/...' (see Workflow, step 2)
annotationFiles/  tx2gene / gene2typesym tables (tracked); Ensembl 113 GRCm39 GTF (not tracked)
singularity_Rstudio_4_6_1.sh   SLURM launcher for RStudio Server in a Singularity container
```

Raw FASTQ files, nf-core output, GSEA run folders, GEO data directories, the
Ensembl GTF and R libraries are **not** tracked (about 165 GB of data);
see `.gitignore`. Unpublished manuscript drafts are also excluded.

## Workflow

1. **Preprocess** (per batch): `nf-core/rnaseq` 3.19.0, STAR + Salmon, Singularity.
   Reference: mouse GRCm39, Ensembl release 113.

   ```bash
   cd 072225_rnaseq && sbatch nf-core_rnaseq_3.19.0.sh
   cd 082526_rnaseq && sbatch nf-core_rnaseq_3.19.0.sh
   ```

   Edit the genome and index paths in the scripts for your environment.
2. **Import and normalise** (R, tximport + DESeq2): the `.Rmd` files in
   `combined/analysis/` (cross-batch) and in each batch's `tximport_mm/`.
   The Rmds read `mm10_ens113_tx2gene.xlsx` and `mm10_ens113_gene2typesym.xlsx`
   from their working directory; the code that builds them from the Ensembl 113
   GTF is commented out near the top of each Rmd. They also `source('Rcode/...')`,
   so make the top-level `Rcode/` visible from the Rmd's directory
   (e.g. `ln -s ../../Rcode Rcode` from `combined/analysis/`).
3. **Pathway analysis:** GSEA preranked (v4.4.0), fgsea on PCA loadings, g:Profiler on DE gene lists.
   Custom signatures are in the `documents/` folders.

R 4.6.1 runs in an RStudio Server Singularity container
(`singularity_Rstudio_4_6_1.sh`, image `bumproo/r4_6_1_singlecell_bulk_rnaseq:2026-10-15`).
Package versions are recorded in `combined/analysis/*_sessionInfo.txt` and the
per-batch `tximport_mm/*_sessionInfo.txt`.

## Data availability

Raw and processed data are deposited at GEO (accession: _TBD_). Submission
metadata is in the `geo_submission/` folders.

## Contact

Charlie Whittaker, MIT Bioinformatics and Computing Core (<charliew@mit.edu>)

## License

See [LICENSE](LICENSE). Third-party code in `Rcode/` is covered by its own license.
