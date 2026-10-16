# Bulk RNA-seq: IL12-MSA vs delta123 — Design Notes & Analysis Constraints

Owner: Charlie Whittaker (MIT)
Status: pre-analysis; two blocking questions unresolved (see §3)
Purpose: T-cell bulk RNA-seq, 5 arms, 14 samples. Primary contrast IL12-MSA vs delta123.

---

## 1. Experimental design

| Condition   | n | Batch | Notes                                  |
|-------------|---|-------|----------------------------------------|
| WT (WT_IP)  | 3 | **1** | role ambiguous — see §3               |
| delta123    | 2 | **1** | labeled `Delta3_IP` / `Δ1+2+3` in plot |
| Control     | 3 | **2** | Control_3 is a possible outlier        |
| Flt3L       | 3 | **2** |                                        |
| IL12-MSA    | 3 | **2** | high within-group variance             |

Batches sequenced **~1 year apart**, same facility.
Sample metadata has a `Dataset` field with values `IntratumorT` and `IP`.

---

## 2. CRITICAL CONSTRAINT: primary contrast is confounded with batch

IL12-MSA (batch 2) and delta123 (batch 1) share no batch.
`design = ~ batch + condition` **will fail** with:

    Error in checkFullRank(modelMatrix) : the model matrix is not full rank

Reason: batch 1 == {WT, delta123}, so batch is a linear combination of the
condition indicators.

**Do not "fix" this by dropping the batch term and reporting the naive
2-group comparison as the primary result.** Every log2FC would be
(true biology + unknown batch offset).

The escape hatch: each batch contains a baseline-like arm (WT in b1,
Control in b2), which makes batch estimable *if* WT is a true baseline (§3).

---

## 3. BLOCKING QUESTIONS — resolve before writing analysis code

### Q1. Is `WT_IP` an untreated baseline, or wild-type IL-12 protein delivered IP?
- **Untreated baseline** → merge WT + Control into a single `baseline` level;
  `~ batch + grp` becomes full rank. Proceed with Option A (§5).
- **Active wild-type IL-12 (parental of delta123)** → **DO NOT MERGE.**
  WT is a treatment arm; merging would absorb real biology into the batch
  coefficient. Instead WT_IP is the correct within-batch comparator for
  delta123, and the contrast becomes mutant-vs-parental (Option B, §5).

### Q2. Does "IP" mean route of administration, or anatomical source of T cells?
If the cells come from different compartments (peritoneal vs intratumoral),
that is a **worse confounder than sequencing batch and is not statistically
correctable.** Plot title "Intratumor T-Cells + IP Combined" needs ruling out.

### Q3 (flag, not blocker). Possible label problem.
If WT_IP is active wild-type IL-12, it sits at far-left PC1 with the
*unresponsive* samples, while the delta123 mutant sits shifted right toward
IL12-MSA. A mutant appearing more active than its parental molecule is either
a major result or a mislabel. Verify delta123 / WT_IP against submission sheets.

---

## 4. PCA findings (14 samples, VST, top-500 default)

PC1 = 44%, PC2 = 24%.

Observed structure:
- Control, WT, Flt3L cluster together at low PC1.
- delta123 and IL12-MSA separate from baselines on PC1.
- **PC2 separates delta123 (+28 to +31) from IL12-MSA (−10 to −21).**
- Control_3 sits at PC1 ≈ +13, displaced toward the treated region.

### Conclusion: the batch effect appears modest, and PC2 is very likely NOT batch

Evidence (WT is delta123's batch-mate, so WT position estimates the batch offset):

| Observation | Interpretation |
|---|---|
| WT (b1) PC2 ≈ −1 to +2.5; delta123 (b1) PC2 ≈ +28 to +31 | batch main effect on PC2 ≈ 0 → **PC2 is condition-driven** |
| Within batch 1, WT PC1 ≈ −20, delta123 PC1 ≈ +19 | PC1 spans nearly full range within one batch → not batch |
| WT (b1) clusters with Control/Flt3L (b2) | cross-batch baseline bridge is real; batch offset on PC1 ≈ 10–13 vs ~50-unit condition span |
| delta123 pair tight (Δ≈4 both PCs) | low within-group variance despite n=2 |

---

## 5. Analysis options (pick per §3 Q1)

### Option A — WT is a true baseline: merge baselines, `~ batch + grp`

```r
library(DESeq2)
coldata$$grp <- as.character(coldata$$condition)
coldata$$grp[coldata$$grp %in% c("WT_IP","Control")] <- "baseline"
coldata$$grp   <- factor(coldata$$grp, levels=c("baseline","Flt3L","IL12_MSA","delta123"))
coldata$$batch <- factor(coldata$$batch)

dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ batch + grp)
keep <- rowSums(counts(dds) >= 10) >= 2      # smallest group = delta123, n=2
dds  <- dds[keep, ]
dds  <- DESeq(dds)

res <- results(dds, contrast=c("grp","IL12_MSA","delta123"), alpha=0.05)
res <- lfcShrink(dds, contrast=c("grp","IL12_MSA","delta123"), res=res, type="ashr")
```

---

## 6. PC2 loading heatmap (top/bottom 10 genes)

Uses `vsd` and `loadings` from the PC1/PC2 loadings section of `combined_processing.Rmd`
(`loadings` has one row per top-variance gene, with `Row.names` = gene_id, `gene_name`,
and the `PC1`/`PC2` loading values). Genes are ordered by PC2 loading, most positive
(delta123 direction) to most negative (IL12-MSA direction); the top 10 and bottom 10
of that ordering are plotted with both dendrograms off, so row order stays exactly the
PC2-loading order and isn't rearranged by clustering.

```r
library(ComplexHeatmap)
library(circlize)

## Order all loaded genes by PC2, most positive to most negative, and take the
## 10 most positive (delta123-direction) and 10 most negative (IL12-MSA-direction)
pc2_ranked   <- loadings[order(loadings$PC2, decreasing = TRUE), ]
pc2_extreme  <- rbind(head(pc2_ranked, 10), tail(pc2_ranked, 10))

## VST expression for those genes, z-scored per gene (row) across samples
expr_mat <- assay(vsd)[pc2_extreme$Row.names, ]
expr_z   <- t(scale(t(expr_mat)))
rownames(expr_z) <- make.unique(ifelse(is.na(pc2_extreme$gene_name),
                                        pc2_extreme$Row.names,
                                        pc2_extreme$gene_name))

## Row split so the top-10/bottom-10 blocks stay visually distinct while row
## order within each block is preserved (not reclustered)
row_grp <- factor(rep(c("PC2 high (delta123)", "PC2 low (IL12-MSA)"), each = 10),
                   levels = c("PC2 high (delta123)", "PC2 low (IL12-MSA)"))

ht <- Heatmap(
  expr_z,
  name             = "z-score",
  col              = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_order        = rownames(expr_z),
  row_split        = row_grp,
  row_title_rot    = 0,
  column_title     = "Top/bottom 10 genes by PC2 loading",
  row_names_gp     = grid::gpar(fontsize = 9)
)

pdf("heatmap_PC2loading_top10bottom10.pdf", width = 6, height = 8)
draw(ht)
dev.off()
ht
```