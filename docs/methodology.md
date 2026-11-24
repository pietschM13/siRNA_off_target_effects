# Methodology notes

Short notes on analysis decisions: the reasoning (why) and how to check it on a new dataset.

---

## Replicates and batch effects: choosing the DESeq2 design

### The question
Each condition has replicates (rep1, rep2, ...). Before running DESeq2, decide whether the
replicates are **independent** or **paired**:

- **Independent:** every sample was prepared separately. "rep1" in siNC has nothing more in
  common with "rep1" in siPARP1-A than with "rep2". The numbering is just a label.
- **Paired (batch):** all rep1 samples were transfected, harvested or sequenced together (same
  day, same batch of cells), and all rep2 samples on another occasion. The replicate number
  then marks a real shared influence, a **batch effect**.

### Why it matters
- With a batch effect, every gene in a batch is shifted in the same direction. Under
  `~ condition` that shift counts as noise, which inflates the variance and hides real
  differences. `~ replicate + condition` estimates the batch shift and removes it, so each
  treatment is compared with its own control *within* the batch. That is more powerful.
- Without a batch effect, adding `replicate` to the design estimates something that isn't
  there. Each extra term costs a degree of freedom, and with n = 2 per condition there are
  very few to spare. The variance estimates get worse, not better.
- So the design must follow how the experiment was actually done, not the sample names.

### How to spot it
1. **Metadata first.** Check the GEO series "Overall design", sample characteristics and
   protocols, and the paper's methods for batch, date, run or lane fields, or wording such as
   "performed on separate days" or "paired". "Biological replicates" alone does not answer it.
2. **Then the data.** Remove each condition's mean from the log-expression. What remains is
   replicate-to-replicate variation. If rep1 samples share a batch, their residuals correlate
   with each other *across conditions*:
   ```python
   lc = np.log2(cpm + 1)                                    # filtered genes x samples
   resid = lc - lc.T.groupby(sample_to_condition).transform("mean").T
   resid.corr()   # compare: rep1 vs rep1 (different conditions) against rep1 vs rep2
   ```
   - Same-number pairs clearly more correlated than different-number pairs: batch effect, so
     use `~ replicate + condition`.
   - Both about zero: no detectable batch, so use `~ condition`.

   A PCA of the samples shows the same thing: if samples split by replicate number rather
   than by condition, there is a batch effect.
3. **Library sizes** (total reads per sample) that differ systematically between rep1 and
   rep2 are another hint, though DESeq2's size factors correct for depth itself.

With n = 2, residuals within one condition are mirror images (+x and −x), so the pattern
only shows up when comparing *across* conditions, which is why the check is set up that way.

### GSE114273 (siPARP1/2 in LNCaP): result
- **Metadata:** GEO overall design says "with biological duplicates". Samples are titled
  `siNC-rep1`, `siNC-rep2`, ... with no batch or date fields. Everything was processed
  identically (72 h siRNA, TopHat2/hg19, HTSeq, edgeR in the original paper). Nothing indicates
  pairing.
- **Data:** residual correlation (genes with ≥ 1 CPM in at least half of the samples,
  n = 13,778), mean over pairs from different conditions:

  | Comparison | Mean residual correlation |
  |---|---|
  | Same replicate number (rep1 vs rep1) | −0.03 |
  | Different replicate number (rep1 vs rep2) | +0.03 |

  Both about zero, so there is no shared replicate batch. Library sizes are even (17.9–21.5 M reads).
- **Decision:** `~ condition`, siNC as reference. Checked 2026-11-24.

---

## Sample QC: reading a PCA

### What it is
Each sample is a list of expression values, one per gene, so it is a point in a space with
one axis per gene (about 12,000 axes here). PCA finds the directions in that space along which
the samples differ most:
- **PC1** is the single biggest difference across all samples, **PC2** the next biggest that is
  independent of PC1, and so on.
- Each PC is a weighted combination of genes that change together across samples. The
  weights are called **loadings**.
- The % on each axis is the share of all sample-to-sample variation that this PC captures.

A PCA plot is therefore a 2D map of the samples that keeps the biggest differences visible.
**Distance means difference.** The axis values have no unit; only relative positions matter.

### How to use it
- **Quality:** replicates should sit close together. A replicate far from its partner is
  either an outlier or a mislabelled or swapped sample.
- **Batch:** if samples group by replicate number or processing day rather than by
  condition, there is a batch effect (see the section above).
- **Biology:** which conditions move together, and which have their own direction.

Run it on log2-CPM of reasonably expressed genes (here ≥ 1 CPM in at least half of the
samples), so read depth and near-zero genes don't dominate. This is descriptive only;
DESeq2 works on raw counts.

### GSE114273 result (2026-11-24)
11,870 genes. Variance explained: PC1 53 %, PC2 21 %, PC3 14 %.

| Axis | Separates |
|---|---|
| PC1 | siPARP2-B from everything else |
| PC2 | siPARP1-A |
| PC3 | siPARP2-A |
| none | siPARP1-B stays close to siNC |

- Replicates nearly overlap (r = 0.996–0.998): good quality, no outliers.
- **siRNAs against the same gene do not group together.** If expression changes were driven
  mainly by losing the target, siPARP2-A and -B (and siPARP1-A and -B) would move in the same
  direction from siNC. Instead each siRNA has its own dominant signature. That is consistent
  with large siRNA-specific effects such as seed-mediated off-targeting, but PCA does not
  show the cause. Other siRNA-specific explanations (knockdown strength, immune response)
  remain possible.
- The very tight replicates also explain the high number of DESeq2-significant genes at
  n = 2: small variance estimates make modest changes significant.

---

## Possible extensions (not implemented)

- **PC loadings as supplementary support.** For each siRNA-specific PC (PC1 = siPARP2-B,
  PC2 = siPARP1-A, PC3 = siPARP2-A), take the genes with the largest loadings in the
  direction of downregulation and test whether they are enriched for that siRNA's seed sites.
  These genes could be compared with the strongly downregulated genes highlighted in the
  volcano plots. Caveat: loadings capture up- and downregulated genes together, so the sign
  must be read relative to where siNC lies on that PC.
