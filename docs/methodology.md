# Methodology notes

Short notes on analysis decisions: the reasoning (why) and how to check it on a new dataset.
Each section starts with a **Context** line: where in the workflow the question came up and
what was being decided, so it can be read on its own later.

---

## Replicates and batch effects: choosing the DESeq2 design

> **Context (2026-11-24, notebook 03 differential expression):** before rebuilding the DESeq2
> step for GSE114273 we had to decide between the design `~ condition` and
> `~ replicate + condition`. The paper only says "two biological replicates", which does not
> answer whether rep1 of every condition was processed together. Outcome: `~ condition`.

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

> **Context (2026-11-24, notebook 03 differential expression):** a sample QC step (correlation
> heatmap and PCA) was added before running DESeq2. The PCA plots were not self-explanatory,
> so this section explains how to read one and what it showed for GSE114273.

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

## Analysis population for the seed CDF

> **Context (2026-11-24, after inspecting GTO's DGE tables):** GTO's precomputed DGE for
> GSE114273 pools siPARP1-A and -B into one "case" group and includes non-coding genes. That
> raised two questions: what pooling does to a per-siRNA seed analysis, and which genes
> should make up the population compared in the CDF. Outcome: our own per-siRNA DESeq2;
> population = genes with a 3'UTR, filtered on control expression; negative group = no site.

### Why per siRNA, never pooled
Each siRNA has its own seed. In a pooled contrast (siRNA-A + siRNA-B vs control), a gene with
a site only for A is repressed in half the samples, so its log2FC, and each seed set's shift,
is roughly **halved**. Pooling also turns each siRNA's own effects into within-group noise.
Effects shared by both siRNAs (the on-target biology) survive pooling; siRNA-specific effects
(the off-target signal) are averaged away. That is why experimenters use two siRNAs, and why
pooled tables (GTO) suit on-target questions, not seed analysis.

### Why the composition of the population matters
The CDF compares log2FC of seed-positive vs seed-negative genes. It only means something if
the two groups differ in seed status and in nothing else that affects log2FC. The population
defines what "seed-negative" means, so it defines the null.
- **Only genes with an annotated 3'UTR.** Seed status comes from TargetScan 3'UTRs (protein-
  coding transcripts). Genes without one (ncRNAs, pseudogenes, ...) would all land in the
  negative group by absence of data, not because they lack sites, and they behave differently
  (expression level, noise). "Protein-coding" is a stand-in; the real requirement is a known
  3'UTR. Genes without one are dropped, not counted as negative.
- **An expression filter** removes genes whose log2FC is mostly noise (next section). It
  depends only on expression, not on seed status, so it cannot bias group membership.

### Seed-negative = no site of any type
Seed-positive groups are defined by `best_site` (8mer, 7mer-m8, 7mer-A1, 6mer). The negative
group is genes with **no site at all** (`best_site == "none"`), not "lacks this site type".
Otherwise the negative group would contain 7mer and 6mer genes, which may also be weakly
repressed, and the difference being tested would shrink. This follows TargetScan-style
analyses. (The old notebook compared 8mer genes against all others.)

### Numbers (GSE114273)
16,638 genes have both a 3'UTR record and a DESeq2 log2FC (19,227 protein-coding genes;
959 without a UTR, 1,987 with zero counts everywhere). `baseMean` is identical in all four
contrasts (one model fit), so one gene population serves all siRNAs.

---

## Expression filter for the seed analysis

> **Context (2026-11-24, planning notebook 04 off-target analysis):** the plan inherited from
> the project handover was to exclude the lowest 5 % and 10 % of genes by expression. On
> inspection these cut-offs were almost meaningless for GSE114273, which led to a discussion of
> what to filter on and how strictly. Outcome: siNC mean ≥ 10 (main), lowest 10 % of siNC
> (sensitivity).

### Why the 5 % / 10 % percentile filter was dropped
In the 16,638-gene population the 5 % and 10 % cut-offs sit at baseMean 0.2 and 0.5, i.e.
below one read. A quarter of the genes have fewer than 10 reads, and their log2FC scatters
2–3 times more than well-expressed genes (SD of log2FC: 1.4–1.75 below 10 reads vs 0.57 above
50). Such noise doesn't create a false shift, but it blurs both CDFs and dilutes a real one.
The percentile idea assumed the low end is small; in a cell line many genes are simply off.

### Filter on the control, never on treated samples
- The hypothesis is that seed-matched genes go **down** after treatment. A filter on treated-
  sample expression preferentially removes genes that went down: selection on the outcome,
  which pushes the result towards "no effect".
- The relevant question is "is the gene expressed enough in untreated cells to be
  repressible?", which only the negative control (siNC) answers.
- A control-only filter is independent of every siRNA's effect and defines one gene set for
  all contrasts, which the specificity control needs.

### Group-wise vs average (baseMean) filters
A baseMean cut-off averages over all samples. A gene at 50 reads in the two siNC samples and 0
in the eight siRNA samples has baseMean 10 and could be cut, although it is clearly expressed
in one condition. Group-wise filters (edgeR `filterByExpr`) instead ask whether the gene
reaches ≥ 10 reads (depth-scaled) in at least as many samples as the smallest group, in any
group. That matters for general DE analysis. For a repression analysis the group that matters
is the control, so filtering on siNC is the group-wise idea applied to the one relevant group.

### Is ≥ 10 counts too strict?
No, it is the field default. With ~19 M reads per siNC sample, 10 reads ≈ 0.5 CPM, which is
barely detected. DESeq2's manual suggests ≥ 10 reads in at least the smallest group size;
edgeR's `filterByExpr` uses 10 by default; the original authors used ≥ 1 CPM (≈ 20 reads).
Removing ~28 % of protein-coding genes is normal: a cell line expresses roughly 11–13k of them,
and many tissue-specific genes are off in LNCaP.

| Filter (mean normalised count in siNC) | Cut-off | Genes kept (of 16,638) |
|---|---|---|
| Lowest 10 % removed | 0.45 | 14,591 |
| Lowest 20 % removed | 2.05 | 13,294 |
| ≥ 5 | 5 | 12,596 |
| **≥ 10 (main)** | 10 | **12,009** |
| ≥ 50 | 50 | 10,605 |

9.9 % of the population has 0 reads in both siNC samples, so "lowest 10 %" removes almost
exactly the undetected genes.

### Sensitivity check: the middle way
- **Why have one:** any threshold is partly arbitrary; a pre-declared check shows the result
  doesn't depend on it, guards against unconsciously tuning the threshold, and is
  informative (a shift that grows with a stricter filter fits noise diluting a real effect).
- **Why only one:** every extra threshold multiplies the outputs and the comparisons, and
  invites reporting whichever looks best. With a field-standard main threshold, more checks
  add little.
- **Rule:** decided before seeing results. One main threshold for all figures and conclusions;
  one sensitivity threshold reported once as a short table of effect sizes ("does the
  conclusion change?").
- **Why "lowest 10 % of siNC" as the sensitivity check:** it is the most permissive principled
  filter (only removes genes not detected in the control), so it tests whether removing the
  noisy low-count genes with ≥ 10 changed the conclusion. A stricter check (≥ 50) would test
  the opposite worry and is a one-line change if wanted.

---

## Possible extensions (not implemented)

> **Context (2026-11-24, after the PCA in notebook 03):** the PCA showed each siRNA has its
> own dominant expression signature. The idea below came up as a possible supplementary
> analysis linking that observation to seed sites; not scheduled.

- **PC loadings as supplementary support.** For each siRNA-specific PC (PC1 = siPARP2-B,
  PC2 = siPARP1-A, PC3 = siPARP2-A), take the genes with the largest loadings in the
  direction of downregulation and test whether they are enriched for that siRNA's seed sites.
  These genes could be compared with the strongly downregulated genes highlighted in the
  volcano plots. Caveat: loadings capture up- and downregulated genes together, so the sign
  must be read relative to where siNC lies on that PC.
