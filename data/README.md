# data/

Everything under `data/` is gitignored except this file. It records **where each file
came from, how it was made, and why**, so the folder can be rebuilt from raw downloads
plus code. Add an entry whenever a file is created, changed or found to be invalid.

- `reference/`: dataset-independent reference data (built once, read by every dataset).
- `datasets/<GSE>/raw/`: downloads, never modified.
- `datasets/<GSE>/processed/`: derived per-dataset files; regenerable.

---

## reference/raw/targetscan/

| File | Source | Notes |
|---|---|---|
| `UTR_Sequences.txt` | TargetScan 7 (exact sub-release 7.0/7.1/7.2 not recorded) | Tab-separated, **has a header row**: `Refseq ID`, `Gene ID`, `Gene Symbol`, `Species ID`, `UTR sequence`. `Refseq ID` actually holds Ensembl transcript IDs (ENST…). Cross-species alignment with `-` gaps. |
| `ORF_Sequences.txt` | TargetScan 7 (exact sub-release not recorded; file dated 2015-11-24) | Tab-separated, **no header**: transcript ID, taxon ID, aligned sequence (lowercase, `-` gaps). |

## reference/processed/

All files are human only (taxon 9606), gaps (`-`) stripped, uppercase DNA (ACGT).

### `human_3UTRs.tsv.gz`
Columns: `transcript_id`, `gene_id`, `gene_symbol`, `utr_sequence`, `utr_length`.
**How:** streamed `UTR_Sequences.txt` row by row, kept `Species ID == 9606`, stripped `-`,
uppercased, U→T. **Why stream:** the raw file is ~7 GB and must not be loaded whole.
28,352 rows. Earliest UTR step; superseded by the clean file.

### `human_3UTRs_clean.tsv.gz` (the UTR file to use)
Columns: as above plus `gene_id_clean` (Ensembl gene ID without the version suffix,
e.g. `ENSG00000004059.6` → `ENSG00000004059`, to match expression data). 28,351 rows,
one per transcript; 19,473 genes, 1–9 UTR transcripts per gene.
- 2026-11-13: created.
- 2026-11-23: dropped column `has_8mer`. **Why:** an old annotation from the siJARID2
  dataset (GSE184929), unrelated to the PARP siRNAs and misleading if reused.
  Written to a temp file, verified (shape and contents unchanged otherwise), then replaced.
- 1 row fewer than `human_3UTRs.tsv.gz`: one gene entry was removed because it had no
  Ensembl gene ID, so it could not be matched to the expression data.

### `human_ORFs_clean.tsv.gz`
Columns: `transcript_id`, `orf_sequence`. 28,351 transcripts (= all 9606 rows in the raw
file; no duplicates or missing values). Created 2026-11-23.
**Why:** to verify which siRNA strand is the sense strand (it should occur in the intended
target's ORF; the guide should not).
- Caveat: TargetScan ORFs are not exact CDSs. Only ~1.4 % have a length divisible by 3.
  E.g. PARP1 `ENST00000366794.5` matches RefSeq NM_001618.4 CDS for 3,044 nt, then drops
  the last stop-codon A and runs ~14 nt into the 3'UTR. Fine for substring searches, not
  for codon/frame work.

### `biotype_check/`
Old one-off check (2026-11-13), via an API call, of how many transcripts in the PARP1/2-related
annotation were protein-coding rather than non-coding or another biotype. Not used downstream.

---

## datasets/GSE114273/

### raw/
| File | Notes |
|---|---|
| `GSE114273_raw_counts_GRCh38.p13_NCBI.tsv.gz` | GEO NCBI-generated raw counts (GRCh38.p13). 18 samples: siRNA arm GSM3138707–16 (used) plus PARP-inhibitor arm GSM3138699–706 (not used). Use for DE. |
| `GSE114273_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` | GEO NCBI-generated TPM. Descriptive use only. |
| `Human.GRCh38.p13.annot.tsv` | NCBI gene annotation accompanying the GEO counts (GeneID, Symbol, EnsemblGeneID, GeneType, …). |
| `All_genes_GSE114273.csv` | GTO GTD0019 "All genes" table (downloaded 2026-11-24). 13,677 genes. **siPARP1-A + B pooled** (GSM3138709–12) vs siNC (GSM3138707–08). Gene symbols only; values rounded to 3 decimals; DE method and gene filter not documented by GTO. Not used for the main analysis (pooling dilutes each seed); kept for a later cross-check notebook. |
| `Differential_genes_GSE114273.csv` | GTO GTD0019 "Differential genes": the 998 genes of the table above with padj < 0.05 (no log2FC cut-off). Never use for a CDF (pre-selected on significance). |

### processed/: seed annotation (2026-11-23)
Made by `notebooks/GSE114273/02_sequence_analysis.ipynb` with the functions in
`src/offtarget/sequences.py`. Built from all human 3'UTRs, **not** filtered to expressed
genes: seed status depends only on the siRNAs and the UTRs, so the gene population
(expression filter) is applied later, when joining to the DESeq2 results.

| File | Contents |
|---|---|
| `GSE114273_sirna_table.tsv` | 4 × 8, index `sirna`. Sense (siPARP1-B corrected), guide, seed nt 2–7, seed nt 2–8, and the four target sites (DNA, UTR orientation, A1 appended). Made by `build_sirna_table`. |
| `GSE114273_seed_gene.tsv` | 19,473 genes × 37, index `gene_id_clean`. `n_utr_transcripts`; per siRNA × site type `_count` and `_present`; per siRNA `_best_site`. Made by `count_sites` (per UTR transcript, overlapping matches counted) then `sites_to_gene`. Load with `index_col=0`. |

- **Isoform rule ("any isoform"):** per gene, each count is the max over its UTR
  transcripts (isoforms share sequence, so summing would count the same site repeatedly).
- **Counts are nested:** an 8mer position also counts as 7mer_m8, 7mer_A1 and 6mer, so
  `_6mer_present` includes every gene with any site. `_best_site` is the single most
  stringent site (8mer > 7mer_m8 > 7mer_A1 > 6mer > none).
- Genes with an 8mer: siPARP1-A 1,292; siPARP1-B 1,091; siPARP2-A 1,088; siPARP2-B 690.

### processed/: differential expression (2026-11-24)
Made by `notebooks/GSE114273/03_differential_expression.ipynb`.

| File | Contents |
|---|---|
| `GSE114273_<siRNA>_DGE_results.tsv` (×4) | 19,227 protein-coding genes × 6, index `EnsemblGeneID`: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`. log2FC = siRNA / siNC (negative = lower after treatment). Load with `index_col=0`. |
| `GSE114273_on_target_DGE.tsv` | PARP1 and PARP2 rows from each contrast, with % remaining expression. |
| `GSE114273_counts_norm.tsv` | 19,227 genes × 10 samples, index `EnsemblGeneID`: DESeq2-normalised counts (raw counts / size factors, from `dds.layers["normed_counts"]`). Used for the siNC expression filter in NB 04 (main: siNC mean ≥ 10). |

- **Input:** GEO raw counts, the 10 siRNA samples, restricted to protein-coding genes with an
  Ensembl ID via `GSE114273_gene_mapping_protein_coding.tsv` (from NB 01).
- **Method:** PyDESeq2 0.5.4, design `~ condition`, siNC as reference. **One fit on all
  10 samples**, so dispersions are estimated across all 5 groups (more stable at n = 2).
  Then one Wald test per siRNA vs siNC. No LFC shrinkage. Cook's outlier refitting is on,
  but inactive with n = 2.
- **Why `~ condition` and not `~ replicate + condition`:** no evidence of paired replicates
  in GEO metadata or in the data (residual-correlation check). See
  `docs/methodology.md`, "Replicates and batch effects".
- **Sample QC:** replicates r = 0.996–0.998; each siRNA forms its own PCA cluster. See
  `docs/methodology.md`, "Sample QC: reading a PCA".
- **Gene counts:** 1,987 genes have zero reads in all samples (log2FC = NaN). The remaining
  17,240 have a log2FC; 12,900–14,600 per contrast have a padj (DESeq2 independent filtering).
  The CDF population is set later by the 5 % / 10 % expression filter, not by padj availability.
- Significant (padj < 0.05): siPARP1-A 6,029; siPARP1-B 3,936; siPARP2-A 5,431; siPARP2-B 7,806.
- Knockdown (log2FC of the target): siPARP1-A −6.52, siPARP1-B −5.69 (PARP1); siPARP2-A
  −3.05, siPARP2-B −4.16 (PARP2). Re-run on 2026-11-24 reproduced the 2026-11-18 results
  exactly (identical log2FC).

### processed/: seed analysis (2026-11-26)
Made by `notebooks/GSE114273/04_off_target_analysis.ipynb`.

| File | Contents |
|---|---|
| `GSE114273_seed_dge_gene.tsv` | 16,638 genes × 11, index `EnsemblGeneID`: `siNC_mean` (mean normalised count in siNC), `in_main` (siNC mean ≥ 10; 12,009 genes), `in_sensitivity` (lowest 10 % of siNC removed; 14,591), and per siRNA `<siRNA>_log2FoldChange` + `<siRNA>_best_site`. Load with `index_col=0`; select a gene set with `df[df["in_main"]]`. |

- Population = genes with a DESeq2 log2FC **and** a 3'UTR record in `seed_gene`; genes
  without a UTR are dropped, not counted as seed-negative.
- Result tables (`seed_stats`, specificity matrices, UTR-adjusted regression) and figures are
  in `results/GSE114273/` (tables tracked by git, figures not).

### Deleted invalid outputs (2026-11-26)
Built before two errors were found, and deleted once the corrected analysis existed:
1. **8-mer bug:** old NB3 built the 8mer as `'A' + 7mer` (A prepended, not appended),
   giving `ATTCATGA`, `AGTGTGAT`, `ATTGACAA`, `ATTCAGCT`.
2. **siPARP1-B typo:** deposited sense `…GUGUGAUA` has its last two nt swapped; corrected
   to `…GUGUGAAU` (19/19 match to PARP1 NM_001618.4 3086–3104). Changes the seed from
   AUCACA to UUCACA. Details in `config/GSE114273.yaml`.

| Deleted file | Affected by | Replaced by |
|---|---|---|
| `processed/GSE114273_expression_data_8mer_A1_annotated.tsv` | 1 and 2 | `GSE114273_seed_gene.tsv` (NB 02) |
| `processed/GSE114273_siPARP*_seed_analysis_input.tsv` (×4) | 1 and 2 | `GSE114273_seed_dge_gene.tsv` (NB 04) |
| `results/GSE114273/figures/GSE114273_CDF_baseline.png`, `…_CDF_filtered.png` | 1 and 2 | new CDF figure from NB 04 |

Not affected by sequence errors: DGE results, expression, knockdown and on-target files.
