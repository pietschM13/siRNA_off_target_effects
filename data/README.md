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
| `GSE114273_raw_counts_GRCh38.p13_NCBI.tsv.gz` | GEO NCBI-generated raw counts (GRCh38.p13). Use for DE. |
| `GSE114273_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz` | GEO NCBI-generated TPM. Descriptive use only. |
| `Human.GRCh38.p13.annot.tsv` | NCBI gene annotation accompanying the GEO counts (GeneID, Symbol, EnsemblGeneID, GeneType, …). |

### processed/: seed annotation (NB3, 2026-11-23)
Made by `notebooks/GSE114273/03_sequence_analysis.ipynb` with the functions in
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

### processed/: known invalid outputs
Built before two errors were found. Do not use; regenerate from the corrected NB3.
1. **8-mer bug:** old NB3 built the 8mer as `'A' + 7mer` (A prepended, not appended),
   giving `ATTCATGA`, `AGTGTGAT`, `ATTGACAA`, `ATTCAGCT`.
2. **siPARP1-B typo:** deposited sense `…GUGUGAUA` has its last two nt swapped; corrected
   to `…GUGUGAAU` (19/19 match to PARP1 NM_001618.4 3086–3104). Changes the seed from
   AUCACA to UUCACA. Details in `config/GSE114273.yaml`.

| File | Affected by |
|---|---|
| `GSE114273_expression_data_8mer_A1_annotated.tsv` | 1 and 2 |
| `GSE114273_siPARP*_seed_analysis_input.tsv` (×4) | Built with the wrong seeds; will be overwritten by the corrected analysis |

Not affected by sequence errors: DGE results, expression, knockdown and on-target files.
