# CLAUDE.md: sirna_shrna_offtargets

## Project purpose

Test whether genes with a 3'UTR seed match to a therapeutic siRNA/shRNA **guide strand** are preferentially downregulated after treatment, using public RNA-seq data. Build a reproducible workflow that can be repeated across datasets.

The hypothesis is **miRNA-like** off-targeting. The seeds are theoretical, derived from the guide strands; they are not established biological regulators, and the mechanism is indifferent to the siRNA's intended target. What is tested is whether a seed sequence is associated with a downregulatory tendency across a population of thousands of genes.

Current priority: **one reliable end-to-end siRNA proof of concept.** No web app, shRNA expansion or machine learning until that works and there is a clear scientific reason.

Full context and open questions: `rna-offtarget-project-handover.md` (Claude project doc, not in this repo). Read its migration notes first.


- Reference data is dataset-independent and expensive to rebuild, so it is kept separate from per-dataset data.
- `data/` holds inputs and regenerable intermediates; `results/` holds output for write-ups.
- Notebooks are for analysis and interpretation. Move reusable logic into `src/offtarget/` as it stabilizes.
- Never overwrite raw data.
- Everything under `data/` is gitignored and reproducible from raw downloads plus code. Record how in `data/README.md`.
- Notebooks live in `notebooks/<GSE>/`, two levels below the project root, so paths use `../../` (e.g. `../../data/...`) and shared code is imported after `sys.path.append("../../src")`. A central `paths.py` may come later.

## Environment

- WSL Ubuntu, VS Code, Python virtual environment at `.venv/`, Git/GitHub.
- Known past issues: kernel/interpreter confusion, Python 3.13 compatibility problems, Git "dubious ownership" error in WSL, large files breaking pushes.
- Reference sequence files are multi-GB. Stream them in a single pass; never load them whole.

## Workflow

Follow **EXPLORE → PLAN → IMPLEMENT → TEST → REVIEW**:

1. **Explore:** read the relevant code, notebooks and data structures before proposing anything.
2. **Plan:** describe the intended change, assumptions and likely failure points. Wait for my approval before editing files.
3. **Implement:** make small, modular changes. Avoid restructuring that wasn't asked for.
4. **Test:** run the relevant tests and check outputs (shapes, counts, value ranges).
5. **Review:** check the biological and statistical logic, not just whether the code runs.

I'm developing my computational skills. Explain **why** an implementation is appropriate, go step by step on unfamiliar concepts, and don't give me too many steps at once.

## Code conventions

- Modular, reusable functions with docstrings that state assumptions.
- Validate inputs. Handle missing values explicitly. Raise informative errors; never silently swallow exceptions.
- Use `pytest` for tests and `ruff` for linting. Priority for tests: sequence and site-construction logic, where an error already cost a full round of results.
- **Keep agreed variable names.** Never silently change a variable's name or definition; flag it if something looks wrong.
- Code that runs is not necessarily correct. Always check:
  - dataframe operations do what is intended (e.g. `isin()` filtering vs `merge()`)
  - sequence orientation (sense vs guide, site vs seed)
  - consistent Ensembl ID cleaning and mapping
  - filters do not introduce bias
  - duplicate UTR records and isoforms are handled by an explicit rule, never silently collapsed

## Naming convention

Dataframes: `<content>_<unit>_<state>`, lowercase snake_case.

- content: `counts`, `annot`, `utr`, `seed`, `dge`
- unit: `tx` (transcript-level), `gene` (gene-level)
- state: `raw`, `tpm`, `pc` (protein-coding), `expr` (mapped to expressed genes), `filtXX` (expression-filtered at XX %)
- A dict holding one dataframe per siRNA takes the suffix `_by_sirna`.
- The dataset ID belongs in file names, not variable names.

Agreed names: `counts_raw`, `counts_tpm`, `annot_raw`, `annot_pc`, `utr_tx`, `utr_tx_expr`, `seed_gene`, `tpm_gene_seed`, `dge_by_sirna`, `seed_dge_by_sirna`, `sirna_sense`, `sirna_guide`, `sirna_seed` (nt 2–7), `sirna_seed_m8` (nt 2–8), `sites_8mer`, `sites_7mer_m8`, `sites_7mer_A1`, `sites_6mer`.

Columns: `<siRNA>_<sitetype>_present`, `<siRNA>_<sitetype>_count`, `<siRNA>_best_site`.

`counts_tpm` is TPM, which is **not** DESeq2 normalisation. Raw counts for differential expression; TPM for descriptive comparisons only.

## Analysis design

### Stage 1: exploratory
Normalized expression → log2FC → seed-positive vs seed-negative → CDF plot → preliminary Wilcoxon rank-sum. Exploratory only.

### Stage 2: main analysis
Raw gene-level counts → DESeq2 with the correct design (paired if samples are paired) → expression filtering → log2FC → map to human 3'UTRs → seed status → CDF + rank-sum test.

### Expression filters
Filter on expression in the **negative control (siNC) only**: mean DESeq2-normalised count of the siNC samples. Never filter on treated samples, because low counts there can be the repression being tested (selection on the outcome).

- **Main:** siNC mean ≥ 10 normalised counts (≈ 0.5 CPM; the DESeq2 and edgeR `filterByExpr` default). All figures and conclusions use this.
- **Sensitivity (one, pre-declared):** remove the lowest 10 % of siNC means (in GSE114273 this removes only genes with 0 reads in siNC). Report as one short table of effect sizes: does the conclusion change?
- The same filtered gene set is used for all siRNA contrasts (needed for the specificity control).
- Percentile filters on all genes (5 % / 10 %) were dropped: in GSE114273 they cut below 1 read and left ~3,900 genes with < 10 reads, whose log2FC is 2–3× noisier.

Reasoning and numbers: `docs/methodology.md`, "Expression filter for the seed analysis". Any change to these thresholds must be agreed with me first.

### Seed definitions

Seed = guide-strand nucleotides **2–7**, counted from the 5' end. Nucleotide 8 is the additional m8 position. The target site is the 3'UTR sequence the seed pairs with, i.e. the reverse complement of the seed, written 5'→3' in UTR orientation.

| Type | Definition (mRNA 5'→3') |
|---|---|
| 8-mer | match to guide nt 2–8 + A opposite nt 1 (A at **3' end** of site) |
| 7mer-m8 | match to guide nt 2–8 |
| 7mer-A1 | match to guide nt 2–7 + A opposite nt 1 |
| 6-mer | match to guide nt 2–7 |

The A1 adenosine is **not encoded in the guide**. It is an unpaired A in the target opposite guide nt 1, which sits at the 3' end of the site in UTR orientation. Append it, never prepend it: `site + 'A'`, not `'A' + site`.

Never use the label "8mer-A1" — the 8mer already includes A1. Keep "seed" for the guide-side sequence and "target site" for the UTR-side sequence.

The four site populations are **nested, not disjoint**: the 8mer site contains both 7mer sites, and each 7mer contains the 6mer. Store all four flags plus a `best_site` column giving the most stringent site present.

Build sites by function from the guide strand. Never hard-code per-siRNA literals.

Keep these definitions distinct in code. The hypothesis that stronger sites give stronger repression (8-mer > 7-mer > 6-mer) is to be **tested**, not assumed.

### Guide strand determination
Guide = reverse complement of the deposited sense strand, **verified**, not assumed. Check against the TargetScan ORF data: a sense strand should occur in the intended target ORF; a guide strand should not. Sequence-based checks cannot establish which strand was actually loaded into RISC.

### Reference sequence data
TargetScan `UTR_Sequences.txt` and `ORF_Sequences.txt` are cross-species alignments, and their layouts differ:
- `UTR_Sequences.txt` **has a header row**: `Refseq ID` (actually Ensembl transcript IDs), `Gene ID`, `Gene Symbol`, `Species ID`, `UTR sequence`.
- `ORF_Sequences.txt` has **no header**. Columns: transcript ID, species taxon ID, aligned sequence.

Filter to human (**taxon 9606**), strip `-` gap characters, and uppercase before matching. Extract each subset once by streaming, save to `data/reference/processed/`, and record the TargetScan version used.

### On-target check
Confirm and document knockdown of the intended target before interpreting off-target effects. Descriptive, on TPM. With n=2 it is a sanity check, not a statistical result.

## Statistical rules

- **Wilcoxon rank-sum (Mann–Whitney U):** seed-positive vs seed-negative genes, which are independent groups. The unit of analysis is the gene, not the sample.
- **Wilcoxon signed-rank / paired tests:** only for paired samples. Never for gene populations.
- Keep these distinctions explicit: raw vs normalized counts, exploratory vs confirmatory analysis, paired vs independent observations, gene-level changes vs seed-associated effects.
- Genes are **not** independent (co-regulation, downstream effects of target knockdown), so rank-sum p-values are anti-conservative. **Always report effect sizes**, not p-values alone.
- With n=2, per-gene log2FC is noisy, which attenuates a true shift towards zero and makes the seed-level test conservative. The power comes from the number of genes, not the number of replicates.
- A distributional shift is an association, not proof that seed matching caused repression.
- Known confounders to be aware of: 3'UTR length (longer UTRs contain more matches by chance), expression level, isoforms and alternative polyadenylation, multiple testing across seed definitions.


### Specificity control
Each siRNA has a different seed, so a given siRNA's seed-positive gene set should shift in its own contrast and **not** in the other three. The gene set is held fixed across contrasts, so this is not confounded by 3'UTR length, GC content or expression level. Implement as a seed-set × contrast matrix of effect sizes. This is the primary control against the long-3'UTR confound; UTR-length-matched comparison is secondary.

## Reproducibility

Every result should record: dataset ID, treatment and control definitions, siRNA and guide sequences, seed definition, expression filter, DE method, gene ID convention, UTR source and version, statistical test, and n genes per group. Per-dataset values live in `config/<GSE>.yaml`.

## Datasets

- **GTO** (Gene Therapy Omnibus) = our own catalogue of RNA therapeutic datasets. Not the same as **GEO**, which is the underlying public repository.
- **Active:** GSE114273 — siPARP1-A/B, siPARP2-A/B against siNC, n=2 per condition.
- **Test only:** GSE164413. Its processed files predate the naming convention; don't let that spread.
- **Historical:** GSE184929 (siJARID2 in hASC). Used for early on-target checks; no longer in use.

### PARP siRNA sequences (GSE114273)

Sense and guide are RNA 5'→3'; target sites are DNA 5'→3', matching the TargetScan UTR sequences. Guide = reverse complement of sense, pending ORF verification.

| siRNA | Sense (as deposited) | Guide | Seed nt 2–7 | Seed nt 2–8 | 8mer | 7mer-m8 | 7mer-A1 | 6mer |
|---|---|---|---|---|---|---|---|---|
| siPARP1-A | CCAUUGAGCACUUCAUGAA | UUCAUGAAGUGCUCAAUGG | UCAUGA | UCAUGAA | TTCATGAA | TTCATGA | TCATGAA | TCATGA |
| siPARP1-B | GAUUUCAUCUGGUGUGAAU (corrected) | AUUCACACCAGAUGAAAUC | UUCACA | UUCACAC | GTGTGAAA | GTGTGAA | TGTGAAA | TGTGAA |
| siPARP2-A | GAAGAAAUUCCUUGACAAA | UUUGUCAAGGAAUUUCUUC | UUGUCA | UUGUCAA | TTGACAAA | TTGACAA | TGACAAA | TGACAA |
| siPARP2-B | AAGUACUAUCUGAUUCAGCUA | UAGCUGAAUCAGAUAGUACUU | AGCUGA | AGCUGAA | TTCAGCTA | TTCAGCT | TCAGCTA | TCAGCT |

siPARP2-B is 21 nt beginning AA, consistent with the older "AA + 19 nt" target-site convention. Its seed is unaffected, but ORF verification must tolerate the leading dinucleotide.

**siPARP1-B correction.** The sense strand as recorded in GTO, `GAUUUCAUCUGGUGUGAUA`, occurs in no human ORF or 3'UTR. Swapping the last two nt (`…GAUA` → `…GAAU`) gives a 19/19 match to PARP1: NM_001618.4 mRNA 3086–3104 (CDS 165..3209, so CDS position 2922), which is the same site as TargetScan ORF `ENST00000366794.5` at 0-based index 2921. NCBI BLAST returns this as the only match. Treated as a transposition typo. Because the swapped nt pair with guide nt 1–2, the seed changes (deposited AUCACA → corrected UUCACA). Any siPARP1-B seed or site results built from the deposited sequence are invalid. Whether the typo is also present in GEO is not yet checked. Both sequences are recorded in `config/GSE114273.yaml`.

Strand verification (TargetScan ORFs, all 28,351 human transcripts): for all four siRNAs (siPARP1-B corrected), the sense strand occurs only in its intended target and the guide occurs in no ORF.

> ⚠️ Earlier NB3 code built the 8-mer as `'A' + target_7mer`, giving the wrong motifs `ATTCATGA`, `AGTGTGAT`, `ATTGACAA`, `ATTCAGCT`. These scan a non-canonical 8-mer with an extra A opposite guide position 9. Any output built on them is invalid. If you find these strings, flag them to me before changing anything.

## Terminology

- "Seed-positive" / "seed-negative" = gene whose 3'UTR does / does not contain the relevant site.
- "Seed" = guide-side sequence. "Target site" = UTR-side sequence. These are reverse complements of each other and must not be conflated in variable names.
EOF
