# CLAUDE.md: sirna_shrna_offtargets

## Project purpose

Test whether genes with a 3'UTR seed match to a therapeutic siRNA/shRNA **guide strand** are preferentially downregulated after treatment, using public RNA-seq data. Build a reproducible workflow that can be repeated across datasets.

Current priority: **one reliable end-to-end siRNA proof of concept.** No web app, shRNA expansion or machine learning until that works and there is a clear scientific reason.

## Repository structure

Expected layout. Verify against the actual repo before making changes, and update this section if it differs.

```
sirna_shrna_offtargets/
├── notebooks/   # exploratory analysis and interpretation, numbered by stage
├── data/        # raw and processed data; large/external files are gitignored
├── src/         # reusable functions and processing logic
├── tests/       # pytest tests
├── README.md
└── .gitignore
```

- Notebooks are for analysis and interpretation. Move reusable logic into `src/` as it stabilizes.
- Keep raw and processed data clearly separated. Never overwrite raw data.
- Never commit large datasets (e.g. TargetScan `UTR_Sequences.txt`, ~7 GB, or raw count matrices).

## Environment

- WSL Ubuntu, VS Code, Python virtual environment at `.venv/`, Git/GitHub.
- **Do not assume the active interpreter or kernel is correct.** Check that it is `.venv/bin/python` before installing packages or running code. Never install into the base environment.
- Known past issues: kernel/interpreter confusion, Python 3.13 compatibility problems, Git "dubious ownership" error in WSL, large files breaking pushes.

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
- Use `pytest` for tests and `ruff` for linting.
- **Keep agreed variable names.** Never silently change a variable's name or definition; flag it if something looks wrong.
- Code that runs is not necessarily correct. Always check:
  - dataframe operations do what is intended (e.g. `isin()` filtering vs `merge()`)
  - sequence orientation (sense vs guide, site vs seed)
  - consistent Ensembl ID cleaning and mapping
  - filters do not introduce bias
  - duplicate UTR records and isoforms are handled by an explicit rule, never silently collapsed

## Analysis design

### Stage 1: exploratory
Normalized expression → log2FC → seed-positive vs seed-negative → CDF plot → preliminary Wilcoxon rank-sum. Exploratory only.

### Stage 2: main analysis
Raw gene-level counts → DESeq2 with the correct design (paired if samples are paired) → expression filtering → log2FC → map to human 3'UTRs → seed status → CDF + rank-sum test.

### Expression filters
- **Strict:** detected in 3/3 replicates in both conditions.
- **Inclusive:** detected in ≥2/3 replicates in both conditions.
- **Sensitivity:** a less stringent threshold, as a robustness check.

Document the exact threshold used.

### Seed definitions
Seed = guide-strand nucleotides 2–8. Site types (mRNA 5'→3'):

| Type | Definition |
|---|---|
| 8-mer | match to guide nt 2–8 + A opposite nt 1 (A at 3' end of site) |
| 7mer-m8 | match to guide nt 2–8 |
| 7mer-A1 | match to guide nt 2–7 + A opposite nt 1 |
| 6-mer | match to guide nt 2–7 |

Keep these definitions distinct in code. The hypothesis that stronger sites give stronger repression (8-mer > 7-mer > 6-mer) is to be **tested**, not assumed.

### UTR reference
TargetScan `UTR_Sequences.txt`, which covers multiple species. Extract **human only** once, efficiently (streaming, not loading the whole file repeatedly), and save the subset to `data/processed/`. Record the TargetScan version used.

### On-target check
Confirm and document knockdown of the intended target before interpreting off-target effects.

## Statistical rules

- **Wilcoxon rank-sum (Mann–Whitney U):** seed-positive vs seed-negative genes, which are independent groups.
- **Wilcoxon signed-rank / paired tests:** only for paired samples. Never for gene populations.
- Keep these distinctions explicit: raw vs normalized counts, exploratory vs confirmatory analysis, paired vs independent observations.
- A distributional shift is an association, not proof that seed matching caused repression.
- Known confounders to address: 3'UTR length (longer UTRs contain more matches by chance), expression level, isoforms and alternative polyadenylation, multiple testing across seed definitions.
- Biological caveats: guide loading is preferential, not absolute; the passenger strand may also cause off-targets (analyze separately, don't assume it's inactive); not every seed match causes repression.

## Reproducibility

Every result should record: dataset ID, treatment and control definitions, siRNA and guide sequences, seed definition, expression filter, DE method, gene ID convention, UTR source and version, statistical test, and n genes per group.

## Datasets

- **GTO** (Gene Therapy Omnibus) = our own catalogue of RNA therapeutic datasets. Not the same as **GEO**, which is the underlying public repository.
- **Candidate:** GTD0003 → GSE106695 (human breast cancer cells, siGRHL2, 4 control / 6 treated). Not yet confirmed: verify replicate structure, paired vs independent design, siRNA sequence and raw-count availability.
- **Earlier test:** GSE184929 (siJARID2 in hASC, day 6 and day 13, 3 paired experiments each). Used for on-target checks and development.

### JARID2 sequences

```python
# siRNA (sense):  GAAGAACGGGUGGUACGUA   (corrected; original was missing a G, BLAST-checked)
# guide (DNA):    TACGTACCACCCGTTCTTC
jarid2_seed_7mer   = "ACGTACC"    # guide nt 2–8
jarid2_target_7mer = "GGTACGT"    # 7mer-m8 site
jarid2_target_8mer = "GGTACGTA"   # 8-mer site
```

> ⚠️ Earlier code may contain `jarid2_target_8mer = "AGGTACGT"`. That string puts the A at the wrong end of the site. If you find it, flag it to me before changing it.

## Terminology

- Use **"seed-sharing miRNAs"** for miRNAs with identical seed sequences, not "miRNA family proxy".
- "Seed-positive" / "seed-negative" = gene whose 3'UTR does / does not contain the relevant site.
