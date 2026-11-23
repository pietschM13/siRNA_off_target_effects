from Bio.Seq import Seq
import numpy as np
import pandas as pd

# Site types, most stringent first (this order decides best_site).
SITE_TYPES = ["8mer", "7mer_m8", "7mer_A1", "6mer"]


def build_sirna_table(sirna_sense):
    """Build guide, seeds and target sites for one or more siRNAs.

    sirna_sense: dict {siRNA name: sense strand, RNA 5'->3'}.
                 For a single siRNA, pass {"siX": "ACGU..."}.
    Assumes guide = reverse complement of sense (verified separately on the ORFs).
    Returns one row per siRNA. Target sites are DNA 5'->3' (UTR orientation).
    """
    rows = []
    for name, sense in sirna_sense.items():
        if set(sense) - set("ACGU"):
            raise ValueError(f"{name}: sense strand must be uppercase RNA (ACGU), got {sense}")

        guide = str(Seq(sense).reverse_complement_rna())
        seed = guide[1:7]       # guide nt 2-7
        seed_m8 = guide[1:8]    # guide nt 2-8

        site_6mer = str(Seq(seed).reverse_complement_rna()).replace("U", "T")
        site_7mer_m8 = str(Seq(seed_m8).reverse_complement_rna()).replace("U", "T")

        rows.append({
            "sirna": name,
            "sirna_sense": sense,
            "sirna_guide": guide,
            "sirna_seed": seed,
            "sirna_seed_m8": seed_m8,
            "site_8mer": site_7mer_m8 + "A",    # A1 appended at the 3' end
            "site_7mer_m8": site_7mer_m8,
            "site_7mer_A1": site_6mer + "A",
            "site_6mer": site_6mer,
        })
    return pd.DataFrame(rows).set_index("sirna")


def count_sites(utr_tx, sirna_table):
    """Count target sites in each 3'UTR transcript, for every siRNA and site type.

    utr_tx: one row per transcript, with columns 'transcript_id', 'gene_id_clean' and
            'utr_sequence' (uppercase DNA, UTR orientation 5'->3').
    sirna_table: output of build_sirna_table (columns 'site_8mer', 'site_7mer_m8', ...).
    Returns one row per transcript: 'transcript_id', 'gene_id_clean' and one
    '<siRNA>_<sitetype>_count' column per siRNA x site type. No sequences are kept.

    Overlapping occurrences are counted.

    Counts are NESTED (cumulative), not exclusive. Each count is the number of positions
    matching AT LEAST that site's definition:
        8mer     seed nt 2-8 match + A1        -> also counted in all three below
        7mer_m8  seed nt 2-8 match             -> includes every 8mer position
        7mer_A1  seed nt 2-7 match + A1        -> includes every 8mer position
        6mer     seed nt 2-7 match             -> includes every position above
    So at one position, an 8mer adds 1 to all four counts; a 7mer_m8 or 7mer_A1 adds 1
    to itself and to 6mer; a bare 6mer adds 1 to 6mer only. Hence, per transcript:
        8mer <= 7mer_m8 <= 6mer   and   8mer <= 7mer_A1 <= 6mer
    The two 7mers are not nested in each other: they differ at one end (nt 8 vs A1).

    Exclusive counts (positions whose best site is exactly that type), if needed:
        8mer only     = 8mer
        7mer_m8 only  = 7mer_m8 - 8mer
        7mer_A1 only  = 7mer_A1 - 8mer
        6mer only     = 6mer - 7mer_m8 - 7mer_A1 + 8mer
    (8mers are added back in the last line because they were subtracted twice.)
    """
    if utr_tx["transcript_id"].duplicated().any():
        raise ValueError("utr_tx has duplicate transcript_id values")
    if utr_tx["utr_sequence"].isna().any():
        raise ValueError("utr_tx has missing utr_sequence values")

    seed_tx = utr_tx[["transcript_id", "gene_id_clean"]].copy()
    for sirna in sirna_table.index:
        for site_type in SITE_TYPES:
            site = sirna_table.loc[sirna, f"site_{site_type}"]
            # (?=...) is a lookahead: it also counts occurrences that overlap each other
            seed_tx[f"{sirna}_{site_type}_count"] = utr_tx["utr_sequence"].str.count(f"(?={site})")
    return seed_tx


def sites_to_gene(seed_tx, sirna_table):
    """Aggregate transcript-level site counts to gene level ("any isoform" rule).

    Rule: per gene, each count is the MAX over its UTR transcripts. Isoforms largely share
    UTR sequence, so summing would count the same site several times; max gives the
    number of distinct sites in the richest isoform (a lower bound if isoforms carry
    different unique sites).

    Returns one row per gene (index 'gene_id_clean') with 'n_utr_transcripts', the
    '<siRNA>_<sitetype>_count' and '<siRNA>_<sitetype>_present' columns, and
    '<siRNA>_best_site' (most stringent site present, or 'none').
    """
    count_cols = [c for c in seed_tx.columns if c.endswith("_count")]
    grouped = seed_tx.groupby("gene_id_clean")

    seed_gene = grouped[count_cols].max()
    seed_gene.insert(0, "n_utr_transcripts", grouped.size())

    for sirna in sirna_table.index:
        present_cols = []
        for site_type in SITE_TYPES:
            col = f"{sirna}_{site_type}_present"
            seed_gene[col] = seed_gene[f"{sirna}_{site_type}_count"] > 0
            present_cols.append(col)
        # np.select picks the first true condition, so SITE_TYPES order = stringency
        seed_gene[f"{sirna}_best_site"] = np.select(
            [seed_gene[c] for c in present_cols], SITE_TYPES, default="none"
        )
    return seed_gene
