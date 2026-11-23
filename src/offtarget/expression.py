
import numpy as np
from scipy.stats import ttest_rel

# Import process:

# import sys
# sys.path.append("../src")

#from expression import (
#    extract_gene_expression,
#    calculate_fold_change,
#    paired_t_test
#)

def extract_gene_expression(df, gene_id, id_col="ENSEMBL ID"):
    """Extract expression values for a specific gene."""
    gene = df.loc[df[id_col] == gene_id]

    if gene.empty:
        raise ValueError(f"Gene ID {gene_id} not found.")

    if len(gene) > 1:
        raise ValueError(f"Multiple rows found for gene ID {gene_id}.")

    return gene.iloc[0].drop(labels=id_col)


def calculate_fold_change(control, treated):
    """Calculate fold change, log2FC and percentage knockdown."""
    control_mean = control.mean()
    treated_mean = treated.mean()

    fold_change = treated_mean / control_mean
    log2fc = np.log2(fold_change)
    knockdown = (1 - fold_change) * 100

    return {
        "control_mean": control_mean,
        "treated_mean": treated_mean,
        "fold_change": fold_change,
        "log2fc": log2fc,
        "knockdown_percent": knockdown
    }


def paired_t_test(control, treated):
    """Perform a paired t-test."""
    return ttest_rel(control, treated)


# Extract and clean human 3'UTR from targetscan UTR seq file

import csv
import gzip
from pathlib import Path

input_path = Path("../data/raw/TargetScan/UTR_Sequences.txt")
output_path = Path("../data/processed/human_3UTRs.tsv.gz")

output_path.parent.mkdir(parents=True, exist_ok=True)

total_records = 0
human_records = 0

with open(input_path, "r", newline="") as infile, \
     gzip.open(output_path, "wt", newline="") as outfile:

    reader = csv.DictReader(infile, delimiter="\t")

    fieldnames = [
        "transcript_id",
        "gene_id",
        "gene_symbol",
        "utr_sequence",
        "utr_length"
    ]

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    for row in reader:
        total_records += 1

        if row["Species ID"].strip() != "9606":
            continue

        sequence = row["UTR sequence"].replace("-", "").strip().upper()
        sequence = sequence.replace("U", "T")

        writer.writerow({
            "transcript_id": row["Refseq ID"].strip(),
            "gene_id": row["Gene ID"].strip(),
            "gene_symbol": row["Gene Symbol"].strip(),
            "utr_sequence": sequence,
            "utr_length": len(sequence)
        })

        human_records += 1

print("Total records processed:", total_records)
print("Human records retained:", human_records)
print("Output file:", output_path)

gene_expression_utr = gene_expression_utr.dropna(
    subset=["n_utr_transcripts"]
).copy()
