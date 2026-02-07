"""
Codon Usage Analyzer

Analyzes codon usage frequency, GC content,
and RSCU values from CDS FASTA files.
"""

from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt


def read_cds_fasta(filename):
    cds_sequences = []

    for record in SeqIO.parse(filename, "fasta"):
        seq = str(record.seq).upper()

        # Basic quality checks
        if len(seq) % 3 != 0:
            continue
        if "N" in seq:
            continue

        cds_sequences.append(seq)

    return cds_sequences


def extract_codons(cds_sequences):
    codons = []

    for seq in cds_sequences:
        for i in range(0, len(seq), 3):
            codons.append(seq[i:i+3])

    return codons


def calculate_gc(cds_sequences):
    g = c = total = 0

    for seq in cds_sequences:
        g += seq.count("G")
        c += seq.count("C")
        total += len(seq)

    return (g + c) / total * 100


def calculate_rscu(codon_counts):
    genetic_code = {
        'TTT':'F','TTC':'F',
        'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I',
        'ATG':'M',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y',
        'CAT':'H','CAC':'H',
        'CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N',
        'AAA':'K','AAG':'K',
        'GAT':'D','GAC':'D',
        'GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C',
        'TGG':'W',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    }

    rscu_values = {}

    for aa in set(genetic_code.values()):
        synonymous_codons = [
            codon for codon, amino_acid in genetic_code.items()
            if amino_acid == aa
        ]

        total = sum(codon_counts.get(c, 0) for c in synonymous_codons)
        if total == 0:
            continue

        expected = total / len(synonymous_codons)

        for codon in synonymous_codons:
            rscu_values[codon] = codon_counts.get(codon, 0) / expected

    return pd.DataFrame.from_dict(
        rscu_values, orient="index", columns=["RSCU"]
    )
    
import os

def main():
    # File paths
    input_fasta = "data/ecoli_cds.fna"
    output_freq = "output/ecoli_codon_frequency.csv"
    output_rscu = "output/ecoli_rscu.csv"
    os.makedirs("output", exist_ok=True)
    # Read CDS
    cds_sequences = read_cds_fasta("input_fasta")
    print("Number of CDS loaded:", len(cds_sequences))

    # Extract codons
    codon_list = extract_codons(cds_sequences)
    codon_counts = Counter(codon_list)

    # Codon frequency table
    codon_df = pd.DataFrame.from_dict(
        codon_counts, orient="index", columns=["Count"]
    ).sort_values(by="Count", ascending=False)

    # Remove stop codons
    stop_codons = ["TAA", "TAG", "TGA"]
    codon_df = codon_df.drop(stop_codons, errors="ignore")

    total_codons = codon_df["Count"].sum()
    codon_df["Frequency_per_1000"] = (
        codon_df["Count"] / total_codons * 1000
    )

    codon_df.to_csv(output_freq)

    # GC content
    gc_percent = calculate_gc(cds_sequences)
    print("GC content (%):", round(gc_percent, 2))

    # RSCU
    rscu_df = calculate_rscu(codon_counts)
    rscu_df.to_csv(output_rscu)

    # Plot top 20 codons
    top_codons = codon_df.head(20)
    plt.figure(figsize=(14, 5))
    plt.bar(top_codons.index, top_codons["Count"])
    plt.xlabel("Codon")
    plt.ylabel("Count")
    plt.title("Top 20 Codons in E. coli")
    plt.tight_layout()
    plt.show()

    print("Codon usage analysis completed.")


if __name__ == "__main__":
    main()
