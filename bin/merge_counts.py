#!/usr/bin/env python3
"""Merge per-sample HTSeq count files (<sample>.counts.txt) into one gene x sample matrix.

Usage: merge_counts.py SAMPLE1.counts.txt SAMPLE2.counts.txt ... > count_matrix.tsv
HTSeq summary rows (__no_feature, __ambiguous, ...) are dropped.
"""
import sys
from pathlib import Path


def read_counts(path):
    counts = {}
    with open(path) as fh:
        for line in fh:
            gene, value = line.rstrip("\n").split("\t")
            if not gene.startswith("__"):
                counts[gene] = int(value)
    return counts


def main(paths):
    if not paths:
        sys.exit("usage: merge_counts.py SAMPLE.counts.txt [...]")
    paths = sorted(paths, key=lambda p: Path(p).name)
    samples = [Path(p).name.removesuffix(".counts.txt") for p in paths]
    tables = [read_counts(p) for p in paths]

    genes = list(tables[0])
    for name, table in zip(samples, tables):
        if list(table) != genes:
            sys.exit(f"Gene list in {name} differs from {samples[0]}; were all samples counted with the same GTF?")

    out = sys.stdout
    out.write("gene_id\t" + "\t".join(samples) + "\n")
    for gene in genes:
        out.write(gene + "\t" + "\t".join(str(t[gene]) for t in tables) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
