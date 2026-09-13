# Tobamovirus sequence analysis

Comparative analysis of tobamovirus genome sequences from wild hosts: multiple sequence alignment, genetic distances and SNP detection.

## What the pipeline does

`pryecti1.py` runs four steps:

1. **Merge the sequences.** Reads the individual `sequence (n).fasta` files with Biopython and writes them into a single multi-FASTA, `secuencias_unicas.fasta`.
2. **Read the alignment.** Loads the Clustal alignment (`.aln`) with `AlignIO` and renders it for inspection.
3. **Genetic distances.** Loads `matriz_distancias.csv`, normalises decimal commas to points, and draws an annotated heatmap of the distance matrix.
4. **SNP detection.** Turns the alignment into a position-by-sequence matrix, flags every column with more than one allele, and writes those positions to `SNPs_detectados.csv`.

## Files

| File | What it is |
|---|---|
| `sequence (n).fasta` | Individual genome sequences downloaded from GenBank |
| `secuencias_unicas.fasta` | All of the above merged into one multi-FASTA |
| `_out...aln` | Clustal multiple sequence alignment |
| `matriz_distancias.csv` | Pairwise genetic distance matrix |
| `SNPs_detectados.csv` | Variable positions found in the alignment |
| `heatmap.png` | Genetic distance heatmap |
| SNP distribution PNGs | SNP distribution and frequency along the genome |

## Requirements

Python 3.9+, with Biopython, pandas, matplotlib and seaborn. Exact versions are pinned in `requirements.txt`:

```
pip install -r requirements.txt
```

## Before running

The script expects the FASTA files, the `.aln` alignment and `matriz_distancias.csv` to sit in the same folder as the script, under exactly those names. The alignment filename is hard-coded near the top.

## Data

Sequences are public genome records; no restricted data is used.
