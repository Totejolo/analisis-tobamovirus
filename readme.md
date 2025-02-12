# Phylogenetic Analysis of Tobamovirus

This repository contains the scripts and datasets used for the phylogenetic analysis and SNP identification of less-studied *Tobamovirus* strains.

## 📂 Contents

- **`proyecto1.py`** → Main Python script for sequence processing and analysis.
- **`secuencias_unicas.fasta`** → Merged FASTA file containing all sequences.
- **`sequence (5-30).fasta`** → Individual FASTA sequence files (5 to 30).
- **`_out.250207034227455MdWkmzDeLfXS1HuBU5RgTlsfnormal.aln`** → Clustal alignment file for the sequences.
- **`SNPs_detectados.csv`** → CSV file containing the detected SNPs.
- **`matriz_distancias.csv`** → Genetic distance matrix generated from the alignment.
- **`requirements.txt`** → List of required Python libraries.

## 🚀 How to Use

1. **Install dependencies** (if not already installed):
    ```bash
    pip install -r requirements.txt
    ```

2. **Run the analysis**:
    Execute the main script using Python:
    ```bash
    python proyecto1.py
    ```

3. **Results**:
   - The results will include the combined sequences file, a CSV of SNPs detected, and a heatmap image of genetic distances.
   - You can explore the generated files such as `SNPs_detectados.csv` and `heatmap.png`.

## 🔗 Link to GitHub Repository
Access the full repository and code at [https://github.com/Totejolo/analisis-tobamovirus](https://github.com/Totejolo/analisis-tobamovirus).
