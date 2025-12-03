# BioStringFlow
BioStringFlow is a fully modular, beginner-friendly, and reproducible DNA-protein analysis pipeline built entirely in R, optimized for life-science students, researchers, and bioinformatics enthusiasts.
It is a modular R pipeline built on top of the Biostrings package for:
- Fetching DNA sequences from NCBI
- Fetching protein sequences from UniProt
- Basic sequence exploration (length, base composition, GC content)
- Motif search
- DNA → protein translation
- Reverse complement and ORF detection
- Codon usage analysis
- Pairwise alignment between translated DNA and UniProt protein
- Simple phylogenetic tree construction

## Features

- Fully modular: each step is in its own R script (`module1`–`module9`).
- Reproducible: same functions can be reused for any gene/species by changing accession IDs.
- Uses standard Bioconductor tools: `Biostrings`, `rentrez`, etc.
- Includes example outputs for the human mitochondrial **COX1** gene.

## Project Structure
modules/
 - module1_fetch_dna.R – fetches DNA from NCBI
 - module2_fetch_protein.R – fetches protein from UniProt
 - module3_gc_content.R – GC content and base composition
 - module4_motifs.R – basic motif search
 - module5_translate.R – translation of DNA and export of protein sequence
 - module6_orfs_reverse_comp.R – reverse complement + ORF detection
 - module7_codon_usage.R – codon usage + heatmap
 - module8_alignment.R – alignment between translated DNA and UniProt protein
 - module9_phylogeny.R – distance matrix + simple phylogenetic tree/
   
main_pipeline.R/

sample output/

Example PDFs and text outputs generated for COX1 (NC_012920.1, P00395)

## Requirements

- R (version ≥ 4.4 recommended)
- Packages:
  - BiocManager
  - Biostrings
  - rentrez
  - seqinr
  - pheatmap
  - ape / pwalign (depending on your alignment/tree modules)

Install dependencies in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("Biostrings", "rentrez", "seqinr", "pheatmap", "ape", "pwalign")
new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) BiocManager::install(new)
