#Before running this script make sure scripts for Module 1 to 9 are saved in same directory in which this script is present...



#Module1 --- Fetches DNA Seq from NCBI ---- Output saved as txt file -----
source("module1.R")
dna <- fetch_dna("NC_012920.1")     #Change only gene accession id
dna
dna$length
dna$dna_result
#Module2 --- Fetches Protein Seq from Uniport ---- Output saved as txt file -----
source("module2.R")
prot <- fetch_protein("P00395")     #change onlt protein id 
prot
prot$length
prot$fasta_path

#Module3 --- Analyze GC Content ---- Output saved as PDF file ----- (make no changes)
source("module3.R")
gc <- analyze_gc(dna$dna, out_pdf = "GC_Content.pdf")
gc
gc$GC_percent
gc$base_freq_df

#Module4 --- Motif Finder ---- Output saved as txt file -----   (nake no changes)
source("module4.R")
motif_output <- find_motifs(
  dna_set = dna$dna,
  motifs = c("ATG", "TATA", "AATAAA", "CG"),
  out_txt = "Motif_Report.txt"
)
motif_output$summary
motif_output$results$ATG$positions

#Module5 --- DNA Seq Translation ---- Output saved as txt file -----  (no changes)
source("module5.R")
translation <- translate_dna(
  dna_set   = dna$dna,
  trim      = TRUE,
  if_fuzzy  = "solve",
  out_fasta = "translated_sequence.fasta"
)
translation$protein
translation$length
translation$protein_string

#Module6 --- Reverse Complement and ORF Finding -----  (no changes)
source("module6.R")
# 1. Reverse complement
rc_seq <- get_reverse_complement(dna$dna)
rc_seq   # prints reverse complement DNAString

# 2. Find ORFs (on both strands)
orfs <- find_orfs(dna$dna,
                  min_aa   = 30,
                  make_plot = TRUE)
head(orfs)


#Module7 --- Codon Usage Analysis -----   (no changes)
source("module7.R")

codon_res <- analyze_codon_usage(dna$dna,
                                 out_pdf = "Codon_Usage.pdf")
codon_res$codon_counts
codon_res$heatmap_matrix

#Module8 --- Alignment ----   (no changes required in case output file name from module 5 and module 2)
source("module8.R")

result_alignment <- align_sequences(
  translated_fasta = "translated_sequence.fasta",
  uniprot_fasta    = "P00395_protein.fasta",
  out_file         = "Module8_Alignment_Report.txt"
)

result_alignment$identity
result_alignment$score
result_alignment$alignment


#Module9 --- Evolution Analysis ---  
source("module9.R")
tree <- build_phylogeny(
  ids    = c("NC_012920.1", "NC_005089.1", "NC_002008.4"),    # change the gene ids as per requirement
  labels = c("Human", "Chimpanzee", "Mouse"),
  seq_type = "DNA",
  out_dir = "COX1_Tree"
)


