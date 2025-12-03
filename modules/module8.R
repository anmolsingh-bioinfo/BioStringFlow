# ==========================================
# MODULE 08: DNA→Protein Alignment
# ==========================================

if (!requireNamespace("Biostrings", quietly = TRUE))
  install.packages("Biostrings")

library(Biostrings)

# -------------------------------------------------------------------
# FUNCTION: align_sequences
# Inputs:
#   translated_fasta: path to translated DNA protein FASTA
#   uniprot_fasta:    path to UniProt protein FASTA
#   out_file:         output text file for alignment summary
# -------------------------------------------------------------------

align_sequences <- function(translated_fasta,
                            uniprot_fasta,
                            out_file = "alignment_output.txt") {
  
  message(">>> MODULE 8: Loading sequences...")
  
  # Load sequences
  dna_prot <- readAAStringSet(translated_fasta)
  uni_prot <- readAAStringSet(uniprot_fasta)
  
  # Extract the first (assuming single-sequence files)
  dna_seq  <- dna_prot[[1]]
  uni_seq  <- uni_prot[[1]]
  
  # Alignment
  message(">>> Performing global alignment...")
  aln <- pairwiseAlignment(dna_seq, uni_seq, type = "global")
  
  # % Identity
  identity <- pid(aln)
  
  # Save report
  lines <- c(
    "===== SEQUENCE ALIGNMENT REPORT =====",
    paste("Translated DNA (AA) file:", translated_fasta),
    paste("UniProt protein file:", uniprot_fasta),
    "",
    paste("Alignment Score:", score(aln)),
    paste("Percent Identity:", round(identity, 2), "%"),
    "",
    "===== ALIGNMENT =====",
    as.character(aln)
  )
  
  writeLines(lines, con = out_file)
  message(paste("Alignment saved to:", out_file))
  
  return(list(
    alignment = aln,
    identity  = identity,
    score     = score(aln),
    out_file  = out_file
  ))
}
