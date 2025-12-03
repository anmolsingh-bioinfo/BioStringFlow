# ===========================================
# MODULE 05: DNA Translation (no width())
# ===========================================

library(Biostrings)

translate_dna <- function(dna_set, 
                          trim = TRUE, 
                          if_fuzzy = "solve",
                          out_fasta = NULL) {
  
  # Validate input
  if (!inherits(dna_set, "DNAStringSet")) {
    stop("Input must be a DNAStringSet object.")
  }
  
  # Extract DNA sequence
  dna <- dna_set[[1]]
  
  # Convert to simple character string
  dna_char <- as.character(dna)
  
  # TRIM TO MULTIPLE OF 3
  if (trim) {
    extra <- nchar(dna_char) %% 3
    if (extra != 0) {
      dna_char <- substr(dna_char, 1, nchar(dna_char) - extra)
      message("Trimmed sequence by ", extra, " base(s) for correct translation.")
    }
  }
  
  # Convert back to DNAString
  dna_trimmed <- DNAString(dna_char)
  
  # Translation using Biostrings
  aa_seq <- translate(
    dna_trimmed,
    if.fuzzy.codon = if_fuzzy
  )
  
  # Convert to AAStringSet
  aa_set <- AAStringSet(aa_seq)
  
  # Optionally save to FASTA
  if (!is.null(out_fasta)) {
    writeXStringSet(aa_set, out_fasta)
    message("Translated protein saved to: ", out_fasta)
  }
  
  return(list(
    trimmed_dna     = DNAStringSet(dna_trimmed),
    protein         = aa_set,
    protein_string  = as.character(aa_set[[1]]),
    length          = nchar(as.character(aa_set[[1]])),
    fuzzy_mode      = if_fuzzy
  ))
}
