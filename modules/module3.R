# ===========================================
# MODULE 03: GC Content Analysis 
# ===========================================

library(Biostrings)

analyze_gc <- function(dna_set, out_pdf = NULL) {
  
  # Validate input
  if (!inherits(dna_set, "DNAStringSet")) {
    stop("Input must be a DNAStringSet object.")
  }
  
  seq <- as.character(dna_set[[1]])  # extract DNA sequence as string
  
  # Base frequencies using Biostrings
  freq <- alphabetFrequency(dna_set, baseOnly = TRUE)
  
  # Extract counts
  A <- freq[1, "A"]
  T <- freq[1, "T"]
  G <- freq[1, "G"]
  C <- freq[1, "C"]
  
  seq_length <- A + T + G + C
  
  # GC %
  gc_percent <- ((G + C) / seq_length) * 100
  
  # Prepare result list
  result <- list(
    length       = seq_length,
    A            = A,
    T            = T,
    G            = G,
    C            = C,
    GC_percent   = gc_percent,
    AT_percent   = ((A + T) / seq_length) * 100,
    base_freq_df = as.data.frame(freq)
  )
  
  # Optionally create GC barplot PDF
  if (!is.null(out_pdf)) {
    pdf(out_pdf)
    barplot(
      c(A, T, G, C),
      names.arg = c("A", "T", "G", "C"),
      col = c("green", "red", "blue", "yellow"),
      main = "Nucleotide Composition",
      ylab = "Count"
    )
    dev.off()
    message("GC PDF saved to: ", out_pdf)
  }
  
  return(result)
}
