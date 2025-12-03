# ===========================================
# MODULE 04: Motif Search (no width())
# ===========================================

library(Biostrings)

find_motifs <- function(dna_set, motifs = c("ATG", "TATA", "AATAAA", "CG"),
                        out_txt = NULL) {
  
  # Validate input
  if (!inherits(dna_set, "DNAStringSet")) {
    stop("Input must be a DNAStringSet object.")
  }
  
  seq <- dna_set[[1]]  # DNAString object
  
  results_list <- list()
  
  # For each motif
  for (m in motifs) {
    matches <- matchPattern(m, seq)
    
    results_list[[m]] <- list(
      count = length(matches),
      positions = start(matches)
    )
  }
  
  # Final summary
  summary_df <- data.frame(
    Motif = motifs,
    Count = sapply(results_list, function(x) x$count)
  )
  
  # Optionally save as txt report
  if (!is.null(out_txt)) {
    cat("=== Motif Search Report ===\n\n", file = out_txt)
    
    for (m in motifs) {
      cat("Motif:", m, "\n", file = out_txt, append = TRUE)
      cat("Count:", results_list[[m]]$count, "\n", file = out_txt, append = TRUE)
      cat("Positions:", 
          paste(results_list[[m]]$positions, collapse = ", "),
          "\n\n", file = out_txt, append = TRUE)
    }
    
    message("Motif report saved to: ", out_txt)
  }
  
  return(list(
    motifs = motifs,
    results = results_list,
    summary = summary_df
  ))
}
