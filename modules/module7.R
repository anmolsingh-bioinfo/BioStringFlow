# ===========================================
# MODULE 07: Codon Usage Analysis
# ===========================================

library(Biostrings)
library(pheatmap)

analyze_codon_usage <- function(dna_set,
                                out_pdf = "Codon_Usage_Heatmap.pdf") {
  
  # 1. Extract the first (and usually only) DNA sequence
  dna_seq <- dna_set[[1]]
  seq_char <- toupper(as.character(dna_seq))
  
  # 2. Trim sequence so length is divisible by 3
  len <- nchar(seq_char)
  trim_len <- len - (len %% 3)
  seq_char <- substr(seq_char, 1, trim_len)
  
  # 3. Break into codons
  codons <- substring(seq_char, seq(1, trim_len, 3), seq(3, trim_len, 3))
  
  # 4. Define all 64 codons
  bases <- c("A", "T", "G", "C")
  codon_list <- as.vector(outer(bases, outer(bases, bases, paste0), paste0))
  
  # 5. Count codons
  codon_table <- table(factor(codons, levels = codon_list))
  codon_df <- as.data.frame(codon_table)
  colnames(codon_df) <- c("Codon", "Count")
  
  # 6. Prepare matrix for heatmap
  # Convert codons into a 16x4 matrix (first two bases group rows)
  heatmap_matrix <- matrix(codon_table, nrow = 16, byrow = TRUE)
  rownames(heatmap_matrix) <- unique(substring(codon_list, 1, 2))
  colnames(heatmap_matrix) <- bases
  
  # 7. Save heatmap
  pdf(out_pdf)
  pheatmap(heatmap_matrix,
           main = "Codon Usage Heatmap",
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  dev.off()
  
  message("Codon usage heatmap saved to: ", out_pdf)
  
  return(list(
    codon_counts = codon_df,
    heatmap_matrix = heatmap_matrix,
    trimmed_length = trim_len,
    total_codons = length(codons)
  ))
}
