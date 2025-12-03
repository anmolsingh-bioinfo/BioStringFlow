# ===========================================
# MODULE 09: Evolutionary Tree Module
# Build phylogeny from multiple sequences
# ===========================================

if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")

library(rentrez)
library(Biostrings)

# -------------------------------------------
# FUNCTION: build_phylogeny()
# -------------------------------------------
build_phylogeny <- function(ids,
                            labels,
                            seq_type = "DNA",   # "DNA" or "Protein"
                            out_dir = "Module9_Tree") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir)
  
  message(">>> Starting Phylogeny Construction...")
  
  # Decide database
  db_name <- if (seq_type == "DNA") "nucleotide" else "protein"
  
  # EMPTY container
  if (seq_type == "DNA") {
    seqs <- DNAStringSet()
  } else {
    seqs <- AAStringSet()
  }
  
  # -------------------------------
  # Fetch each sequence
  # -------------------------------
  for (i in seq_along(ids)) {
    id <- ids[i]
    sp <- labels[i]
    
    message(paste("Fetching:", sp, "| Accession:", id))
    
    tryCatch({
      # Fetch FASTA
      fasta_txt <- rentrez::entrez_fetch(
        db = db_name,
        id = id,
        rettype = "fasta",
        retmode = "text"
      )
      
      if (!grepl("^>", fasta_txt)) stop("Invalid FASTA returned.")
      
      # Write to temp
      tf <- tempfile(fileext = ".fasta")
      writeLines(fasta_txt, tf)
      
      # Read
      s <- if (seq_type == "DNA") readDNAStringSet(tf) else readAAStringSet(tf)
      
      # Rename
      names(s) <- sp
      
      # Add
      seqs <- c(seqs, s)
      
      Sys.sleep(0.4) # prevent NCBI overload
      
    }, error = function(e) {
      message("  FAILED for ", sp)
      message("  Error: ", conditionMessage(e))
    })
  }
  
  if (length(seqs) < 2) stop("Not enough valid sequences to build tree.")
  
  # Save combined FASTA
  fasta_out <- file.path(out_dir, paste0("Combined_", seq_type, ".fasta"))
  writeXStringSet(seqs, fasta_out)
  
  # -------------------------------
  # Compute distance matrix
  # -------------------------------
  message("Computing distance matrix...")
  dist_mat <- stringDist(seqs, method = "levenshtein")
  
  # -------------------------------
  # Build tree
  # -------------------------------
  message("Building tree using hierarchical clustering...")
  tree <- hclust(dist_mat, method = "average")
  
  # -------------------------------
  # Save Tree PDF
  # -------------------------------
  pdf_path <- file.path(out_dir, paste0("Phylogeny_", seq_type, ".pdf"))
  pdf(pdf_path)
  plot(tree,
       main = paste("Evolutionary Tree -", seq_type),
       xlab = "Species",
       sub = "Distance: Levenshtein, Clustering: UPGMA")
  dev.off()
  
  message(">>> Phylogeny completed. Output saved in: ", out_dir)
  
  return(list(
    sequences = seqs,
    dist_matrix = dist_mat,
    tree = tree,
    fasta = fasta_out,
    pdf = pdf_path
  ))
}
