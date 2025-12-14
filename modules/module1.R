# ===========================================
# MODULE 01: Fetch DNA from NCBI
# ===========================================

library(rentrez)
library(Biostrings)

fetch_dna <- function(ncbi_id, out_fasta = NULL) {
  
  if (is.null(out_fasta)) {
    out_fasta <- paste0(ncbi_id, "_dna.fasta")
  }
  
  message(">>> Fetching DNA for: ", ncbi_id)
  
  # Try POST + FETCH method
  fasta_txt <- tryCatch({
    wh <- rentrez::entrez_post(db="nucleotide", id=ncbi_id)
    rentrez::entrez_fetch(db="nucleotide", web_history=wh,
                          rettype="fasta", retmode="text")
  }, error = function(e) {
    message("NCBI Fetch failed. Attempting local fallback...")
    return(NA_character_)
  })
  
  # If NCBI failed, use local file IF exists
  if (is.na(fasta_txt)) {
    if (!file.exists(out_fasta)) {
      stop("NCBI fetch failed and no local FASTA found: ", out_fasta)
    }
    message("Using local FASTA instead: ", out_fasta)
    dna_set <- Biostrings::readDNAStringSet(out_fasta)
    seq_len <- nchar(as.character(dna_set[[1]]))
    
    return(list(
      dna = dna_set,
      fasta_path = out_fasta,
      length = seq_len,
      source = "local"
    ))
  }
  
  # Save FASTA
  writeLines(fasta_txt, out_fasta)
  message("Saved FASTA to: ", out_fasta)
  
  # Load sequence
  dna_set <- Biostrings::readDNAStringSet(out_fasta)
  seq_len <- nchar(as.character(dna_set[[1]]))
  
  return(list(
    dna = dna_set,
    fasta_path = out_fasta,
    length = seq_len,
    raw_fasta = fasta_txt,
    source = "NCBI"
  ))
}

