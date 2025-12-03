# ===========================================
# MODULE 02: Fetch Protein Sequence from UniProt (NO WIDTH())
# ===========================================

library(Biostrings)

fetch_protein <- function(uniprot_id, out_fasta = NULL) {
  
  if (is.null(out_fasta)) {
    out_fasta <- paste0(uniprot_id, "_protein.fasta")
  }
  
  message(">>> Fetching protein for UniProt ID: ", uniprot_id)
  
  # UniProt URL
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  
  ok <- tryCatch({
    download.file(url, destfile = out_fasta, quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  
  if (!ok) {
    if (!file.exists(out_fasta)) {
      stop("UniProt download failed and no local file found: ", out_fasta)
    }
    message("Using local FASTA file.")
  }
  
  # Load protein FASTA
  aa_set <- Biostrings::readAAStringSet(out_fasta)
  
  # Instead of width(), use nchar()
  seq_len <- nchar(as.character(aa_set[[1]]))
  
  return(list(
    protein = aa_set,
    fasta_path = out_fasta,
    length = seq_len,
    source = ifelse(ok, "UniProt", "local")
  ))
}
