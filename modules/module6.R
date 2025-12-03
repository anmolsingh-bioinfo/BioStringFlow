## =========================
## MODULE 6: RC + ORF DETECTION
## =========================

# Helper: ensure we have a plain DNAString
.get_dna_string <- function(dna_input) {
  if (methods::is(dna_input, "DNAStringSet")) {
    if (length(dna_input) < 1L) stop("DNAStringSet is empty.")
    return(dna_input[[1]])
  } else if (methods::is(dna_input, "DNAString")) {
    return(dna_input)
  } else {
    # try to coerce from character
    dna_char <- as.character(dna_input)[1]
    return(Biostrings::DNAString(dna_char))
  }
}

# ---------- 6A. Reverse complement ----------

get_reverse_complement <- function(dna_input) {
  dna <- .get_dna_string(dna_input)
  rc  <- Biostrings::reverseComplement(dna)
  return(rc)
}

# ---------- 6B. ORF finder (6 frames) ----------

find_orfs <- function(dna_input,
                      min_aa   = 30,
                      start_codon = "ATG",
                      stop_codons = c("TAA", "TAG", "TGA"),
                      make_plot   = TRUE) {
  # Get a single DNAString
  dna <- .get_dna_string(dna_input)
  
  # Convert to character for manual scanning (no width())
  seq_plus  <- as.character(dna)
  seq_minus <- as.character(Biostrings::reverseComplement(dna))
  
  seq_len_plus  <- nchar(seq_plus)
  seq_len_minus <- nchar(seq_minus)
  
  orf_list <- list()
  idx <- 1
  
  scan_one_strand <- function(seq_char, strand_label) {
    seq_len <- nchar(seq_char)
    res <- list()
    idx_local <- 1
    
    for (frame in 0:2) {
      in_orf <- FALSE
      orf_start <- NA
      
      pos <- frame + 1
      while (pos + 2 <= seq_len) {
        codon <- substr(seq_char, pos, pos + 2)
        
        if (!in_orf && codon == start_codon) {
          in_orf <- TRUE
          orf_start <- pos
        } else if (in_orf && codon %in% stop_codons) {
          orf_end <- pos + 2
          aa_len  <- (orf_end - orf_start + 1) / 3
          
          if (aa_len >= min_aa) {
            res[[idx_local]] <- data.frame(
              strand   = strand_label,
              frame    = frame,          # frame 0,1,2 for that strand
              start_nt = orf_start,
              end_nt   = orf_end,
              length_aa = aa_len,
              stringsAsFactors = FALSE
            )
            idx_local <- idx_local + 1
          }
          in_orf <- FALSE
          orf_start <- NA
        }
        
        pos <- pos + 3
      }
    }
    
    if (length(res) == 0L) {
      return(NULL)
    } else {
      return(do.call(rbind, res))
    }
  }
  
  # Scan + strand
  plus_orfs  <- scan_one_strand(seq_plus,  strand_label = "+")
  # Scan - strand (on reverse complement)
  minus_orfs <- scan_one_strand(seq_minus, strand_label = "-")
  
  if (!is.null(plus_orfs)) {
    orf_list[[idx]] <- plus_orfs
    idx <- idx + 1
  }
  if (!is.null(minus_orfs)) {
    orf_list[[idx]] <- minus_orfs
  }
  
  if (length(orf_list) == 0L) {
    message("No ORFs found with minimum length >= ", min_aa, " aa.")
    return(invisible(NULL))
  }
  
  orf_df <- do.call(rbind, orf_list)
  
  #  simple length distribution plot
  if (make_plot) {
    graphics::hist(
      orf_df$length_aa,
      breaks = 20,
      main   = "ORF Length Distribution",
      xlab   = "ORF length (aa)"
    )
  }
  
  return(orf_df)
}
