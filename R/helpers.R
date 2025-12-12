#' Internal helpers
#' @keywords internal
#' @noRd

# helpers.R
# 
# This script contains small utility functions used internally across the package.
# These helpers provide generic, reusable functionality (e.g., data reshaping,
# validation checks, and formatting) to keep main plotting functions concise.

# -------------------------------------------------------------------------
# melt(): reshape a matrix or data.frame into long format
# -------------------------------------------------------------------------
melt <- function(mat, value_name = "value") {
  # Check if input is a matrix
  if (!is.matrix(mat)) {
    message("Input is not a matrix - coercing to matrix before melting.")
    mat <- as.matrix(mat)
  }

  # Check if input is empty
  if (nrow(mat) == 0) {
    message("Input has 0 rows.")
    return(data.frame())
  }
  
  # Convert to long format
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c("Var1", "Var2", value_name)
  
  out
}


# -------------------------------------------------------------------------
# cor_p(): Compute correlation and p-value matrix
# -------------------------------------------------------------------------
cor_p <- function(data, use = "complete.obs", method = "pearson") {
  # Keep numeric variables
  num <- data[, sapply(data, is.numeric), drop = FALSE]
  # Variable count
  p   <- ncol(num)
  if (p < 2L) stop("Need at least two numeric columns.")

  # Correlation
  R <- cor(num, use = use, method = method)

  # Sample sizes per pair
  if (use == "complete.obs") {
    # All pairs have same sample size
    n_all <- sum(complete.cases(num))
    N <- matrix(n_all, nrow = p, ncol = p,
                dimnames = list(colnames(num), colnames(num)))
  } else {
    # Calculate pairwise sample size
    M <- !is.na(num)
    N <- t(M) %*% M
    N <- as.matrix(N)
    dimnames(N) <- list(colnames(num), colnames(num))
  }

  # Ensures all correlation (or similar) values stay safely inside (−1, 1)
  # by replacing tiny rounding overflows/underflows with stable limits.
  # e.g. due to floating-point rounding, it might contain slightly 
  # out-of-bounds numbers like 1.0000000002 or -1.0000000002
  R_clipped <- pmax(pmin(R,  1 - 1e-15), -1 + 1e-15)
  # Subtracts 2 degrees of freedom for each correlation pair because a 
  # correlation effectively fits a 2-parameter regression (intercept + slope)
  df_mat    <- N - 2
  T         <- R_clipped * sqrt((df_mat) / pmax(1e-15, 1 - R_clipped^2))
  P         <- 2 * pt(abs(T), df = pmax(1, df_mat), lower.tail = FALSE)
  # Diagonal has r=1, set p=0
  diag(P) <- 0
  list(r = R, p = P)
}


# -------------------------------------------------------------------------
# parse_clustal(): Parse Clustal alignment file
# -------------------------------------------------------------------------
parse_clustal <- function(file) {
  if (!file.exists(file)) {
    stop(sprintf("File not found: %s", file))
  }
  
  lines <- readLines(file, warn = FALSE)
  
  # Skip header and empty lines
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) == 0) {
    stop("Empty alignment file")
  }
  
  # Skip first line if it's Clustal header
  if (grepl("^CLUSTAL", lines[1], ignore.case = TRUE)) {
    lines <- lines[-1]
  }
  
  # Filter out conservation lines (*, :, .)
  lines <- lines[!grepl("^\\s*[*:. ]+$", lines)]
  
  # Parse sequences
  seq_list <- list()
  for (line in lines) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) >= 2) {
      seq_name <- parts[1]
      seq_chunk <- parts[2]
      
      if (!seq_name %in% names(seq_list)) {
        seq_list[[seq_name]] <- ""
      }
      seq_list[[seq_name]] <- paste0(seq_list[[seq_name]], seq_chunk)
    }
  }
  
  if (length(seq_list) == 0) {
    stop("No sequences found in alignment file")
  }
  
  # Replace gaps with x
  seq_list <- lapply(seq_list, function(s) gsub("-", "x", s, fixed = TRUE))
  
  # Convert to data frame
  data.frame(
    name = names(seq_list),
    sequence = unlist(seq_list, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------------------
# needleman_wunsch(): Global sequence alignment
# -------------------------------------------------------------------------
#' Needleman-Wunsch global sequence alignment
#'
#' @param seq1 Character string of first sequence
#' @param seq2 Character string of second sequence  
#' @param match Numeric score for matching characters (default 1)
#' @param mismatch Numeric score for mismatching characters (default -1)
#' @param gap Numeric penalty for gaps (default -3)
#'
#' @return List with aligned1, aligned2 (character strings with gaps as "-"),
#'   and score (numeric alignment score)
#'
#' @keywords internal
#' @noRd
needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -3) {
  # Input validation
  if (!is.character(seq1) || length(seq1) != 1 || nchar(seq1) == 0) {
    stop("'seq1' must be a single non-empty character string")
  }
  if (!is.character(seq2) || length(seq2) != 1 || nchar(seq2) == 0) {
    stop("'seq2' must be a single non-empty character string")
  }
  if (!is.numeric(match) || length(match) != 1) {
    stop("'match' must be a single numeric value")
  }
  if (!is.numeric(mismatch) || length(mismatch) != 1) {
    stop("'mismatch' must be a single numeric value")
  }
  if (!is.numeric(gap) || length(gap) != 1) {
    stop("'gap' must be a single numeric value")
  }
  
  # Split sequences
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  n <- length(s1)
  m <- length(s2)
  
  # Initialize matrices
  score <- matrix(0, nrow = n + 1, ncol = m + 1)
  trace <- matrix("", nrow = n + 1, ncol = m + 1)
  
  # Initialize first column and row
  for (i in 2:(n + 1)) {
    score[i, 1] <- score[i - 1, 1] + gap
    trace[i, 1] <- "U"
  }
  for (j in 2:(m + 1)) {
    score[1, j] <- score[1, j - 1] + gap
    trace[1, j] <- "L"
  }
  
  # Fill score and trace matrices
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      diag_score <- score[i - 1, j - 1] + ifelse(s1[i - 1] == s2[j - 1], match, mismatch)
      up_score <- score[i - 1, j] + gap
      left_score <- score[i, j - 1] + gap
      
      best <- max(diag_score, up_score, left_score)
      score[i, j] <- best
      
      trace[i, j] <- if (best == diag_score) {
        "D"
      } else if (best == up_score) {
        "U"
      } else {
        "L"
      }
    }
  }
  
  # Traceback to construct alignment
  i <- n + 1
  j <- m + 1
  align1 <- character(0)
  align2 <- character(0)
  
  while (i > 1 || j > 1) {
    step <- trace[i, j]
    if (step == "D") {
      align1 <- c(s1[i - 1], align1)
      align2 <- c(s2[j - 1], align2)
      i <- i - 1
      j <- j - 1
    } else if (step == "U") {
      align1 <- c(s1[i - 1], align1)
      align2 <- c("-", align2)
      i <- i - 1
    } else {
      align1 <- c("-", align1)
      align2 <- c(s2[j - 1], align2)
      j <- j - 1
    }
  }
  
  list(
    aligned1 = paste(align1, collapse = ""),
    aligned2 = paste(align2, collapse = ""),
    score = score[n + 1, m + 1]
  )
}


# -------------------------------------------------------------------------
# mut_to_seq(): Apply mutations to protein sequence
# -------------------------------------------------------------------------
#' Apply mutations to protein sequence
#'
#' @param ref Character string of the reference protein sequence
#' @param mut Character vector describing mutations. Each element can be a
#'   single mutation (e.g. "A20S") or multiple mutations separated by `sep`
#'   (e.g. "A20S / N100Q / A20- / -20P / -121HHHHHH")
#' @param sep Character string separator used to split multiple mutations
#'   within elements of `mut` (default " / ")
#'
#' @return Character string of the mutated protein sequence
#'
#' @details
#' Supported mutation formats (1-based positions relative to reference):
#'
#' - Substitution: \code{"A20S"} replaces residue at position 20 (must be A
#'   in reference) with S
#'
#' - Deletion: \code{"A20-"} deletes residue at position 20 (must be A in
#'   reference)
#'
#' - Insertion: \code{"-50P"} or \code{"-50HHHHHH"} inserts residue(s) before
#'   position 50 in the reference. Position 1 inserts before the first residue.
#'   Position \code{length(ref) + 1} appends at the end
#'
#' All mutation positions are interpreted relative to the original reference
#' sequence, so mutation order does not matter. Conflicting mutations at the
#' same position will raise an error.
#'
#' @keywords internal
#' @noRd
mut_to_seq <- function(ref, mut, sep = " / ") {
  # Input validation
  if (!is.character(ref) || length(ref) != 1 || nchar(ref) == 0) {
    stop("'ref' must be a single non-empty character string")
  }

  if (missing(mut) || is.null(mut)) {
    return(ref)
  }

  if (!is.character(mut)) {
    stop("'mut' must be a character vector")
  }

  if (!is.character(sep) || length(sep) != 1 || nchar(sep) == 0) {
    stop("'sep' must be a single non-empty character string")
  }

  # Expand and clean mutation codes
  mut_codes <- unlist(strsplit(mut, sep, fixed = TRUE), use.names = FALSE)
  mut_codes <- trimws(mut_codes)
  mut_codes <- mut_codes[nchar(mut_codes) > 0]

  if (length(mut_codes) == 0) {
    return(ref)
  }

  # Prepare reference
  ref_vec <- strsplit(ref, "")[[1]]
  ref_len <- length(ref_vec)

  # Define mutation patterns
  mut_pattern_sub <- "^([A-Za-z\\*])([0-9]+)([A-Za-z\\*])$"
  mut_pattern_del <- "^([A-Za-z\\*])([0-9]+)-$"
  mut_pattern_ins <- "^-([0-9]+)([A-Za-z\\*]+)$"

  # Storage for mutations by position
  substitutions <- list()  # pos -> new_aa
  deletions <- integer()    # vector of positions to delete
  insertions <- list()      # pos -> vector of aa to insert before pos

  # Parse all mutations
  for (code in mut_codes) {
    # Substitution: A20S
    if (grepl(mut_pattern_sub, code)) {
      m <- regexec(mut_pattern_sub, code)
      parts <- regmatches(code, m)[[1]]
      
      if (length(parts) != 4) {
        stop(sprintf("Invalid mutation format: '%s'", code))
      }
      
      from_aa <- toupper(parts[2])
      pos <- as.integer(parts[3])
      to_aa <- toupper(parts[4])
      
      if (is.na(pos) || pos < 1 || pos > ref_len) {
        stop(sprintf("Position out of range in '%s' (must be 1-%d)", code, ref_len))
      }
      
      if (toupper(ref_vec[pos]) != from_aa) {
        stop(sprintf("Residue mismatch in '%s': expected %s at position %d, found %s",
                     code, from_aa, pos, ref_vec[pos]))
      }
      
      if (!is.null(substitutions[[as.character(pos)]])) {
        stop(sprintf("Multiple substitutions at position %d", pos))
      }
      
      substitutions[[as.character(pos)]] <- to_aa
      
    } else if (grepl(mut_pattern_del, code)) {
      # Deletion: A20-
      m <- regexec(mut_pattern_del, code)
      parts <- regmatches(code, m)[[1]]
      
      if (length(parts) != 3) {
        stop(sprintf("Invalid deletion format: '%s'", code))
      }
      
      from_aa <- toupper(parts[2])
      pos <- as.integer(parts[3])
      
      if (is.na(pos) || pos < 1 || pos > ref_len) {
        stop(sprintf("Position out of range in '%s' (must be 1-%d)", code, ref_len))
      }
      
      if (toupper(ref_vec[pos]) != from_aa) {
        stop(sprintf("Residue mismatch in '%s': expected %s at position %d, found %s",
                     code, from_aa, pos, ref_vec[pos]))
      }
      
      if (pos %in% deletions) {
        stop(sprintf("Multiple deletions at position %d", pos))
      }
      
      deletions <- c(deletions, pos)
      
    } else if (grepl(mut_pattern_ins, code)) {
      # Insertion: -50P or -50HHHHHH
      m <- regexec(mut_pattern_ins, code)
      parts <- regmatches(code, m)[[1]]
      
      if (length(parts) != 3) {
        stop(sprintf("Invalid insertion format: '%s'", code))
      }
      
      pos <- as.integer(parts[2])
      ins_str <- toupper(parts[3])
      ins_vec <- strsplit(ins_str, "")[[1]]
      
      if (is.na(pos) || pos < 1 || pos > (ref_len + 1)) {
        stop(sprintf("Position out of range in '%s' (must be 1-%d)", code, ref_len + 1))
      }
      
      if (!is.null(insertions[[as.character(pos)]])) {
        stop(sprintf("Multiple insertions at position %d", pos))
      }
      
      insertions[[as.character(pos)]] <- ins_vec
      
    } else {
      stop(sprintf("Invalid mutation format: '%s'", code))
    }
  }

  # Check for conflicts between deletions and substitutions
  for (pos in deletions) {
    if (!is.null(substitutions[[as.character(pos)]])) {
      stop(sprintf("Conflicting mutations at position %d (both deletion and substitution)", pos))
    }
  }

  # Build result sequence
  result <- character(0)
  
  for (i in seq_len(ref_len + 1)) {
    # Check for insertion before position i
    if (!is.null(insertions[[as.character(i)]])) {
      result <- c(result, insertions[[as.character(i)]])
    }
    
    # If we've gone past the reference, stop
    if (i > ref_len) break
    
    # Check if position i is deleted
    if (i %in% deletions) {
      next
    }
    
    # Check if position i is substituted
    if (!is.null(substitutions[[as.character(i)]])) {
      result <- c(result, substitutions[[as.character(i)]])
    } else {
      result <- c(result, ref_vec[i])
    }
  }

  paste(result, collapse = "")
}


# -------------------------------------------------------------------------
# seq_to_mut(): Mutant sequence back to mutation notation
# -------------------------------------------------------------------------
#' Convert mutant sequence back to mutation notation
#'
#' @param ref Character string of the reference protein sequence
#' @param mutant Character vector of mutant sequence(s) to analyze
#' @param sep Character string separator to use between multiple mutations
#'   (default " / ")
#' @param match Numeric score for matching characters in alignment (default 2)
#' @param mismatch Numeric score for mismatching characters in alignment
#'   (default -2)
#' @param gap Numeric penalty for gaps in alignment (default -1)
#'
#' @return Character vector of mutation notation string(s), one per input
#'   mutant sequence. Empty string ("") if mutant equals reference.
#'
#' @details
#' This function performs the reverse operation of \code{mut_to_seq}. It aligns
#' each mutant sequence to the reference using Needleman-Wunsch global alignment,
#' then extracts mutations in the standard notation:
#'
#' - Substitution: \code{"A20S"} (reference A at position 20 changed to S)
#' - Deletion: \code{"A50-"} (reference A at position 50 deleted)
#' - Insertion: \code{"-50P"} or \code{"-50HHHHHH"} (residues inserted before
#'   position 50 in reference)
#'
#' All positions are 1-based and refer to the reference sequence. Multiple
#' mutations are separated by \code{sep}.
#'
#' The alignment scoring parameters can be adjusted. Default values favor gaps
#' over mismatches, which produces biologically intuitive alignments for typical
#' protein engineering scenarios.
#'
#' \strong{Important limitation:} The reverse mapping from sequence to mutation
#' notation is not always unique. Multiple mutation descriptions can produce
#' identical mutant sequences, particularly with:
#' \itemize{
#'   \item Deletions in repeated regions (e.g., \code{"AAAA"} → \code{"AA"}
#'     could be \code{"A1- / A2-"} or \code{"A3- / A4-"})
#'   \item Insertions and deletions near each other
#'   \item Substitutions in homopolymeric regions
#' }
#'
#' The roundtrip guarantee is: if \code{mut_to_seq(ref, mut1)} produces the
#' same sequence as \code{mut_to_seq(ref, mut2)}, then \code{mut1} and
#' \code{mut2} are functionally equivalent, even if their notation differs.
#' The alignment algorithm will return one valid representation, but it may
#' not match the original notation used to generate the sequence.
#'
#' @keywords internal
#' @noRd
seq_to_mut <- function(ref, mutant, sep = " / ",
                       match = 2, mismatch = -2, gap = -1) {
  # Input validation
  if (!is.character(ref) || length(ref) != 1 || nchar(ref) == 0) {
    stop("'ref' must be a single non-empty character string")
  }
  
  if (missing(mutant) || is.null(mutant)) {
    stop("'mutant' is required")
  }
  
  if (!is.character(mutant)) {
    stop("'mutant' must be a character vector")
  }
  
  if (any(nchar(mutant) == 0)) {
    stop("'mutant' contains empty strings")
  }
  
  if (!is.character(sep) || length(sep) != 1 || nchar(sep) == 0) {
    stop("'sep' must be a single non-empty character string")
  }
  
  if (!is.numeric(match) || length(match) != 1) {
    stop("'match' must be a single numeric value")
  }
  
  if (!is.numeric(mismatch) || length(mismatch) != 1) {
    stop("'mismatch' must be a single numeric value")
  }
  
  if (!is.numeric(gap) || length(gap) != 1) {
    stop("'gap' must be a single numeric value")
  }
  
  # Process each mutant sequence
  result <- lapply(mutant, function(mut_seq) {
    # If identical to reference, return empty string
    if (mut_seq == ref) {
      return("")
    }
    
    # Align sequences
    aln <- needleman_wunsch(ref, mut_seq, match = match, 
                           mismatch = mismatch, gap = gap)
    ref_aln <- strsplit(aln$aligned1, "")[[1]]
    mut_aln <- strsplit(aln$aligned2, "")[[1]]
    
    # Track mutations
    mutations <- character(0)
    ref_pos <- 0
    i <- 1
    
    while (i <= length(ref_aln)) {
      ref_char <- ref_aln[i]
      mut_char <- mut_aln[i]
      
      if (ref_char != "-") {
        ref_pos <- ref_pos + 1
      }
      
      if (ref_char == "-") {
        # Insertion: collect consecutive insertions
        ins_chars <- character(0)
        j <- i
        while (j <= length(ref_aln) && ref_aln[j] == "-") {
          ins_chars <- c(ins_chars, mut_aln[j])
          j <- j + 1
        }
        
        # Insert before position ref_pos + 1
        insert_pos <- ref_pos + 1
        ins_str <- paste(toupper(ins_chars), collapse = "")
        mutations <- c(mutations, paste0("-", insert_pos, ins_str))
        
        i <- j
        
      } else if (mut_char == "-") {
        # Deletion at position ref_pos
        mutations <- c(mutations, paste0(toupper(ref_char), ref_pos, "-"))
        i <- i + 1
        
      } else if (toupper(ref_char) != toupper(mut_char)) {
        # Substitution at position ref_pos
        mutations <- c(mutations, paste0(toupper(ref_char), ref_pos, toupper(mut_char)))
        i <- i + 1
        
      } else {
        # Match - do nothing
        i <- i + 1
      }
    }
    
    paste(mutations, collapse = sep)
  })
  
  unlist(result)
}


# -------------------------------------------------------------------------
# parse_fasta(): Parse FASTA file into named list
# -------------------------------------------------------------------------
parse_fasta <- function(file) {
  if (!file.exists(file)) {
    stop(sprintf("File not found: %s", file))
  }
  
  lines <- readLines(file, warn = FALSE)
  
  if (length(lines) == 0) {
    stop("Empty fasta file")
  }
  
  seq_list <- list()
  current_name <- NULL
  current_seq <- character(0)
  
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0) next
    
    if (grepl("^>", line)) {
      # Save previous sequence
      if (!is.null(current_name)) {
        seq_list[[current_name]] <- paste(current_seq, collapse = "")
      }
      
      # Start new sequence
      current_name <- sub("^>\\s*", "", line)
      current_seq <- character(0)
    } else {
      current_seq <- c(current_seq, line)
    }
  }
  
  # Save last sequence
  if (!is.null(current_name)) {
    seq_list[[current_name]] <- paste(current_seq, collapse = "")
  }
  
  if (length(seq_list) == 0) {
    stop("No sequences found in fasta file")
  }
  
  seq_list
}


# -------------------------------------------------------------------------
# seq_to_mat(): Convert sequences to mutation matrix
# -------------------------------------------------------------------------
#' Convert sequences to mutation matrix
#'
#' @param ref Character string of the reference protein sequence
#' @param seq Either a named list of sequences or a fasta file path
#' @param align Logical indicating whether to perform alignment. If FALSE,
#'   assumes all sequences are the same length as reference with no indels
#'   (faster). Default is TRUE.
#' @param sep Character string separator used in mutation notation when
#'   align = TRUE (default " / ")
#' @param prepend_pos Logical indicating whether to prepend position names
#'   to mutated residues. If TRUE, values like "S" in column "A20" become
#'   "A20S". Reference positions (empty strings) remain empty. Default is FALSE.
#'
#' @return Data frame with sequence names as rownames, mutation positions as
#'   columns (e.g., "A20"), and mutated residues as values. Positions matching
#'   reference are empty strings. Only positions that are mutated in at least
#'   one sequence are included. If \code{prepend_pos = TRUE}, non-empty values
#'   are formatted as "A20S" instead of "S".
#'
#' @keywords internal
#' @noRd
seq_to_mat <- function(ref, seq, align = FALSE, sep = " / ", prepend_pos = FALSE) {
  # Input validation
  if (!is.character(ref) || length(ref) != 1 || nchar(ref) == 0) {
    stop("'ref' must be a single non-empty character string")
  }
  
  if (!is.logical(align) || length(align) != 1) {
    stop("'align' must be TRUE or FALSE")
  }
  
  if (!is.logical(prepend_pos) || length(prepend_pos) != 1) {
    stop("'prepend_pos' must be TRUE or FALSE")
  }
  
  # Handle list vs file input
  if (is.list(seq)) {
    if (is.null(names(seq)) || any(names(seq) == "")) {
      stop("'seq' must be a named list")
    }
    seq_list <- seq
  } else if (is.character(seq) && length(seq) == 1) {
    seq_list <- parse_fasta(seq)
  } else {
    stop("'seq' must be a named list or a fasta file path")
  }
  
  seq_names <- names(seq_list)
  seq_vectors <- unlist(seq_list, use.names = FALSE)
  
  ref_vec <- strsplit(ref, "")[[1]]
  ref_len <- length(ref_vec)
  
  # Fast path: no alignment needed
  if (!align) {
    # Check all sequences are same length as reference
    seq_lengths <- nchar(seq_vectors)
    if (any(seq_lengths != ref_len)) {
      bad_idx <- which(seq_lengths != ref_len)
      stop(sprintf(
        "When align = FALSE, all sequences must be same length as reference (%d). Mismatched: %s",
        ref_len,
        paste(seq_names[bad_idx], collapse = ", ")
      ))
    }
    
    # Find mutated positions
    mutated_positions <- integer(0)
    for (i in seq_along(seq_vectors)) {
      seq_vec <- strsplit(seq_vectors[i], "")[[1]]
      diffs <- which(toupper(seq_vec) != toupper(ref_vec))
      mutated_positions <- unique(c(mutated_positions, diffs))
    }
    
    if (length(mutated_positions) == 0) {
      return(data.frame(row.names = seq_names))
    }
    
    mutated_positions <- sort(mutated_positions)
    col_names <- paste0(toupper(ref_vec[mutated_positions]), mutated_positions)
    
    # Build matrix
    result_mat <- matrix(
      "",
      nrow = length(seq_names),
      ncol = length(mutated_positions),
      dimnames = list(seq_names, col_names)
    )
    
    for (i in seq_along(seq_names)) {
      seq_vec <- strsplit(seq_vectors[i], "")[[1]]
      for (j in seq_along(mutated_positions)) {
        pos <- mutated_positions[j]
        if (toupper(seq_vec[pos]) != toupper(ref_vec[pos])) {
          result_mat[i, j] <- toupper(seq_vec[pos])
        }
      }
    }
    
    # Prepend column names to values if requested
    if (prepend_pos) {
      for (col in colnames(result_mat)) {
        non_empty <- result_mat[, col] != ""
        result_mat[non_empty, col] <- paste0(col, result_mat[non_empty, col])
      }
    }
    
    return(as.data.frame(result_mat, stringsAsFactors = FALSE))
  }
  
  # Slow path: alignment-based
  mutations_list <- seq_to_mut(ref, seq_vectors, sep = sep)
  names(mutations_list) <- seq_names
  
  # Collect all mutated positions
  mutated_positions <- list()
  
  for (i in seq_along(mutations_list)) {
    muts <- mutations_list[i]
    if (muts == "") next
    
    mut_codes <- unlist(strsplit(muts, sep, fixed = TRUE), use.names = FALSE)
    mut_codes <- trimws(mut_codes)
    mut_codes <- mut_codes[nchar(mut_codes) > 0]
    
    for (code in mut_codes) {
      # Skip insertions
      if (grepl("^-", code)) next
      
      # Parse position
      m <- regexec("^([A-Za-z\\*])([0-9]+)", code)
      parts <- regmatches(code, m)[[1]]
      
      if (length(parts) >= 3) {
        ref_aa <- toupper(parts[2])
        pos <- as.integer(parts[3])
        pos_name <- paste0(ref_aa, pos)
        
        if (!pos_name %in% names(mutated_positions)) {
          mutated_positions[[pos_name]] <- list(pos = pos, ref_aa = ref_aa)
        }
      }
    }
  }
  
  # Return empty dataframe if no mutations
  if (length(mutated_positions) == 0) {
    return(data.frame(row.names = seq_names))
  }
  
  # Sort positions numerically
  pos_order <- order(sapply(mutated_positions, function(x) x$pos))
  mutated_positions <- mutated_positions[pos_order]
  col_names <- names(mutated_positions)
  
  # Initialize matrix with empty strings
  result_mat <- matrix(
    "",
    nrow = length(seq_names),
    ncol = length(col_names),
    dimnames = list(seq_names, col_names)
  )
  
  # Fill matrix
  for (i in seq_along(seq_names)) {
    seq_name <- seq_names[i]
    muts <- mutations_list[seq_name]
    
    if (muts != "") {
      mut_codes <- unlist(strsplit(muts, sep, fixed = TRUE), use.names = FALSE)
      mut_codes <- trimws(mut_codes)
      mut_codes <- mut_codes[nchar(mut_codes) > 0]
      
      for (code in mut_codes) {
        if (grepl("^-", code)) next
        
        m <- regexec("^([A-Za-z\\*])([0-9]+)(.+)$", code)
        parts <- regmatches(code, m)[[1]]
        
        if (length(parts) == 4) {
          ref_aa <- toupper(parts[2])
          pos <- as.integer(parts[3])
          to_aa <- toupper(parts[4])
          
          pos_name <- paste0(ref_aa, pos)
          if (pos_name %in% col_names) {
            result_mat[seq_name, pos_name] <- to_aa
          }
        }
      }
    }
  }
  
  # Prepend column names to values if requested
  if (prepend_pos) {
    for (col in colnames(result_mat)) {
      non_empty <- result_mat[, col] != ""
      result_mat[non_empty, col] <- paste0(col, result_mat[non_empty, col])
    }
  }
  
  as.data.frame(result_mat, stringsAsFactors = FALSE)
}