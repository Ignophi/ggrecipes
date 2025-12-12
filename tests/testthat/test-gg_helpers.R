# ------------------------------------------------------------------------------
# needleman_wunsch(): Global sequence alignment – tests
# ------------------------------------------------------------------------------

test_that("needleman_wunsch() validates inputs", {
  # seq1 validation
  expect_error(
    needleman_wunsch(1, "ACGT"),
    "'seq1' must be a single non-empty character string"
  )
  expect_error(
    needleman_wunsch(c("A", "C"), "ACGT"),
    "'seq1' must be a single non-empty character string"
  )
  expect_error(
    needleman_wunsch("", "ACGT"),
    "'seq1' must be a single non-empty character string"
  )

  # seq2 validation
  expect_error(
    needleman_wunsch("ACGT", 1),
    "'seq2' must be a single non-empty character string"
  )
  expect_error(
    needleman_wunsch("ACGT", c("A", "C")),
    "'seq2' must be a single non-empty character string"
  )
  expect_error(
    needleman_wunsch("ACGT", ""),
    "'seq2' must be a single non-empty character string"
  )

  # scoring parameter validation
  expect_error(
    needleman_wunsch("ACGT", "ACGT", match = c(1, 2)),
    "'match' must be a single numeric value"
  )
  expect_error(
    needleman_wunsch("ACGT", "ACGT", mismatch = "x"),
    "'mismatch' must be a single numeric value"
  )
  expect_error(
    needleman_wunsch("ACGT", "ACGT", gap = c(-1, -2)),
    "'gap' must be a single numeric value"
  )
})

# Basic alignment behaviour

test_that("needleman_wunsch() aligns identical sequences without gaps", {
  res <- needleman_wunsch("GATTACA", "GATTACA")

  expect_equal(res$aligned1, "GATTACA")
  expect_equal(res$aligned2, "GATTACA")
  # 7 matches * 1
  expect_equal(res$score, 7)
})

test_that("needleman_wunsch() handles all mismatches without gaps when gap penalty is harsher", {
  res <- needleman_wunsch("AAAA", "TTTT")  # match = 1, mismatch = -1, gap = -2

  expect_equal(res$aligned1, "AAAA")
  expect_equal(res$aligned2, "TTTT")
  # 4 mismatches * -1
  expect_equal(res$score, -4)
})

test_that("needleman_wunsch() works for single-character sequences", {
  res <- needleman_wunsch("A", "G")  # default scoring

  expect_equal(res$aligned1, "A")
  expect_equal(res$aligned2, "G")
  expect_equal(res$score, -1)       # 1 mismatch
})

# Gaps and differing sequence lengths

test_that("needleman_wunsch() inserts terminal gap for shorter sequence", {
  res <- needleman_wunsch("ACGT", "ACG")

  expect_equal(res$aligned1, "ACGT")
  expect_equal(res$aligned2, "ACG-")
  # matches: 3 * 1, gap: -3
  expect_equal(res$score, 0)
})

test_that("needleman_wunsch() inserts leading gap when beneficial", {
  res <- needleman_wunsch("CGT", "ACGT")

  expect_equal(res$aligned1, "-CGT")
  expect_equal(res$aligned2, "ACGT")
  # matches: C,C and G,G and T,T => 3 * 1; one gap: -3
  expect_equal(res$score, 0)
})

test_that("needleman_wunsch() allows zero gap penalty and uses gaps instead of mismatches", {
  res <- needleman_wunsch("AA", "TTAA", match = 1, mismatch = -1, gap = 0)

  # With no gap penalty, it should skip the leading 'TT' via gaps
  expect_equal(res$aligned1, "--AA")
  expect_equal(res$aligned2, "TTAA")
  # 2 matches * 1, 2 gaps * 0
  expect_equal(res$score, 2)
})

# Classic example and parameterisation

test_that("needleman_wunsch() reproduces classic textbook alignment with custom gap", {
  # seq1 and seq2 from a common Needleman-Wunsch example
  res <- needleman_wunsch("GATTACA", "GCATGCU", match = 1, mismatch = -1, gap = -1)

  expect_equal(res$aligned1, "G-ATTACA")
  expect_equal(res$aligned2, "GCA-TGCU")
  # Score: +1, -1, +1, -1, +1, -1, +1, -1 = 0
  expect_equal(res$score, 0)
})

# Symmetry and internal consistency

test_that("needleman_wunsch() yields symmetric alignments under symmetric scoring", {
  seq1 <- "GATTACA"
  seq2 <- "GCATGCU"

  res1 <- needleman_wunsch(seq1, seq2, match = 1, mismatch = -1, gap = -1)
  res2 <- needleman_wunsch(seq2, seq1, match = 1, mismatch = -1, gap = -1)

  # Score should be identical in both directions
  expect_equal(res1$score, res2$score)

  # Alignments should swap when sequences are swapped
  expect_equal(res1$aligned1, res2$aligned2)
  expect_equal(res1$aligned2, res2$aligned1)
})

test_that("needleman_wunsch() score matches sum over returned alignment", {
  # Helper to recompute score from aligned strings
  score_alignment <- function(aligned1, aligned2, match, mismatch, gap) {
    s1 <- strsplit(aligned1, "")[[1]]
    s2 <- strsplit(aligned2, "")[[1]]

    stopifnot(length(s1) == length(s2))

    pair_scores <- vapply(
      seq_along(s1),
      function(i) {
        if (s1[i] == "-" || s2[i] == "-") {
          gap
        } else if (s1[i] == s2[i]) {
          match
        } else {
          mismatch
        }
      },
      numeric(1)
    )

    sum(pair_scores)
  }

  # Use a couple of non-trivial examples
  res1 <- needleman_wunsch("GATTACA", "GCATGCU", match = 1, mismatch = -1, gap = -1)
  res2 <- needleman_wunsch("ACGT", "ACG", match = 2, mismatch = -1, gap = -2)

  expect_equal(
    res1$score,
    score_alignment(res1$aligned1, res1$aligned2, match = 1, mismatch = -1, gap = -1)
  )
  expect_equal(
    res2$score,
    score_alignment(res2$aligned1, res2$aligned2, match = 2, mismatch = -1, gap = -2)
  )
})


# ------------------------------------------------------------------------------
# mut_to_seq(): Apply mutations to protein sequence â€" tests
# ------------------------------------------------------------------------------

test_that("mut_to_seq() validates inputs", {
  # ref validation
  expect_error(
    mut_to_seq(1, "M1A"),
    "'ref' must be a single non-empty character string"
  )
  expect_error(
    mut_to_seq(c("A", "C"), "M1A"),
    "'ref' must be a single non-empty character string"
  )
  expect_error(
    mut_to_seq("", "M1A"),
    "'ref' must be a single non-empty character string"
  )
  
  # mut validation
  expect_error(
    mut_to_seq("ACGT", 123),
    "'mut' must be a character vector"
  )
  
  # sep validation
  expect_error(
    mut_to_seq("ACGT", "A1C", sep = c("/", ",")),
    "'sep' must be a single non-empty character string"
  )
  expect_error(
    mut_to_seq("ACGT", "A1C", sep = ""),
    "'sep' must be a single non-empty character string"
  )
})

test_that("mut_to_seq() returns reference when mut is NULL or missing", {
  ref <- "MKTAYIAK"
  
  # Purpose: ensure no mutations returns original sequence
  expect_equal(mut_to_seq(ref), ref)
  expect_equal(mut_to_seq(ref, NULL), ref)
  expect_equal(mut_to_seq(ref, ""), ref)
  expect_equal(mut_to_seq(ref, "   "), ref)
})

test_that("mut_to_seq() handles single substitutions", {
  ref <- "MKTAYIAK"
  
  # Purpose: single substitutions should work
  expect_equal(mut_to_seq(ref, "M1A"), "AKTAYIAK")
  expect_equal(mut_to_seq(ref, "K2R"), "MRTAYIAK")
  expect_equal(mut_to_seq(ref, "K8R"), "MKTAYIAR")
})

test_that("mut_to_seq() handles single deletions", {
  ref <- "MKTAYIAK"
  
  # Purpose: single deletions should work
  expect_equal(mut_to_seq(ref, "M1-"), "KTAYIAK")
  expect_equal(mut_to_seq(ref, "K2-"), "MTAYIAK")
  expect_equal(mut_to_seq(ref, "K8-"), "MKTAYIA")
})

test_that("mut_to_seq() handles single insertions", {
  ref <- "MKTAYIAK"
  
  # Purpose: single insertions should work
  expect_equal(mut_to_seq(ref, "-1H"), "HMKTAYIAK")
  expect_equal(mut_to_seq(ref, "-1HHH"), "HHHMKTAYIAK")
  expect_equal(mut_to_seq(ref, "-2P"), "MPKTAYIAK")
  expect_equal(mut_to_seq(ref, "-9HHHHHH"), "MKTAYIAKHHHHHH")
})

test_that("mut_to_seq() handles multiple substitutions", {
  ref <- "MKTAYIAK"
  
  # Purpose: multiple substitutions at different positions
  expect_equal(mut_to_seq(ref, "M1A / K2R"), "ARTAYIAK")
  expect_equal(mut_to_seq(ref, "M1A / K2R / K8L"), "ARTAYIAL")
  expect_equal(mut_to_seq(ref, "K2R / M1A"), "ARTAYIAK")  # order independent
})

test_that("mut_to_seq() handles multiple deletions", {
  ref <- "MKTAYIAK"
  
  # Purpose: multiple deletions at different positions
  expect_equal(mut_to_seq(ref, "M1- / K2-"), "TAYIAK")
  expect_equal(mut_to_seq(ref, "M1- / T3-"), "KAYIAK")
  expect_equal(mut_to_seq(ref, "K2- / M1-"), "TAYIAK")  # order independent
})

test_that("mut_to_seq() handles multiple insertions", {
  ref <- "MKTAYIAK"
  
  # Purpose: multiple insertions at different positions
  expect_equal(mut_to_seq(ref, "-1H / -2P"), "HMPKTAYIAK")
  expect_equal(mut_to_seq(ref, "-1HH / -9GGG"), "HHMKTAYIAKGGG")
  # order independent
  expect_equal(mut_to_seq(ref, "-9GGG / -1HH"), "HHMKTAYIAKGGG")
})

test_that("mut_to_seq() handles mixed mutation types", {
  ref <- "MKTAYIAK"
  
  # Purpose: combinations of substitutions, deletions, and insertions
  expect_equal(mut_to_seq(ref, "M1A / K2-"), "ATAYIAK")
  expect_equal(mut_to_seq(ref, "M1- / -1G"), "GKTAYIAK")
  expect_equal(mut_to_seq(ref, "M1A / -1HHH"), "HHHAKTAYIAK")
  expect_equal(mut_to_seq(ref, "M1A / K2R / T3- / -5P"), "ARAPYIAK")
})

test_that("mut_to_seq() validates mutation format", {
  ref <- "MKTAYIAK"
  
  # Purpose: invalid formats should raise errors
  expect_error(mut_to_seq(ref, "INVALID"), "Invalid mutation format")
  expect_error(mut_to_seq(ref, "M1"), "Invalid mutation format")
  expect_error(mut_to_seq(ref, "1A"), "Invalid mutation format")
  expect_error(mut_to_seq(ref, "M1A2"), "Invalid mutation format")
  expect_error(mut_to_seq(ref, "--1A"), "Invalid mutation format")
})

test_that("mut_to_seq() validates position ranges", {
  ref <- "MKTAYIAK"
  
  # Purpose: positions must be within valid range
  expect_error(mut_to_seq(ref, "M0A"), "Position out of range")
  expect_error(mut_to_seq(ref, "M9A"), "Position out of range")
  expect_error(mut_to_seq(ref, "M999S"), "Position out of range")
  expect_error(mut_to_seq(ref, "A0-"), "Position out of range")
  expect_error(mut_to_seq(ref, "A9-"), "Position out of range")
  expect_error(mut_to_seq(ref, "-0P"), "Position out of range")
  expect_error(mut_to_seq(ref, "-10P"), "Position out of range")
})

test_that("mut_to_seq() validates reference residue matches", {
  ref <- "MKTAYIAK"
  
  # Purpose: mutations must specify correct reference residue
  expect_error(mut_to_seq(ref, "A1K"), "Residue mismatch.*expected A.*found M")
  expect_error(mut_to_seq(ref, "R2A"), "Residue mismatch.*expected R.*found K")
  expect_error(mut_to_seq(ref, "X8-"), "Residue mismatch.*expected X.*found K")
})

test_that("mut_to_seq() detects conflicting mutations", {
  ref <- "MKTAYIAK"
  
  # Purpose: multiple mutations at same position should error
  expect_error(mut_to_seq(ref, "M1A / M1C"), "Multiple substitutions at position 1")
  expect_error(mut_to_seq(ref, "M1- / M1-"), "Multiple deletions at position 1")
  expect_error(mut_to_seq(ref, "M1A / M1-"), "Conflicting mutations at position 1")
  expect_error(mut_to_seq(ref, "-2P / -2G"), "Multiple insertions at position 2")
})

test_that("mut_to_seq() handles edge positions correctly", {
  ref <- "MKTAYIAK"
  
  # Purpose: first and last positions should work correctly
  expect_equal(mut_to_seq(ref, "M1A"), "AKTAYIAK")
  expect_equal(mut_to_seq(ref, "K8R"), "MKTAYIAR")
  expect_equal(mut_to_seq(ref, "M1-"), "KTAYIAK")
  expect_equal(mut_to_seq(ref, "K8-"), "MKTAYIA")
  expect_equal(mut_to_seq(ref, "-1HHH"), "HHHMKTAYIAK")
  expect_equal(mut_to_seq(ref, "-9HHHHHH"), "MKTAYIAKHHHHHH")
})

test_that("mut_to_seq() handles custom separators", {
  ref <- "MKTAYIAK"
  
  # Purpose: different separators should work
  expect_equal(mut_to_seq(ref, "M1A,K2R", sep = ","), "ARTAYIAK")
  expect_equal(mut_to_seq(ref, "M1A;K2R;T3S", sep = ";"), "ARSAYIAK")
  expect_equal(mut_to_seq(ref, "M1A | K2R", sep = " | "), "ARTAYIAK")
})

test_that("mut_to_seq() is case-insensitive for amino acids", {
  ref <- "MKTAYIAK"
  
  # Purpose: lowercase mutations should work
  expect_equal(mut_to_seq(ref, "m1a"), "AKTAYIAK")
  expect_equal(mut_to_seq(ref, "M1a"), "AKTAYIAK")
  expect_equal(mut_to_seq(ref, "m1A"), "AKTAYIAK")
})

test_that("mut_to_seq() handles very short sequences", {
  # Purpose: single amino acid sequences should work
  expect_equal(mut_to_seq("M", "M1A"), "A")
  expect_equal(mut_to_seq("M", "M1-"), "")
  expect_equal(mut_to_seq("M", "-1H"), "HM")
  expect_equal(mut_to_seq("M", "-2H"), "MH")
})

test_that("mut_to_seq() handles complex real-world scenarios", {
  ref <- "MKTAYIAKQRQISFVKSHFSRQL"
  
  # Purpose: realistic protein engineering mutations
  
  # Add N-terminal His-tag and substitute lysines
  result1 <- mut_to_seq(ref, "-1HHHHHH / K2R / K16R")
  expect_equal(result1, "HHHHHHMRTAYIAKQRQISFVRSHFSRQL")
  
  # Delete first residue and add C-terminal His-tag
  result2 <- mut_to_seq(ref, "M1- / -24HHHHHH")
  expect_equal(result2, "KTAYIAKQRQISFVKSHFSRQLHHHHHH")
  
  # Multiple point mutations for stability
  result3 <- mut_to_seq(ref, "M1L / A4G / K8R / Q9E")
  expect_equal(result3, "LKTGYIARERQISFVKSHFSRQL")
  
  # Complex combination
  result4 <- mut_to_seq(ref, "M1- / -1G / K2R / -24GGS")
  expect_equal(result4, "GRTAYIAKQRQISFVKSHFSRQLGGS")
})

test_that("mut_to_seq() mutation order independence", {
  ref <- "MKTAYIAK"
  
  # Purpose: different orderings should give identical results
  mut_set1 <- c("M1A", "K2R", "T3S", "-1HH", "K8-")
  mut_set2 <- c("K8-", "T3S", "-1HH", "M1A", "K2R")
  mut_set3 <- c("K2R", "K8-", "M1A", "T3S", "-1HH")
  
  result1 <- mut_to_seq(ref, paste(mut_set1, collapse = " / "))
  result2 <- mut_to_seq(ref, paste(mut_set2, collapse = " / "))
  result3 <- mut_to_seq(ref, paste(mut_set3, collapse = " / "))
  
  expect_equal(result1, result2)
  expect_equal(result2, result3)
  expect_equal(result1, result3)
})

test_that("mut_to_seq() handles special amino acid characters", {
  ref <- "M*TAYIAK"
  
  # Purpose: stop codon (*) should be handled
  expect_equal(mut_to_seq(ref, "*2A"), "MATAYIAK")
  expect_equal(mut_to_seq(ref, "*2-"), "MTAYIAK")
  expect_equal(mut_to_seq(ref, "M1*"), "**TAYIAK")
})

test_that("mut_to_seq() handles whitespace in mutation strings", {
  ref <- "MKTAYIAK"
  
  # Purpose: extra whitespace should be handled gracefully
  expect_equal(mut_to_seq(ref, "  M1A  /  K2R  "), "ARTAYIAK")
  expect_equal(mut_to_seq(ref, "M1A / K2R / "), "ARTAYIAK")
  expect_equal(mut_to_seq(ref, " / M1A / K2R"), "ARTAYIAK")
})

test_that("mut_to_seq() handles vector of mutation strings", {
  ref <- "MKTAYIAK"
  
  # Purpose: mutation vector with multiple elements should work
  result <- mut_to_seq(ref, c("M1A / K2R", "T3S / -1HH"))
  expect_equal(result, "HHARSAYIAK")
})

test_that("mut_to_seq() handles all deletions", {
  ref <- "ACDE"
  
  # Purpose: deleting all residues should give empty string
  result <- mut_to_seq(ref, "A1- / C2- / D3- / E4-")
  expect_equal(result, "")
})

test_that("mut_to_seq() handles insertions creating long sequences", {
  ref <- "MK"
  
  # Purpose: large insertions should work
  result <- mut_to_seq(ref, "-1HHHHHH / -3GSGSGS")
  expect_equal(result, "HHHHHHMKGSGSGS")
})

test_that("mut_to_seq() positions refer to reference throughout", {
  ref <- "MKTAYIAK"
  
  # Purpose: critical test that positions always refer to original reference
  
  # Even though M1 is deleted, K2R still refers to K at position 2 in reference
  result1 <- mut_to_seq(ref, "M1- / K2R")
  expect_equal(result1, "RTAYIAK")
  
  # Even though HHH is inserted at start, M1A still refers to M at position 1
  result2 <- mut_to_seq(ref, "-1HHH / M1A")
  expect_equal(result2, "HHHAKTAYIAK")
  
  # Complex case: delete, insert, substitute all at adjacent positions
  result3 <- mut_to_seq(ref, "M1- / -2PP / K2R")
  expect_equal(result3, "PPRTAYIAK")
})


# ------------------------------------------------------------------------------
# seq_to_mut(): Convert mutant sequence to mutation notation – tests
# ------------------------------------------------------------------------------

test_that("seq_to_mut() validates inputs", {
  # ref validation
  expect_error(
    seq_to_mut(1, "ACGT"),
    "'ref' must be a single non-empty character string"
  )
  expect_error(
    seq_to_mut(c("A", "C"), "ACGT"),
    "'ref' must be a single non-empty character string"
  )
  expect_error(
    seq_to_mut("", "ACGT"),
    "'ref' must be a single non-empty character string"
  )
  
  # mutant validation
  expect_error(
    seq_to_mut("ACGT"),
    "'mutant' is required"
  )
  expect_error(
    seq_to_mut("ACGT", NULL),
    "'mutant' is required"
  )
  expect_error(
    seq_to_mut("ACGT", 123),
    "'mutant' must be a character vector"
  )
  expect_error(
    seq_to_mut("ACGT", ""),
    "'mutant' contains empty strings"
  )
  expect_error(
    seq_to_mut("ACGT", c("ACGT", "")),
    "'mutant' contains empty strings"
  )
  
  # sep validation
  expect_error(
    seq_to_mut("ACGT", "ACGT", sep = c("/", ",")),
    "'sep' must be a single non-empty character string"
  )
  expect_error(
    seq_to_mut("ACGT", "ACGT", sep = ""),
    "'sep' must be a single non-empty character string"
  )
  
  # scoring parameter validation
  expect_error(
    seq_to_mut("ACGT", "ACGT", match = c(1, 2)),
    "'match' must be a single numeric value"
  )
  expect_error(
    seq_to_mut("ACGT", "ACGT", mismatch = "x"),
    "'mismatch' must be a single numeric value"
  )
  expect_error(
    seq_to_mut("ACGT", "ACGT", gap = c(-1, -2)),
    "'gap' must be a single numeric value"
  )
})

test_that("seq_to_mut() returns empty string for identical sequences", {
  ref <- "MKTAYIAK"
  
  # Purpose: identical sequence should return empty string
  expect_equal(seq_to_mut(ref, ref), "")
})

test_that("seq_to_mut() detects single substitutions", {
  ref <- "MKTAYIAK"
  
  # Purpose: single substitutions should be detected
  expect_equal(seq_to_mut(ref, "AKTAYIAK"), "M1A")
  expect_equal(seq_to_mut(ref, "MRTAYIAK"), "K2R")
  expect_equal(seq_to_mut(ref, "MKTAYIAR"), "K8R")
})

test_that("seq_to_mut() detects single deletions", {
  ref <- "MKTAYIAK"
  
  # Purpose: single deletions should be detected
  expect_equal(seq_to_mut(ref, "KTAYIAK"), "M1-")
  expect_equal(seq_to_mut(ref, "MTAYIAK"), "K2-")
  expect_equal(seq_to_mut(ref, "MKTAYIA"), "K8-")
})

test_that("seq_to_mut() detects single insertions", {
  ref <- "MKTAYIAK"
  
  # Purpose: single insertions should be detected
  expect_equal(seq_to_mut(ref, "HMKTAYIAK"), "-1H")
  expect_equal(seq_to_mut(ref, "HHHMKTAYIAK"), "-1HHH")
  expect_equal(seq_to_mut(ref, "MPKTAYIAK"), "-2P")
  expect_equal(seq_to_mut(ref, "MKTAYIAKH"), "-9H")
  expect_equal(seq_to_mut(ref, "MKTAYIAKHHHHHH"), "-9HHHHHH")
})

test_that("seq_to_mut() handles vector of mutants", {
  ref <- "MKTAYIAK"
  
  # Purpose: multiple mutants should return vector of mutation strings
  result <- seq_to_mut(ref, c("AKTAYIAK", "MKTAYIAR", "KTAYIAK"))
  expect_equal(result, c("M1A", "K8R", "M1-"))
  expect_length(result, 3)
})

test_that("seq_to_mut() handles custom separators", {
  ref <- "MKTAYIAK"
  mutant <- "ARTAYIAK"
  
  # Purpose: custom separator should be used
  result <- seq_to_mut(ref, mutant, sep = ",")
  expect_true(grepl(",", result))
  
  result2 <- seq_to_mut(ref, mutant, sep = ";")
  expect_true(grepl(";", result2))
})

# Helper function to compare mutation strings regardless of order
compare_mutations <- function(mut1, mut2, sep = " / ") {
  if (mut1 == "" && mut2 == "") return(TRUE)
  if (mut1 == "" || mut2 == "") return(FALSE)
  
  m1 <- sort(trimws(strsplit(mut1, sep, fixed = TRUE)[[1]]))
  m2 <- sort(trimws(strsplit(mut2, sep, fixed = TRUE)[[1]]))
  
  identical(m1, m2)
}

test_that("seq_to_mut() roundtrip with single substitution", {
  ref <- "MKTAYIAK"
  original_mut <- "M1A"
  
  # Purpose: roundtrip should preserve mutation
  mutant <- mut_to_seq(ref, original_mut)
  result_mut <- seq_to_mut(ref, mutant)
  
  expect_true(compare_mutations(original_mut, result_mut))
})

test_that("seq_to_mut() roundtrip with single deletion", {
  ref <- "MKTAYIAK"
  original_mut <- "M1-"
  
  # Purpose: roundtrip should preserve mutation
  mutant <- mut_to_seq(ref, original_mut)
  result_mut <- seq_to_mut(ref, mutant)
  
  expect_true(compare_mutations(original_mut, result_mut))
})

test_that("seq_to_mut() roundtrip with single insertion", {
  ref <- "MKTAYIAK"
  
  # Purpose: roundtrip should preserve insertions
  for (original_mut in c("-1H", "-1HHH", "-2P", "-9HHHHHH")) {
    mutant <- mut_to_seq(ref, original_mut)
    result_mut <- seq_to_mut(ref, mutant)
    expect_true(compare_mutations(original_mut, result_mut),
                info = sprintf("Failed for mutation: %s", original_mut))
  }
})

test_that("seq_to_mut() roundtrip with multiple substitutions", {
  ref <- "MKTAYIAK"
  original_mut <- "M1A / K2R / K8L"
  
  # Purpose: roundtrip should preserve all substitutions
  mutant <- mut_to_seq(ref, original_mut)
  result_mut <- seq_to_mut(ref, mutant)
  
  expect_true(compare_mutations(original_mut, result_mut))
})

test_that("seq_to_mut() roundtrip with multiple deletions", {
  ref <- "MKTAYIAK"
  original_mut <- "M1- / K2-"
  
  # Purpose: roundtrip should preserve all deletions
  mutant <- mut_to_seq(ref, original_mut)
  result_mut <- seq_to_mut(ref, mutant)
  
  expect_true(compare_mutations(original_mut, result_mut))
})

test_that("seq_to_mut() roundtrip with multiple insertions", {
  ref <- "MKTAYIAK"
  original_mut <- "-1HH / -9GGG"
  
  # Purpose: roundtrip should preserve all insertions
  mutant <- mut_to_seq(ref, original_mut)
  result_mut <- seq_to_mut(ref, mutant)
  
  expect_true(compare_mutations(original_mut, result_mut))
})

test_that("seq_to_mut() roundtrip with mixed mutation types", {
  ref <- "MKTAYIAK"
  
  # Purpose: roundtrip should handle complex combinations
  test_cases <- c(
    "M1A / -1HHH",
    "-1HHH / M1A / K2R / K8-"
  )
  
  original_mut <- test_cases[1]

  for (original_mut in test_cases) {
    mutant <- mut_to_seq(ref, original_mut)
    result_mut <- seq_to_mut(ref, mutant)
    expect_true(compare_mutations(original_mut, result_mut),
                info = sprintf("Failed for mutation: %s", original_mut))
  }
})

test_that("seq_to_mut() roundtrip with edge positions", {
  ref <- "MKTAYIAK"
  
  # Purpose: first and last positions should roundtrip correctly
  test_cases <- c(
    "M1A",
    "K8R",
    "M1-",
    "K8-",
    "-1HHH",
    "-9HHHHHH"
  )
  
  for (original_mut in test_cases) {
    mutant <- mut_to_seq(ref, original_mut)
    result_mut <- seq_to_mut(ref, mutant)
    expect_true(compare_mutations(original_mut, result_mut),
                info = sprintf("Failed for mutation: %s", original_mut))
  }
})

test_that("seq_to_mut() roundtrip with very short sequences", {
  ref <- "MK"
  
  # Purpose: minimal sequences should roundtrip
  test_cases <- c(
    "M1A",
    "M1-",
    "-1H",
    "-3H",
    "M1A / K2R",
    "-1HH / K2-"
  )
  
  for (original_mut in test_cases) {
    mutant <- mut_to_seq(ref, original_mut)
    result_mut <- seq_to_mut(ref, mutant)
    expect_true(compare_mutations(original_mut, result_mut),
                info = sprintf("Failed for mutation: %s", original_mut))
  }
})

test_that("seq_to_mut() roundtrip with complex real-world scenarios", {
  ref <- "MKTAYIAKQRQISFVKSHFSRQL"
  
  # Purpose: realistic protein engineering mutations should roundtrip
  test_cases <- c(
    "-1HHHHHH / K2R / K16R",
    "M1- / -24HHHHHH",
    "M1L / A4G / K8R / Q9E"
  )
  
  for (original_mut in test_cases) {
    mutant <- mut_to_seq(ref, original_mut)
    result_mut <- seq_to_mut(ref, mutant)
    expect_true(compare_mutations(original_mut, result_mut),
                info = sprintf("Failed for mutation: %s", original_mut))
  }
})

test_that("seq_to_mut() roundtrip is mutation order independent", {
  ref <- "MKTAYIAK"
  
  # Purpose: different orderings should produce equivalent results after roundtrip
  mut_set1 <- "M1A / K2R / T3S / -1HH / K8-"
  mut_set2 <- "K8- / T3S / -1HH / M1A / K2R"
  mut_set3 <- "K2R / K8- / M1A / T3S / -1HH"
  
  mutant1 <- mut_to_seq(ref, mut_set1)
  mutant2 <- mut_to_seq(ref, mut_set2)
  mutant3 <- mut_to_seq(ref, mut_set3)
  
  # All should produce same mutant sequence
  expect_equal(mutant1, mutant2)
  expect_equal(mutant2, mutant3)
  
  # Converting back should give equivalent mutation sets
  result1 <- seq_to_mut(ref, mutant1)
  result2 <- seq_to_mut(ref, mutant2)
  result3 <- seq_to_mut(ref, mutant3)
  
  expect_true(compare_mutations(result1, result2))
  expect_true(compare_mutations(result2, result3))
  expect_true(compare_mutations(mut_set1, result1))
})

test_that("seq_to_mut() roundtrip with special amino acid characters", {
  ref <- "M*TAYIAK"
  
  # Purpose: stop codon should roundtrip
  test_cases <- c(
    "*2A",
    "*2-",
    "M1*"
  )
  
  for (original_mut in test_cases) {
    mutant <- mut_to_seq(ref, original_mut)
    result_mut <- seq_to_mut(ref, mutant)
    expect_true(compare_mutations(original_mut, result_mut),
                info = sprintf("Failed for mutation: %s", original_mut))
  }
})

test_that("seq_to_mut() roundtrip with vector of mutations", {
  ref <- "MKTAYIAK"
  
  # Purpose: vector input should roundtrip correctly
  original_muts <- c("M1A / K2R", "T3S / -1HH", "K8-")
  
  mutants <- vapply(original_muts, function(m) mut_to_seq(ref, m), character(1))
  result_muts <- seq_to_mut(ref, mutants)
  
  expect_length(result_muts, 3)
  for (i in seq_along(original_muts)) {
    expect_true(compare_mutations(original_muts[i], result_muts[i]),
                info = sprintf("Failed for mutation %d: %s", i, original_muts[i]))
  }
})

test_that("seq_to_mut() handles consecutive insertions correctly", {
  ref <- "MKTAYIAK"
  
  # Purpose: multiple consecutive insertions should be grouped
  mutant <- "HHHPPPGGGMKTAYIAK"
  result <- seq_to_mut(ref, mutant)
  
  # Should detect as single insertion block
  expect_equal(result, "-1HHHPPPGGG")
  
  # Roundtrip verification
  expect_equal(mut_to_seq(ref, result), mutant)
})

test_that("seq_to_mut() case insensitivity matches mut_to_seq", {
  ref <- "MKTAYIAK"
  
  # Purpose: lowercase mutations should roundtrip with uppercase result
  original_mut <- "m1a / k2r"
  mutant <- mut_to_seq(ref, original_mut)
  result_mut <- seq_to_mut(ref, mutant)
  
  # Result should be uppercase
  expect_true(compare_mutations(toupper(original_mut), result_mut))
})

test_that("seq_to_mut() double roundtrip consistency", {
  ref <- "MKTAYIAK"
  original_mut <- "M1A / K2R / -1HHH / K8-"
  
  # Purpose: applying seq_to_mut -> mut_to_seq -> seq_to_mut should be stable
  mutant1 <- mut_to_seq(ref, original_mut)
  detected_mut1 <- seq_to_mut(ref, mutant1)
  mutant2 <- mut_to_seq(ref, detected_mut1)
  detected_mut2 <- seq_to_mut(ref, mutant2)
  
  expect_equal(mutant1, mutant2)
  expect_true(compare_mutations(detected_mut1, detected_mut2))
})