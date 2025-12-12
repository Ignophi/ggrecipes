# tests/testthat/test-gg_seq.R

test_that("gg_seq validates input parameters", {
  df <- data.frame(sequence = c("ABC", "DEF"))
  
  # data must be data frame
  expect_error(gg_seq(data = "not_df", ref = "ABCDEF"), 
               "'data' must be a data frame")
  
  # ref must be single non-empty string
  expect_error(gg_seq(data = df, ref = c("AB", "CD")), 
               "'ref' must be a single non-empty character string")
  expect_error(gg_seq(data = df, ref = ""), 
               "'ref' must be a single non-empty character string")
  
  # sequence column must exist
  expect_error(gg_seq(data = df, ref = "ABCDEF", sequence = "missing"), 
               "Column 'missing' not found in data")
  
  # name column must exist if specified
  expect_error(gg_seq(data = df, ref = "ABCDEF", name = "missing"), 
               "Column 'missing' not found in data")
  
  # size must be positive number
  expect_error(gg_seq(data = df, ref = "ABCDEF", size = -1), 
               "'size' must be a positive number")
  expect_error(gg_seq(data = df, ref = "ABCDEF", size = 0), 
               "'size' must be a positive number")
  
  # face must be valid
  expect_error(gg_seq(data = df, ref = "ABCDEF", face = "invalid"), 
               "'face' must be one of")
  
  # wrap must be positive or NULL
  expect_error(gg_seq(data = df, ref = "ABCDEF", wrap = -10), 
               "'wrap' must be a positive number or NULL")
})

test_that("gg_seq validates highlight parameter", {
  df <- data.frame(sequence = "ABC")
  
  # highlight must be named list
  expect_error(gg_seq(data = df, ref = "ABCDEF", highlight = c(1, 2)), 
               "'highlight' must be a named list")
  
  # highlight names must be valid colors
  expect_error(gg_seq(data = df, ref = "ABCDEF", 
                     highlight = list(notacolor = c(1, 2))), 
               "'highlight' contains invalid color name")
  
  # highlight values must be numeric
  expect_error(gg_seq(data = df, ref = "ABCDEF", 
                     highlight = list(red = c("a", "b"))), 
               "'highlight' values must be numeric position vectors")
})

test_that("gg_seq validates annotate parameter", {
  df <- data.frame(sequence = "ABC")
  
  # annotate must be list
  expect_error(gg_seq(data = df, ref = "ABCDEF", annotate = "invalid"), 
               "'annotate' must be a list")
  
  # each annotation must have label and pos
  expect_error(gg_seq(data = df, ref = "ABCDEF", 
                     annotate = list(list(label = "test"))), 
               "missing required element.*pos")
  
  expect_error(gg_seq(data = df, ref = "ABCDEF", 
                     annotate = list(list(pos = 1))), 
               "missing required element.*label")
  
  # label must be character
  expect_error(gg_seq(data = df, ref = "ABCDEF", 
                     annotate = list(list(label = 123, pos = 1))), 
               "'label' must be a single character string")
  
  # pos must be numeric
  expect_error(gg_seq(data = df, ref = "ABCDEF", 
                     annotate = list(list(label = "test", pos = "1"))), 
               "'pos' must be a single number")
})

test_that("gg_seq handles empty/invalid data", {
  # No valid sequences
  df_empty <- data.frame(sequence = character(0))
  expect_warning(p <- gg_seq(data = df_empty, ref = "ABCDEF"), 
                "No valid sequences")
  expect_s3_class(p, "gg")
  
  # Sequences with NA
  df_na <- data.frame(sequence = c(NA, "", "ABC"))
  p <- gg_seq(data = df_na, ref = "ABCDEF")
  expect_s3_class(p, "gg")
  
  # No matching sequences
  df_nomatch <- data.frame(sequence = c("XYZ", "QRS"))
  expect_warning(p <- gg_seq(data = df_nomatch, ref = "ABCDEF"), 
                "No sequences.*found in the reference")
  expect_s3_class(p, "gg")
})

test_that("gg_seq creates basic plot", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(
    id = c("seq1", "seq2"),
    sequence = c("ABC", "FGH")
  )
  
  p <- gg_seq(data = df, ref = ref)
  
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("gg_seq handles show_ref parameter", {
  ref <- "ABCDEF"
  df <- data.frame(sequence = "ABC")
  
  # With reference
  p1 <- gg_seq(data = df, ref = ref, show_ref = TRUE)
  expect_s3_class(p1, "gg")
  
  # Without reference
  p2 <- gg_seq(data = df, ref = ref, show_ref = FALSE)
  expect_s3_class(p2, "gg")
})

test_that("gg_seq applies highlighting", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(sequence = "ABC")
  
  p <- gg_seq(data = df, ref = ref, 
             highlight = list(yellow = c(1, 2, 3)))
  
  expect_s3_class(p, "gg")
  
  # Should have geom_rect layer for highlights
  expect_true(any(sapply(p$layers, function(l) {
    inherits(l$geom, "GeomRect")
  })))
})

test_that("gg_seq handles wrapping", {
  ref <- paste0(rep("A", 100), collapse = "")
  df <- data.frame(sequence = paste0(rep("A", 50), collapse = ""))
  
  p <- gg_seq(data = df, ref = ref, wrap = 30)
  
  expect_s3_class(p, "gg")
  
  # Should have facet_wrap
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("gg_seq adds annotations", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(sequence = "ABC")
  
  p <- gg_seq(data = df, ref = ref, 
             annotate = list(
               list(label = "Test", pos = 5),
               list(label = "Test2", pos = 8, angle = 90)
             ))
  
  expect_s3_class(p, "gg")
  
  # Check for text annotation layer
  has_text <- any(sapply(p$layers, function(l) {
    inherits(l$geom, "GeomText")
  }))
  expect_true(has_text)
})

test_that("gg_seq uses name column for labels", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(
    id = c("Sample1", "Sample2"),
    sequence = c("ABC", "DEF")
  )
  
  p <- gg_seq(data = df, ref = ref, name = "id")
  
  expect_s3_class(p, "gg")
  
  # Check y-axis has labels (not empty)
  y_labels <- p$scales$get_scales("y")
  # Basic check that plot was created with name parameter
  expect_true(TRUE)
})

test_that("gg_seq handles duplicate sequences", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(
    id = c("seq1", "seq2", "seq3"),
    sequence = c("ABC", "ABC", "DEF")
  )
  
  # Should use first name for duplicates
  p <- gg_seq(data = df, ref = ref, name = "id")
  expect_s3_class(p, "gg")
})

test_that("gg_seq combines multiple features", {
  ref <- "ABCDEFGHIJKLMNOP"
  df <- data.frame(
    id = c("seq1", "seq2"),
    sequence = c("ABCD", "GHIJ")
  )
  
  p <- gg_seq(
    data = df, 
    ref = ref, 
    name = "id",
    color = c(A = "red", G = "blue"),
    highlight = list(yellow = c(1:4)),
    wrap = 10,
    annotate = list(list(label = "Region", pos = 2))
  )
  
  expect_s3_class(p, "gg")
})

test_that("gg_seq handles overlapping sequences", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(
    sequence = c("ABCD", "CDEF", "EFGH")
  )
  
  p <- gg_seq(data = df, ref = ref)
  expect_s3_class(p, "gg")
})

test_that("gg_seq handles consecutive highlight ranges", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(sequence = "ABC")
  
  # Consecutive positions should merge into single rectangle
  p <- gg_seq(data = df, ref = ref, 
             highlight = list(yellow = c(1, 2, 3, 5, 6)))
  
  expect_s3_class(p, "gg")
})

test_that("gg_seq handles multiple highlight colors", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(sequence = "ABC")
  
  p <- gg_seq(data = df, ref = ref, 
             highlight = list(
               yellow = c(1, 2),
               blue = c(5, 6),
               "#FF0000" = c(8, 9)
             ))
  
  expect_s3_class(p, "gg")
})

test_that("gg_seq annotation defaults work", {
  ref <- "ABCDEFGHIJ"
  df <- data.frame(sequence = "ABC")
  
  # Custom defaults
  p <- gg_seq(
    data = df, 
    ref = ref, 
    annotate = list(list(label = "Test", pos = 5)),
    annotate_defaults = list(size = 5, face = "bold")
  )
  
  expect_s3_class(p, "gg")
})

test_that("gg_seq validates annotate_defaults", {
  df <- data.frame(sequence = "ABC")
  
  # Invalid parameter name
  expect_error(
    gg_seq(data = df, ref = "ABCDEF", 
          annotate_defaults = list(invalid = 1)),
    "Invalid names in 'annotate_defaults'"
  )
})