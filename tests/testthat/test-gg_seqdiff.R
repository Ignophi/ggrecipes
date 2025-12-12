test_that("gg_seqdiff basic functionality works", {
  ref <- "ABCDEFGH"
  df <- data.frame(
    seq = c(ref, "ABXDEFGH", "ABCDEFXH"),
    id = c("ref", "mut1", "mut2")
  )
  
  p <- gg_seqdiff(data = df, ref = ref, sequence = "seq", name = "id")
  
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff handles sequences with no differences", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = c(ref, ref), id = c("s1", "s2"))
  
  p <- gg_seqdiff(data = df, ref = ref, sequence = "seq")
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff validates sequence length", {
  ref <- "ABCDEFGH"
  df <- data.frame(
    seq = c("ABCDEFGH", "ABC"),  # Second sequence too short
    id = c("s1", "s2")
  )
  
  expect_error(
    gg_seqdiff(data = df, ref = ref, sequence = "seq", name = "id"),
    "All sequences must be the same length"
  )
})

test_that("gg_seqdiff validates required parameters", {
  expect_error(gg_seqdiff(data = "not_df", ref = "ABC"), "'data' must be a data frame")
  expect_error(gg_seqdiff(data = data.frame(x = 1, sequence = 1), ref = NULL), "'ref' must be")
  expect_error(gg_seqdiff(data = data.frame(x = 1, sequence = 1), ref = ""), "'ref' must be")
})

test_that("gg_seqdiff validates column names", {
  df <- data.frame(other = "ABC")
  
  expect_error(
    gg_seqdiff(data = df, ref = "ABC", sequence = "missing"),
    "Column 'missing' not found"
  )
  
  expect_error(
    gg_seqdiff(data = df, ref = "ABC", sequence = "other", name = "missing"),
    "Column 'missing' not found"
  )
})

test_that("gg_seqdiff handles color parameter", {
  ref <- "ABCKEFGH"
  df <- data.frame(seq = c(ref, "ABXKEFGH"))
  
  p <- gg_seqdiff(
    data = df, 
    ref = ref, 
    sequence = "seq",
    color = c(K = "blue", X = "red")
  )
  
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff handles highlight parameter", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = c(ref, "ABXDEFGH"))
  
  p <- gg_seqdiff(
    data = df, 
    ref = ref, 
    sequence = "seq",
    highlight = list("yellow" = c(1, 2, 3))
  )
  
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff validates highlight colors", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = ref)
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = ref, 
      sequence = "seq",
      highlight = list("invalid_color" = c(1, 2))
    ),
    "invalid color"
  )
})

test_that("gg_seqdiff validates highlight format", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = ref)
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = ref, 
      sequence = "seq",
      highlight = c("blue" = 1)  # Not a list
    ),
    "'highlight' must be a named list"
  )
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = ref, 
      sequence = "seq",
      highlight = list("blue")  # Unnamed
    ),
    "'highlight' must be a named list"
  )
})

test_that("gg_seqdiff handles wrap parameter", {
  ref <- "ABCDEFGHIJKLMNOP"
  df <- data.frame(seq = c(ref, "ABCDEFGHXJKLMNOP"))
  
  p <- gg_seqdiff(data = df, ref = ref, sequence = "seq", wrap = 8)
  
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff validates wrap parameter", {
  df <- data.frame(seq = "ABC")
  
  expect_error(
    gg_seqdiff(data = df, ref = "ABC", sequence = "seq", wrap = -1),
    "'wrap' must be a positive number"
  )
  
  expect_error(
    gg_seqdiff(data = df, ref = "ABC", sequence = "seq", wrap = "not_numeric"),
    "'wrap' must be a positive number"
  )
})

test_that("gg_seqdiff validates size parameter", {
  df <- data.frame(seq = "ABC")
  
  expect_error(
    gg_seqdiff(data = df, ref = "ABC", sequence = "seq", size = -1),
    "'size' must be a positive number"
  )
  
  expect_error(
    gg_seqdiff(data = df, ref = "ABC", sequence = "seq", size = "large"),
    "'size' must be a positive number"
  )
})

test_that("gg_seqdiff handles annotate parameter", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = c(ref, "ABXDEFGH"))
  
  p <- gg_seqdiff(
    data = df, 
    ref = ref, 
    sequence = "seq",
    annotate = list(
      list(label = "Site 1", pos = 3),
      list(label = "Site 2", pos = 5, angle = 90)
    )
  )
  
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff validates annotate structure", {
  df <- data.frame(seq = "ABC")
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = "ABC", 
      sequence = "seq",
      annotate = "not_list"
    ),
    "'annotate' must be a list"
  )
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = "ABC", 
      sequence = "seq",
      annotate = list("not_inner_list")
    ),
    "Annotation .* must be a list"
  )
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = "ABC", 
      sequence = "seq",
      annotate = list(list(label = "test"))  # Missing pos
    ),
    "missing required element"
  )
})

test_that("gg_seqdiff validates annotate_defaults", {
  df <- data.frame(seq = "ABC")
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = "ABC", 
      sequence = "seq",
      annotate_defaults = "not_list"
    ),
    "'annotate_defaults' must be a list"
  )
  
  expect_error(
    gg_seqdiff(
      data = df, 
      ref = "ABC", 
      sequence = "seq",
      annotate_defaults = list(invalid_param = 1)
    ),
    "Invalid names in 'annotate_defaults'"
  )
})

test_that("gg_seqdiff handles empty sequences", {
  df <- data.frame(seq = c(NA, "", "ABC"))
  
  p <- gg_seqdiff(data = df, ref = "ABC", sequence = "seq")
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff warns when no valid sequences", {
  df <- data.frame(seq = c(NA, ""))
  
  expect_warning(
    gg_seqdiff(data = df, ref = "ABC", sequence = "seq"),
    "No valid sequences"
  )
})

test_that("gg_seqdiff handles show_ref parameter", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = c(ref, "ABXDEFGH"))
  
  p1 <- gg_seqdiff(data = df, ref = ref, sequence = "seq", show_ref = TRUE)
  p2 <- gg_seqdiff(data = df, ref = ref, sequence = "seq", show_ref = FALSE)
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("gg_seqdiff handles clustal input", {
  # Create temp clustal file
  clustal_content <- c(
    "CLUSTAL W",
    "",
    "seq1    ABCDEFGH",
    "seq2    ABXDEFGH",
    "        ** ****"
  )
  
  clustal_file <- tempfile(fileext = ".aln")
  writeLines(clustal_content, clustal_file)
  on.exit(unlink(clustal_file))
  
  p <- gg_seqdiff(clustal = clustal_file, ref = "ABCDEFGH")
  
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff validates clustal parameter", {
  expect_error(
    gg_seqdiff(clustal = c("file1", "file2"), ref = "ABC"),
    "'clustal' must be a single file path"
  )
  
  expect_error(
    gg_seqdiff(clustal = 123, ref = "ABC"),
    "'clustal' must be a single file path"
  )
})

test_that("gg_seqdiff handles gaps (x) in sequences", {
  ref <- "ABCDEFGH"
  df <- data.frame(seq = c("ABCDxFGH", "xBCDEFGH"))  # gaps shown as x
  
  p <- gg_seqdiff(data = df, ref = ref, sequence = "seq")
  expect_s3_class(p, "ggplot")
})

test_that("gg_seqdiff integration: complex scenario", {
  ref <- "ABCDEFGHIJKLMNOP"
  df <- data.frame(
    seq = c(
      ref,
      "ABXDEFGHIJKLMNOP",
      "ABCDEFGHxJKLMNOP"
    ),
    id = c("WT", "Mut1", "Mut2")
  )
  
  p <- gg_seqdiff(
    data = df,
    ref = ref,
    sequence = "seq",
    name = "id",
    color = c(X = "red", K = "blue"),
    highlight = list("#FFE0B2" = c(5:8), "lightblue" = c(12:15)),
    wrap = 8,
    annotate = list(
      list(label = "Region1", pos = 6),
      list(label = "Region2", pos = 13, angle = 90)
    ),
    show_ref = TRUE
  )
  
  expect_s3_class(p, "ggplot")
})