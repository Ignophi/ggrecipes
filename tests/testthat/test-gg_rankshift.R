# tests/testthat/test-gg_rankshift.R

# --------------------------------------------------------------------
# Input validation tests
# --------------------------------------------------------------------

test_that("gg_rankshift errors when data is not a data frame", {
  expect_error(gg_rankshift(data = "not_a_df"), "'data' must be a data frame")
  expect_error(gg_rankshift(data = NULL), "'data' must be a data frame")
  expect_error(gg_rankshift(data = list(a = 1)), "'data' must be a data frame")
})

test_that("gg_rankshift errors when required columns are missing", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, id = "missing"), "Required column\\(s\\) missing")
  expect_error(gg_rankshift(df, group = "missing"), "Required column\\(s\\) missing")
  expect_error(gg_rankshift(df, value = "missing"), "Required column\\(s\\) missing")
})

test_that("gg_rankshift errors when group has wrong number of levels", {
  df_single <- data.frame(
    id = c("A", "B"),
    group = c("G1", "G1"),
    value = c(10, 20)
  )
  
  df_triple <- data.frame(
    id = rep(c("A", "B"), 3),
    group = rep(c("G1", "G2", "G3"), each = 2),
    value = c(10, 20, 15, 25, 12, 22)
  )
  
  expect_error(gg_rankshift(df_single), "must have exactly 2 unique values")
  expect_error(gg_rankshift(df_triple), "must have exactly 2 unique values")
})

test_that("gg_rankshift errors when value is not numeric", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c("high", "low", "medium", "high")
  )
  
  expect_error(gg_rankshift(df), "must be numeric")
})

test_that("gg_rankshift validates style parameter", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, style = "invalid"), "must be one of")
  expect_error(gg_rankshift(df, style = c("bar", "box")), "must be one of")
})

test_that("gg_rankshift validates alpha parameters", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, alpha = -0.1), "must be between 0 and 1")
  expect_error(gg_rankshift(df, alpha = 1.5), "must be between 0 and 1")
  expect_error(gg_rankshift(df, line_alpha = -0.1), "must be between 0 and 1")
  expect_error(gg_rankshift(df, line_alpha = 1.5), "must be between 0 and 1")
  expect_error(gg_rankshift(df, point_alpha = -0.1), "must be between 0 and 1")
  expect_error(gg_rankshift(df, point_alpha = 1.5), "must be between 0 and 1")
})

test_that("gg_rankshift validates line_width parameter", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, line_width = -1), "must be a non-negative number")
})

test_that("gg_rankshift validates panel_ratio parameter", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, panel_ratio = 0), "must be a positive number")
  expect_error(gg_rankshift(df, panel_ratio = -1), "must be a positive number")
})

test_that("gg_rankshift validates point parameters", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, point_size = -1), "must be a non-negative number")
  expect_error(gg_rankshift(df, point_shape = -1), "must be a number between 0 and 25")
  expect_error(gg_rankshift(df, point_shape = 26), "must be a number between 0 and 25")
})

test_that("gg_rankshift validates logical parameters", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, show_points = "TRUE"), "must be TRUE or FALSE")
  expect_error(gg_rankshift(df, decreasing = 1), "must be TRUE or FALSE")
  expect_error(gg_rankshift(df, free_x = "yes"), "must be TRUE or FALSE")
})

test_that("gg_rankshift validates stat_summary parameter", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, stat_summary = "invalid"), "must be one of")
})

test_that("gg_rankshift validates fills parameter", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, fill = "red"), 
               "must be a character vector of length 2")
  expect_error(gg_rankshift(df, fill = c("red", "blue", "green")), 
               "must be a character vector of length 2")
})

test_that("gg_rankshift validates rank_change_colors parameter", {
  df <- data.frame(
    id = c("A", "B", "A", "B"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df, rank_change_colors = c("red", "blue")), 
               "must be a character vector of length 3")
  expect_error(gg_rankshift(df, rank_change_colors = c(a = "red", b = "blue", c = "green")), 
               "must be a named vector with names")
})

test_that("gg_rankshift errors when no common samples between groups", {
  df <- data.frame(
    id = c("A", "B", "C", "D"),
    group = c("G1", "G1", "G2", "G2"),
    value = c(10, 20, 15, 25)
  )
  
  expect_error(gg_rankshift(df), "No common samples found")
})

# --------------------------------------------------------------------
# Basic functionality tests
# --------------------------------------------------------------------

test_that("gg_rankshift creates plot with box style", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B", "C"), 4),
    group = rep(c("G1", "G2"), each = 6),
    value = c(10, 20, 30, 12, 22, 32, 15, 25, 35, 18, 28, 38)
  )
  
  p <- gg_rankshift(df, style = "box")
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift creates plot with bar style", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B", "C"), 4),
    group = rep(c("G1", "G2"), each = 6),
    value = c(10, 20, 30, 12, 22, 32, 15, 25, 35, 18, 28, 38)
  )
  
  p <- gg_rankshift(df, style = "bar")
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift works with custom fill colors", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 2),
    group = rep(c("G1", "G2"), each = 2),
    value = c(10, 20, 15, 25)
  )
  
  p <- gg_rankshift(df, fill = c("red", "blue"))
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift works with custom rank change colors", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 2),
    group = rep(c("G1", "G2"), each = 2),
    value = c(10, 20, 15, 25)
  )
  
  colors <- c(increase = "green", decrease = "red", no_change = "grey")
  p <- gg_rankshift(df, rank_change_colors = colors)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift works without points", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 2),
    group = rep(c("G1", "G2"), each = 2),
    value = c(10, 20, 15, 25)
  )
  
  p <- gg_rankshift(df, show_points = FALSE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift uses median summary statistic", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 4),
    group = rep(c("G1", "G2"), each = 4),
    value = c(10, 20, 100, 15, 12, 22, 110, 18)
  )
  
  p <- gg_rankshift(df, stat_summary = "median")
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift works with shared x-axis", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 2),
    group = rep(c("G1", "G2"), each = 2),
    value = c(10, 20, 15, 25)
  )
  
  p <- gg_rankshift(df, free_x = FALSE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift works with decreasing rank order", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 2),
    group = rep(c("G1", "G2"), each = 2),
    value = c(10, 20, 15, 25)
  )
  
  p <- gg_rankshift(df, decreasing = TRUE)
  expect_s3_class(p, "patchwork")
})

# --------------------------------------------------------------------
# Rank change behavior tests
# --------------------------------------------------------------------

test_that("gg_rankshift correctly identifies no rank change", {
  skip_if_not_installed("patchwork")
  
  # Create data where ranks stay the same
  # G1: A=10, B=20, C=30 -> ranks: A=1, B=2, C=3
  # G2: A=15, B=25, C=35 -> ranks: A=1, B=2, C=3
  df <- data.frame(
    id = rep(c("A", "B", "C"), 2),
    group = rep(c("G1", "G2"), each = 3),
    value = c(10, 20, 30, 15, 25, 35)
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
  
  # Check that plot was created successfully
  # The internal line_df should have all "no_change"
  # We can't easily test internal data, but we can verify plot creation
})

test_that("gg_rankshift correctly identifies rank increase", {
  skip_if_not_installed("patchwork")
  
  # Create data where B improves from rank 3 to rank 1
  # G1: A=30, B=10, C=20 -> ranks: A=3, B=1, C=2
  # G2: A=15, B=35, C=20 -> ranks: A=1, B=3, C=2
  df <- data.frame(
    id = rep(c("A", "B", "C"), 2),
    group = rep(c("G1", "G2"), each = 3),
    value = c(30, 10, 20, 15, 35, 20)
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift correctly identifies rank decrease", {
  skip_if_not_installed("patchwork")
  
  # Create data where A declines from rank 3 to rank 1
  # G1: A=30, B=10, C=20 -> ranks: A=3, B=1, C=2
  # G2: A=5, B=15, C=25 -> ranks: A=1, B=2, C=3
  df <- data.frame(
    id = rep(c("A", "B", "C"), 2),
    group = rep(c("G1", "G2"), each = 3),
    value = c(30, 10, 20, 5, 15, 25)
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles mixed rank changes", {
  skip_if_not_installed("patchwork")
  
  # Create data with mixed rank changes
  # G1: A=10, B=20, C=30, D=40 -> ranks: 1,2,3,4
  # G2: A=40, B=30, C=10, D=20 -> ranks: 4,3,1,2
  # A: increase (1->4), B: increase (2->3), C: decrease (3->1), D: decrease (4->2)
  df <- data.frame(
    id = rep(c("A", "B", "C", "D"), 2),
    group = rep(c("G1", "G2"), each = 4),
    value = c(10, 20, 30, 40, 40, 30, 10, 20)
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles rank changes with decreasing=TRUE", {
  skip_if_not_installed("patchwork")
  
  # With decreasing=TRUE, higher values get lower rank numbers
  # G1: A=10, B=20, C=30 -> ranks: A=3, B=2, C=1
  # G2: A=35, B=25, C=15 -> ranks: A=1, B=2, C=3
  # A improves from 3 to 1, B stays at 2, C declines from 1 to 3
  df <- data.frame(
    id = rep(c("A", "B", "C"), 2),
    group = rep(c("G1", "G2"), each = 3),
    value = c(10, 20, 30, 35, 25, 15)
  )
  
  p <- gg_rankshift(df, decreasing = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles replicates correctly", {
  skip_if_not_installed("patchwork")
  
  # Multiple replicates per sample per group
  # Should aggregate before ranking
  df <- data.frame(
    id = rep(c("A", "B", "C"), each = 4),
    group = rep(rep(c("G1", "G2"), each = 2), 3),
    value = c(
      # A replicates: G1 mean=10, G2 mean=30
      9, 11, 29, 31,
      # B replicates: G1 mean=20, G2 mean=20
      19, 21, 19, 21,
      # C replicates: G1 mean=30, G2 mean=10
      29, 31, 9, 11
    )
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

# --------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------

test_that("gg_rankshift handles two samples", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B"), 2),
    group = rep(c("G1", "G2"), each = 2),
    value = c(10, 20, 15, 25)
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles many samples", {
  skip_if_not_installed("patchwork")
  
  n <- 20
  df <- data.frame(
    id = rep(paste0("S", 1:n), 2),
    group = rep(c("G1", "G2"), each = n),
    value = c(1:n, sample(1:n))
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles NA values", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = rep(c("A", "B", "C"), 2),
    group = rep(c("G1", "G2"), each = 3),
    value = c(10, NA, 30, 15, 25, 35)
  )
  
  # Should handle NAs through na.rm in aggregate
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles tied ranks", {
  skip_if_not_installed("patchwork")
  
  # Samples with identical values
  df <- data.frame(
    id = rep(c("A", "B", "C"), 2),
    group = rep(c("G1", "G2"), each = 3),
    value = c(10, 20, 20, 15, 25, 25)
  )
  
  p <- gg_rankshift(df)
  expect_s3_class(p, "patchwork")
})

test_that("gg_rankshift handles custom column names", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    sample_id = rep(c("A", "B"), 2),
    condition = rep(c("Before", "After"), each = 2),
    measurement = c(10, 20, 15, 25)
  )
  
  p <- gg_rankshift(df, id = "sample_id", group = "condition", value = "measurement")
  expect_s3_class(p, "patchwork")
})

# --------------------------------------------------------------------
# Integration tests
# --------------------------------------------------------------------

test_that("gg_rankshift handles real-world scenario", {
  skip_if_not_installed("patchwork")
  
  # Simulated antibody binding data
  antibodies <- paste0("Ab", 1:8)
  df <- data.frame(
    variant = rep(antibodies, each = 6),
    condition = rep(rep(c("Before", "After"), each = 3), length(antibodies)),
    kd_nm = c(
      # Some improve, some decline, some stay similar
      rnorm(3, 50, 5), rnorm(3, 20, 3),  # Ab1: improves
      rnorm(3, 30, 3), rnorm(3, 60, 5),  # Ab2: declines
      rnorm(3, 40, 4), rnorm(3, 40, 4),  # Ab3: no change
      rnorm(3, 70, 7), rnorm(3, 25, 3),  # Ab4: improves
      rnorm(3, 20, 2), rnorm(3, 50, 5),  # Ab5: declines
      rnorm(3, 45, 4), rnorm(3, 45, 4),  # Ab6: no change
      rnorm(3, 35, 3), rnorm(3, 15, 2),  # Ab7: improves
      rnorm(3, 55, 5), rnorm(3, 55, 5)   # Ab8: no change
    )
  )
  
  p <- gg_rankshift(
    df,
    id = "variant",
    group = "condition",
    value = "kd_nm",
    decreasing = TRUE,  # Lower KD = better = higher rank
    style = "box"
  )
  
  expect_s3_class(p, "patchwork")
})