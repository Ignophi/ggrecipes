# tests/testthat/test-gg_kdmap.R

test_that("gg_kdmap handles a very tight KD range (similar affinities)", {
  # Data: all KD values ~100 nM
  test_tight <- data.frame(
    id = paste0("Var", 1:5),
    ka = c(1.2e5, 1.3e5, 1.1e5, 1.25e5, 1.15e5),
    kd = c(0.012, 0.013, 0.011, 0.0125, 0.0115)
  )

  p <- gg_kdmap(test_tight, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure tight-range data plots successfully as a ggplot object
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap handles a very wide KD range (pM to uM)", {
  # Data: KD spans orders of magnitude from 0.1 nM to 10,000 nM
  test_wide <- data.frame(
    id = paste0("Var", 1:6),
    ka = c(1e7, 1e6, 1e5, 1e4, 1e5, 1e6),
    kd = c(0.001, 0.01, 0.1, 1, 10, 100)
  )

  p <- gg_kdmap(test_wide, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure wide-range data (several orders of magnitude) is handled
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap handles ultra-high affinity binders (sub-nM)", {
  # Data: KD in the picomolar / sub-nM range
  test_high_aff <- data.frame(
    id = paste0("Binder", 1:4),
    ka = c(5e6, 8e6, 1e7, 6e6),
    kd = c(0.001, 0.0005, 0.0008, 0.0012)
  )

  p <- gg_kdmap(test_high_aff, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure very small KD values plot correctly
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap handles low-affinity binders (high nM to uM)", {
  # Data: KD in the ~500–2000 nM range
  test_low_aff <- data.frame(
    id = paste0("Weak", 1:4),
    ka = c(1e4, 5e3, 8e3, 1.2e4),
    kd = c(5, 10, 8, 6)
  )

  p <- gg_kdmap(test_low_aff, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure low-affinity interactions (large KD) are plotted
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap works for a single data point", {
  # Data: single binding variant
  test_single <- data.frame(
    id = "Solo",
    ka = 1e5,
    kd = 0.01
  )

  p <- gg_kdmap(test_single, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure function gracefully handles a single row input
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap works for two points in the same order of magnitude", {
  # Data: two similar KD values
  test_two <- data.frame(
    id = c("A", "B"),
    ka = c(1e5, 2e5),
    kd = c(0.01, 0.02)
  )

  p <- gg_kdmap(test_two, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure small datasets with similar KD are supported
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap handles extremely fast on-rates", {
  # Data: very fast association rate constants
  test_extreme_on <- data.frame(
    id = paste0("Fast", 1:3),
    ka = c(1e8, 5e7, 8e7),
    kd = c(0.1, 0.05, 0.08)
  )

  p <- gg_kdmap(test_extreme_on, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure large ka values do not break the plot
  expect_s3_class(p, "ggplot")
})

test_that("gg_kdmap handles extremely fast off-rates", {
  # Data: very fast dissociation rate constants
  test_extreme_off <- data.frame(
    id = paste0("Unstable", 1:3),
    ka = c(1e4, 5e4, 2e4),
    kd = c(50, 100, 75)
  )

  p <- gg_kdmap(test_extreme_off, iso_n = 15, show_anno = TRUE)

  # Purpose: ensure very large kd values do not break the plot
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Additional feature tests
# --------------------------------------------------------------------

test_that("gg_kdmap respects show_anno argument", {
  # Reuse a moderate dataset
  test_wide <- data.frame(
    id = paste0("Var", 1:6),
    ka = c(1e7, 1e6, 1e5, 1e4, 1e5, 1e6),
    kd = c(0.001, 0.01, 0.1, 1, 10, 100)
  )

  p_anno    <- gg_kdmap(test_wide, iso_n = 15, show_anno = TRUE)
  p_no_anno <- gg_kdmap(test_wide, iso_n = 15, show_anno = FALSE)

  # Purpose: both variants should produce valid ggplot objects
  expect_s3_class(p_anno, "ggplot")
  expect_s3_class(p_no_anno, "ggplot")
})

test_that("gg_kdmap accepts different iso_n values without error", {
  # Reuse a representative dataset
  test_tight <- data.frame(
    id = paste0("Var", 1:5),
    ka = c(1.2e5, 1.3e5, 1.1e5, 1.25e5, 1.15e5),
    kd = c(0.012, 0.013, 0.011, 0.0125, 0.0115)
  )

  # Purpose: ensure iso_n parameter is flexible and does not cause failures
  expect_s3_class(gg_kdmap(test_tight, iso_n = 5,  show_anno = TRUE), "ggplot")
  expect_s3_class(gg_kdmap(test_tight, iso_n = 10, show_anno = TRUE), "ggplot")
  expect_s3_class(gg_kdmap(test_tight, iso_n = 20, show_anno = TRUE), "ggplot")
})

test_that("gg_kdmap errors when required columns are missing", {
  # Start from a valid dataset
  valid <- data.frame(
    id = c("A", "B"),
    ka = c(1e5, 2e5),
    kd = c(0.01, 0.02)
  )

  # Purpose: missing kd column should produce an error
  expect_error(
    gg_kdmap(valid[, c("id", "ka")]),
    regexp = ""
  )

  # Purpose: missing ka column should also produce an error
  expect_error(
    gg_kdmap(valid[, c("id", "kd")]),
    regexp = ""
  )
})

test_that("gg_kdmap works when id is factor or character", {
  base_df <- data.frame(
    id = c("A", "B", "C"),
    ka = c(1e5, 2e5, 1.5e5),
    kd = c(0.01, 0.02, 0.015)
  )

  df_char <- base_df
  df_factor <- base_df
  df_factor$id <- factor(df_factor$id)

  # Purpose: support both character and factor id types
  p_char   <- gg_kdmap(df_char,   iso_n = 15, show_anno = TRUE)
  p_factor <- gg_kdmap(df_factor, iso_n = 15, show_anno = TRUE)

  expect_s3_class(p_char, "ggplot")
  expect_s3_class(p_factor, "ggplot")
})
