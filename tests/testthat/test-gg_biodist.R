# tests/testthat/test-gg_biodist.R

# --------------------------------------------------------------------
# Input validation tests
# --------------------------------------------------------------------

test_that("gg_biodist errors when data is not a data frame", {
  # Purpose: ensure non-data.frame inputs are rejected
  expect_error(gg_biodist(data = "not_a_df"), "'data' must be a data frame")
  expect_error(gg_biodist(data = NULL), "'data' must be a data frame")
  expect_error(gg_biodist(data = list(a = 1)), "'data' must be a data frame")
})

test_that("gg_biodist errors when required columns are missing", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: verify column existence checks work
  expect_error(gg_biodist(df, id = "missing_col"), "Required column\\(s\\) missing")
  expect_error(gg_biodist(df, value = "missing_col"), "No columns in .* match the pattern")
})

test_that("gg_biodist coerces non-numeric value column", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c("2", "2.5", "1.8", "15", "14", "16")
  )
  
  # Run quietly
  p <- gg_biodist(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
  
  df$value <- letters[1:6]
  expect_error(gg_biodist(df, quiet = TRUE), "all values are NA after as.numeric")
})

test_that("gg_biodist validates group parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16),
    group = rep(c("A", "B", "A"), 2)
  )

  # Purpose: ensure group must be single character string
  expect_error(gg_biodist(df, group = c("a", "b")), "must be a single character string")
  expect_error(gg_biodist(df, group = 123), "must be a single character string")
  expect_error(gg_biodist(df, group = "nonexistent"), "not found in data")
})

test_that("gg_biodist validates separate parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure separate must be character vector with valid organ names
  expect_error(gg_biodist(df, separate = 123), "must be a character vector")
  expect_error(gg_biodist(df, separate = "Nonexistent"), "not found in data")
})

test_that("gg_biodist validates stat_summary parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure only valid summary stats are accepted
  expect_error(gg_biodist(df, stat_summary = "invalid"), "must be one of")
})

test_that("gg_biodist validates bar_alpha parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure alpha is between 0 and 1
  expect_error(gg_biodist(df, bar_alpha = -0.1), "must be between 0 and 1")
  expect_error(gg_biodist(df, bar_alpha = 1.5), "must be between 0 and 1")
  expect_error(gg_biodist(df, bar_alpha = "0.5"), "must be between 0 and 1")
})

test_that("gg_biodist validates point_size parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure point size is non-negative
  expect_error(gg_biodist(df, point_size = -1), "must be a non-negative number")
  expect_error(gg_biodist(df, point_size = "large"), "must be a non-negative number")
})

test_that("gg_biodist validates error_bars parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure error_bars is logical
  expect_error(gg_biodist(df, error_bars = "TRUE"), "must be TRUE or FALSE")
  expect_error(gg_biodist(df, error_bars = 1), "must be TRUE or FALSE")
})

test_that("gg_biodist validates fill_colors parameter", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16),
    group = rep(c("A", "B", "A"), 2)
  )

  # Purpose: ensure fill_colors has enough colors for groups
  expect_error(gg_biodist(df, group = "group", fill_colors = "red"),
               "must have at least as many colors")
})

# --------------------------------------------------------------------
# Wide to long conversion tests
# --------------------------------------------------------------------

test_that("gg_biodist converts wide data using regex pattern", {
  df <- data.frame(
    sample = paste0("S", 1:3),
    Blood_val = c(2, 2.5, 1.8),
    Liver_val = c(15, 14, 16),
    Kidney_val = c(8, 7.5, 9)
  )

  # Purpose: ensure wide format is properly converted to long
  p <- gg_biodist(df, value = "_val$")
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist errors when no columns match wide pattern", {
  df <- data.frame(
    sample = paste0("S", 1:3),
    Blood_val = c(2, 2.5, 1.8),
    Liver_val = c(15, 14, 16)
  )

  # Purpose: ensure meaningful error when pattern doesn't match
  expect_error(gg_biodist(df, value = "nomatch"), "No columns.*match the pattern")
})

test_that("gg_biodist preserves non-measurement columns in wide conversion", {
  df <- data.frame(
    sample = paste0("S", 1:3),
    group = c("Control", "Treated", "Control"),
    extra = 1:3,
    Blood_val = c(2, 2.5, 1.8),
    Liver_val = c(15, 14, 16)
  )

  # Purpose: ensure metadata columns are retained after reshaping
  p <- gg_biodist(df, value = "_val$", group = "group")
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Basic functionality tests
# --------------------------------------------------------------------

test_that("gg_biodist creates plot from long format data", {
  df <- data.frame(
    id = rep(c("Blood", "Liver", "Kidney"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16, 8, 7.5, 9)
  )

  # Purpose: ensure basic long format plotting works
  p <- gg_biodist(df)
  expect_s3_class(p, "ggplot")
  expect_true("GeomBar" %in% class(p$layers[[1]]$geom))
})

test_that("gg_biodist creates plot from wide format data", {
  df <- data.frame(
    sample = paste0("S", 1:3),
    Blood_val = c(2, 2.5, 1.8),
    Liver_val = c(15, 14, 16)
  )

  # Purpose: ensure basic wide format plotting works
  p <- gg_biodist(df, value = "_val$")
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist handles grouped data", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 4),
    value = c(2, 2.5, 1.8, 2.2, 15, 14, 16, 15.5),
    group = rep(c("Control", "Treated"), 4)
  )

  # Purpose: ensure grouping creates dodged bars
  p <- gg_biodist(df, group = "group")
  expect_s3_class(p, "ggplot")
  expect_true("fill" %in% names(p$labels))
})

test_that("gg_biodist creates separate facets", {
  df <- data.frame(
    id = rep(c("Blood", "Liver", "Kidney"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16, 8, 7.5, 9)
  )

  # Purpose: ensure separate parameter creates faceted plot
  p <- gg_biodist(df, separate = "Blood")
  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("gg_biodist applies custom colors without grouping", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure single custom color is applied
  p <- gg_biodist(df, fill_colors = "purple")
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist applies custom colors with grouping", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 4),
    value = c(2, 2.5, 1.8, 2.2, 15, 14, 16, 15.5),
    group = rep(c("Control", "Treated"), 4)
  )

  # Purpose: ensure multiple custom colors map to groups
  p <- gg_biodist(df, group = "group", fill_colors = c("red", "blue"))
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist can disable points", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure points layer is not added when point_size = 0
  p <- gg_biodist(df, point_size = 0)
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 1)
})

test_that("gg_biodist adds error bars without grouping", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure error bars are added
  p <- gg_biodist(df, error_bars = TRUE)
  expect_s3_class(p, "ggplot")
  expect_true(any(sapply(p$layers, function(l) "GeomErrorbar" %in% class(l$geom))))
})

test_that("gg_biodist adds error bars with grouping", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 4),
    value = c(2, 2.5, 1.8, 2.2, 15, 14, 16, 15.5),
    group = rep(c("Control", "Treated"), 4)
  )

  # Purpose: ensure error bars work with grouped data
  p <- gg_biodist(df, group = "group", error_bars = TRUE)
  expect_s3_class(p, "ggplot")
  expect_true(any(sapply(p$layers, function(l) "GeomErrorbar" %in% class(l$geom))))
})

test_that("gg_biodist uses median summary statistic", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure median can be used instead of mean
  p <- gg_biodist(df, stat_summary = "median")
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist applies custom y-axis label", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure y_label parameter works
  p <- gg_biodist(df, y_label = "Custom Label")
  expect_equal(p$labels$y, "Custom Label")
})

# --------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------

test_that("gg_biodist handles NA values", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(NA, 2.5, 1.8, 15, 14, 16)
  )

  # Purpose: ensure NA values don't break plotting
  p <- gg_biodist(df)
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist handles single organ", {
  df <- data.frame(id = "Blood", value = 1:3)

  # Purpose: ensure single organ data is plotted
  p <- gg_biodist(df)
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist can separate all organs", {
  df <- data.frame(
    id = rep(c("Blood", "Liver", "Kidney"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16, 8, 7.5, 9)
  )
  organs <- unique(df$id)

  # Purpose: ensure all organs can be in separate facets
  p <- gg_biodist(df, separate = organs)
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist handles factor id column", {
  df <- data.frame(
    id = rep(c("Blood", "Liver"), each = 3),
    value = c(2, 2.5, 1.8, 15, 14, 16)
  )
  df$id <- factor(df$id)

  # Purpose: ensure factor id column is handled correctly
  p <- gg_biodist(df)
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Integration tests
# --------------------------------------------------------------------

test_that("gg_biodist works with all features combined", {
  df <- data.frame(
    id = rep(c("Blood", "Liver", "Kidney"), each = 4),
    value = c(2, 2.5, 1.8, 2.2, 15, 14, 16, 15.5, 8, 7.5, 9, 8.5),
    group = rep(c("Control", "Treated"), 6)
  )

  # Purpose: ensure all parameters work together
  p <- gg_biodist(
    df,
    group = "group",
    separate = "Blood",
    fill_colors = c("lightblue", "lightcoral"),
    bar_alpha = 0.5,
    point_size = 2,
    stat_summary = "median",
    error_bars = TRUE,
    y_label = "%ID/g"
  )
  expect_s3_class(p, "ggplot")
})

test_that("gg_biodist works with wide data and all features", {
  df <- data.frame(
    sample = paste0("S", 1:4),
    group = rep(c("Control", "Treated"), 2),
    Blood_val = c(2, 2.5, 1.8, 2.2),
    Liver_val = c(15, 14, 16, 15.5),
    Kidney_val = c(8, 7.5, 9, 8.5)
  )

  # Purpose: ensure wide format works with multiple features
  p <- gg_biodist(
    df,
    value = "_val$",
    group = "group",
    separate = "Blood",
    fill_colors = c("red", "blue"),
    error_bars = TRUE
  )
  expect_s3_class(p, "ggplot")
})