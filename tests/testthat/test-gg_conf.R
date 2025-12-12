# tests/testthat/test-gg_conf.R

# --------------------------------------------------------------------
# Input validation tests
# --------------------------------------------------------------------

test_that("gg_conf errors when data is not a data frame", {
  # Purpose: ensure non-data.frame inputs are rejected
  expect_error(gg_conf(data = "not_a_df", x = "a", y = "b"), "'data' must be a data frame")
  expect_error(gg_conf(data = NULL, x = "a", y = "b"), "'data' must be a data frame")
  expect_error(gg_conf(data = list(a = 1), x = "a", y = "b"), "'data' must be a data frame")
})

test_that("gg_conf errors when x or y are not single character strings", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure x and y must be single character strings
  expect_error(gg_conf(df, x = c("cat1", "cat2"), y = "cat2"), "must be a single character string")
  expect_error(gg_conf(df, x = "cat1", y = c("cat1", "cat2")), "must be a single character string")
  expect_error(gg_conf(df, x = 123, y = "cat2"), "must be a single character string")
  expect_error(gg_conf(df, x = "cat1", y = 456), "must be a single character string")
})

test_that("gg_conf errors when required columns are missing", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: verify column existence checks work
  expect_error(gg_conf(df, x = "missing_col", y = "cat2"), "Required column\\(s\\) missing")
  expect_error(gg_conf(df, x = "cat1", y = "missing_col"), "Required column\\(s\\) missing")
})

test_that("gg_conf errors when facet columns are missing", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10),
    group = rep(c("G1", "G2"), 5)
  )
  
  # Purpose: ensure facet column validation works
  expect_error(gg_conf(df, x = "cat1", y = "cat2", facet_x = "missing"), 
               "not found in data")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", facet_y = "missing"), 
               "not found in data")
})

test_that("gg_conf validates text_size parameter", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure text_size is positive numeric
  expect_error(gg_conf(df, x = "cat1", y = "cat2", text_size = -1), 
               "must be a positive number")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", text_size = 0), 
               "must be a positive number")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", text_size = "large"), 
               "must be a positive number")
})

test_that("gg_conf validates point_size_range parameter", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure point_size_range is valid numeric vector of length 2
  expect_error(gg_conf(df, x = "cat1", y = "cat2", point_size_range = 5), 
               "must be a numeric vector of length 2")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", point_size_range = c(5, 10, 15)), 
               "must be a numeric vector of length 2")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", point_size_range = c(-1, 10)), 
               "must have positive values")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", point_size_range = c(10, 5)), 
               "must have positive values with \\[2\\] > \\[1\\]")
})

test_that("gg_conf validates show_grid parameter", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure show_grid is logical
  expect_error(gg_conf(df, x = "cat1", y = "cat2", show_grid = "TRUE"), 
               "must be TRUE or FALSE")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", show_grid = 1), 
               "must be TRUE or FALSE")
})

test_that("gg_conf validates expand parameter", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure expand is non-negative numeric
  expect_error(gg_conf(df, x = "cat1", y = "cat2", expand = -0.1), 
               "must be a non-negative number")
  expect_error(gg_conf(df, x = "cat1", y = "cat2", expand = "small"), 
               "must be a non-negative number")
})

# --------------------------------------------------------------------
# Basic functionality tests
# --------------------------------------------------------------------

test_that("gg_conf creates plot from simple categorical data", {
  df <- data.frame(
    cat1 = rep(c("A", "B", "C"), each = 4),
    cat2 = rep(c("X", "Y"), 6)
  )
  
  # Purpose: ensure basic plotting works
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% class(p$layers[[1]]$geom))
})

test_that("gg_conf handles factor variables", {
  df <- data.frame(
    cat1 = factor(rep(c("A", "B", "C"), each = 4)),
    cat2 = factor(rep(c("X", "Y"), 6))
  )
  
  # Purpose: ensure factor variables are handled correctly
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf handles character variables", {
  df <- data.frame(
    cat1 = rep(c("A", "B", "C"), each = 4),
    cat2 = rep(c("X", "Y"), 6),
    stringsAsFactors = FALSE
  )
  
  # Purpose: ensure character variables are handled correctly
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf applies custom fill color", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure custom fill color is applied
  p <- gg_conf(df, x = "cat1", y = "cat2", fill = "lightcoral")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf applies custom text size", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure custom text size is applied
  p <- gg_conf(df, x = "cat1", y = "cat2", text_size = 6)
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf applies custom point size range", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure custom point size range is applied
  p <- gg_conf(df, x = "cat1", y = "cat2", point_size_range = c(5, 20))
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf hides grid lines when show_grid = FALSE", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure grid lines can be hidden
  p <- gg_conf(df, x = "cat1", y = "cat2", show_grid = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf shows grid lines when show_grid = TRUE", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure grid lines are shown by default
  p <- gg_conf(df, x = "cat1", y = "cat2", show_grid = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf applies custom expand value", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 5),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure custom expand value is applied
  p <- gg_conf(df, x = "cat1", y = "cat2", expand = 0.3)
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf applies horizontal faceting", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 10),
    cat2 = rep(c("X", "Y"), 10),
    group = rep(c("G1", "G2"), each = 10)
  )
  
  # Purpose: ensure horizontal faceting works
  p <- gg_conf(df, x = "cat1", y = "cat2", facet_x = "group")
  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("gg_conf applies vertical faceting", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 10),
    cat2 = rep(c("X", "Y"), 10),
    group = rep(c("G1", "G2"), each = 10)
  )
  
  # Purpose: ensure vertical faceting works
  p <- gg_conf(df, x = "cat1", y = "cat2", facet_y = "group")
  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("gg_conf applies grid faceting", {
  df <- data.frame(
    cat1 = rep(c("A", "B"), 20),
    cat2 = rep(c("X", "Y"), 20),
    group_x = rep(c("G1", "G2"), each = 20),
    group_y = rep(c("H1", "H2"), 20)
  )
  
  # Purpose: ensure grid faceting works
  p <- gg_conf(df, x = "cat1", y = "cat2", facet_x = "group_x", facet_y = "group_y")
  expect_s3_class(p, "ggplot")
  expect_true("FacetGrid" %in% class(p$facet))
})

# --------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------

test_that("gg_conf handles single category in x", {
  df <- data.frame(
    cat1 = rep("A", 10),
    cat2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure single-category data is handled
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf handles single category in y", {
  df <- data.frame(
    cat1 = rep(c("A", "B", "C"), length.out = 10),
    cat2 = rep("X", 10)
  )
  
  # Purpose: ensure single-category data is handled
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf handles unbalanced categories", {
  df <- data.frame(
    cat1 = c(rep("A", 8), rep("B", 2)),
    cat2 = rep(c("X", "Y"), 5)
  )
  
  # Purpose: ensure unbalanced counts are handled
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf handles many categories", {
  df <- data.frame(
    cat1 = rep(LETTERS[1:10], each = 10),
    cat2 = rep(letters[1:10], 10)
  )
  
  # Purpose: ensure large number of categories is handled
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf handles NA values", {
  df <- data.frame(
    cat1 = c(rep("A", 5), rep("B", 5), NA, NA),
    cat2 = c(rep("X", 5), rep("Y", 5), "Z", "Z")
  )
  
  # Purpose: ensure NA values don't break the plot
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf handles minimum dataset (2x2)", {
  df <- data.frame(
    cat1 = c("A", "B"),
    cat2 = c("X", "Y")
  )
  
  # Purpose: ensure minimum viable data creates plot
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Integration tests
# --------------------------------------------------------------------

test_that("gg_conf works with all features combined", {
  df <- data.frame(
    cat1 = rep(c("A", "B", "C"), 20),
    cat2 = rep(c("X", "Y"), 30),
    group_x = rep(c("G1", "G2"), each = 30),
    group_y = rep(c("H1", "H2"), 30)
  )
  
  # Purpose: ensure all parameters work together
  p <- gg_conf(
    df,
    x = "cat1",
    y = "cat2",
    fill = "lightgreen",
    text_size = 5,
    text_color = "darkblue",
    point_size_range = c(4, 18),
    border_color = "darkgreen",
    show_grid = FALSE,
    expand = 0.2,
    facet_x = "group_x",
    facet_y = "group_y"
  )
  expect_s3_class(p, "ggplot")
})

test_that("gg_conf produces correct axis labels", {
  df <- data.frame(
    category_1 = rep(c("A", "B"), 5),
    category_2 = rep(c("X", "Y", "Z"), length.out = 10)
  )
  
  # Purpose: ensure axis labels match column names
  p <- gg_conf(df, x = "category_1", y = "category_2")
  expect_equal(p$labels$x, "category_1")
  expect_equal(p$labels$y, "category_2")
})

test_that("gg_conf works with numeric categories", {
  df <- data.frame(
    cat1 = rep(c(1, 2, 3), length.out = 12),
    cat2 = rep(c(10, 20), 6)
  )
  
  # Purpose: ensure numeric categories are treated as categorical
  p <- gg_conf(df, x = "cat1", y = "cat2")
  expect_s3_class(p, "ggplot")
})