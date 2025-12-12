# tests/testthat/test-gg_splitcorr.R

# --------------------------------------------------------------------
# Input validation tests
# --------------------------------------------------------------------

test_that("gg_splitcorr errors when data is not a data frame", {
  # Purpose: ensure non-data.frame inputs are rejected
  expect_error(gg_splitcorr(data = "not_a_df", split = "x"), "'data' must be a data frame")
  expect_error(gg_splitcorr(data = NULL, split = "x"), "'data' must be a data frame")
  expect_error(gg_splitcorr(data = list(a = 1), split = "x"), "'data' must be a data frame")
})

test_that("gg_splitcorr errors when split column is missing", {
  df <- data.frame(x = 1:10, y = 1:10, z = 1:10)
  
  # Purpose: verify split column existence check works
  expect_error(gg_splitcorr(df, split = "missing"), "not found in data")
})

test_that("gg_splitcorr errors when split variable is not binary", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c("A", "B", "C"), length.out = 10)
  )
  
  # Purpose: ensure split must have exactly 2 unique values
  expect_error(gg_splitcorr(df, split = "group"), "must have exactly 2 unique values")
  
  df_single <- data.frame(x = 1:10, y = 1:10, group = "A")
  expect_error(gg_splitcorr(df_single, split = "group"), "must have exactly 2 unique values")
})

test_that("gg_splitcorr validates style parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure only valid styles are accepted
  expect_error(gg_splitcorr(df, split = "group", style = "invalid"), "must be one of")
})

test_that("gg_splitcorr validates method parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure only valid correlation methods are accepted
  expect_error(gg_splitcorr(df, split = "group", method = "invalid"), "must be one of")
})

test_that("gg_splitcorr validates linetype parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure only valid linetypes are accepted
  expect_error(gg_splitcorr(df, split = "group", linetype = "invalid"), "must be one of")
})

test_that("gg_splitcorr validates text_size parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure text_size is positive number
  expect_error(gg_splitcorr(df, split = "group", text_size = -1), "must be a positive number")
  expect_error(gg_splitcorr(df, split = "group", text_size = 0), "must be a positive number")
  expect_error(gg_splitcorr(df, split = "group", text_size = "large"), "must be a positive number")
})

test_that("gg_splitcorr validates linealpha parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure linealpha is between 0 and 1
  expect_error(gg_splitcorr(df, split = "group", linealpha = -0.1), "must be between 0 and 1")
  expect_error(gg_splitcorr(df, split = "group", linealpha = 1.5), "must be between 0 and 1")
})

test_that("gg_splitcorr validates colors parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure colors has exactly 3 elements
  expect_error(gg_splitcorr(df, split = "group", colors = c("red", "white")), 
               "must be a vector of length 3")
  expect_error(gg_splitcorr(df, split = "group", colors = c("red", "white", "blue", "green")), 
               "must be a vector of length 3")
})

test_that("gg_splitcorr validates text_colors parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure text_colors has exactly 2 elements
  expect_error(gg_splitcorr(df, split = "group", text_colors = "black"), 
               "must be a vector of length 2")
  expect_error(gg_splitcorr(df, split = "group", text_colors = c("black", "white", "red")), 
               "must be a vector of length 2")
})

test_that("gg_splitcorr validates offset parameter", {
  df <- data.frame(
    x = 1:10, 
    y = 1:10, 
    z = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure offset is non-negative
  expect_error(gg_splitcorr(df, split = "group", offset = -1), "must be a non-negative number")
})

test_that("gg_splitcorr errors with insufficient numeric columns", {
  df <- data.frame(
    x = 1:10,
    group = rep(c(0, 1), each = 5)
  )
  
  # Purpose: ensure at least 2 numeric columns after removing split
  expect_error(gg_splitcorr(df, split = "group"), "Need at least 2 numeric columns")
})

# --------------------------------------------------------------------
# Basic functionality tests
# --------------------------------------------------------------------

test_that("gg_splitcorr creates plot with tile style", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure basic tile plotting works
  p <- gg_splitcorr(df, split = "group", style = "tile")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
})

test_that("gg_splitcorr creates plot with point style", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure basic point plotting works
  p <- gg_splitcorr(df, split = "group", style = "point")
  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% class(p$layers[[1]]$geom))
})

test_that("gg_splitcorr uses pearson correlation", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure pearson method works
  p <- gg_splitcorr(df, split = "group", method = "pearson")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr uses spearman correlation", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure spearman method works
  p <- gg_splitcorr(df, split = "group", method = "spearman")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr applies custom colors", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure custom color palette is applied
  p <- gg_splitcorr(df, split = "group", colors = c("purple", "yellow", "green"))
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr applies custom prefix", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure custom prefix is used
  p <- gg_splitcorr(df, split = "group", prefix = "Custom: ")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr handles different linetype options", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure all valid linetypes work
  expect_s3_class(gg_splitcorr(df, split = "group", linetype = "solid"), "ggplot")
  expect_s3_class(gg_splitcorr(df, split = "group", linetype = "dashed"), "ggplot")
  expect_s3_class(gg_splitcorr(df, split = "group", linetype = "dotdash"), "ggplot")
})

# --------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------

test_that("gg_splitcorr handles NA values", {
  df <- data.frame(
    x = c(1:18, NA, NA),
    y = c(NA, 2:19, NA),
    z = 21:40,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure NA values don't break plotting
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr handles minimum number of variables", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure minimum viable input (2 numeric + split) works
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr handles many numeric columns", {
  df <- data.frame(
    v1 = 1:20,
    v2 = 21:40,
    v3 = 41:60,
    v4 = 61:80,
    v5 = 81:100,
    v6 = 101:120,
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure larger correlation matrices work
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr handles character split variable", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = rep(c("A", "B"), each = 10)
  )
  
  # Purpose: ensure character split variable works
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr handles factor split variable", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = factor(rep(c("Control", "Treated"), each = 10))
  )
  
  # Purpose: ensure factor split variable works
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr handles unbalanced groups", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    group = c(rep(0, 5), rep(1, 15))
  )
  
  # Purpose: ensure unequal group sizes don't break plotting
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

test_that("gg_splitcorr ignores non-numeric columns", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    char_col = letters[1:20],
    group = rep(c(0, 1), each = 10)
  )
  
  # Purpose: ensure non-numeric columns are automatically excluded
  p <- gg_splitcorr(df, split = "group")
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Integration tests
# --------------------------------------------------------------------

test_that("gg_splitcorr works with all parameters combined", {
  df <- data.frame(
    x = 1:20,
    y = 21:40,
    z = 41:60,
    w = 61:80,
    group = rep(c("Control", "Treated"), each = 10)
  )
  
  # Purpose: ensure all parameters work together
  p <- gg_splitcorr(
    df,
    split = "group",
    style = "point",
    method = "spearman",
    padjust = "bonferroni",
    colors = c("darkblue", "grey90", "darkred"),
    text_colors = c("yellow", "purple"),
    text_size = 4,
    border_color = "grey50",
    prefix = "Treatment: ",
    linetype = "solid",
    linealpha = 0.8,
    offset = 1
  )
  expect_s3_class(p, "ggplot")
})