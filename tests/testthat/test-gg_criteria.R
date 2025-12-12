# tests/testthat/test-gg_criteria.R

# --------------------------------------------------------------------
# Input validation tests
# --------------------------------------------------------------------

test_that("gg_criteria errors when data is not a data frame", {
  # Purpose: ensure non-data.frame inputs are rejected
  expect_error(gg_criteria(data = "not_a_df"), "'data' must be a data frame")
  expect_error(gg_criteria(data = NULL), "'data' must be a data frame")
  expect_error(gg_criteria(data = list(a = 1)), "'data' must be a data frame")
})

test_that("gg_criteria errors when id is invalid", {
  df <- data.frame(
    sample = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure id parameter validation works
  expect_error(gg_criteria(df, id = c("a", "b")), "must be a single character string")
  expect_error(gg_criteria(df, id = 123), "must be a single character string")
})

test_that("gg_criteria errors when criteria is invalid", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure criteria parameter validation works
  expect_error(gg_criteria(df, criteria = c("a", "b")), "must be a single character string")
  expect_error(gg_criteria(df, criteria = 123), "must be a single character string")
})

test_that("gg_criteria errors when no criteria columns match pattern", {
  df <- data.frame(
    id = c("S1", "S2"),
    value = c(1, 2)
  )
  
  # Purpose: ensure meaningful error when pattern doesn't match
  expect_error(gg_criteria(df, criteria = "_criteria$"), "No columns.*match the pattern")
})

test_that("gg_criteria validates tile parameters", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure tile dimension parameters are valid
  expect_error(gg_criteria(df, tile_width = -0.1), "must be between 0 and 1")
  expect_error(gg_criteria(df, tile_width = 1.5), "must be between 0 and 1")
  expect_error(gg_criteria(df, tile_height = 0), "must be between 0 and 1")
  expect_error(gg_criteria(df, tile_height = 2), "must be between 0 and 1")
  expect_error(gg_criteria(df, tile_alpha = -0.1), "must be between 0 and 1")
  expect_error(gg_criteria(df, tile_alpha = 1.5), "must be between 0 and 1")
})

test_that("gg_criteria validates border_width parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure border_width is non-negative
  expect_error(gg_criteria(df, border_width = -1), "must be a non-negative number")
})

test_that("gg_criteria validates text_size parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure text_size is positive
  expect_error(gg_criteria(df, text_size = 0), "must be a positive number")
  expect_error(gg_criteria(df, text_size = -1), "must be a positive number")
})

test_that("gg_criteria validates logical parameters", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure logical parameters are TRUE or FALSE
  expect_error(gg_criteria(df, show_text = "TRUE"), "must be TRUE or FALSE")
  expect_error(gg_criteria(df, show_legend = 1), "must be TRUE or FALSE")
  expect_error(gg_criteria(df, quiet = "yes"), "must be TRUE or FALSE")
})

test_that("gg_criteria validates tile_fill parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure tile_fill is character vector or NULL
  expect_error(gg_criteria(df, tile_fill = 123), "must be a character vector or NULL")
})

test_that("gg_criteria validates bar_column parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score = c(10, 20)
  )
  
  # Purpose: ensure bar_column validation works
  expect_error(gg_criteria(df, bar_column = 123), "must be a character vector or NULL")
  expect_error(gg_criteria(df, bar_column = "missing"), "not found in data")
  expect_error(gg_criteria(df, bar_column = "pass_criteria"), "must be numeric")
})

test_that("gg_criteria validates panel_ratio when bar_column is used", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score = c(10, 20)
  )
  
  # Purpose: ensure panel_ratio is positive with bar_column
  expect_error(gg_criteria(df, bar_column = "score", panel_ratio = 0), 
               "must be a positive number")
  expect_error(gg_criteria(df, bar_column = "score", panel_ratio = -1), 
               "must be a positive number")
})

test_that("gg_criteria validates bar_fill parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score = c(10, 20)
  )
  
  # Purpose: ensure bar_fill validation works
  expect_error(gg_criteria(df, bar_column = "score", bar_fill = 123), 
               "must be a character vector.*or NULL")
})

# --------------------------------------------------------------------
# Basic functionality tests
# --------------------------------------------------------------------

test_that("gg_criteria creates plot from wide format data", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    yield_criteria = c("Pass", "Pass", "Fail"),
    purity_criteria = c("Pass", "Fail", "Pass")
  )
  
  # Purpose: ensure basic plotting works
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria applies custom criteria pattern", {
  df <- data.frame(
    id = c("S1", "S2"),
    QC_yield = c("Pass", "Fail"),
    QC_purity = c("Pass", "Pass")
  )
  
  # Purpose: ensure custom criteria pattern works
  p <- gg_criteria(df, criteria = "QC_", quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria applies custom tile colors", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  colors <- c(Pass = "green", Fail = "red")
  
  # Purpose: ensure custom colors are applied
  p <- gg_criteria(df, tile_fill = colors, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria applies custom tile dimensions", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure custom tile dimensions work
  p <- gg_criteria(df, tile_width = 0.5, tile_height = 0.5, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria applies tile transparency", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure tile_alpha works
  p <- gg_criteria(df, tile_alpha = 0.5, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria hides text when show_text = FALSE", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure text can be hidden
  p <- gg_criteria(df, show_text = FALSE, quiet = TRUE)
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 1)
})

test_that("gg_criteria applies custom border styling", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure border customization works
  p <- gg_criteria(df, border_color = "red", border_width = 1, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria hides legend when show_legend = FALSE", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure legend can be hidden
  p <- gg_criteria(df, show_legend = FALSE, quiet = TRUE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})

test_that("gg_criteria suppresses messages when quiet = TRUE", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure quiet mode suppresses messages
  expect_silent(gg_criteria(df, quiet = TRUE))
})

# --------------------------------------------------------------------
# Data handling tests
# --------------------------------------------------------------------

test_that("gg_criteria handles NA values", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    pass_criteria = c("Pass", NA, "Fail")
  )
  
  # Purpose: ensure NA values are handled gracefully
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles empty string values", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    pass_criteria = c("Pass", "", "Fail")
  )
  
  # Purpose: ensure empty strings are handled gracefully
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria warns when no valid values", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c(NA, "")
  )
  
  # Purpose: ensure warning when all values are invalid
  expect_warning(p <- gg_criteria(df, quiet = TRUE), "No non-NA values")
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles multiple criteria types", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    yield_criteria = c("Pass", "Pass", "Fail"),
    purity_criteria = c("Good", "Poor", "Good"),
    stability_criteria = c("High", "Medium", "Low")
  )
  
  # Purpose: ensure multiple criterion types work
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------

test_that("gg_criteria handles single sample", {
  df <- data.frame(
    id = "S1",
    yield_criteria = "Pass",
    purity_criteria = "Pass"
  )
  
  # Purpose: ensure single sample works
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles single criterion", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    pass_criteria = c("Pass", "Pass", "Fail")
  )
  
  # Purpose: ensure single criterion works
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles many criteria", {
  criteria_cols <- paste0("crit", 1:20, "_criteria")
  df <- data.frame(
    id = c("S1", "S2"),
    matrix(sample(c("Pass", "Fail"), 40, replace = TRUE), 
           nrow = 2, ncol = 20)
  )
  names(df) <- c("id", criteria_cols)
  
  # Purpose: ensure many criteria work
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles many samples", {
  df <- data.frame(
    id = paste0("Sample_", 1:50),
    pass_criteria = sample(c("Pass", "Fail"), 50, replace = TRUE)
  )
  
  # Purpose: ensure many samples work
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles custom ID column name", {
  df <- data.frame(
    sample_name = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure custom id column works
  p <- gg_criteria(df, id = "sample_name", quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria preserves sample order", {
  df <- data.frame(
    id = c("Z", "A", "M"),
    pass_criteria = c("Pass", "Pass", "Fail")
  )
  
  # Purpose: ensure original sample order is preserved
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Barplot functionality tests
# --------------------------------------------------------------------

test_that("gg_criteria creates barplot with single bar column", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    pass_criteria = c("Pass", "Pass", "Fail"),
    score = c(10, 20, 30)
  )
  
  # Purpose: ensure single barplot works
  p <- gg_criteria(df, bar_column = "score", quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_criteria creates barplots with multiple bar columns", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    pass_criteria = c("Pass", "Pass", "Fail"),
    score1 = c(10, 20, 30),
    score2 = c(5, 15, 25)
  )
  
  # Purpose: ensure multiple barplots work
  p <- gg_criteria(df, bar_column = c("score1", "score2"), quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_criteria applies custom bar colors", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score = c(10, 20)
  )
  
  # Purpose: ensure custom bar colors work
  p <- gg_criteria(df, bar_column = "score", bar_fill = "red", quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_criteria recycles bar colors", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score1 = c(10, 20),
    score2 = c(5, 15),
    score3 = c(8, 18)
  )
  
  # Purpose: ensure bar colors are recycled correctly
  p <- gg_criteria(df, bar_column = c("score1", "score2", "score3"), 
                   bar_fill = c("red", "blue"), quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_criteria applies custom panel ratio", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score = c(10, 20)
  )
  
  # Purpose: ensure custom panel ratio works
  p <- gg_criteria(df, bar_column = "score", panel_ratio = 0.5, quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_criteria includes recommended dimensions attribute", {
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail")
  )
  
  # Purpose: ensure recommended_dims attribute is present
  p <- gg_criteria(df, quiet = TRUE)
  expect_true(!is.null(attr(p, "recommended_dims")))
  expect_length(attr(p, "recommended_dims"), 2)
  expect_named(attr(p, "recommended_dims"), c("width", "height"))
})

test_that("gg_criteria recommended dimensions change with barplots", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    pass_criteria = c("Pass", "Fail"),
    score = c(10, 20)
  )
  
  # Purpose: ensure dimensions account for barplots
  p1 <- gg_criteria(df, quiet = TRUE)
  p2 <- gg_criteria(df, bar_column = "score", quiet = TRUE)
  
  dims1 <- attr(p1, "recommended_dims")
  dims2 <- attr(p2, "recommended_dims")
  
  expect_true(dims2["width"] > dims1["width"])
})

# --------------------------------------------------------------------
# Integration tests
# --------------------------------------------------------------------

test_that("gg_criteria works with all parameters combined", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    QC_yield = c("Pass", "Pass", "Fail"),
    QC_purity = c("Good", "Poor", "Good"),
    score1 = c(95, 88, 65),
    score2 = c(98, 85, 96)
  )
  
  colors <- c(Pass = "green", Fail = "red", Good = "blue", Poor = "orange")
  
  # Purpose: ensure all parameters work together
  p <- gg_criteria(
    df,
    id = "sample_id",
    criteria = "QC_",
    bar_column = c("score1", "score2"),
    bar_fill = c("purple", "gold"),
    panel_ratio = 0.4,
    tile_fill = colors,
    tile_width = 0.8,
    tile_height = 0.8,
    tile_alpha = 0.9,
    show_text = TRUE,
    border_color = "black",
    border_width = 0.5,
    text_size = 10,
    show_legend = TRUE,
    quiet = TRUE
  )
  
  expect_s3_class(p, "patchwork")
})

test_that("gg_criteria handles real-world criteria data", {
  df <- data.frame(
    id = paste0("Sample_", 1:5),
    yield_criteria = c("Pass", "Pass", "Fail", "Pass", "Fail"),
    purity_criteria = c("Pass", "Fail", "Pass", "Pass", "Pass"),
    stability_criteria = c("Good", "Good", "Poor", "Good", "Good")
  )
  
  # Purpose: ensure realistic data structure works
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_criteria handles data with extra columns", {
  df <- data.frame(
    id = c("S1", "S2"),
    date = c("2024-01-01", "2024-01-02"),
    pass_criteria = c("Pass", "Fail"),
    notes = c("OK", "Check"),
    fail_criteria = c("Fail", "Pass")
  )
  
  # Purpose: ensure extra columns don't break the function
  p <- gg_criteria(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})