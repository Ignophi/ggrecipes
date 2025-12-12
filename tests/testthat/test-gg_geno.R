# tests/testthat/test-gg_geno.R

# --------------------------------------------------------------------
# Input validation tests
# --------------------------------------------------------------------

test_that("gg_geno errors when data is not a data frame", {
  # Purpose: ensure non-data.frame inputs are rejected
  expect_error(gg_geno(data = "not_a_df"), "'data' must be a data frame")
  expect_error(gg_geno(data = NULL), "'data' must be a data frame")
  expect_error(gg_geno(data = list(a = 1)), "'data' must be a data frame")
})

test_that("gg_geno errors when id is invalid", {
  df <- data.frame(
    sample = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure id parameter validation works
  expect_error(gg_geno(df, id = c("a", "b")), "must be a single character string")
  expect_error(gg_geno(df, id = 123), "must be a single character string")
  expect_error(gg_geno(df, id = "missing"), "not found in data")
})

test_that("gg_geno errors when geno is invalid", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure geno parameter validation works
  expect_error(gg_geno(df, geno = c("a", "b")), "must be a single character string")
  expect_error(gg_geno(df, geno = 123), "must be a single character string")
})

test_that("gg_geno errors when no genotype columns match pattern", {
  df <- data.frame(
    id = c("S1", "S2"),
    value = c(1, 2)
  )
  
  # Purpose: ensure meaningful error when pattern doesn't match
  expect_error(gg_geno(df, geno = "_geno$"), "No columns match pattern")
})

test_that("gg_geno validates tile parameters", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure tile dimension parameters are valid
  expect_error(gg_geno(df, tile_width = -0.1), "must be between 0 and 1")
  expect_error(gg_geno(df, tile_width = 1.5), "must be between 0 and 1")
  expect_error(gg_geno(df, tile_height = 0), "must be between 0 and 1")
  expect_error(gg_geno(df, tile_height = 2), "must be between 0 and 1")
})

test_that("gg_geno validates border_width parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure border_width is non-negative
  expect_error(gg_geno(df, border_width = -1), "must be a non-negative number")
})

test_that("gg_geno validates text_size parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure text_size is positive
  expect_error(gg_geno(df, text_size = 0), "must be a positive number")
  expect_error(gg_geno(df, text_size = -1), "must be a positive number")
})

test_that("gg_geno validates logical parameters", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure logical parameters are TRUE or FALSE
  expect_error(gg_geno(df, show_legend = "TRUE"), "must be TRUE or FALSE")
  expect_error(gg_geno(df, quiet = 1), "must be TRUE or FALSE")
})

test_that("gg_geno validates tile_fill parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure tile_fill is character vector or NULL
  expect_error(gg_geno(df, tile_fill = 123), "must be a named character vector or NULL")
})

test_that("gg_geno validates bar_column parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1"),
    score = c(10, 20)
  )
  
  # Purpose: ensure bar_column validation works
  expect_error(gg_geno(df, bar_column = 123), "must be a character vector or NULL")
  expect_error(gg_geno(df, bar_column = "missing"), "not found in data")
  expect_error(gg_geno(df, bar_column = "rs123_geno"), "must be numeric")
})

test_that("gg_geno validates panel_ratio when bar_column is used", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1"),
    score = c(10, 20)
  )
  
  # Purpose: ensure panel_ratio is positive with bar_column
  expect_error(gg_geno(df, bar_column = "score", panel_ratio = 0), 
               "must be a positive number")
  expect_error(gg_geno(df, bar_column = "score", panel_ratio = -1), 
               "must be a positive number")
})

test_that("gg_geno validates bar_fill parameter", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1"),
    score = c(10, 20)
  )
  
  # Purpose: ensure bar_fill has enough colors
  expect_error(gg_geno(df, bar_column = "score", bar_fill = 123), 
               "must be a character vector or NULL")
  expect_error(gg_geno(df, bar_column = "score", bar_fill = character(0)), 
               "must have at least 1 colors")
})

# --------------------------------------------------------------------
# Basic functionality tests
# --------------------------------------------------------------------

test_that("gg_geno creates plot from valid unphased genotype data", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "0/1", "1/1"),
    rs456_geno = c("0/1", "1/1", "0/0")
  )
  
  # Purpose: ensure basic plotting works with unphased genotypes
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno creates plot from valid phased genotype data", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0|0", "0|1", "1|0"),
    rs456_geno = c("0|1", "1|1", "0|0")
  )
  
  # Purpose: ensure basic plotting works with phased genotypes
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles mixed phased and unphased genotypes", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "0|1", "1/1"),
    rs456_geno = c("0|1", "1/1", "0/0")
  )
  
  # Purpose: ensure mixed phasing is handled correctly
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno applies custom genotype pattern", {
  df <- data.frame(
    id = c("S1", "S2"),
    SNP_rs123 = c("0/0", "0/1"),
    SNP_rs456 = c("1/1", "0/1")
  )
  
  # Purpose: ensure custom geno pattern works
  p <- gg_geno(df, geno = "SNP_", quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno applies custom tile colors", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  colors <- c("0" = "red", "1" = "blue")
  
  # Purpose: ensure custom colors are applied
  p <- gg_geno(df, tile_fill = colors, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno applies custom tile dimensions", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure custom tile dimensions work
  p <- gg_geno(df, tile_width = 0.5, tile_height = 0.5, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno applies custom border width", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure custom border width works
  p <- gg_geno(df, border_width = 1, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno hides legend when show_legend = FALSE", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure legend can be hidden
  p <- gg_geno(df, show_legend = FALSE, quiet = TRUE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})

test_that("gg_geno suppresses messages when quiet = TRUE", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure quiet mode suppresses messages
  expect_silent(gg_geno(df, quiet = TRUE))
})

# --------------------------------------------------------------------
# Genotype parsing tests
# --------------------------------------------------------------------

test_that("gg_geno handles multi-allelic genotypes", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("A/T", "G/C")
  )
  
  # Purpose: ensure non-numeric alleles work
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno skips invalid genotype formats", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "invalid", "1/1")
  )
  
  # Purpose: ensure invalid genotypes are skipped without error
  expect_warning(p <- gg_geno(df, quiet = TRUE), NA)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles NA genotypes", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", NA, "1/1")
  )
  
  # Purpose: ensure NA genotypes are handled gracefully
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles empty genotypes", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "", "1/1")
  )
  
  # Purpose: ensure empty strings are handled gracefully
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno warns when no valid genotypes", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c(NA, "invalid")
  )
  
  # Purpose: ensure warning when all genotypes are invalid
  expect_warning(p <- gg_geno(df, quiet = TRUE), "No valid genotypes")
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Edge cases
# --------------------------------------------------------------------

test_that("gg_geno handles single sample", {
  df <- data.frame(
    id = "S1",
    rs123_geno = "0/0",
    rs456_geno = "0/1"
  )
  
  # Purpose: ensure single sample data works
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles single SNP", {
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "0/1", "1/1")
  )
  
  # Purpose: ensure single SNP data works
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles many SNPs", {
  snps <- paste0("rs", 1:20, "_geno")
  df <- data.frame(
    id = c("S1", "S2"),
    matrix(sample(c("0/0", "0/1", "1/1"), 40, replace = TRUE), 
           nrow = 2, ncol = 20)
  )
  names(df) <- c("id", snps)
  
  # Purpose: ensure many SNPs work
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles many samples", {
  df <- data.frame(
    id = paste0("Sample_", 1:50),
    rs123_geno = sample(c("0/0", "0/1", "1/1"), 50, replace = TRUE)
  )
  
  # Purpose: ensure many samples work
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles homozygous reference", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/0")
  )
  
  # Purpose: ensure all-reference genotypes work
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles homozygous alternate", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("1/1", "1/1")
  )
  
  # Purpose: ensure all-alternate genotypes work
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles custom ID column name", {
  df <- data.frame(
    sample_name = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure custom id column works
  p <- gg_geno(df, id = "sample_name", quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

# --------------------------------------------------------------------
# Barplot functionality tests
# --------------------------------------------------------------------

test_that("gg_geno creates barplot with single bar column", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "0/1", "1/1"),
    score = c(10, 20, 30)
  )
  
  # Purpose: ensure single barplot works
  p <- gg_geno(df, bar_column = "score", quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_geno creates barplots with multiple bar columns", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2", "S3"),
    rs123_geno = c("0/0", "0/1", "1/1"),
    score1 = c(10, 20, 30),
    score2 = c(5, 15, 25)
  )
  
  # Purpose: ensure multiple barplots work
  p <- gg_geno(df, bar_column = c("score1", "score2"), quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_geno applies custom bar colors", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1"),
    score = c(10, 20)
  )
  
  # Purpose: ensure custom bar colors work
  p <- gg_geno(df, bar_column = "score", bar_fill = "red", quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_geno applies custom panel ratio", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1"),
    score = c(10, 20)
  )
  
  # Purpose: ensure custom panel ratio works
  p <- gg_geno(df, bar_column = "score", panel_ratio = 0.5, quiet = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("gg_geno includes recommended dimensions attribute", {
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1")
  )
  
  # Purpose: ensure recommended_dims attribute is present
  p <- gg_geno(df, quiet = TRUE)
  expect_true(!is.null(attr(p, "recommended_dims")))
  expect_length(attr(p, "recommended_dims"), 2)
  expect_named(attr(p, "recommended_dims"), c("width", "height"))
})

test_that("gg_geno recommended dimensions change with barplots", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    id = c("S1", "S2"),
    rs123_geno = c("0/0", "0/1"),
    score = c(10, 20)
  )
  
  # Purpose: ensure dimensions account for barplots
  p1 <- gg_geno(df, quiet = TRUE)
  p2 <- gg_geno(df, bar_column = "score", quiet = TRUE)
  
  dims1 <- attr(p1, "recommended_dims")
  dims2 <- attr(p2, "recommended_dims")
  
  expect_true(dims2["width"] > dims1["width"])
})

# --------------------------------------------------------------------
# Integration tests
# --------------------------------------------------------------------

test_that("gg_geno works with all parameters combined", {
  skip_if_not_installed("patchwork")
  
  df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    SNP_rs123 = c("0|0", "0|1", "1|0"),
    SNP_rs456 = c("A/T", "T/T", "A/A"),
    score1 = c(10, 20, 30),
    score2 = c(5, 15, 25)
  )
  
  colors <- c("0" = "red", "1" = "blue", "A" = "green", "T" = "yellow")
  
  # Purpose: ensure all parameters work together
  p <- gg_geno(
    df,
    id = "sample_id",
    geno = "SNP_",
    bar_column = c("score1", "score2"),
    bar_fill = c("purple", "orange"),
    panel_ratio = 0.4,
    tile_fill = colors,
    tile_width = 0.6,
    tile_height = 0.6,
    border_width = 0.8,
    text_size = 10,
    show_legend = TRUE,
    quiet = TRUE
  )
  
  expect_s3_class(p, "patchwork")
})

test_that("gg_geno handles real-world genotype data structure", {
  df <- data.frame(
    id = paste0("Sample_", 1:5),
    rs1234_geno = c("0/0", "0/1", "1/1", "0|1", "1|0"),
    rs5678_geno = c("A/A", "A/G", "G/G", "A|G", "G|A"),
    rs9012_geno = c("C/T", "T/T", "C/C", "C|T", "T|C")
  )
  
  # Purpose: ensure realistic data structure works
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("gg_geno handles data with extra columns", {
  df <- data.frame(
    id = c("S1", "S2"),
    age = c(25, 30),
    rs123_geno = c("0/0", "0/1"),
    gender = c("M", "F"),
    rs456_geno = c("1/1", "0/1")
  )
  
  # Purpose: ensure extra columns don't break the function
  p <- gg_geno(df, quiet = TRUE)
  expect_s3_class(p, "ggplot")
})