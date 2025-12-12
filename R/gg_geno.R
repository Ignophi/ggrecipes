#' Genotype Heatmap with Optional Barplots
#'
#' @description
#' Creates a heatmap visualization of biallelic genotypes (e.g., SNPs, variants)
#' with split-tile representation showing phased or unphased alleles. Each tile
#' is divided diagonally to display both alleles, with border color indicating
#' phasing status (black for phased, white for unphased). Optional barplots can
#' be added to display associated numeric data.
#'
#' @param geno Character string specifying the regular expression pattern that
#'   identifies genotype columns (e.g., \code{"_geno$"} to match columns like
#'   \code{rs123_geno}, \code{rs456_geno}). The pattern is removed to create
#'   SNP labels. Default is "_geno".
#' @inheritParams gg_criteria
#'
#' @return A ggplot2 object (or patchwork object if \code{bar_column} is
#'   specified). The plot displays:
#'   \itemize{
#'     \item Split tiles representing biallelic genotypes
#'     \item Top-left triangle: first allele
#'     \item Bottom-right triangle: second allele
#'     \item Black borders: phased genotypes (separator: |)
#'     \item White borders: unphased genotypes (separator: /)
#'     \item Optional barplots showing numeric data for each sample
#'     \item Color legend mapping alleles to colors
#'   }
#'   The returned object has an attribute \code{"recommended_dims"} containing
#'   suggested plot width and height in inches.
#'
#' @details
#' The function expects genotype data in wide format with:
#' \itemize{
#'   \item One column for sample IDs (specified by \code{id})
#'   \item Multiple columns matching the \code{geno} pattern, each representing
#'     a SNP/variant
#'   \item Genotypes encoded as "allele1/allele2" (unphased) or "allele1|allele2"
#'     (phased)
#' }
#'
#' Genotype format examples:
#' \itemize{
#'   \item Unphased: \code{"0/0"}, \code{"0/1"}, \code{"1/1"}
#'   \item Phased: \code{"0|0"}, \code{"0|1"}, \code{"1|0"}
#' }
#'
#' Each genotype tile is split diagonally with the first allele in the top-left
#' triangle and the second allele in the bottom-right triangle. The border color
#' indicates phasing: black borders for phased genotypes (|) and white borders
#' for unphased genotypes (/).
#'
#' When \code{bar_column} is specified, the function requires the patchwork
#' package to combine the heatmap with barplots. Barplots are added to the
#' right of the heatmap and display numeric values with text labels.
#'
#' WARNING: Alignment between heatmap and barplots depends on plot dimensions.
#' The function provides recommended dimensions (accessible via
#' \code{attr(plot, "recommended_dims")}) that ensure proper alignment. You
#' can adjust these dimensions to improve appearance (e.g., reduce width to
#' tighten spacing, or scale proportionally for size) while maintaining
#' alignment.
#'
#' @import ggplot2
#' @importFrom stats setNames
#'
#' @examples
#' # Create example SNP and phenotype data
#' set.seed(123)
#' 
#' snp_data <- data.frame(
#'   id = paste0("P", sprintf("%03d", 1:12)),
#'   # SNP columns
#'   rs1234_geno = sample(c(c("0/0", "0/1", "1/1"), NA),
#'                        12, replace = TRUE, 
#'                        prob = c(0.4, 0.4, 0.15, 0.05)),
#'   rs5678_geno = sample(c("0/0", "0/1", "0/2", "1/1", "1/2", "2/2", NA),
#'                        12, replace = TRUE, 
#'                        prob = c(0.25, 0.25, 0.1, 0.15, 0.15, 0.05, 0.05)),
#'   rs9012_geno = sample(c(c("0|0", "0|1", "1|1", "0/1", "1/2"), NA),
#'                        12, replace = TRUE,  
#'                        prob = c(0.2, 0.2, 0.15, 0.2, 0.15, 0.1)),
#'   rs3456_geno = sample(c(c("0/0", "0/1", "1/1"), NA),
#'                        12, replace = TRUE,  
#'                        prob = c(0.45, 0.35, 0.15, 0.05)),
#'   rs7890_geno = sample(c("0/0", "0/1", "0/2", "1/3", "2/2", NA),
#'                        12, replace = TRUE,  
#'                        prob = c(0.3, 0.25, 0.15, 0.1, 0.15, 0.05)),
#'   rs2468_geno = sample(c("0|0", "0|1", "1|1", "1|2", NA),
#'                        12, replace = TRUE, 
#'                        prob = c(0.3, 0.35, 0.2, 0.1, 0.05)),
#'   
#'   # Phenotype columns for bar plots
#'   Age = sample(25:75, 12, replace = TRUE),
#'   BMI = round(rnorm(12, mean = 26, sd = 4), 1),
#'   Insulin = round(rnorm(12, mean = 12, sd = 3), 1)
#' )
#' 
#' # Base genotype plot
#' gg_geno(
#'   data = snp_data,
#'   id = "id",
#'   geno = "_geno$"
#' )
#' 
#' # Show optional barplots
#' gg_geno(
#'   data = snp_data,
#'   id = "id",
#'   geno = "_geno$",
#'   show_legend = TRUE,
#'   panel_ratio = 1,
#'   bar_column = c("Age", "BMI", "Insulin"),
#'   bar_fill = c("#c77d77", "#e0b46e", "#c7bc77"),
#'   text_size = 10
#' )
#' @export
gg_geno <- function(data,
                    id = "id",
                    geno = "_geno",
                    bar_column = NULL,
                    bar_fill = NULL,
                    panel_ratio = 0.3,
                    tile_fill = NULL,
                    tile_width = 0.7,
                    tile_height = 0.7,
                    border_width = 0.5,
                    text_size = 10,
                    show_legend = TRUE,
                    quiet = FALSE) {
  
  # Check packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }
  
  if (!is.null(bar_column) && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for bar_column but not installed")
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!is.character(id) || length(id) != 1) {
    stop("'id' must be a single character string")
  }
  
  if (!id %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", id))
  }
  
  if (!is.character(geno) || length(geno) != 1) {
    stop("'geno' must be a single character string")
  }
  
  if (!is.null(tile_fill) && !is.character(tile_fill)) {
    stop("'tile_fill' must be a named character vector or NULL")
  }
  
  if (!is.numeric(tile_width) || tile_width <= 0 || tile_width > 1) {
    stop("'tile_width' must be between 0 and 1")
  }
  
  if (!is.numeric(tile_height) || tile_height <= 0 || tile_height > 1) {
    stop("'tile_height' must be between 0 and 1")
  }
  
  if (!is.numeric(border_width) || border_width < 0) {
    stop("'border_width' must be a non-negative number")
  }
  
  if (!is.numeric(text_size) || text_size <= 0) {
    stop("'text_size' must be a positive number")
  }
  
  if (!is.logical(show_legend) || !is.logical(quiet)) {
    stop("'show_legend' and 'quiet' must be TRUE or FALSE")
  }
  
  # Validate bar_column
  if (!is.null(bar_column)) {
    if (!is.character(bar_column) || length(bar_column) < 1) {
      stop("'bar_column' must be a character vector or NULL")
    }
    
    missing_cols <- setdiff(bar_column, colnames(data))
    if (length(missing_cols) > 0) {
      stop(sprintf("Column(s) not found in data: %s",
                   paste(missing_cols, collapse = ", ")))
    }
    
    non_numeric <- bar_column[!sapply(data[bar_column], is.numeric)]
    if (length(non_numeric) > 0) {
      stop(sprintf("Column(s) must be numeric: %s",
                   paste(non_numeric, collapse = ", ")))
    }
    
    if (!is.numeric(panel_ratio) || panel_ratio <= 0) {
      stop("'panel_ratio' must be a positive number")
    }
  }
  
  # Validate bar_fill
  if (!is.null(bar_fill)) {
    if (!is.character(bar_fill)) {
      stop("'bar_fill' must be a character vector or NULL")
    }
    if (!is.null(bar_column) && length(bar_fill) < length(bar_column)) {
      stop(sprintf("'bar_fill' must have at least %d colors for the number of bar columns", 
                   length(bar_column)))
    }
  }

  # Define paired colors
  paired_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                     "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
                     "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  
  # Find genotype columns
  geno_cols <- grep(geno, names(data), value = TRUE)
  
  if (length(geno_cols) == 0) {
    stop(sprintf("No columns match pattern '%s'", geno))
  }
  
  # Store original order
  id_levels <- unique(data[[id]])
  
  # Reshape to long format
  base_cols <- setdiff(names(data), geno_cols)
  long_list <- lapply(geno_cols, function(col_name) {
    tmp <- data[base_cols]
    tmp$snp <- sub(geno, "", col_name)
    tmp$genotype <- data[[col_name]]
    tmp
  })
  
  df_plot <- do.call(rbind, long_list)
  rownames(df_plot) <- NULL
  
  snp_levels <- sub(geno, "", geno_cols)
  df_plot$snp <- factor(df_plot$snp, levels = snp_levels)
  
  # Parse genotypes
  parse_geno <- function(gt) {
    if (is.na(gt) || gt == "") return(NULL)
    
    if (grepl("\\|", gt)) {
      list(alleles = strsplit(gt, "\\|")[[1]], phased = TRUE)
    } else if (grepl("/", gt)) {
      list(alleles = strsplit(gt, "/")[[1]], phased = FALSE)
    } else {
      NULL
    }
  }
  
  parsed <- lapply(df_plot$genotype, parse_geno)
  df_plot$allele1 <- sapply(parsed, function(x) if (is.null(x) || length(x$alleles) != 2) NA else x$alleles[1])
  df_plot$allele2 <- sapply(parsed, function(x) if (is.null(x) || length(x$alleles) != 2) NA else x$alleles[2])
  df_plot$phased <- sapply(parsed, function(x) if (is.null(x)) NA else x$phased)
  
  # Remove invalid rows
  df_plot <- df_plot[complete.cases(df_plot[c("allele1", "allele2")]), ]
  
  if (nrow(df_plot) == 0) {
    warning("No valid genotypes to plot")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  
  # Set factor levels
  df_plot[[id]] <- factor(df_plot[[id]], levels = id_levels)
  
  # Setup colors
  all_alleles <- unique(c(df_plot$allele1, df_plot$allele2))
  
  if (is.null(tile_fill)) {
    tile_fill <- setNames(paired_colors[seq_along(all_alleles)], sort(all_alleles))
  }
  
  # Create polygon data
  w <- tile_width / 2
  h <- tile_height / 2
  
  tile_data <- do.call(rbind, lapply(seq_len(nrow(df_plot)), function(i) {
    row <- df_plot[i, ]
    x <- as.numeric(row$snp)
    y <- as.numeric(row[[id]])
    border_col <- if (row$phased) "black" else "white"
    
    # Top-left and bottom-right triangles
    rbind(
      data.frame(
        x = c(x - w, x + w, x - w, x - w),
        y = c(y + h, y + h, y - h, y + h),
        group = paste(i, "tri1", sep = "_"),
        allele = row$allele1,
        border_color = border_col
      ),
      data.frame(
        x = c(x + w, x + w, x - w, x + w),
        y = c(y + h, y - h, y - h, y + h),
        group = paste(i, "tri2", sep = "_"),
        allele = row$allele2,
        border_color = border_col
      )
    )
  }))
  
  # Build heatmap
  p <- ggplot2::ggplot(tile_data, ggplot2::aes(
    x = .data$x,
    y = .data$y,
    group = .data$group,
    fill = .data$allele,
    color = .data$border_color
  )) +
    ggplot2::geom_polygon(linewidth = border_width) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_manual(values = tile_fill, name = "Allele") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(snp_levels),
      labels = snp_levels,
      expand = c(0, 0),
      position = "top"
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(id_levels),
      labels = id_levels,
      limits = c(0.5, length(id_levels) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, colour = "black"),
      axis.text.y = ggplot2::element_text(colour = "black", vjust = 0.5, size = text_size),
      axis.text.x = ggplot2::element_text(size = text_size),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "black", linewidth = 0.1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0, unit = "pt"),
      legend.position = if (show_legend) "left" else "none"
    )
  
  # Calculate dimensions
  n_rows <- length(id_levels)
  n_cols <- length(snp_levels)
  heatmap_height <- n_rows * 0.4 + 1.5
  heatmap_width <- n_cols * 0.4 + 2
  
  # Add barplots if requested
  if (!is.null(bar_column)) {
    # Get bar colors
    if (!is.null(bar_fill)) {
      bar_colors <- bar_fill
    } else {
      bar_colors <- paired_colors
    }
    
    # Create barplot function
    make_barplot <- function(col, color) {
      bar_data <- df_plot[!duplicated(df_plot[[id]]), c(id, col)]
      bar_data[[id]] <- factor(bar_data[[id]], levels = id_levels)
      
      ggplot2::ggplot(bar_data, ggplot2::aes(x = .data[[col]], y = .data[[id]])) +
        ggplot2::geom_bar(stat = "identity", fill = color, width = tile_height,
                         color = "black", linewidth = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = round(.data[[col]], 1)),
                          hjust = -0.1, size = 3, color = "black") +
        ggplot2::scale_x_continuous(expand = c(0, 0, 0.15, 0), position = "top") +
        ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = 0, add = c(0.5, 0.5)), drop = FALSE) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.text.x.top = ggplot2::element_text(vjust = 0, colour = "black", size = text_size),
          axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
          axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
          panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          plot.background = ggplot2::element_rect(fill = "white", color = NA)
        ) +
        ggplot2::labs(x = col, y = NULL)
    }
    
    bar_plots <- mapply(make_barplot, bar_column, 
                       bar_colors[((seq_along(bar_column) - 1) %% length(bar_colors)) + 1],
                       SIMPLIFY = FALSE)
    
    widths <- c(1, rep(panel_ratio / length(bar_column), length(bar_column)))
    combined <- Reduce(`+`, c(list(p), bar_plots)) + 
      patchwork::plot_layout(widths = widths)
    
    total_width <- heatmap_width + length(bar_column) * 2
    attr(combined, "recommended_dims") <- c(width = total_width, height = heatmap_height)
    
    if (!quiet) {
      message(sprintf("Recommended dimensions: %.1f x %.1f inches",
                     total_width, heatmap_height))
    }
    
    return(combined)
  }
  
  # Return heatmap only
  attr(p, "recommended_dims") <- c(width = heatmap_width, height = heatmap_height)
  
  if (!quiet) {
    message(sprintf("Recommended dimensions: %.1f x %.1f inches",
                   heatmap_width, heatmap_height))
  }
  
  p
}