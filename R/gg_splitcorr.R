#' Split-Correlation Heatmap
#'
#' Creates a split-correlation plot where the upper triangle shows correlations 
#' for one subgroup and the lower triangle for another, based on a binary 
#' splitting variable. This allows quick visual comparison of correlation 
#' structures between two groups. Significant correlations (after multiple 
#' testing adjustment) can be labeled directly in the plot.
#'
#' @param data A data frame containing numeric variables to correlate and the 
#'   variable to split by.
#' @param split A character string specifying the name of the binary variable 
#'   in \code{data} used to split the dataset.
#' @param style Type of visualization; either \code{"tile"} (default) for a 
#'   heatmap or \code{"point"} for a bubble-style plot.
#' @param method Correlation method to use; either \code{"pearson"} (default) 
#'   or \code{"spearman"}.
#' @param padjust Method for p-value adjustment, passed to 
#'   \code{\link[stats]{p.adjust}} (default \code{"BH"}).
#' @param use Handling of missing values, passed to 
#'   \code{\link[stats]{cor}} (default \code{"complete.obs"}).
#' @param colors A vector of three colors for the low, mid, and high values 
#'   of the correlation scale (default \code{c("blue", "white", "red")}).
#' @param text_colors A vector of two colors for the text labels, used for 
#'   negative and positive correlations (default \code{c("white", "black")}).
#' @param text_size Numeric value giving the size of correlation text labels 
#'   (default \code{3.5}).
#' @param border_color Color for tile or point borders (default \code{"black"}).
#' @param prefix Character string prefix for group labels. If NULL (default),
#'   uses "split_variable = " format. Default is NULL.
#' @param linetype Type of diagonal line separating upper and lower triangles; 
#'   one of \code{"solid"}, \code{"dashed"}, or \code{"dotdash"} (default \code{"dashed"}).
#' @param linealpha Alpha transparency for the diagonal line (default \code{0.5}).
#' @param offset Numeric offset for the position of the group labels (default \code{0.75}).
#'
#' @return A ggplot2 object showing the split-correlation heatmap. The plot displays:
#'   \itemize{
#'     \item Upper triangle: correlations for the first level of the split variable
#'     \item Lower triangle: correlations for the second level of the split variable
#'     \item Diagonal line separating the two triangles
#'     \item Group labels indicating which split level is shown in each triangle
#'     \item Correlation values displayed only for significant pairs (p ≤ 0.05 after adjustment)
#'     \item Color gradient representing correlation strength (-1 to 1)
#'     \item Optional point size (if \code{style = "point"}) indicating absolute correlation strength
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Splits the dataset into two groups using the variable specified by \code{split}.
#'   \item Computes pairwise correlations and p-values for each group (via a helper \code{cor_p()}).
#'   \item Combines the upper triangle from one group and lower triangle from the other.
#'   \item Adjusts p-values using the selected method and annotates significant cells.
#' }
#' The result is a heatmap (or bubble plot) showing both groups' correlation patterns
#' in a single compact visualization.
#'
#' @examples
#' # Compare correlations between V-shaped vs straight engines
#' data(mtcars)
#' 
#' gg_splitcorr(
#'   data = mtcars,
#'   split = "vs",
#'   prefix = "Engine Type: "
#' )
#' 
#' # Alternative style "point"
#' gg_splitcorr(
#'   data = mtcars,
#'   split = "vs",
#'   style = "point",
#'   method = "spearman",
#'   prefix = "Engine Type: "
#' )
#' @seealso
#' \code{\link[stats]{cor}} for correlation computation,
#' \code{\link[stats]{p.adjust}} for multiple testing correction,
#' \code{\link[ggplot2]{geom_tile}},
#' \code{\link[ggplot2]{geom_point}}
#'
#' @import ggplot2
#' @importFrom stats p.adjust cor pt complete.cases
#' @export
gg_splitcorr <- function(data, 
                         split,
                         style = "tile",
                         method = "pearson",
                         padjust = "BH",
                         use = "complete.obs",
                         colors = c("blue", "white", "red"),
                         text_colors = c("white", "black"),
                         text_size = 3.5,
                         border_color = "black",
                         prefix = NULL,
                         linetype = "dashed",
                         linealpha = 0.5,
                         offset = 0.75) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!split %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in data", split))
  }

  # Check binary split variable
  split_vals <- unique(data[[split]][!is.na(data[[split]])])
  if (length(split_vals) != 2) {
    stop(sprintf("'%s' must have exactly 2 unique values, found %d", 
                 split, length(split_vals)))
  }
  
  valid_styles <- c("tile", "point")
  if (!style %in% valid_styles) {
    stop(sprintf("'style' must be one of %s", 
                 paste(shQuote(valid_styles), collapse = ", ")))
  }
  
  valid_methods <- c("pearson", "spearman")
  if (!method %in% valid_methods) {
    stop(sprintf("'method' must be one of %s", 
                 paste(shQuote(valid_methods), collapse = ", ")))
  }
  
  valid_linetypes <- c("solid", "dashed", "dotdash")
  if (!linetype %in% valid_linetypes) {
    stop(sprintf("'linetype' must be one of %s", 
                 paste(shQuote(valid_linetypes), collapse = ", ")))
  }
  
  if (!is.numeric(text_size) || text_size <= 0) {
    stop("'text_size' must be a positive number")
  }
  
  if (!is.numeric(linealpha) || linealpha < 0 || linealpha > 1) {
    stop("'linealpha' must be between 0 and 1")
  }
  
  if (length(colors) != 3) {
    stop("'colors' must be a vector of length 3")
  }
  
  if (length(text_colors) != 2) {
    stop("'text_colors' must be a vector of length 2")
  }

  if (!is.numeric(offset) || offset < 0) {
    stop("'offset' must be a non-negative number")
  }

  # Split dataframe based on variable
  split_df <- split(data, data[[split]])

  # Calculate pairwise correlation
  cor_df <- lapply(split_df, function(df) {
    # drop split variable
    df <- df[ , !(names(df) == split), drop = FALSE]
    # keep numeric only
    df <- df[ , sapply(df, is.numeric), drop = FALSE]
    if (ncol(df) < 2) {
      stop("Need at least 2 numeric columns after removing split variable")
    }
    cor_p(df, use = use, method = method)
  })

  # Set diagonal to 1 in order not to display its text
  r_mat <- cor_df[[1]]$r
  r_mat[lower.tri(r_mat)] <- cor_df[[2]]$r[lower.tri(cor_df[[2]]$r)]
  diag(r_mat) <- 0

  p_mat <- cor_df[[1]]$p
  p_mat[lower.tri(p_mat)] <- cor_df[[2]]$p[lower.tri(cor_df[[2]]$p)]
  diag(p_mat) <- 1

  # Compute upper.tri/lower.tri
  ut <- upper.tri(p_mat); lt <- lower.tri(p_mat)
  p_mat[ut] <- p.adjust(p_mat[ut], method = padjust)
  p_mat[lt] <- p.adjust(p_mat[lt], method = padjust)

  # Convert the correlation matrix to a long format
  plot_df <- melt(r_mat)
  # Add p-values to the melted data
  plot_df$p_value <- melt(p_mat)$value

  n_vars <- ncol(r_mat)
  prefix <- if (is.null(prefix)) paste0(split, " = ") else prefix
  lab_ut <- paste0(prefix, names(cor_df)[1])
  lab_lt <- paste0(prefix, names(cor_df)[2])

  # Plot heatmap using ggplot2
  style_layer <- if (style == "point") {
    list(
      ggplot2::geom_point(ggplot2::aes(size = abs(.data$value)), shape = 21, color = border_color, stroke = 0.5),
      ggplot2::scale_size_area(max_size = 12, limits = c(0, 1), guide = "none")
    )
  } else {
    ggplot2::geom_tile(color = border_color, linewidth = 0.5)
  }
  
  p <- ggplot2::ggplot(data = plot_df,
                       ggplot2::aes(.data$Var1, .data$Var2,
                                    fill = .data$value)) +
    style_layer +
    ggplot2::geom_text(ggplot2::aes(
      label = ifelse(.data$p_value <= 0.05, round(.data$value, 1), ""),
      color = ifelse(.data$value < 0, text_colors[1], text_colors[2])),
      fontface = "bold", size = text_size, show.legend = FALSE) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3],
                                   midpoint = 0, limits = c(-1, 1),
                                   name = "Correlation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = ggplot2::margin(t = 20, r = 20, b = 5, l = 5, unit = "pt")
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::annotate("text", x = 0.5, y = n_vars + offset,
                      label = lab_ut, size = 6, fontface = "bold", color = "black",
                      hjust = 0, vjust = 0) +
    ggplot2::annotate("text", x = n_vars + offset, y = 0.5,
                      label = lab_lt, size = 6, fontface = "bold", color = "black",
                      angle = 90, hjust = 0, vjust = 1) +
    ggplot2::geom_segment(x = 0.5, y = 0.5,
                          xend = n_vars + 0.5, yend = n_vars + 0.5,
                          linewidth = 0.02, color = "black",
                          linetype = linetype, alpha = linealpha) +
    ggplot2::coord_fixed(clip = "off", expand = FALSE)
  return(p)
}