#' Rank Shift Plot
#'
#' @description
#' Creates a three-panel visualization comparing ranks of samples between two
#' groups. Left and right panels show distributions (bars or boxplots) for each
#' group ordered by rank. The center panel displays connecting lines colored by
#' rank change direction, revealing which samples improved, declined, or
#' maintained their relative position.
#'
#' @param data A data frame containing the data to visualize.
#' @param id Character string specifying the column name in \code{data} that
#'   contains sample identifiers. Default is "id".
#' @param group Character string specifying the column name in \code{data} that
#'   defines the two groups to compare. Must contain exactly 2 unique values.
#'   Default is "group".
#' @param value Character string specifying the column name in \code{data} that
#'   contains the numeric values used for ranking. Default is "value".
#' @param style Character string specifying visualization type: \code{"bar"} for
#'   barplots or \code{"box"} for boxplots. Default is "box".
#' @param fill Character vector of length 2 specifying fill colors for
#'   the two groups. If NULL (default), uses \code{c("#92b9de", "#e8927c")}.
#' @param alpha Numeric value (0-1) specifying transparency of bars/boxes.
#'   Default is 0.7.
#' @param line_alpha Numeric value (0-1) specifying transparency of connecting
#'   lines. Default is 0.5.
#' @param line_width Numeric value specifying width of connecting lines.
#'   Default is 0.3.
#' @param panel_ratio Numeric value specifying the relative width of the center
#'   line panel compared to side panels (which have width 1). Default is 1.
#' @param text_size Numeric value specifying the base size for axis text.
#'   Default is 11.
#' @param show_points Logical indicating whether to overlay individual data
#'   points on bars/boxes. Default is TRUE.
#' @param point_size Numeric value specifying size of points when
#'   \code{show_points = TRUE}. Default is 1.
#' @param point_shape Numeric value (0-25) specifying point shape when
#'   \code{show_points = TRUE}. Default is 21 (filled circle).
#' @param point_alpha Numeric value (0-1) specifying transparency of points when
#'   \code{show_points = TRUE}. Default is 0.7.
#' @param stat_summary Character string specifying the summary statistic used
#'   for ranking: \code{"mean"} or \code{"median"}. Default is "mean".
#' @param decreasing Logical indicating rank order direction. If TRUE, higher
#'   values receive lower rank numbers (rank 1 = highest value). If FALSE,
#'   lower values receive lower rank numbers. Default is FALSE.
#' @param free_x Logical indicating whether x-axes are independent between
#'   panels. If FALSE, both panels share the same x-axis limits. Default is TRUE.
#' @param rank_change_colors Named character vector of length 3 specifying
#'   colors for rank changes. Must contain names \code{"increase"},
#'   \code{"decrease"}, and \code{"no_change"}. Default is
#'   \code{c(increase = "#d73027", decrease = "#4575b4", no_change = "grey50")}.
#'
#' @return A patchwork object combining three plots: left panel showing the first
#'   group's distribution, center panel showing rank change lines, and right
#'   panel showing the second group's distribution. Line colors indicate whether
#'   samples increased in rank (moved up), decreased in rank (moved down), or
#'   maintained the same rank between groups.
#'
#' @details
#' The function identifies samples common to both groups, calculates summary
#' statistics (mean or median) for each sample within each group, and assigns
#' ranks based on these values. Samples are then displayed in rank order with
#' connecting lines showing how ranks changed between groups.
#'
#' When \code{style = "box"} and \code{show_points = TRUE}, boxplots display
#' the distribution while individual points show all replicates. This is useful
#' for visualizing measurement variability alongside rank changes.
#'
#' The \code{decreasing} parameter controls ranking direction:
#' \itemize{
#'   \item \code{FALSE} (default): Rank 1 = lowest value, higher ranks =
#'     higher values
#'   \item \code{TRUE}: Rank 1 = highest value, lower ranks = lower values
#' }
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stats reorder
#'
#' @examples
#' # Example data: bacterial strain growth rate control vs antibiotic treated
#' growth_data <- data.frame(
#'   strain = rep(paste0("Strain", 1:13), each = 6),
#'   condition = rep(c("Control", "Treated"), each = 3, times = 13),
#'   growth_rate = c(
#'     rnorm(39, mean = 0.85, sd = 0.12),  # Control
#'     rnorm(39, mean = 0.45, sd = 0.10)   # Treated (reduced growth)
#'   )
#' )
#'
#' # Basic rank shift plot
#' gg_rankshift(
#'   data = growth_data,
#'   id = "strain",
#'   group = "condition",
#'   value = "growth_rate",
#' )
#'
#' # With barplots instead of boxplots
#' gg_rankshift(
#'   data = growth_data,
#'   id = "strain",
#'   group = "condition",
#'   value = "growth_rate",
#'   style = "bar"
#' )
#'
#' # Custom styling
#' gg_rankshift(
#'   data = growth_data,
#'   id = "strain",
#'   group = "condition",
#'   value = "growth_rate",
#'   fill = c("#e41a1c", "#377eb8"),
#'   rank_change_colors = c(
#'     increase = "#1b9e77",
#'     decrease = "#d95f02",
#'     no_change = "#7570b3"
#'   ),
#'   panel_ratio = 0.5,
#'   point_size = 2.5,
#'   line_width = 1,
#'   decreasing = TRUE
#' )
#' @export
gg_rankshift <- function(data,
                         id = "id",
                         group = "group",
                         value = "value",
                         style = "box",
                         fill = NULL,
                         alpha = 0.7,
                         line_alpha = 0.5,
                         line_width = 0.5,
                         panel_ratio = 0.5,
                         text_size = 11,
                         show_points = TRUE,
                         point_size = 1,
                         point_shape = 21,
                         point_alpha = 0.7,
                         stat_summary = "mean",
                         decreasing = FALSE,
                         free_x = TRUE,
                         rank_change_colors = c(
                           increase = "#d73027",
                           decrease = "#4575b4",
                           no_change = "grey50"
                         )) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  valid_styles <- c("bar", "box")
  if (length(style) != 1 || is.na(style) || !style %in% valid_styles) {
    stop(sprintf("'style' must be one of %s",
                paste(shQuote(valid_styles), collapse = ", ")))
  }

  if (!is.numeric(panel_ratio) || panel_ratio <= 0) {
    stop("'line_panel' must be a positive number")
  }

  if (!is.logical(show_points)) {
    stop("'show_points' must be TRUE or FALSE")
  }

  if (!is.numeric(point_size) || point_size < 0) {
    stop("'point_size' must be a non-negative number")
  }

  if (!is.numeric(point_shape) || !point_shape %in% 0:25) {
    stop("'point_shape' must be a number between 0 and 25")
  }

  if (!is.numeric(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("'point_alpha' must be between 0 and 1")
  }

  if (!is.logical(decreasing)) {
    stop("'decreasing' must be TRUE or FALSE")
  }

  if (!is.logical(free_x)) {
    stop("'free_x' must be TRUE or FALSE")
  }

  valid_stats <- c("mean", "median")
  if (!stat_summary %in% valid_stats) {
    stop(sprintf("'stat_summary' must be one of %s",
                paste(shQuote(valid_stats), collapse = ", ")))
  }
    
  required_cols <- c(id, group, value)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Required column(s) missing from data frame: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  group_levels <- unique(data[[group]][!is.na(data[[group]])])
  if (length(group_levels) != 2) {
    stop(sprintf(
      "'%s' must have exactly 2 unique values, found %d",
      group, length(group_levels)
    ))
  }
  
  if (!is.numeric(data[[value]])) {
    stop(sprintf("Column '%s' must be numeric", value))
  }
  
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be between 0 and 1")
  }
  
  if (!is.numeric(line_alpha) || line_alpha < 0 || line_alpha > 1) {
    stop("'line_alpha' must be between 0 and 1")
  }
  
  if (!is.numeric(line_width) || line_width < 0) {
    stop("'line_width' must be a non-negative number")
  }
  
  if (!is.null(fill)) {
    if (!is.character(fill) || length(fill) != 2) {
      stop("'fill' must be a character vector of length 2")
    }
  } else {
    fill <- c("#92b9de", "#e8927c")
  }

  # Define variables

  
  # Validate rank_change_colors
  required_names <- c("increase", "decrease", "no_change")
  if (!is.character(rank_change_colors) || length(rank_change_colors) != 3) {
    stop("'rank_change_colors' must be a character vector of length 3")
  }
  if (!all(required_names %in% names(rank_change_colors))) {
    stop(sprintf(
      "'rank_change_colors' must be a named vector with names: %s",
      paste(required_names, collapse = ", ")
    ))
  }
  
  # Find common IDs
  group_data <- split(data, data[[group]])
  common_ids <- Reduce(intersect, lapply(group_data, function(x) x[[id]]))
  
  if (length(common_ids) == 0) {
    stop("No common samples found between the two groups")
  }
  
  # Filter and rank
  data_common <- data[data[[id]] %in% common_ids, ]
  group_data <- split(data_common, data_common[[group]])
  
  # Aggregate and rank
  stat_fun <- if (stat_summary == "mean") mean else median

  df_1 <- aggregate(
    group_data[[1]][[value]],
    by = list(id = group_data[[1]][[id]]),
    FUN = stat_fun,
    na.rm = TRUE
  )
  names(df_1) <- c(id, value)
  df_1 <- df_1[order(df_1[[value]], decreasing = decreasing), ]
  df_1$rank <- seq_len(nrow(df_1))

  df_2 <- aggregate(
    group_data[[2]][[value]],
    by = list(id = group_data[[2]][[id]]),
    FUN = stat_fun,
    na.rm = TRUE
  )
  names(df_2) <- c(id, value)
  df_2 <- df_2[order(df_2[[value]], decreasing = decreasing), ]
  df_2$rank <- seq_len(nrow(df_2))

  # Prepare raw data for points with correct rank ordering
  raw_1 <- merge(group_data[[1]], df_1[, c(id, "rank")], by = id)
  raw_2 <- merge(group_data[[2]], df_2[, c(id, "rank")], by = id)

  # Calculate shared x-axis limits if needed
  if (!free_x) {
    x_min <- min(data_common[[value]], na.rm = TRUE)
    x_max <- max(data_common[[value]], na.rm = TRUE)
    x_range <- x_max - x_min
    
    # Bars need to start from 0, boxes don't
    if (style == "bar") {
      x_lim <- c(0, x_max + 0.1 * x_range)
    } else {
      x_lim <- c(x_min - 0.1 * x_range, x_max + 0.1 * x_range)
    }
  }
  
  # Line data
  line_df <- merge(
    df_1[, c(id, "rank")],
    df_2[, c(id, "rank")],
    by = id,
    suffixes = c("_1", "_2")
  )
  
  # Calculate rank change direction
  line_df$rank_change <- ifelse(
    line_df$rank_2 < line_df$rank_1, "increase",
    ifelse(line_df$rank_2 > line_df$rank_1, "decrease", "no_change")
  )
  line_df$rank_change <- factor(
    line_df$rank_change,
    levels = c("increase", "decrease", "no_change")
  )

  n_samples <- nrow(df_1)
  
  # Calculate maximum label width for dynamic middle panel sizing
  max_label_chars <- max(nchar(as.character(line_df[[id]])))
  
  # Left barplot
  p1 <- ggplot2::ggplot(df_1, ggplot2::aes(y = reorder(.data[[id]], .data$rank), x = .data[[value]])) +
    {
      if (style == "bar") {
        ggplot2::geom_bar(stat = "identity", fill = fill[1], alpha = alpha, 
                        color = "black", linewidth = 0.2)
      } else {
        ggplot2::geom_boxplot(
          data = raw_1,
          ggplot2::aes(y = reorder(.data[[id]], .data$rank), x = .data[[value]]),
          fill = fill[1],
          alpha = alpha,
          color = "black",
          linewidth = 0.2,
          # Hide outliers since we show points
          outlier.shape = NA
        )
      }
    } +
    {
      if (show_points) {
        ggplot2::geom_point(
          data = raw_1,
          ggplot2::aes(y = reorder(.data[[id]], .data$rank), x = .data[[value]]),
          color = "black",
          fill = fill[1],
          size = point_size,
          shape = point_shape,
          alpha = point_alpha
        )
      }
    } +
    {
      if (free_x) {
        ggplot2::scale_x_reverse(expand = ggplot2::expansion(mult = c(0.1, 0)))
      } else {
        ggplot2::scale_x_reverse(limits = rev(x_lim), expand = c(0, 0))
      }
    } +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = 0, add = c(0.6, 0.5)),
      position = "right",
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = group_levels[1]) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = element_text(hjust = 1, color = "black"),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.2),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.2),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.2),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Center line plot with colored lines
  p2 <- ggplot2::ggplot(line_df, ggplot2::aes(
    x = 0, y = .data$rank_1, xend = 1, yend = .data$rank_2,
    color = .data$rank_change
  )) +
    ggplot2::geom_segment(alpha = line_alpha, linewidth = line_width) +
    ggplot2::scale_color_manual(
      values = rank_change_colors,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(n_samples + 0.6, 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  
  # Right barplot
  p3 <- ggplot2::ggplot(df_2, ggplot2::aes(y = reorder(.data[[id]], .data$rank), x = .data[[value]])) +
    {
      if (style == "bar") {
        ggplot2::geom_bar(stat = "identity", fill = fill[2], alpha = alpha,
                        color = "black", linewidth = 0.2)
      } else {
        ggplot2::geom_boxplot(
          data = raw_2,
          ggplot2::aes(y = reorder(.data[[id]], .data$rank), x = .data[[value]]),
          fill = fill[2],
          alpha = alpha,
          color = "black",
          linewidth = 0.2,
          outlier.shape = NA
        )
      }
    } +
    {
      if (show_points) {
        ggplot2::geom_point(
          data = raw_2,
          ggplot2::aes(y = reorder(.data[[id]], .data$rank), x = .data[[value]]),
          color = "black",
          fill = fill[2],
          size = point_size,
          shape = point_shape,
          alpha = point_alpha
        )
      }
    } +
   {
      if (free_x) {
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
      } else {
        ggplot2::scale_x_continuous(limits = x_lim, expand = c(0, 0))
      }
    } +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = 0, add = c(0.6, 0.5)),
      position = "left",
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = group_levels[2]) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = element_text(hjust = 0, color = "black"),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.2),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.2),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.2),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Combine with dynamic middle width
  p1 + p2 + p3 + patchwork::plot_layout(widths = c(1, panel_ratio, 1))
}