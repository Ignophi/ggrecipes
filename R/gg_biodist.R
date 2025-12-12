#' Biodistribution Barplot with Optional Free-Scale Faceting
#'
#' @description
#' Creates a barplot visualization of biodistribution data (e.g., %ID/g across
#' organs) with optional separation of specific organs onto free y-scales to
#' prevent squishing of lower values. Points are overlaid on bars and all
#' facets are displayed in a single row.
#'
#' The function accepts:
#' \itemize{
#'   \item Long data: \code{value} is the name of the numeric column containing
#'     the measurements (original behaviour).
#'   \item Wide data: \code{value} is a regular expression pattern that matches
#'     the measurement columns (e.g. \code{"_val$"} for \code{Blood_val},
#'     \code{Liver_val}, etc.). In this case, the data are internally converted
#'     to long format, creating a column named by \code{id} for organ/tissue and
#'     a column named by \code{value} for the measurements.
#' }
#'
#' @param data A data frame containing biodistribution measurements.
#' @param id Character string specifying the column name in \code{data} that
#'   contains organ/tissue identifiers (in long format) or the name that will
#'   be used for the organ/tissue column created from wide data. Default is
#'   "id".
#' @param value Character string specifying either: (1) the column name in
#'   \code{data} that contains the measurement values (long data), or (2) a
#'   regular expression pattern that identifies the measurement columns (wide
#'   data, e.g., \code{"_val$"} to match \code{Blood_val}, \code{Liver_val},
#'   etc.). Default is "value".
#' @param y_label Character string specifying the label to use for the y-axis.
#'   Default is the same as \code{value}.
#' @param group Character string specifying the column name in \code{data} that
#'   defines groups (e.g., treated vs control). Default is NULL (no grouping).
#'   When provided, bars and points are dodged and colored by group.
#' @param separate Character vector of organ/tissue names (matching values in
#'   the \code{id} column) to display on separate y-axis scales. This prevents
#'   high-uptake organs from compressing visualization of lower-uptake organs.
#'   Default is NULL (no separation).
#' @param fill_colors Character vector of colors used to fill bars and points.
#'   If \code{group} is NULL, the first color is used (default \code{"#92b9de"}
#'   if \code{fill_colors} is NULL). If \code{group} is provided, one color per
#'   group level is used (in the order of factor levels). If \code{fill_colors}
#'   is NULL in the grouped case, a default brewer palette (\code{"Paired"}) is
#'   used.
#' @param bar_alpha Numeric value (0-1) specifying transparency of bars.
#'   Default is 0.7.
#' @param point_size Numeric value specifying size of points. Default is 1.5.
#'   If set to 0, points are not drawn.
#' @param stat_summary Character string specifying summary statistic for bars.
#'   One of \code{"mean"} or \code{"median"}. Default is \code{"mean"}. Points
#'   show individual values.
#' @param error_bars Logical indicating whether to add error bars (mean ± SD).
#'   Default is FALSE.
#' @param quiet Logical indicating whether to suppress messages during data
#'   processing (e.g., when coercing non-numeric columns to numeric). Default
#'   is FALSE.
#'
#' @return A ggplot2 object showing the biodistribution as barplot with points.
#'   The plot displays:
#'   \itemize{
#'     \item Bars representing summary statistics (mean or median) for each organ
#'     \item Individual data points overlaid on bars
#'     \item Optional grouping with dodged bars and colored by group
#'     \item Optional separation of high-uptake organs onto independent y-scales
#'     \item Optional error bars showing standard deviation
#'   }
#'
#' @details
#' The function automatically handles both long and wide format data:
#' \itemize{
#'   \item Long format: Each row represents one measurement with columns for
#'     organ ID and value
#'   \item Wide format: Each row represents one sample/replicate with separate
#'     columns for each organ (detected via regex pattern)
#' }
#'
#' When \code{separate} is specified, high-uptake organs are displayed in
#' separate facets with independent y-axis scales, preventing compression of
#' lower values in the main plot.
#'
#' @import ggplot2
#' @importFrom stats aggregate median sd
#'
#' @examples
#' bio_data <- data.frame(
#'   id        = paste0("sample_", 1:6),
#'   condition = rep(c("Control", "Treated"), each = 3),
#'   replicate = rep(1:3, times = 2),
#' 
#'   Blood_val  = c(4.8, 5.2, 4.5, 4.1, 4.3, 4.0),
#'   Heart_val  = c(1.9, 2.1, 2.0, 1.6, 1.8, 1.7),
#'   Lung_val   = c(3.5, 3.8, 3.2, 3.0, 3.1, 2.9),
#'   Liver_val  = c(14.2, 15.1, 13.8, 11.5, 12.0, 11.2),
#'   Spleen_val = c(9.1, 8.7, 9.4, 7.2, 7.5, 7.0),
#'   Kidney_val = c(125.0, 112.8, 121.9, 111.1, 102.4, 103.0),
#'   Tumor_val  = c(22.5, 24.1, 23.3, 28.2, 29.5, 27.8),
#'   Muscle_val = c(0.7, 0.6, 0.8, 0.5, 0.4, 0.6),
#'   Bone_val   = c(1.4, 1.6, 1.5, 1.1, 1.2, 1.0)
#' )
#' 
#' # Base biodist plot
#' gg_biodist(bio_data, id = "organ",
#'            value = "_val", group = "condition",
#'            point_size = 1.25,
#'            y_label = "%ID/g")
#' 
#' # Separate high uptake organs on separate axis
#' gg_biodist(bio_data, id = "organ",
#'            value = "_val", group = "condition",
#'            point_size = 1.25,
#'            y_label = "%ID/g",
#'            separate = c("Tumor", "Kidney"))
#'
#' # Customization
#' gg_biodist(bio_data, id = "organ",
#'            value = "_val", group = "condition",
#'            point_size = 0, error_bars = TRUE,
#'            fill_colors = c("#e41a1c", "#377eb8"),
#'            y_label = "%ID/g",
#'            separate = c("Tumor", "Kidney"))
#' @export
gg_biodist <- function(data,
                       id = "id",
                       value = "value",
                       y_label = value,
                       group = NULL,
                       separate = NULL,
                       fill_colors = NULL,
                       bar_alpha = 0.7,
                       point_size = 1.5,
                       stat_summary = "mean",
                       error_bars = FALSE,
                       quiet = FALSE) {
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  # Reshape wide to long if needed
  if (!value %in% names(data)) {
    measure_cols <- grep(value, names(data), value = TRUE)
    
    if (length(measure_cols) == 0) {
      stop(sprintf(
        "No columns in 'data' match the pattern specified in 'value': %s",
        shQuote(value)
      ))
    }
    
    base_cols <- setdiff(names(data), measure_cols)
    long_list <- vector("list", length(measure_cols))
    
    for (i in seq_along(measure_cols)) {
      col_name <- measure_cols[i]
      tmp <- data[base_cols]
      organ_name <- sub(value, "", col_name)
      tmp[[id]] <- organ_name
      tmp[[value]] <- data[[col_name]]
      long_list[[i]] <- tmp
    }
    
    data <- do.call(rbind, long_list)
    rownames(data) <- NULL
    
    organ_levels <- sub(value, "", measure_cols)
    data[[id]] <- factor(data[[id]], levels = unique(organ_levels))
  }
  
  # Validate required columns
  required_cols <- c(id, value)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Required column(s) missing from data frame: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  if (!is.numeric(data[[value]])) {
    if (!quiet) {  # Add this condition
      message(sprintf(
        "Column '%s' is not numeric; attempting to coerce to numeric with as.numeric().",
        value
      ))
    }

    original <- data[[value]]
    
    if (is.factor(original) || is.character(original)) {
      data[[value]] <- suppressWarnings(as.numeric(as.character(original)))
    } else {
      data[[value]] <- suppressWarnings(as.numeric(original))
    }
    
    if (all(is.na(data[[value]]))) {
      stop(sprintf(
        "Coercion of column '%s' to numeric failed: all values are NA after as.numeric().",
        value
      ))
    }
    
    n_na <- sum(is.na(data[[value]]))
    if (n_na > sum(is.na(original))) {
      warning(sprintf(
        "Coercion of column '%s' to numeric introduced %d additional NA value(s).",
        value,
        n_na - sum(is.na(original))
      ))
    }
  }
  
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1) {
      stop("'group' must be a single character string naming a column in 'data'")
    }
    if (!group %in% names(data)) {
      stop(sprintf("Column '%s' specified in 'group' not found in data", group))
    }
  }
  
  if (!is.null(separate)) {
    if (!is.character(separate)) {
      stop("'separate' must be a character vector or NULL")
    }
    missing_organs <- setdiff(separate, data[[id]])
    if (length(missing_organs) > 0) {
      stop(sprintf(
        "Organ(s) specified in 'separate' not found in data: %s",
        paste(missing_organs, collapse = ", ")
      ))
    }
  }
  
  valid_stats <- c("mean", "median")
  if (!stat_summary %in% valid_stats) {
    stop(sprintf("'stat_summary' must be one of %s",
                 paste(shQuote(valid_stats), collapse = ", ")))
  }
  
  if (!is.numeric(bar_alpha) || bar_alpha < 0 || bar_alpha > 1) {
    stop("'bar_alpha' must be between 0 and 1")
  }
  
  if (!is.numeric(point_size) || point_size < 0) {
    stop("'point_size' must be a non-negative number")
  }
  
  if (!is.logical(error_bars)) {
    stop("'error_bars' must be TRUE or FALSE")
  }
  
  # Prepare faceting variable
  data$facet_group <- ifelse(
    !is.null(separate) & data[[id]] %in% separate,
    as.character(data[[id]]),
    "Main"
  )
  
  # Prepare grouping variable
  if (!is.null(group)) {
    data$.group_var <- as.factor(data[[group]])
    group_levels <- levels(data$.group_var)
    n_groups <- length(group_levels)
    
    if (!is.null(fill_colors)) {
      if (!is.character(fill_colors)) {
        stop("'fill_colors' must be a character vector of colors")
      }
      if (length(fill_colors) < n_groups) {
        stop(sprintf(
          "'fill_colors' must have at least as many colors as there are groups (%d)",
          n_groups
        ))
      }
    }
  } else {
    if (!is.null(fill_colors) && length(fill_colors) < 1) {
      stop("'fill_colors' must have at least one color when 'group' is NULL")
    }
  }
  
  # Compute summary statistics
  stat_fun <- if (stat_summary == "mean") mean else median
  
  if (is.null(group)) {
    summary_df <- aggregate(
      data[[value]],
      by = list(id = data[[id]], facet_group = data$facet_group),
      FUN = stat_fun,
      na.rm = TRUE
    )
    names(summary_df)[3] <- "stat_value"
    
    if (error_bars) {
      sd_df <- aggregate(
        data[[value]],
        by = list(id = data[[id]], facet_group = data$facet_group),
        FUN = sd,
        na.rm = TRUE
      )
      names(sd_df)[3] <- "sd_value"
      summary_df <- merge(summary_df, sd_df, by = c("id", "facet_group"))
    }
  } else {
    summary_df <- aggregate(
      data[[value]],
      by = list(
        id = data[[id]],
        facet_group = data$facet_group,
        .group_var = data$.group_var
      ),
      FUN = stat_fun,
      na.rm = TRUE
    )
    names(summary_df)[4] <- "stat_value"
    
    if (error_bars) {
      sd_df <- aggregate(
        data[[value]],
        by = list(
          id = data[[id]],
          facet_group = data$facet_group,
          .group_var = data$.group_var
        ),
        FUN = sd,
        na.rm = TRUE
      )
      names(sd_df)[4] <- "sd_value"
      summary_df <- merge(
        summary_df,
        sd_df,
        by = c("id", "facet_group", ".group_var")
      )
    }
  }
  
  # Set facet levels
  facet_levels <- c("Main", separate)
  data$facet_group <- factor(data$facet_group, levels = facet_levels)
  summary_df$facet_group <- factor(summary_df$facet_group, levels = facet_levels)
  
  # Build plot
  p <- ggplot2::ggplot()
  
  if (is.null(group)) {
    fill_col <- if (is.null(fill_colors)) "#92b9de" else fill_colors[1]
    
    p <- p +
      ggplot2::geom_bar(
        data = summary_df,
        ggplot2::aes(x = .data$id, y = .data$stat_value),
        stat = "identity",
        fill = fill_col,
        color = "black",
        alpha = bar_alpha
      )
    
    if (point_size > 0) {
      p <- p +
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(x = .data[[id]], y = .data[[value]]),
          color = "black",
          fill = fill_col,
          size = point_size,
          shape = 21,
          alpha = 0.7
        )
    }
  } else {
    dodge_width <- 0.9
    
    p <- p +
      ggplot2::geom_bar(
        data = summary_df,
        ggplot2::aes(
          x = .data$id,
          y = .data$stat_value,
          fill = .data$.group_var
        ),
        stat = "identity",
        color = "black",
        alpha = bar_alpha,
        position = ggplot2::position_dodge(width = dodge_width)
      )
    
    if (point_size > 0) {
      p <- p +
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(
            x = .data[[id]],
            y = .data[[value]],
            fill = .data$.group_var
          ),
          color = "black",
          size = point_size,
          shape = 21,
          alpha = 0.7,
          position = ggplot2::position_jitterdodge(
            jitter.width = 0.1,
            jitter.height = 0,
            dodge.width = dodge_width
          )
        )
    }
    
    p <- p + ggplot2::labs(fill = group)
  }
  
  # Add error bars
  if (error_bars) {
    if (is.null(group)) {
      p <- p + ggplot2::geom_errorbar(
        data = summary_df,
        ggplot2::aes(
          x = .data$id,
          ymin = .data$stat_value - .data$sd_value,
          ymax = .data$stat_value + .data$sd_value
        ),
        width = 0.3
      )
    } else {
      dodge_width <- 0.9
      p <- p + ggplot2::geom_errorbar(
        data = summary_df,
        ggplot2::aes(
          x = .data$id,
          ymin = .data$stat_value - .data$sd_value,
          ymax = .data$stat_value + .data$sd_value,
          group = .data$.group_var
        ),
        width = 0.3,
        position = ggplot2::position_dodge(width = dodge_width)
      )
    }
  }
  
  # Apply faceting and theme
  p <- p +
    ggplot2::facet_wrap(
      ~ facet_group,
      nrow = 1,
      scales = "free",
      space = "free_x",
      strip.position = "top"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1, colour = "black"
      ),
      axis.text.y = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5)
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_label
    )
  
  # Apply fill scales
  if (!is.null(group)) {
    if (is.null(fill_colors)) {
      p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
    } else {
      values <- fill_colors[seq_along(group_levels)]
      names(values) <- group_levels
      p <- p + ggplot2::scale_fill_manual(values = values)
    }
  }
  
  return(p)
}