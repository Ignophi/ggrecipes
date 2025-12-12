#' Confusion/Contingency Table Bubble Plot
#'
#' @description
#' The function automatically computes frequency counts for each unique 
#' combination of the x and y categorical variables using \code{table()}. 
#' Bubble sizes are scaled proportionally to represent counts, with the 
#' range controlled by \code{point_size_range}. Useful for visualizing 
#' cross-tabulations, confusion matrices, or any bivariate categorical data.
#'
#' @param data A data frame containing the categorical variables.
#' @param x Character string specifying the column name in \code{data} for
#'   the x-axis categorical variable.
#' @param y Character string specifying the column name in \code{data} for
#'   the y-axis categorical variable.
#' @param fill Character string specifying the fill color for bubbles.
#'   Default is "skyblue".
#' @param text_size Numeric value specifying the size of count labels.
#'   Default is 4.
#' @param text_color Character string specifying the color of count labels.
#'   Default is "black".
#' @param point_size_range Numeric vector of length 2 specifying the minimum
#'   and maximum bubble sizes. Default is \code{c(3, 15)}.
#' @param border_color Character string specifying the color of bubble borders.
#'   Default is "black".
#' @param show_grid Logical indicating whether to show major grid lines.
#'   Default is TRUE.
#' @param expand Numeric value specifying the expansion multiplier for both
#'   axes. Default is 0.15.
#' @param facet_x Character string specifying an optional column name in
#'   \code{data} for horizontal faceting. Default is NULL (no faceting).
#' @param facet_y Character string specifying an optional column name in
#'   \code{data} for vertical faceting. Default is NULL (no faceting).
#'
#' @return A ggplot2 object showing the confusion table as a bubble plot.
#'   The plot displays:
#'   \itemize{
#'     \item Bubbles at each x-y combination, sized by frequency count
#'     \item Count labels displayed in the center of each bubble
#'     \item Optional faceting for additional categorical variables
#'     \item Customizable colors, sizes, and grid visibility
#'   }
#'
#' @import ggplot2
#' @importFrom stats aggregate reformulate
#'
#' @examples
#' data(mtcars)
#' # Create synthetic categorical variables
#' # Bin horsepower & fuel efficiency into ordered categories
#' mtcars$horsepower <- 
#'   cut(mtcars$hp, breaks = 5, 
#'       labels = c("Very Low", "Low", "Medium", "High", "Very High"))
#' mtcars$`miles per gallon` <- 
#'   cut(mtcars$mpg, breaks = 5,
#'       labels = c("Very Low", "Low", "Medium", "High", "Very High"))
#' 
#' # Base plot
#' gg_conf(data = mtcars, x = "horsepower", y = "miles per gallon")
#' 
#' # Custom styling
#' gg_conf(data = mtcars, x = "horsepower", y = "miles per gallon",
#'         fill = "lightcoral", point_size_range = c(5, 20),
#'         show_grid = FALSE)
#'
#' # With faceting by "vs" column
#' gg_conf(data = mtcars, x = "horsepower", y = "miles per gallon",
#'         fill = "lightcoral", point_size_range = c(5, 20),
#'         facet_x = "vs",
#'         show_grid = FALSE)
#' @export
gg_conf <- function(data,
                    x,
                    y,
                    fill = "skyblue",
                    text_size = 4,
                    text_color = "black",
                    point_size_range = c(3, 15),
                    border_color = "black",
                    show_grid = TRUE,
                    expand = 0.15,
                    facet_x = NULL,
                    facet_y = NULL) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!is.character(x) || length(x) != 1) {
    stop("'x' must be a single character string naming a column in 'data'")
  }
  
  if (!is.character(y) || length(y) != 1) {
    stop("'y' must be a single character string naming a column in 'data'")
  }
  
  required_cols <- c(x, y)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Required column(s) missing from data frame: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  if (!is.null(facet_x) && !facet_x %in% colnames(data)) {
    stop(sprintf("Column '%s' specified in 'facet_x' not found in data", facet_x))
  }
  
  if (!is.null(facet_y) && !facet_y %in% colnames(data)) {
    stop(sprintf("Column '%s' specified in 'facet_y' not found in data", facet_y))
  }
  
  if (!is.numeric(text_size) || text_size <= 0) {
    stop("'text_size' must be a positive number")
  }
  
  if (!is.numeric(point_size_range) || length(point_size_range) != 2) {
    stop("'point_size_range' must be a numeric vector of length 2")
  }
  
  if (point_size_range[1] < 0 || point_size_range[2] <= point_size_range[1]) {
    stop("'point_size_range' must have positive values with [2] > [1]")
  }
  
  if (!is.logical(show_grid)) {
    stop("'show_grid' must be TRUE or FALSE")
  }
  
  if (!is.numeric(expand) || expand < 0) {
    stop("'expand' must be a non-negative number")
  }
  
  # Compute counts
  group_vars <- c(x, y)
  if (!is.null(facet_x)) group_vars <- c(group_vars, facet_x)
  if (!is.null(facet_y)) group_vars <- c(group_vars, facet_y)
  
  count_df <- aggregate(
    rep(1, nrow(data)),
    by = data[group_vars],
    FUN = length
  )
  names(count_df)[ncol(count_df)] <- "count"
  
  # Build plot
  p <- ggplot2::ggplot(count_df, ggplot2::aes(
    x = .data[[x]],
    y = .data[[y]]
  )) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$count),
      shape = 21,
      fill = fill,
      color = border_color
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$count),
      size = text_size,
      color = text_color,
      vjust = 0.5
    ) +
    ggplot2::scale_size(
      range = point_size_range,
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(mult = c(expand, expand))
    ) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = c(expand, expand))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.title.x = ggplot2::element_text(color = "black"),
      axis.title.y = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA)
    ) +
    ggplot2::labs(x = x, y = y)
  
  # Apply grid visibility
  if (!show_grid) {
    p <- p + ggplot2::theme(
      panel.grid.major = ggplot2::element_blank()
    )
  }
  
  # Apply faceting if specified
  if (!is.null(facet_x) && !is.null(facet_y)) {
    p <- p + ggplot2::facet_grid(
      reformulate(facet_x, facet_y)
    )
  } else if (!is.null(facet_x)) {
    p <- p + ggplot2::facet_wrap(
      reformulate(facet_x),
      nrow = 1
    )
  } else if (!is.null(facet_y)) {
    p <- p + ggplot2::facet_wrap(
      reformulate(facet_y),
      ncol = 1
    )
  }
  
  return(p)
}