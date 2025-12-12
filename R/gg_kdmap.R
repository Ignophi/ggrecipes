#' Kinetic Rate Map (Association/Dissociation plot)
#'
#' @description
#' Generates a log-log plot of association rate (ka) vs dissociation rate (kd)
#' with iso-affinity (KD) lines. Useful for visualizing kinetic binding data
#' from surface plasmon resonance (SPR), biolayer interferometry (BLI), or
#' other biophysical assays.
#'
#' @param data A data frame containing kinetic data. Must include columns
#'   specified by the \code{id}, \code{ka}, and \code{kd} parameters.
#' @param id Character string specifying the column name in \code{data} that
#'   contains identifiers for each measurement. Used to group replicates and
#'   identify reference points. Default is "id".
#' @param ka Character string specifying the column name in \code{data} that
#'   contains association rate constants. Values must be in M^-1s^-1 units.
#'   Default is "ka".
#' @param kd Character string specifying the column name in \code{data} that
#'   contains dissociation rate constants. Values must be in s^-1 units.
#'   Default is "kd".
#' @param labels Character string specifying a column name in \code{data}
#'   containing text labels to display next to points. If NULL (default), no
#'   labels are shown.
#' @param size Numeric value for point size, or character string specifying
#'   a column name in \code{data} to map to point size. Default is 3.
#' @param shape Numeric value (0-25) for point shape, or character string
#'   specifying a column name in \code{data} to map to point shape. Default
#'   is 21 (filled circle).
#' @param fill Character string specifying a color for point fill, or a column
#'   name in \code{data} to map to fill aesthetic. Default is "grey".
#' @param color Character string specifying a color for point border, or a
#'   column name in \code{data} to map to color aesthetic. Default is "black".
#' @param ref_id Character string specifying the ID of a reference point to
#'   highlight with different aesthetics. Must match a value in the column
#'   specified by \code{id}. Default is NULL (no reference point).
#' @param ref_shape Numeric value (0-25) specifying the shape for the reference
#'   point. Default is 21.
#' @param ref_color Character string specifying the border color for the
#'   reference point. Default is "black".
#' @param ref_fill Character string specifying the fill color for the reference
#'   point. Default is "white".
#' @param ref_size Numeric value specifying the size for the reference point.
#'   Default is the value of the \code{size} parameter.
#' @param rep_lines Logical indicating whether to connect replicate points
#'   (points with the same ID) with lines. Default is TRUE.
#' @param iso_color Character string specifying the color of iso-KD lines.
#'   Default is "black".
#' @param iso_alpha Numeric value (0-1) specifying the transparency of iso-KD
#'   lines. Default is 1.
#' @param iso_width Numeric value specifying the line width of iso-KD lines.
#'   Default is 0.3.
#' @param iso_type Character string specifying the line type of iso-KD lines.
#'   Must be one of "solid", "dashed", "dotted", "dotdash", "longdash", or
#'   "twodash". Default is "dashed".
#' @param iso_n Numeric value specifying the number of iso-KD lines to draw,
#'   spaced evenly in log10(KD) across the current plot range. Default is 8.
#' @param text_padding Numeric value specifying the padding around text labels
#'   when \code{labels} is used. Passed to \code{ggrepel::geom_text_repel}.
#'   Default is 1.
#' @param title Character string for the plot title. If NULL (default), no
#'   title is displayed.
#' @param show_anno Logical indicating whether to show corner annotations
#'   indicating fast on-rate (top left) and fast off-rate (bottom right)
#'   regions. Default is FALSE.
#'
#' @return A ggplot2 object showing the kinetic rate map. The plot displays:
#'   \itemize{
#'     \item X-axis: dissociation rate (kd) in s^-1 on log scale
#'     \item Y-axis: association rate (ka) in M^-1s^-1 on log scale
#'     \item Points representing individual kinetic measurements
#'     \item Diagonal iso-affinity lines representing constant KD values
#'     \item Secondary axes (top and right) showing KD values in appropriate units (pM, nM, µM, or mM)
#'     \item Optional connecting lines between replicates (same ID)
#'     \item Optional highlighted reference point
#'     \item Optional text labels for points
#'     \item Optional corner annotations indicating rate directions
#'   }
#'   KD values are automatically calculated from ka and kd rates and displayed
#'   in the most appropriate units based on the data range.
#'
#' @details
#' The function creates a log-log plot with:
#' \itemize{
#'   \item X-axis: dissociation rate (kd, in s^-1)
#'   \item Y-axis: association rate (ka, in M^-1s^-1)
#'   \item Diagonal lines representing iso-affinity contours (constant KD values)
#' }
#'
#' \strong{Required units:}
#' \itemize{
#'   \item ka: M^-1s^-1 (association rate constant)
#'   \item kd: s^-1 (dissociation rate constant)
#' }
#'
#' The function calculates equilibrium dissociation constant as KD = kd/ka and
#' automatically converts to the most appropriate units (pM, nM, µM, or mM) for
#' display on secondary axes.
#'
#' Iso-KD lines are generated so that there are always \code{iso_n} lines spanning
#' the visible plot window. They are equally spaced in log10(KD) across the
#' KD range implied by the current x/y limits. The top (x) and right (y)
#' secondary axes are labeled with the KD values corresponding to these lines.
#'
#' @seealso \code{\link[ggplot2]{geom_point}} for point customization,
#'   \code{\link[ggrepel]{geom_text_repel}} for label placement.
#'
#' @examples
#' # Basic example: 5 variants with single measurements
#' kinetic_data <- data.frame(
#'   id = c("WT", "Mut1", "Mut2", "Mut3", "Mut4"),
#'   ka = c(1.2e5, 2.5e5, 2e5, 8.0e4, 1.8e5),
#'   kd = c(1.5e-3, 2.0e-3, 1.5e-3, 1.2e-3, 1.8e-3)
#' )
#' 
#' gg_kdmap(data = kinetic_data, show_anno = TRUE)
#' 
#' # With replicates: lines connect points with same ID
#' kinetic_rep <- data.frame(
#'   id = c("WT", "WT", "WT", "Mut1", "Mut1", "Mut2", "Mut3", "Mut4"),
#'   ka = c(1.2e5, 1.5e5, 1.1e5, 2.5e5, 2.4e5, 2e5, 8.0e4, 1.8e5),
#'   kd = c(1.5e-3, 1.6e-3, 1.4e-3, 2.0e-3, 1.9e-3, 1.5e-3, 1.2e-3, 1.8e-3)
#' )
#' 
#' gg_kdmap(data = kinetic_rep, show_anno = TRUE, fill = "id")
#' 
#' # Add labels and highlight reference
#' gg_kdmap(
#'   data = kinetic_rep,
#'   labels = "id",
#'   ref_id = "WT",
#'   ref_fill = "white",
#'   ref_color = "red",
#'   fill = "id"
#' )
#' 
#' # Customize iso-KD lines
#' gg_kdmap(
#'   data = kinetic_rep,
#'   iso_n = 12,
#'   iso_color = "#7192ad",
#'   iso_type = "solid"
#' )
#' 
#' # Turn off replicate lines
#' gg_kdmap(data = kinetic_rep, rep_lines = FALSE)
#' 
#' # Custom column names
#' custom_data <- data.frame(
#'   sample = c("WT", "Mut1", "Mut2"),
#'   kon = c(1.2e5, 2.5e5, 5.0e4),
#'   koff = c(1.5e-3, 2.0e-3, 5.0e-4)
#' )
#' 
#' gg_kdmap(data = custom_data, id = "sample", ka = "kon", kd = "koff")
#' @import ggplot2
#' @importFrom grDevices col2rgb
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid textGrob gpar unit arrow
#'
#' @export
gg_kdmap <- function(data,
                     id = "id",
                     ka = "ka", 
                     kd = "kd",
                     labels = NULL,
                     size = 4,
                     shape = 21,
                     fill = "grey",
                     color = "black",
                     ref_id = NULL,
                     ref_shape = 21,
                     ref_color = "black",
                     ref_fill = "white",
                     ref_size = size,
                     rep_lines = TRUE,
                     iso_color = "black",
                     iso_alpha = 1,
                     iso_width = 0.3,
                     iso_type = "dashed",
                     iso_n = 8,
                     text_padding = 1,
                     title = NULL,
                     show_anno = FALSE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  required_cols <- c(id, ka, kd)
  missing_required <- setdiff(required_cols, colnames(data))
  if (length(missing_required) > 0) {
    stop(sprintf(
      "Required column(s) missing from data frame: %s",
      paste(missing_required, collapse = ", ")
    ))
  }
  
  if (!is.null(labels) && is.character(labels) && !labels %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in data frame", labels))
  }
  
  if (!is.null(ref_id) && !ref_id %in% data[[id]]) {
    stop(sprintf(
      "Reference ID '%s' not found in data. Available IDs: %s",
      ref_id, paste(unique(data[[id]]), collapse = ", ")
    ))
  }
  
  if (!is.numeric(size) && !size %in% colnames(data)) {
    stop("'size' must be a number or column name in data")
  }
  
  if (is.numeric(shape)) {
    valid_shapes <- c(0:25)
    if (!shape %in% valid_shapes) {
      stop(sprintf("'shape' must be between 0 and 25, got %s", shape))
    }
  } else if (!shape %in% colnames(data)) {
    stop("'shape' must be a valid shape number or column name in data")
  }
  
  if (is.numeric(ref_shape)) {
    valid_shapes <- c(0:25)
    if (!ref_shape %in% valid_shapes) {
      stop(sprintf("'ref_shape' must be between 0 and 25, got %s", ref_shape))
    }
  }
  
  aes_params <- list(fill = fill, color = color)
  for (aes_name in names(aes_params)) {
    aes_val <- aes_params[[aes_name]]
    if (is.character(aes_val) && length(aes_val) == 1) {
      is_valid_color <- tryCatch({
        grDevices::col2rgb(aes_val)
        TRUE
      }, error = function(e) FALSE)
      
      if (!is_valid_color && !aes_val %in% colnames(data)) {
        stop(sprintf("'%s' must be a valid color or column name in data", aes_name))
      }
    }
  }
  
  if (!is.logical(rep_lines)) {
    stop("'rep_lines' must be TRUE or FALSE")
  }
  
  if (!is.numeric(iso_alpha) || iso_alpha < 0 || iso_alpha > 1) {
    stop("'iso_alpha' must be between 0 and 1")
  }
  
  if (!is.numeric(iso_width) || iso_width <= 0) {
    stop("'iso_width' must be a positive number")
  }
  
  valid_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  if (!iso_type %in% valid_linetypes) {
    stop(sprintf("'iso_type' must be one of %s", 
                 paste(shQuote(valid_linetypes), collapse = ", ")))
  }
  
  if (!is.numeric(iso_n) || length(iso_n) != 1 || iso_n < 1) {
    stop("'iso_n' must be a positive number")
  }
  iso_n <- as.integer(iso_n)
  
  if (!is.numeric(text_padding) || text_padding < 0) {
    stop("'text_padding' must be a non-negative number")
  }
  
  if (!is.logical(show_anno)) {
    stop("'show_anno' must be TRUE or FALSE")
  }
  
  point_data <- data
  if (!is.null(ref_id)) {
    point_data <- data[data[[id]] != ref_id, , drop = FALSE]
  }
  
  y_vals <- data[[ka]] / 1e6
  x_range <- range(data[[kd]], na.rm = TRUE)
  y_range <- range(y_vals, na.rm = TRUE)
  
  x_range[1] <- max(x_range[1], .Machine$double.eps)
  y_range[1] <- max(y_range[1], .Machine$double.eps)
  
  x_lim <- x_range * c(0.8, 1.2)
  y_lim <- y_range * c(0.8, 1.2)
  x_lim[1] <- max(x_lim[1], .Machine$double.eps)
  y_lim[1] <- max(y_lim[1], .Machine$double.eps)
  
  KD_min_plot <- (x_lim[1] / y_lim[2]) * 1e3
  KD_max_plot <- (x_lim[2] / y_lim[1]) * 1e3
  KD_min_plot <- max(KD_min_plot, .Machine$double.eps)
  KD_max_plot <- max(KD_max_plot, KD_min_plot * 1.01)
  
  KD_iso <- 10 ^ seq(
    from = log10(KD_min_plot),
    to   = log10(KD_max_plot),
    length.out = iso_n
  )
  
  format_kd_label <- function(k_nM) {
    if (k_nM < 1) {
      paste0(signif(k_nM * 1e3, 3), " pM")
    } else if (k_nM < 1e3) {
      paste0(signif(k_nM, 3), " nM")
    } else if (k_nM < 1e6) {
      paste0(signif(k_nM / 1e3, 3), " \u00b5M")
    } else {
      paste0(signif(k_nM / 1e6, 3), " mM")
    }
  }
  
  KD_iso_labels <- vapply(KD_iso, format_kd_label, character(1))
  
  iso_df <- do.call(rbind, lapply(KD_iso, function(k) {
    data.frame(
      KD_nM = k,
      x = x_lim,
      y = (x_lim * 1e3) / k
    )
  }))
  
  g <- ggplot(point_data, aes(
    x = .data[[kd]],
    y = .data[[ka]] / 1e6
  )) +
    
    geom_line(
      data = iso_df,
      mapping = aes(x = .data[["x"]], y = .data[["y"]], group = factor(.data[["KD_nM"]])),
      inherit.aes = FALSE,
      color = iso_color,
      linetype = iso_type,
      alpha = iso_alpha,
      linewidth = iso_width
    ) +
    
    {
      if (rep_lines) {
        geom_line(
          data = data,
          mapping = aes(
            x = .data[[kd]],
            y = .data[[ka]] / 1e6,
            group = factor(.data[[id]])
          ),
          linewidth = 1,
          inherit.aes = FALSE
        )
      }
    } +

    # Points on top of lines
    {
      point_aes <- aes(
        size = if (!is.null(size) && size %in% colnames(point_data)) .data[[size]] else NULL,
        shape = if (!is.null(shape) && shape %in% colnames(point_data)) as.factor(.data[[shape]]) else NULL,
        fill = if (!is.null(fill) && fill %in% colnames(point_data)) .data[[fill]] else NULL,
        color = if (!is.null(color) && color %in% colnames(point_data)) .data[[color]] else NULL
      )
      
      point_params <- list()
      if (!is.null(size) && !size %in% colnames(point_data)) point_params$size <- size
      if (!is.null(shape) && !shape %in% colnames(point_data)) point_params$shape <- shape
      if (!is.null(fill) && !fill %in% colnames(point_data)) point_params$fill <- fill
      if (!is.null(color) && !color %in% colnames(point_data)) {
        point_params$color <- color
      }
      
      do.call(geom_point, c(list(mapping = point_aes), point_params))
    } +
    
    scale_x_log10(
      expand = c(0, 0),
      sec.axis = sec_axis(
        ~ .,
        breaks = (y_lim[2] * KD_iso) / 1e3,
        labels = KD_iso_labels,
        name = ""
      )
    ) +
    scale_y_log10(
      sec.axis = sec_axis(
        ~ .,
        breaks = (x_lim[2] * 1e3) / KD_iso,
        labels = KD_iso_labels,
        name = ""
      )
    ) +
    
    {
      if (!is.null(shape) && shape %in% colnames(point_data)) {
        scale_shape_manual(
          values = rep(
            c(21, 22, 23, 24, 25),
            length.out = length(unique(point_data[[shape]]))
          )
        )
      }
    } +
    
    guides(
      fill = guide_legend(override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(override.aes = list(fill = "grey"))
    ) +
    
    {
      if (!is.null(labels)) {
        ggrepel::geom_text_repel(
          data = data,
          aes(
            x = .data[[kd]],
            y = .data[[ka]] / 1e6,
            label = .data[[labels]]
          ),
          size = 4,
          box.padding = text_padding,
          point.padding = 1,
          #min.segment.length = 0.2,
          #force = 1,
          max.overlaps = Inf
        )
      }
    } +
    
    coord_cartesian(
      xlim = x_lim,
      ylim = y_lim
    ) +
    
    {
      if (!is.null(ref_id) && ref_id %in% data[[id]]) {
        geom_point(
          data = data[data[[id]] == ref_id, , drop = FALSE],
          mapping = aes(
            x = .data[[kd]],
            y = .data[[ka]] / 1e6
          ),
          size = ref_size,
          shape = ref_shape,
          color = ref_color,
          fill = ref_fill,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
      }
    } +

    {
      if (show_anno) {
        pad_frac_x <- 0.01
        pad_frac_y <- 0.01

        lx_min <- log10(x_lim[1])
        lx_max <- log10(x_lim[2])
        ly_min <- log10(y_lim[1])
        ly_max <- log10(y_lim[2])

        x_left   <- 10^(lx_min + pad_frac_x * (lx_max - lx_min))
        x_right  <- 10^(lx_max - pad_frac_x * (lx_max - lx_min))
        y_bottom <- 10^(ly_min + pad_frac_y * (ly_max - ly_min))
        y_top    <- 10^(ly_max - pad_frac_y * (ly_max - ly_min))

        list(
          annotate(
            "text",
            x = x_right,
            y = y_bottom,
            label = "Fast\nOff-rate",
            hjust = 1, vjust = 0,
            fontface = "bold",
            family = "mono",
            size = 5
          ),
          annotate(
            "text",
            x = x_left,
            y = y_top,
            label = "Fast\nOn-rate",
            hjust = 0, vjust = 1,
            fontface = "bold",
            family = "mono",
            size = 5
          )
        )
      }
    } +
    
    labs(
      x = expression("Kd (s"^-1*")"),
      y = expression("Ka (M"^-1*"s"^-1*")"),
      title = title,
      fill = if (!is.null(fill) && fill %in% colnames(point_data)) fill else NULL,
      color = if (!is.null(color) && color %in% colnames(point_data)) color else NULL,
      shape = if (!is.null(shape) && shape %in% colnames(point_data)) shape else NULL,
      size = if (!is.null(size) && size %in% colnames(point_data)) size else NULL
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      plot.margin = margin(10, 15, 10, 10),
      axis.line.x.bottom = ggplot2::element_line(
        color = "black",
        arrow = grid::arrow(length = grid::unit(0.3, "cm"), type = "closed")
      ),
      axis.line.y.left = ggplot2::element_line(
        color = "black",
        arrow = grid::arrow(length = grid::unit(0.3, "cm"), type = "closed")
      ),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.ticks.x.top = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  return(g)
}