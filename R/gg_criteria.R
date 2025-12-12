#' Criteria Heatmap with Optional Barplots
#'
#' @description
#' Creates a tile-based heatmap visualization where samples are displayed on
#' the y-axis and criteria are shown on the x-axis. Each tile shows a
#' criterion value with color coding. Optional barplots can be added to the
#' right side to display numeric metrics for each sample.
#'
#' @param data A data frame containing the data to visualize.
#' @param id Character string specifying the column name in \code{data} that
#'   contains sample identifiers. Default is "id".
#' @param criteria Character string specifying either: (1) the column name in
#'   \code{data} that contains criterion names (long data, not yet supported),
#'   or (2) a regular expression pattern that identifies the criterion columns
#'   (wide data, e.g., \code{"_criteria$"} to match \code{Pass_criteria},
#'   \code{Fail_criteria}, etc.). Default is "_criteria".
#' @param tile_fill Named character vector of colors for criterion values
#'   (e.g., \code{c(Pass = "green", Fail = "red")}). If NULL (default), uses
#'   the Brewer "Paired" palette.
#' @param bar_column Character vector specifying column name(s) in \code{data}
#'   to display as horizontal barplot(s) to the right of the heatmap. Columns
#'   must be numeric. Default is NULL (no barplots).
#' @param bar_fill Character vector specifying the fill color(s) for bars.
#'   When \code{bar_column} contains multiple elements, colors are recycled
#'   if necessary to match the number of bars. If NULL (default), uses colors
#'   from the Brewer "Paired" palette.
#' @param panel_ratio Numeric value specifying the total relative width of the
#'   barplot panel(s) compared to the heatmap panel (which has a reference
#'   width of 1). When multiple \code{bar_column} values are provided, this
#'   width is divided equally among them. For example, \code{panel_ratio = 0.3}
#'   creates a 1:0.3 width ratio between heatmap and barplot panels. Only used
#'   when \code{bar_column} is not NULL. Default is 0.3.
#' @param tile_width Numeric value (0-1) specifying the width of tiles as a
#'   proportion of available space. Default is 0.7.
#' @param tile_height Numeric value (0-1) specifying the height of tiles as a
#'   proportion of available space. Default is 0.7.
#' @param tile_alpha Numeric value (0-1) specifying transparency of tiles.
#'   Default is 1 (fully opaque).
#' @param show_text Logical indicating whether to show criterion values as
#'   text labels on tiles. Default is TRUE.
#' @param border_color Character string specifying the color of tile borders.
#'   Default is "black".
#' @param border_width Numeric value specifying the width of tile borders.
#'   Default is 0.5.
#' @param text_size Numeric value specifying the size of axis text. Default is 8.
#' @param show_legend Logical indicating whether to show the fill legend.
#'   Default is TRUE.
#' @param quiet Logical indicating whether to suppress messages.
#'   Default is FALSE.
#'
#' @return A ggplot2 object (or patchwork object if \code{bar_column} is used)
#'   showing the criteria heatmap. The plot displays:
#'   \itemize{
#'     \item Tiles at each sample-criterion intersection
#'     \item Tile colors indicating criterion values
#'     \item Criterion values as text labels within tiles
#'     \item Optional horizontal barplots on the right side
#'     \item A "recommended_dims" attribute with suggested width and height in inches
#'   }
#'
#' @details
#' The function currently supports wide format data where each criterion is a
#' separate column. The column names should match the pattern specified in
#' \code{criteria} (e.g., columns ending in "_crit"). The pattern is removed
#' from column names to create criterion labels.
#'
#' When \code{bar_column} is specified, horizontal barplots are added to the
#' right side of the heatmap. Multiple barplots can be created by providing a
#' vector of column names. Each barplot shows the numeric values with text
#' labels positioned outside the bars.
#'
#' The function calculates and suggests plot dimensions based on the number of
#' samples and criteria. Access these via \code{attr(plot, "recommended_dims")}.
#'
#' WARNING: Alignment between heatmap and barplots depends on plot dimensions.
#' The function provides recommended dimensions (accessible via
#' \code{attr(plot, "recommended_dims")}) that ensure proper alignment. You
#' can adjust these dimensions to improve appearance (e.g., reduce width to
#' tighten spacing, or scale proportionally for size) while maintaining
#' alignment.
#'
#' @import ggplot2
#' @importFrom patchwork plot_layout
#'
#' @examples
#' # Example: Gene prioritization criteria
#' gene_data <- data.frame(
#'   gene = c("BRCA1", "TP53", "EGFR", "KRAS", "MYC", 
#'            "PTEN", "APC", "CDKN2A", "RB1", "VHL"),
#'   `Missense Variant_crit` = c("Yes", "Yes", "Yes", NA, "Yes", 
#'                                "Yes", NA, "Yes", NA, "Yes"),
#'   `eQTL_crit` = c("Yes", "Yes", NA, "Yes", "Yes", 
#'                   "Yes", "Yes", "Yes", "Yes", NA),
#'   `pQTL_crit` = c("Yes", NA, "Yes", "Yes", NA, 
#'                   "Yes", "Yes", NA, "Yes", "Yes"),
#'   `GWAS Hit_crit` = c("Yes", "Yes", "Yes", "Yes", NA, 
#'                       "Yes", "Yes", "Yes", NA, NA),
#'   `Loss of Function_crit` = c(NA, "Yes", NA, NA, "Yes", 
#'                               "Yes", "Yes", NA, "Yes", NA),
#'   `High Conservation_crit` = c("Yes", "Yes", "Yes", "Yes", "Yes", 
#'                                "Yes", "Yes", "Yes", "Yes", "Yes"),
#'   `mRNA DE_crit` = c("Yes", NA, "Yes", NA, NA,
#'                      "Yes", "Yes", "Yes", NA, "Yes"),
#'   `Prot DE_crit` = c(NA, "Yes", NA, NA, NA,
#'                      "Yes", NA, NA, NA, "Yes"),
#'   check.names = FALSE
#' )
#' 
#' # Calculate total criteria met
#' crit_cols <- grep("_crit$", names(gene_data), value = TRUE)
#' gene_data$`Total` <- rowSums(gene_data[crit_cols] == "Yes", na.rm = TRUE)
#' 
#' # Base criteria plot
#' # **Warn**: space between heatmap and barplots depends on plot width
#' gg_criteria(
#'   data = gene_data,
#'   id = "gene",
#'   criteria = "_crit$",
#'   bar_column = "Total",
#'   show_text = FALSE,
#'   tile_fill = c(Yes = "#A6CEE3", No = "white"),
#'   bar_fill = "#A6CEE3",
#'   panel_ratio = 1
#' )
#' 
#' # Example: VHH Variant Analysis
#' # Define amino acid chemistry colors
#' aa_colors <- c(
#'   "D" = "#E60A0A", "E" = "#E60A0A", # Acidic (red)
#'   "K" = "#145AFF", "R" = "#145AFF", # Basic (blue)
#'   "H" = "#8282D2",                  # Histidine (purple)
#'   "S" = "#FA9600", "T" = "#FA9600", # Polar uncharged (orange)
#'   "N" = "#00DCDC", "Q" = "#00DCDC", # Polar amides (cyan)
#'   "C" = "#E6E600",                  # Cysteine (yellow)
#'   "G" = "#EBEBEB",                  # Glycine (light gray)
#'   "P" = "#DC9682",                  # Proline (tan)
#'   "A" = "#C8C8C8",                  # Alanine (gray)
#'   "V" = "#0F820F", "I" = "#0F820F", # Hydrophobic (green)
#'   "L" = "#0F820F", "M" = "#0F820F",
#'   "F" = "#3232AA", "W" = "#B45AB4", # Aromatic (dark blue/purple)
#'   "Y" = "#3232AA"
#' )
#' 
#' vhh_variants <- data.frame(
#'   variant = c("WT", "Mut1", "Mut2", "Mut3", "Mut4", "Mut7", "Mut5",
#'               "Mut6", "Mut8", "Mut9", "Mut10", "Mut11"),
#'   Q5_mut = c(NA, "H", NA, NA, NA, NA, "H", "H", "H", "D", NA, "H"),
#'   S55_mut = c(NA, NA, "P", NA, NA, "P", "P", NA, "P", "P", NA, NA),
#'   N73_mut = c(NA, NA, NA, "E", NA, NA, NA, "E", "E", NA, "E", NA),
#'   K80_mut = c(NA, "L", NA, NA, "S", "V", NA, NA, NA, "L", "S", NA),
#'   F99_mut = c(NA, NA, "L", NA, NA, NA, NA, NA, NA, "W", NA, "W"),
#'   KD_nM = c(45, 18, 5.2, 38, 42, 20, 3.8, 15, 3.2, 4.5, 40, 22),
#'   yield_mg_L = c(12, 11.8, 10, 13, 11, 10, 10, 12, 10, 7.8, 12.5, 8.5),
#'   Tm_C = c(68.5, 67.8, 68, 72.3, 35, 66, 67.5, 70, 70.5, 72, 38, 74)
#' )
#' 
#' # Create the systematic variant heatmap
#' gg_criteria(
#'   data = vhh_variants,
#'   id = "variant",
#'   criteria = "_mut$",
#'   tile_fill = aa_colors,
#'   bar_column = c("KD_nM", "yield_mg_L", "Tm_C"),
#'   panel_ratio = 2,
#'   tile_width = 0.70,
#'   tile_height = 0.70,
#'   show_text = TRUE,
#'   border_color = "grey40",
#'   border_width = 0.4,
#'   text_size = 10,
#'   show_legend = FALSE
#' )
#' @export
gg_criteria <- function(data,
                        id = "id",
                        criteria = "_criteria",
                        bar_column = NULL,
                        bar_fill = NULL,
                        panel_ratio = 0.3,
                        tile_width = 0.7,
                        tile_height = 0.7,
                        tile_alpha = 1,
                        tile_fill = NULL,
                        show_text = TRUE,
                        border_color = "black",
                        border_width = 0.5,
                        text_size = 8,
                        show_legend = TRUE,
                        quiet = FALSE) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.null(tile_fill) && !is.character(tile_fill)) {
    stop("'tile_fill' must be a character vector or NULL")
  }
  
  if (!is.logical(show_legend)) {
    stop("'show_legend' must be TRUE or FALSE")
  }
  
  if (!is.logical(quiet)) {
    stop("'quiet' must be TRUE or FALSE")
  }
  
  if (!is.character(id) || length(id) != 1) {
    stop("'id' must be a single character string naming a column in 'data'")
  }
  
  if (!is.character(criteria) || length(criteria) != 1) {
    stop("'criteria' must be a single character string")
  }
  
  if (!is.null(bar_fill)) {
    if (!is.character(bar_fill) || length(bar_fill) < 1) {
      stop("'bar_fill' must be a character vector with at least one color or NULL")
    }
  }
  
  if (!is.numeric(tile_width) || tile_width <= 0 || tile_width > 1) {
    stop("'tile_width' must be between 0 and 1")
  }
  
  if (!is.numeric(tile_height) || tile_height <= 0 || tile_height > 1) {
    stop("'tile_height' must be between 0 and 1")
  }
  
  if (!is.numeric(tile_alpha) || tile_alpha < 0 || tile_alpha > 1) {
    stop("'tile_alpha' must be between 0 and 1")
  }

  if (!is.logical(show_text)) {
    stop("'show_text' must be TRUE or FALSE")
  }
  
  if (!is.numeric(border_width) || border_width < 0) {
    stop("'border_width' must be a non-negative number")
  }
  
  if (!is.numeric(text_size) || text_size <= 0) {
    stop("'text_size' must be a positive number")
  }
  
  # Store original sample order
  original_id_order <- if (id %in% names(data)) unique(data[[id]]) else NULL
  
  # Store original data for barplot
  original_data <- data
  
  # Reshape wide to long if needed
  if (!criteria %in% names(data)) {
    measure_cols <- grep(criteria, names(data), value = TRUE)
    
    if (length(measure_cols) == 0) {
      stop(sprintf(
        "No columns in 'data' match the pattern specified in 'criteria': %s",
        shQuote(criteria)
      ))
    }
    
    base_cols <- setdiff(names(data), measure_cols)
    long_list <- vector("list", length(measure_cols))
    
    for (i in seq_along(measure_cols)) {
      col_name <- measure_cols[i]
      tmp <- data[base_cols]
      criterion_name <- sub(criteria, "", col_name)
      tmp$.criterion <- criterion_name
      tmp$.value <- data[[col_name]]
      long_list[[i]] <- tmp
    }
    
    data <- do.call(rbind, long_list)
    rownames(data) <- NULL
    
    criterion_levels <- sub(criteria, "", measure_cols)
    data$.criterion <- factor(data$.criterion, levels = criterion_levels)
    criteria_col <- ".criterion"
    value_col <- ".value"
  } else {
    stop("Long format not yet implemented - use wide format with regex pattern")
  }
  
  # Validate bar_column
  if (!is.null(bar_column)) {
    if (!is.character(bar_column) || length(bar_column) < 1) {
      stop("'bar_column' must be a character vector or NULL")
    }
    
    missing_cols <- setdiff(bar_column, colnames(original_data))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "Column(s) specified in 'bar_column' not found in data: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }
    
    non_numeric <- bar_column[!sapply(original_data[bar_column], is.numeric)]
    if (length(non_numeric) > 0) {
      stop(sprintf(
        "Column(s) in 'bar_column' must be numeric: %s",
        paste(non_numeric, collapse = ", ")
      ))
    }
    
    if (!is.numeric(panel_ratio) || panel_ratio <= 0) {
      stop("'panel_ratio' must be a positive number")
    }
  }
  
  # Remove NA/empty values
  data_plot <- data[!is.na(data[[value_col]]) & data[[value_col]] != "", ]
  
  if (nrow(data_plot) == 0) {
    warning("No non-NA values to plot")
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
            ggplot2::labs(title = "No data to display"))
  }
  
  # Set factor levels and add numeric positions
  data_plot[[value_col]] <- as.factor(data_plot[[value_col]])
  data_plot[[id]] <- factor(data_plot[[id]], levels = rev(original_id_order))
  data_plot$.y_pos <- as.numeric(data_plot[[id]])
  
  # Build heatmap plot with numeric y-axis
  p <- ggplot2::ggplot(data_plot, ggplot2::aes(
    x = .data[[criteria_col]],
    y = .data[[".y_pos"]],
    fill = .data[[value_col]]
  ))
  
  p <- p + 
    ggplot2::geom_tile(
      width = tile_width,
      height = tile_height,
      color = border_color,
      linewidth = border_width,
      alpha = tile_alpha
    )

  if (show_text) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[[value_col]]),
      size = 3,
      color = "black"
    )
  }

  p <- p +
    ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
    ggplot2::scale_y_continuous(
      breaks = seq_along(rev(original_id_order)),
      labels = rev(original_id_order),
      limits = c(0.5, length(rev(original_id_order)) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 0, vjust = 0, colour = "black"
      ),
      axis.text.y = ggplot2::element_text(colour = "black", vjust = 0.5),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(size = text_size),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "black", linewidth = 0.1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "left"
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL)
    
  # Apply fill colors
  if (!is.null(tile_fill)) {
    p <- p + ggplot2::scale_fill_manual(values = tile_fill)
  } else {
    p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
  }
  
  # Hide legend if requested
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  # Calculate recommended dimensions (inches)
  n_rows <- length(unique(data_plot[[id]]))
  n_cols <- length(unique(data_plot[[criteria_col]]))
  
  heatmap_height <- n_rows * 0.4 + 1.5
  heatmap_width <- n_cols * 0.4 + 2
  
  # Create barplot(s) if requested
  if (!is.null(bar_column)) {
    # Define default colors for barplots
    paired_colors <- c(
      "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
      "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
    )

    bar_plots <- lapply(seq_along(bar_column), function(i) {
      col <- bar_column[i]
      
      # Use bar_fill if provided, otherwise use Paired palette
      if (is.null(bar_fill)) {
        bar_color <- paired_colors[((i - 1) %% length(paired_colors)) + 1]
      } else {
        bar_color <- bar_fill[((i - 1) %% length(bar_fill)) + 1]
      }
      
      bar_data <- original_data[, c(id, col), drop = FALSE]
      bar_data[[id]] <- factor(bar_data[[id]], levels = rev(original_id_order))
      
      p_bar <- ggplot2::ggplot(bar_data, ggplot2::aes(
        x = .data[[col]],
        y = .data[[id]]
      )) +
        ggplot2::geom_bar(
          stat = "identity",
          fill = bar_color,
          width = tile_height,
          color = "black",
          linewidth = 0.3
        ) +
        ggplot2::geom_text(
          ggplot2::aes(label = round(.data[[col]], 1)),
          hjust = -0.1,
          size = 3,
          color = "black"
        ) +
        ggplot2::scale_x_continuous(
          expand = c(0, 0, 0.15, 0),
          position = "top"
        ) +
        ggplot2::scale_y_discrete(
          expand = ggplot2::expansion(mult = 0, add = c(0.5, 0.5)),
          drop = FALSE
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(colour = "black"),
          axis.text = ggplot2::element_text(size = text_size),
          axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.25),
          axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.25),
          panel.grid.major = ggplot2::element_line(
            color = "grey90", linewidth = 0.3
          ),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          plot.background = ggplot2::element_rect(fill = "white", color = NA),
          axis.text.x.top = ggplot2::element_text(vjust = 0, colour = "black"),
          axis.title.x.top = ggplot2::element_text(),
          axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25)
        ) +
        ggplot2::labs(x = col, y = NULL)
      
      p_bar
    })
    
    # Calculate widths - heatmap gets 1, bars share panel_ratio
    individual_panel_ratio <- panel_ratio / length(bar_column)
    widths <- c(1, rep(individual_panel_ratio, length(bar_column)))
    
    # Combine plots
    combined <- Reduce(`+`, c(list(p), bar_plots)) + 
      patchwork::plot_layout(widths = widths)
    
    # Calculate total width
    bar_total_width <- length(bar_column) * 2
    total_width <- heatmap_width + bar_total_width
    
    attr(combined, "recommended_dims") <- c(
      width = total_width,
      height = heatmap_height
    )
    
    if (!quiet) {
      message(sprintf(
        "Recommended dimensions: %.1f x %.1f inches",
        total_width,
        heatmap_height
      ))
    }
    
    return(combined)
  }
  
  # Return heatmap only
  attr(p, "recommended_dims") <- c(width = heatmap_width, height = heatmap_height)
  p <- p + 
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0, unit = "pt")
    )
  
  if (!quiet) {
    message(sprintf(
      "Recommended dimensions: %.1f x %.1f inches",
      heatmap_width,
      heatmap_height
    ))
  }
  
  return(p)
}