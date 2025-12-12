#' Sequence Coverage Plot
#'
#' @description
#' Plots sequences that are substrings of a reference sequence, with each
#' unique sequence shown as a row at its aligned position. Useful for
#' visualizing peptide mapping coverage, proteomics experiments, or any
#' analysis where you need to show which parts of a reference sequence
#' are covered by shorter sequences. Supports character coloring and
#' region highlighting (e.g., CDRs, tags, binding sites).
#'
#' @param data A data frame containing sequences.
#' @param ref Character string of the reference sequence against which sequences
#'   will be aligned.
#' @param sequence Character string specifying the column name in \code{data} that
#'   contains the sequences. Default is "sequence".
#' @param name Character string specifying the column name in \code{data} that
#'   contains names/IDs for the sequences. If NULL (default), no names are shown on
#'   the y-axis. If specified, these names will be displayed instead of blank labels.
#'   Note: If the same sequence appears multiple times in \code{data} with different
#'   names, only the first occurrence's name will be used in the plot.
#' @param color Named character vector of colors for specific characters, or NULL
#'   (default). If NULL, all characters are displayed in black. Characters specified
#'   here are automatically displayed in bold. Example: \code{c(K = "blue", R = "red")}.
#' @param highlight Named list where names are valid R colors and values are numeric
#'   vectors of positions to highlight with vertical shading. The shading spans the
#'   full height of the plot with semi-transparent color (alpha = 0.3) and black
#'   borders. Consecutive positions are automatically merged into continuous bands.
#'   To specify multiple regions with the same color, provide them as a single vector:
#'   \code{list("#FFE0B2" = c(1:10, 50:60))}. Default is NULL (no highlighting).
#' @param wrap Numeric value specifying the maximum number of positions to display
#'   per row before wrapping to the next panel. If NULL (default), no wrapping is
#'   applied and the entire sequence is displayed in a single row. When specified,
#'   the plot is split into multiple vertically stacked panels, each showing up to
#'   \code{wrap} positions. Useful for displaying long sequences in a more compact
#'   format.
#' @param size Numeric value for the size of sequence characters. Default is 3.
#' @param face Character string specifying font face for characters. Must be one of
#'   "plain", "bold", "italic", or "bold.italic". Default is "plain". Note: Characters
#'   specified in \code{color} are always displayed in bold regardless of this setting.
#' @param show_ref Logical indicating whether to show the reference sequence
#'   in the plot. Default is TRUE.
#' @param margin_t Numeric value for the top margin in points. Increases space above
#'   the plot to prevent annotation clipping. Adjust this value upward if rotated
#'   annotations are cut off. Default is 30.
#' @param annotate List of annotation specifications to add text labels above the plot.
#'   Each element must be a list with required elements \code{label} (character string)
#'   and \code{pos} (numeric position). Optional elements: \code{angle} (rotation angle),
#'   \code{hjust} (horizontal justification, defaults to 0 for 90-degree angles and 0.5
#'   for other angles), \code{vjust} (vertical justification), \code{size} (text size),
#'   \code{face} (font face), \code{color} (text color). Default is NULL (no annotations).
#' @param annotate_defaults List of default values for annotation styling. Valid
#'   elements: \code{angle} (default 0), \code{vjust} (default 0.5), \code{size}
#'   (default 3), \code{face} (default "plain"), \code{color} (default "black").
#'   These defaults are used when individual annotations don't specify these parameters.
#'   Default is \code{list(angle = 0, vjust = 0.5, size = 3, face = "plain", color = "black")}.
#'
#' @return A ggplot2 object showing the aligned sequences. The plot displays:
#'   \itemize{
#'     \item Character sequences aligned horizontally by their position in the reference
#'     \item Each sequence on a separate row, ordered by starting position
#'     \item Optional reference sequence at the top (if \code{show_ref = TRUE})
#'     \item Specified characters highlighted in bold with custom colors
#'     \item Optional vertical shading bands at specified positions with black borders
#'     \item Optional text annotations above the plot
#'   }
#'   If no sequences from \code{data} match the reference sequence, returns an empty
#'   ggplot2 object with a warning message. Sequences that don't match the reference
#'   are silently excluded from the plot without warning.
#'
#' @seealso \code{\link{gg_seqdiff}} for visualizing sequence differences relative
#'   to a reference.
#'
#' @import ggplot2
#' @import scales
#' @importFrom grDevices col2rgb
#' @importFrom utils modifyList
#'
#' @examples
#' # Create synthetic example of peptide mapping data
#' # Reference sequence
#' ref_seq <- paste0(
#'   "QVQLVESGGGLVQAGGSLRLSCAASGFTFSSYAMGWFRQAPGKEREFVAAINSGGST",
#'   "YYPDSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADLRGTTVNNYWGQGTQV",
#'   "TVSSEQKLISEEDL"
#' )
#' 
#' # Peptides with RT and intensity
#' df_peptides <- data.frame(
#'   id = c("Pep_1004", "Pep_1010", "Pep_1007", "Pep_1011", 
#'          "Pep_1009", "Pep_1005", "Pep_1013", "Pep_1003", 
#'          "Pep_1001", "Pep_1012", "Pep_1006", "Pep_1008", 
#'          "Pep_1002"),
#'   sequence = c(
#'     "QAPGKER",
#'     "GRFTISR",
#'     "GTTVNNYWGQGTQVTVSSEQKLISEEDL",
#'     "GRFTISRDNAKNTVYLQMNSLK",
#'     "EREFVAAINSGGSTYYPDSVK",
#'     "QAPGKEREFVAAINSGGSTYYPDSVKGR",
#'     "NTVYLQMNSLKPEDTAVYYCAADLR",
#'     "LSCAASGFTFSSYAMGWFRQAPGKER",
#'     "QVQLVESGGGLVQAGGSLR",
#'     "PEDTAVYYCAADLRGTTVNNYWGQGTQVTVSSEQKLISEEDL",
#'     "FTISRDNAKNTVYLQMNSLKPEDTAVYYCAADLR",
#'     "LSCAASGFTFSSYAMGWFRQAPGK",
#'     "LSCAASGFTFSSYAMGWFR"
#'   ),
#'   rt_min = c(10, 28.5, 34.4, 34.4, 36, 36.5, 40.8, 
#'              42.5, 42.8, 43.3, 44.1, 44.8, 46.7),
#'   intensity = c(2769840, 2248170, 2172370, 1698280, 2202810, 
#'                 983267, 659246, 1064906, 1988932, 1438544, 
#'                 639990, 1017811, 1112824),
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Base coverage map
#' gg_seq(data = df_peptides, ref = ref_seq, wrap = 70)
#' 
#' # With peptide IDs and residue coloring
#' gg_seq(
#'   data = df_peptides, 
#'   ref = ref_seq, 
#'   name = "id",
#'   color = c(C = "red", K = "blue", R = "#468c2d"), 
#'   highlight = list(
#'     "#ffb4b4" = c(27:33, 51:57, 96:107),
#'     "#70bcfa" = c(1, 43, 64, 75, 86)
#'   ),
#'   wrap = 70
#' )
#' 
#' # With annotations
#' gg_seq(
#'   data = df_peptides, 
#'   ref = ref_seq, 
#'   name = "id",
#'   color = c(C = "red", K = "blue", R = "#468c2d"),
#'   highlight = list(
#'     "#ffb4b4" = c(27:33, 51:57, 96:107),  # CDR regions
#'     "#70bcfa" = c(1, 43, 64, 75, 86),     # Lysines
#'     "#d68718" = c(105:106),               # Liability site
#'     "#94d104" = c(119:128)                # c-Myc tag
#'   ),
#'   annotate = list(
#'     list(label = "CDR1", pos = 30),
#'     list(label = "CDR2", pos = 54),
#'     list(label = "CDR3", pos = 101),
#'     list(label = "N-term", pos = 1, angle = 90, vjust = 1),
#'     list(label = "K43", pos = 43, angle = 90),
#'     list(label = "K64", pos = 64, angle = 90),
#'     list(label = "K75", pos = 75, angle = 90),
#'     list(label = "K86", pos = 86, angle = 90),
#'     list(label = "Liability", pos = 106, angle = 90),
#'     list(label = "c-Myc tag", pos = 124)
#'   ),
#'   annotate_defaults = list(face = "bold"),
#'   wrap = 80
#' )
#' @export
gg_seq <- function(data,
                   ref,
                   sequence = "sequence",
                   name = NULL,
                   color = NULL,
                   highlight = NULL,
                   wrap = NULL,
                   size = 3,
                   face = "plain",
                   show_ref = TRUE,
                   margin_t = 30,
                   annotate = NULL,
                   annotate_defaults = 
                    list(angle = 0,
                         vjust = 0.5,
                         size = 3,
                         face = "plain",
                         color = "black")
                  ) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!sequence %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in data", sequence))
  }
  
  if (!is.null(name) && !name %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in data", name))
  }

  if (!is.character(ref) || length(ref) != 1 || nchar(ref) == 0) {
    stop("'ref' must be a single non-empty character string")
  }

  if (!is.null(wrap) && (!is.numeric(wrap) || wrap <= 0)) {
    stop("'wrap' must be a positive number or NULL")
  }
  
  if (!is.numeric(size) || size <= 0) {
    stop("'size' must be a positive number")
  }

  valid_faces <- c("plain", "bold", "italic", "bold.italic")
  if (!face %in% valid_faces) {
    stop(sprintf("'face' must be one of %s", paste(shQuote(valid_faces), collapse = ", ")))
  }

  if (!is.null(highlight)) {
    if (!is.list(highlight) || is.null(names(highlight))) {
      stop("'highlight' must be a named list with color names as keys and position vectors as values.\n  Example: highlight = list('lightblue' = c(5, 10), '#FF0000' = c(15, 20))")
    }
    
    # Check that all names are valid colors
    invalid_colors <- character(0)
    for (hl_color in names(highlight)) {
      # Try to convert color name to RGB - will fail for invalid colors
      tryCatch({
        grDevices::col2rgb(hl_color)
      }, error = function(e) {
        invalid_colors <<- c(invalid_colors, hl_color)
      })
    }
    
    if (length(invalid_colors) > 0) {
      stop(sprintf(
        "'highlight' contains invalid color name(s): %s\nUse valid R color names (e.g., 'yellow', 'lightblue', '#FF0000') or check colors() for available names",
        paste(shQuote(invalid_colors), collapse = ", ")
      ))
    }
    
    # Check that all values are numeric vectors
    non_numeric <- names(highlight)[!sapply(highlight, is.numeric)]
    if (length(non_numeric) > 0) {
      stop(sprintf(
        "'highlight' values must be numeric position vectors. Non-numeric found for: %s",
        paste(shQuote(non_numeric), collapse = ", ")
      ))
    }
  }

  # Validate annotate_defaults
  if (!is.list(annotate_defaults)) {
    stop("'annotate_defaults' must be a list")
  }

  # Merge user defaults with built-in defaults
  default_ann <- list(angle = 0, vjust = -0.5, size = 3, face = "plain", color = "black")
  annotate_defaults <- utils::modifyList(default_ann, annotate_defaults)
  valid_default_names <- c("angle", "vjust", "size", "face", "color")
  invalid_defaults <- setdiff(names(annotate_defaults), valid_default_names)
  if (length(invalid_defaults) > 0) {
    stop(sprintf("Invalid names in 'annotate_defaults': %s",
                paste(invalid_defaults, collapse = ", ")))
  }

  # Validate annotate parameter
  if (!is.null(annotate)) {
    if (!is.list(annotate)) {
      stop("'annotate' must be a list of annotation specifications")
    }
    
    for (i in seq_along(annotate)) {
      ann <- annotate[[i]]
      
      if (!is.list(ann)) {
        stop(sprintf("Annotation %d must be a list", i))
      }
      
      required_ann <- c("label", "pos")
      missing_ann <- setdiff(required_ann, names(ann))
      if (length(missing_ann) > 0) {
        stop(sprintf("Annotation %d missing required element(s): %s",
                    i, paste(missing_ann, collapse = ", ")))
      }
      
      if (!is.character(ann$label) || length(ann$label) != 1) {
        stop(sprintf("Annotation %d: 'label' must be a single character string", i))
      }
      
      if (!is.numeric(ann$pos) || length(ann$pos) != 1) {
        stop(sprintf("Annotation %d: 'pos' must be a single number", i))
      }
      
      # Validate optional parameters if present
      if (!is.null(ann$angle) && (!is.numeric(ann$angle) || length(ann$angle) != 1)) {
        stop(sprintf("Annotation %d: 'angle' must be a single number", i))
      }
      
      if (!is.null(ann$vjust) && (!is.numeric(ann$vjust) || length(ann$vjust) != 1)) {
        stop(sprintf("Annotation %d: 'vjust' must be a single number", i))
      }
      
      if (!is.null(ann$size) && (!is.numeric(ann$size) || ann$size <= 0)) {
        stop(sprintf("Annotation %d: 'size' must be a positive number", i))
      }
      
      if (!is.null(ann$face) && !ann$face %in% valid_faces) {
        stop(sprintf("Annotation %d: 'face' must be one of %s", 
                    i, paste(shQuote(valid_faces), collapse = ", ")))
      }
      
      if (!is.null(ann$color) && is.character(ann$color)) {
        is_valid <- tryCatch({
          grDevices::col2rgb(ann$color)
          TRUE
        }, error = function(e) FALSE)
        
        if (!is_valid) {
          stop(sprintf("Annotation %d: '%s' is not a valid color", i, ann$color))
        }
      }
    }
  }

  # Get unique sequences and filter out NA/empty
  all_seqs <- data[[sequence]]
  unique_seqs <- unique(all_seqs[!is.na(all_seqs) & nchar(all_seqs) > 0])

  if (length(unique_seqs) == 0) {
    warning("No valid sequences in data")
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
            ggplot2::labs(title = "No valid sequences"))
  }

  # Find positions using sapply (vectorized over patterns)
  positions <- sapply(unique_seqs, function(seq) {
    regexpr(seq, ref, fixed = TRUE)[1]
  })
  idx_found <- positions > 0

  if (!any(idx_found)) {
    warning("No sequences from 'data' were found in the reference sequence")
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
            ggplot2::labs(title = "No matching sequences found"))
  }

  # Filter to found sequences
  unique_seqs <- unique_seqs[idx_found]
  starts <- positions[idx_found]
  ends <- starts + nchar(unique_seqs) - 1

  # Get names (match back to original data)
  seq_names <- if (!is.null(name)) {
    sapply(unique_seqs, function(s) {
      data[[name]][match(s, data[[sequence]])]
    })
  } else {
    unique_seqs
  }

  # Create seq_data
  seq_data <- data.frame(
    sequence = unique_seqs,
    seq_name = seq_names,
    start = starts,
    end = ends,
    stringsAsFactors = FALSE
  )
  
  # Sort by start position
  seq_data <- seq_data[order(seq_data$start), ]
  seq_data$row <- seq_len(nrow(seq_data))
  
  # Build plot data efficiently
  plot_list <- list()

  # Add reference if requested
  if (show_ref) {
    ref_chars <- strsplit(ref, "")[[1]]
    plot_list[[1]] <- data.frame(
      pos = seq_along(ref_chars),
      char = ref_chars,
      seq_name = "Reference",
      row = 0,
      stringsAsFactors = FALSE
    )
  }

  # Add sequences using lapply (more efficient than loop with rbind)
  plot_list <- c(plot_list, lapply(seq_len(nrow(seq_data)), function(i) {
    chars <- strsplit(seq_data$sequence[i], "")[[1]]
    data.frame(
      pos = seq(seq_data$start[i], seq_data$end[i]),
      char = chars,
      seq_name = seq_data$seq_name[i],
      row = seq_data$row[i],
      stringsAsFactors = FALSE
    )
  }))

  df_plot <- do.call(rbind, plot_list)
  
  # Add color and fontface columns
  df_plot$color <- "black"
  df_plot$face <- face
  
  for (char in names(color)) {
    df_plot$color[df_plot$char == char] <- color[char]
    df_plot$face[df_plot$char == char] <- "bold"
  }

  # Apply wrapping if specified
  if (!is.null(wrap)) {
    df_plot$wrap_block <- ceiling(df_plot$pos / wrap)
    
    # Add invisible boundary points to align all blocks to same width
    blocks <- unique(df_plot$wrap_block)
    dummy_rows <- lapply(blocks, function(b) {
      data.frame(
        pos = c((b - 1) * wrap + 1, b * wrap),
        char = "",
        seq_name = "",
        row = if (show_ref) 0 else 1,
        color = "white",
        face = "plain",
        wrap_block = b,
        stringsAsFactors = FALSE
      )
    })
    df_plot <- rbind(df_plot, do.call(rbind, dummy_rows))
  }

  # Determine y-axis levels
  y_levels <- rev(seq(if (show_ref) 0 else 1, max(df_plot$row)))
  
  # Create y-axis labels
  y_labels <- if (!is.null(name)) {
    # Create labels vector with unique seq_names per row
    label_vec <- sapply(y_levels, function(r) {
      unique(df_plot$seq_name[df_plot$row == r])[1]
    })
    label_vec
  } else {
    rep("", length(y_levels))
  }

  # Apply wrapping if specified
  if (!is.null(wrap)) {
    df_plot$wrap_block <- ceiling(df_plot$pos / wrap)
  }

  # Create plot
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(
    x = .data$pos, 
    y = factor(row, levels = y_levels),
    label = char
  ))

  # Add highlight rectangles if specified (spans full plot height)
  if (!is.null(highlight)) {
    for (hl_color in names(highlight)) {
      pos_hl <- highlight[[hl_color]]
      if (is.numeric(pos_hl) && length(pos_hl) > 0) {
        # Sort positions
        pos_hl <- sort(unique(pos_hl))
        
        # Find consecutive ranges
        ranges <- list()
        start <- pos_hl[1]
        end <- pos_hl[1]
        
        for (i in seq_along(pos_hl)[-1]) {
          if (pos_hl[i] == end + 1) {
            end <- pos_hl[i]
          } else {
            ranges[[length(ranges) + 1]] <- c(start, end)
            start <- pos_hl[i]
            end <- pos_hl[i]
          }
        }
        ranges[[length(ranges) + 1]] <- c(start, end)
        
        # Create rectangles for each range
        df_hl <- do.call(rbind, lapply(ranges, function(r) {
          df_range <- data.frame(
            xmin = r[1] - 0.5,
            xmax = r[2] + 0.5,
            ymin = 0.5,
            ymax = length(y_levels) + 1,
            fill = hl_color
          )
          
          # Add wrap_block if wrapping is enabled
          if (!is.null(wrap)) {
            # Assign block based on start position
            df_range$wrap_block <- ceiling(r[1] / wrap)
          }
          
          df_range
        }))
        
        p <- p + ggplot2::geom_rect(
          data = df_hl,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$ymin, ymax = .data$ymax,
                       fill = .data$fill),
          alpha = 0.3,
          color = "black",
          linewidth = 0.1,
          inherit.aes = FALSE
        )
      }
    }
    p <- p + ggplot2::scale_fill_identity()
  }

  # Add text annotations if specified
  if (!is.null(annotate)) {
    ann_df <- do.call(rbind, lapply(annotate, function(ann) {
      angle_val <- if (!is.null(ann$angle)) ann$angle else annotate_defaults$angle
      data.frame(
        x = ann$pos,
        label = ann$label,
        angle = angle_val,
        hjust = if (!is.null(ann$hjust)) {
          ann$hjust
        } else if (angle_val == 90) {
          0
        } else {
          0.5
        },
        vjust = if (!is.null(ann$vjust)) ann$vjust else annotate_defaults$vjust,
        size = if (!is.null(ann$size)) ann$size else annotate_defaults$size,
        face = if (!is.null(ann$face)) ann$face else annotate_defaults$face,
        color = if (!is.null(ann$color)) ann$color else annotate_defaults$color,
        stringsAsFactors = FALSE
      )
    }))
    
    # Add wrap_block if wrapping is enabled
    if (!is.null(wrap)) {
      ann_df$wrap_block <- ceiling(ann_df$x / wrap)
    }
    
    # Calculate padding needed based on annotation angles and sizes
    max_angle <- max(abs(ann_df$angle))
    max_size <- max(ann_df$size)
    
    y_padding <- if (max_angle >= 75) {
      max_size * 0.8
    } else if (max_angle >= 30) {
      max_size * 0.5
    } else {
      max_size * 0.5
    }
    
    p <- p + ggplot2::geom_text(
      data = ann_df,
      ggplot2::aes(x = .data$x, label = .data$label),
      y = Inf,
      angle = ann_df$angle,
      hjust = ann_df$hjust,
      vjust = 0.5,
      size = ann_df$size,
      fontface = ann_df$face,
      color = ann_df$color,
      inherit.aes = FALSE
    )
  } else {
    y_padding <- 0
  }

  p <- p + ggplot2::geom_text(
      ggplot2::aes(color = .data$color, fontface = .data$face),
      size = size, 
      family = "mono"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_y_discrete(
      labels = y_labels,
      expand = ggplot2::expansion(add = c(0, y_padding))
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0.5),
      breaks = function(lims) {
        breaks <- scales::pretty_breaks()(lims)
        breaks[breaks > 0]
      }
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Position", y = "", title = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.ticks.x = ggplot2::element_line(color = "black"),  # Show x-axis ticks
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = if (!is.null(name)) {
        ggplot2::element_text(color = "black", size = 8)
      } else {
        ggplot2::element_blank()
      },
      plot.margin = ggplot2::margin(t = margin_t, r = 5, b = 5, l = 5),
      panel.spacing.y = ggplot2::unit(2, "lines"),
      strip.text = ggplot2::element_blank()
    ) +
    {
      if (!is.null(wrap)) {
        ggplot2::facet_wrap(~wrap_block, ncol = 1, scales = "free_x", 
                            strip.position = "right")
      }
    }
  return(p)
}