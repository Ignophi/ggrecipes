#' Sequence Difference Plot
#'
#' @description
#' Visualizes mutations across aligned sequences by displaying only positions
#' that differ from a reference sequence. Accepts either a data frame with
#' sequences (which are aligned internally using Needleman-Wunsch global
#' alignment) or a pre-aligned Clustal file. Supports character coloring and region
#' highlighting (e.g., CDRs, tags, binding sites). Useful for focusing on
#' differences in sequence comparisons. Inspired by Jalview's behavior when
#' setting a reference sequence and hiding matching positions.
#'
#' @inheritParams gg_seq
#' @param clustal Character string specifying path to a Clustal alignment file.
#'   If provided, overrides \code{data}, \code{ref}, \code{sequence}, and \code{name}.
#'   Default is NULL.
#'
#' @return A ggplot2 object showing the aligned sequences with differences highlighted.
#'   The plot displays:
#'   \itemize{
#'     \item Each sequence on a separate row
#'     \item Actual characters for positions that differ from the reference
#'     \item Optional reference sequence at the top (if \code{show_ref = TRUE})
#'     \item Specified characters highlighted in bold with custom colors
#'     \item Optional vertical shading bands at specified positions
#'     \item Optional text annotations above the plot
#'   }
#'   If no valid sequences are found, returns an empty ggplot2 object with a message.
#'
#' @details
#' The function aligns each sequence to the reference using Needleman-Wunsch
#' global alignment, then displays only differences. Positions matching the
#' reference are not displayed to reduce visual clutter and emphasize mutations.
#' Gaps from alignment are also shown as dashes.
#'
#' @import ggplot2
#' @import scales
#' @importFrom grDevices col2rgb
#' @importFrom utils modifyList capture.output
#'
#' @examples
#' # -----------------------------------------------------------------------
#' # Example with Clustal alignment file
#' # -----------------------------------------------------------------------
#' # Create a temporary Clustal file
#' clustal_file <- tempfile(fileext = ".aln")
#' writeLines(c(
#'   "CLUSTAL W (1.83) multiple sequence alignment",
#'   "",
#'   "WT              EQKLISEEDLMKTAYIAKQRQISFVKSHFSRQLERIEKKIEAHFDDLHP",
#'   "Mutant1         EQKLISEEDLMKTAYIAKQRQISFVKSHFSRQLERIEKKIEAHFDDLHP",
#'   "Mutant2         EQKLISEEDLMKTAYIAKQRQRSFVKSHFSRQLERIEKKWEAHFDDLHP",
#'   "Mutant3         EQKLISEEDLMKTAYIAKQRQISFVKSHFSRQLER----IEAHFDDLHP",
#'   "Mutant4         EQKLISEEDLMKTAYIAKQRQISFVKSHFSRQAERIEKKIEAHFDDLHP",
#'   "Mutant5         EQKLISEEDLAKTAYIAKQRQISFVKSHFSRQLERIEKKIEAHFDDRHP",
#'   "Mutant6         EQKLISEEDLMKTAYIAKQRQISFVKSHFSRQLERIEKKIEAHFDDLHP",
#'   "                *********** ***************** * ******* *******:**",
#'   "",
#'   "WT              DIVALSGHTFGKTHGAGKQSSHHHHHH",
#'   "Mutant1         DIVALSGHTFGKTHGAGKQSSHHHHHH",
#'   "Mutant2         DIVALSGHTFGKTHGAGKQSSHHHHHH",
#'   "Mutant3         DIVALSGHTFGKTHGAGKQSS------",
#'   "Mutant4         DIVALSGHTFGKTHGAGKQSSHHHHHH",
#'   "Mutant5         DIVALSGHTFGKTHGAGKQSSHHHHHH",
#'   "Mutant6         DRVALSGHTFAKTHGAGKQSS------",
#'   "                * ******** **********      "
#' ), clustal_file)
#' 
#' # Plot Clustal alignment
#' gg_seqdiff(
#'   clustal = clustal_file,
#'   ref = paste0("EQKLISEEDLMKTAYIAKQRQISFVKSHFSRQLERIEKKIEAHFDDLHP",
#'                "DIVALSGHTFGKTHGAGKQSSHHHHHH"),
#'   color = c(K = "#285bb8", R = "#285bb8",    # Basic
#'             E = "#a12b20", D = "#a12b20",    # Acidic
#'             W = "#9b59b6", F = "#9b59b6",    # Aromatic
#'             H = "#f39c12"),                  # Histidine
#'   highlight = list(
#'     "#94d104" = 1:10,      # N-terminal c-Myc tag
#'     "#FFE0B2" = 30:45,     # Active site
#'     "#94d104" = 72:77      # C-terminal His-tag
#'   ),
#'   annotate = list(
#'     list(label = "c-Myc", pos = 5),
#'     list(label = "Active site", pos = 37),
#'     list(label = "6xHis", pos = 74)
#'   ),
#'   wrap = 60
#' )
#' 
#' # Clean up
#' unlink(clustal_file)
#'
#' # -----------------------------------------------------------------------
#' # Example with DNA sequences - gene structure with regulatory elements
#' # -----------------------------------------------------------------------
#' dna_ref <- paste0(
#'   "TATAAA",                       # TATA box (promoter)
#'   "ATGCGATCGATCGATCGTAGCTAGCT",   # Exon 1
#'   "GTAAGTATCGATCGAT",             # Intron 1 (splice sites: GT...AG)
#'   "ACGTACGTACGTAGCTAGCTAGCTAC",   # Exon 2
#'   "GTACGTACGTACGTAC",             # Intron 2
#'   "GTACGTACGTAGCTAGCTAGCTACGT",   # Exon 3
#'   "ACGTACGTAAATAA"                # 3'UTR with poly-A signal
#' )
#' 
#' dna_df <- data.frame(
#'   sequence = c(
#'     dna_ref,                         
#'     sub("TATAAA", "TATATA", dna_ref),
#'     gsub("GTAAGT", "ATAAGT", dna_ref),
#'     gsub("CGATAG", "CGATAA", dna_ref),
#'     sub("ATG", "AAG", dna_ref),
#'     gsub("AATAA$", "AACAA", dna_ref),
#'     sub("GCGATCGATCGATCG", "GCGATCAATCGATCG", dna_ref),
#'     gsub("ACGTACGTACGTAG", "ACGTACATACGTAG", dna_ref)
#'   ),
#'   id = c("WT", "Promoter_mut", "Splice_donor",
#'          "Splice_acceptor", "Start_codon", "PolyA_mut",
#'          "Exon1_missense", "Exon2_frameshift")
#' )
#' 
#' # Highlight gene structure elements
#' gg_seqdiff(
#'   data = dna_df, 
#'   ref = dna_ref, 
#'   name = "id",
#'   color = c(G = "#4e8fb5", C = "#845cab"),
#'   highlight = list(
#'     "#FFE0B2" = 1:6,                     # TATA box (promoter)
#'     "#C8E6C9" = c(7:32, 49:74, 91:116),  # Exons
#'     "#FFCCBC" = 117:130                  # 3'UTR with poly-A
#'   ),
#'   annotate = list(
#'     list(label = "TATA", pos = 1, angle = 90),
#'     list(label = "ATG", pos = 7, angle = 90, color = "red"),
#'     list(label = "Exon1", pos = 19),
#'     list(label = "GT", pos = 33, angle = 90, size = 2.5),
#'     list(label = "GA", pos = 46, angle = 90, size = 2.5),
#'     list(label = "Exon2", pos = 61),
#'     list(label = "GT", pos = 75, angle = 90, size = 2.5),
#'     list(label = "AC", pos = 89, angle = 90, size = 2.5),
#'     list(label = "Exon3", pos = 103),
#'     list(label = "AATAAA", pos = 125, angle = 90, color = "blue")
#'   ),
#'   wrap = 80
#' )
#'
#' # -----------------------------------------------------------------------
#' # Example with antibody sequences with CDR mutations
#' # -----------------------------------------------------------------------
#' ref_seq <- paste0(
#'   "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAINSGGSTYYP",
#'   "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADLRGTTVKDYWGQGTQVTVSSEQKLISEEDL"
#' )
#' 
#' # All sequences must be same length as reference
#' mutant_df <- data.frame(
#'   sequence = c(
#'     ref_seq,  # Wild-type
#'     # CDR1 mutations (27-33)
#'     sub("GRTFSSYAMG", "GRTASSYAMG", ref_seq),
#'     # CDR2 mutations (51-57)
#'     sub("AINSGGS", "AINSAGS", ref_seq),
#'     # CDR3 mutations (96-107)
#'     sub("AADLRGTTVKDY", "AADLRGTTAKDY", ref_seq),
#'     # Framework mutations
#'     sub("QVQLVES", "EVQLVAS", ref_seq),
#'     # Multiple CDR mutations
#'     sub("AADLRGTTVKDY", "AADWRGTTVKDY", 
#'         sub("GRTFSSYAMG", "GYTASSAAMG", ref_seq))
#'   ),
#'   id = c("WT", "CDR1_F30A", "CDR2_G54A", "CDR3_V104A", 
#'          "FR1_E5A", "CDR1+3_multi")
#' )
#' 
#' # Highlight CDRs and tags, color key residues
#' gg_seqdiff(
#'   data = mutant_df, 
#'   ref = ref_seq, 
#'   name = "id",
#'   color = c(R = "#285bb8", K = "#285bb8",  # positive
#'             E = "#a12b20", D = "#a12b20",  # negative
#'             W = "#9b59b6", F = "#9b59b6"), # aromatic
#'   highlight = list(
#'     "#70bcfa" = 1,                       # N-terminal
#'     "#ffb4b4" = c(27:33, 51:57, 96:107), # CDRs
#'     "#94d104" = 119:128                  # c-Myc tag
#'   ),
#'   annotate = list(
#'     list(label = "N-term", pos = 1, angle = 90),
#'     list(label = "CDR1", pos = 30),
#'     list(label = "CDR2", pos = 54),
#'     list(label = "CDR3", pos = 102),
#'     list(label = "c-Myc", pos = 123)
#'   ),
#'   wrap = 66
#' )
#' @export
gg_seqdiff <- function(data, 
                       ref,
                       sequence = "sequence",
                       name = NULL,
                       clustal = NULL,
                       color = NULL,
                       highlight = NULL,
                       wrap = NULL,
                       size = 3,
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
  # Handle Clustal input
  if (!is.null(clustal)) {
    if (!is.character(clustal) || length(clustal) != 1) {
      stop("'clustal' must be a single file path")
    }
    
    # Parse file
    aln_data <- parse_clustal(clustal)
    
    # Set data
    data <- aln_data
    name <- "name"
    sequence <- "sequence"
  }

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
  
  if (!is.numeric(size) || size <= 0) {
    stop("'size' must be a positive number")
  }

  if (!is.null(wrap) && (!is.numeric(wrap) || wrap <= 0)) {
    stop("'wrap' must be a positive number or NULL")
  }

  # Validate annotate_defaults
  if (!is.list(annotate_defaults)) {
    stop("'annotate_defaults' must be a list")
  }

  # Merge user defaults with built-in defaults
  default_ann <- list(angle = 0, vjust = -0.5, size = 3, 
                    face = "plain", color = "black")
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
      
      valid_faces <- c("plain", "bold", "italic", "bold.italic")
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

  if (!is.null(highlight)) {
    if (!is.list(highlight) || is.null(names(highlight))) {
      stop("'highlight' must be a named list with color names as keys and position vectors as values.\n  Example: highlight = list('lightblue' = c(5, 10), '#FF0000' = c(15, 20))")
    }
    invalid_colors <- character(0)
    for (hl_color in names(highlight)) {
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
    non_numeric <- names(highlight)[!sapply(highlight, is.numeric)]
    if (length(non_numeric) > 0) {
      stop(sprintf(
        "'highlight' values must be numeric position vectors. Non-numeric found for: %s",
        paste(shQuote(non_numeric), collapse = ", ")
      ))
    }
  }

  # Rest of validation...
  if (is.null(data) || is.null(ref)) {
    stop("Must provide either 'data' and 'ref', or 'clustal' and 'ref'")
  }

  # Split reference into characters
  ref_chars <- strsplit(ref, "")[[1]]
  ref_length <- length(ref_chars)
  
  # Get unique sequences and filter out NA/empty
  all_seqs <- data[[sequence]]
  unique_seqs <- unique(all_seqs[!is.na(all_seqs) & nchar(all_seqs) > 0])

  if (length(unique_seqs) == 0) {
    warning("No valid sequences in data")
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
            ggplot2::labs(title = "No valid sequences"))
  }

  # Validate that all sequences are the same length as reference
  seq_lengths <- nchar(unique_seqs)
  idx_invalid <- seq_lengths != ref_length
  
  if (any(idx_invalid)) {
    invalid_seqs <- unique_seqs[idx_invalid]
    invalid_lengths <- seq_lengths[idx_invalid]
    
    # Build error message with details
    error_details <- data.frame(
      sequence_length = invalid_lengths,
      expected_length = ref_length
    )
    
    # Add names if available
    if (!is.null(name)) {
      seq_names <- sapply(invalid_seqs, function(s) {
        data[[name]][match(s, data[[sequence]])]
      })
      error_details <- cbind(name = seq_names, error_details)
    }
    
    error_msg <- sprintf(
      "All sequences must be the same length as the reference sequence.\nReference length: %d\n\nInvalid sequences found:\n%s",
      ref_length,
      paste(capture.output(print(error_details, row.names = FALSE)), collapse = "\n")
    )
    
    stop(error_msg)
  }

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
    stringsAsFactors = FALSE
  )
  
  # Assign row numbers
  seq_data$row <- seq_len(nrow(seq_data))
  
  # Build plot data efficiently
  plot_list <- list()

  # Add reference if requested
  if (show_ref) {
    plot_list[[1]] <- data.frame(
      pos = seq_along(ref_chars),
      char = ref_chars,
      seq_name = "Reference",
      row = 0,
      stringsAsFactors = FALSE
    )
  }

  # Add sequences using lapply - show only differences
  plot_list <- c(plot_list, lapply(seq_len(nrow(seq_data)), function(i) {
    seq_chars <- strsplit(seq_data$sequence[i], "")[[1]]
    
    # Compare each character to reference - show dash if match, character if different
    display_chars <- ifelse(seq_chars == ref_chars, "-", seq_chars)
    
    data.frame(
      pos = seq_along(seq_chars),
      char = display_chars,
      seq_name = seq_data$seq_name[i],
      row = seq_data$row[i],
      stringsAsFactors = FALSE
    )
  }))

  df_plot <- do.call(rbind, plot_list)
  
  # Add color and fontface columns
  df_plot$color <- ifelse(df_plot$char == "-", "grey70", "black")
  df_plot$face <- ifelse(df_plot$char == "-", "plain", "bold")
  
  # Apply highlighting only to non-dash characters
  for (char in names(color)) {
    match_idx <- df_plot$char == char & df_plot$char != "-"
    df_plot$color[match_idx] <- color[char]
  }
  # Replace x to "-"
  df_plot$char <- gsub("x", "-", df_plot$char, fixed = TRUE)

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
      vjust = ann_df$vjust,
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
      axis.ticks.x = ggplot2::element_line(color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = if (!is.null(name)) {
        ggplot2::element_text(color = "black")
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