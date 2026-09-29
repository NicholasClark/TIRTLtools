#' Create a PDF report of quality control plots for a TIRTL-seq dataset
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function runs the main quality control (QC) functions of the package on a
#' dataset and writes the resulting tables and plots to a single PDF file.
#'
#' The report contains:
#' * a summary table of reads, unique chains, and pairs for each sample
#'   (\code{\link{summarize_data}()} and \code{\link{get_pair_stats}()})
#' * the number of reads per sample (\code{\link{plot_n_reads}()})
#' * the number of pairs called by each algorithm (\code{\link{plot_paired}()})
#' * the fraction of chains paired by read fraction range (\code{\link{plot_paired_by_read_fraction_range}()})
#' * the number of partners of each chain (\code{\link{plot_num_partners}()} and \code{\link{plot_n_paired_step}()})
#' * clonotype rank plots (\code{\link{plot_ranks}()})
#' * the fraction of reads from the most frequent clonotypes (\code{\link{plot_clonotype_indices}()})
#' * a heatmap of overlap between samples (\code{\link{plot_sample_overlap}()}, requires
#'   the ComplexHeatmap and dendsort packages and more than one sample)
#' * (optional) per-sample plots of pairing status for the most frequent chains
#'   (\code{\link{plot_paired_vs_rank}()}, \code{\link{plot_read_fraction_vs_pair_status}()},
#'   and \code{\link{plot_pairs_with_eachother}()})
#'
#' Two page layouts are available:
#' * `layout = "one_per_page"` (default): one plot per page on landscape pages
#'   (11 x 8.5 inches), good for viewing on screen.
#' * `layout = "multi_per_page"`: several plots per page on standard portrait
#'   pages (8.5 x 11 inches), good for printing. Dataset-level plots are placed two per
#'   page and the per-sample plots for each sample are placed together on one page.
#'
#' In both layouts, text sizes of the plots are set so that axis labels and titles
#' are readable when the PDF is printed at its page size.
#'
#' If any individual plot fails, the error message is written in its place and a
#' warning is issued, so that the rest of the report is still created.
#'
#' @param data a TIRTLseqDataSet object created by \code{\link{load_tirtlseq}()}
#' @param file the path of the output PDF file (default is "TIRTLseq_QC_report.pdf")
#' @param layout the page layout, either "one_per_page" (one plot per landscape page)
#' or "multi_per_page" (multiple plots per portrait page). Default is "one_per_page".
#' @param samples (optional) the samples to include in the report (default is all samples)
#' @param per_sample whether to include per-sample plots of pairing status for the most
#' frequent single-chains (default is TRUE)
#' @param n_max the number of most frequent single-chains to use for the per-sample plots (default is 50)
#' @param title (optional) a title for the first page of the report
#' @param width (optional) the width of each PDF page in inches (default is 11 for
#' "one_per_page" and 8.5 for "multi_per_page")
#' @param height (optional) the height of each PDF page in inches (default is 8.5 for
#' "one_per_page" and 11 for "multi_per_page")
#' @param verbose whether to print progress messages (default is TRUE)
#' @family qc
#'
#' @returns invisibly returns the path to the PDF file
#' @export
#'
#' @examples
#' \dontrun{
#' data = load_tirtlseq("path/to/data")
#' ## one plot per landscape page
#' qc_report(data, file = "my_QC_report.pdf")
#' ## multiple plots per portrait (8.5 x 11) page, for printing
#' qc_report(data, file = "my_QC_report_print.pdf", layout = "multi_per_page")
#' }
qc_report = function(
    data,
    file = "TIRTLseq_QC_report.pdf",
    layout = c("multi_per_page", "one_per_page"),
    samples = NULL,
    per_sample = TRUE,
    n_max = 50,
    title = NULL,
    width = NULL,
    height = NULL,
    verbose = TRUE
) {
  if(!inherits(data, "TIRTLseqDataSet")) {
    stop("'data' must be a TIRTLseqDataSet object created by load_tirtlseq()")
  }
  layout = match.arg(layout)
  checkmate::assert_string(file)
  checkmate::assert_flag(per_sample)
  checkmate::assert_count(n_max, positive = TRUE)
  checkmate::assert_number(width, lower = 1, null.ok = TRUE)
  checkmate::assert_number(height, lower = 1, null.ok = TRUE)
  checkmate::assert_flag(verbose)

  if(!is.null(samples)) {
    ind = if(is.numeric(samples)) samples else match(samples, names(data$data))
    if(any(is.na(ind))) {
      stop(paste("Samples not found in data:", paste(samples[is.na(ind)], collapse = ", ")))
    }
    data$data = data$data[ind]
    data$meta = data$meta[ind, , drop = FALSE]
  }
  sample_names = names(data$data)
  n_samples = length(sample_names)
  if(is.null(title)) title = "TIRTL-seq QC report"

  multi = layout == "multi_per_page"
  if(is.null(width)) width = if(multi) 8.5 else 11
  if(is.null(height)) height = if(multi) 11 else 8.5
  ## base font size (pt) for plot text -- pages in the "multi" layout hold smaller plots
  base_size = if(multi) 10 else 14

  ## open the device and make sure it is closed even if something fails
  grDevices::pdf(file, width = width, height = height, onefile = TRUE)
  dev_num = grDevices::dev.cur()
  on.exit(if(dev_num %in% grDevices::dev.list()) grDevices::dev.off(dev_num), add = TRUE)

  .qc_msg = function(...) if(verbose) message(...)

  ## title page and summary table
  .qc_msg("Summarizing data...")
  # summary_df = .qc_try("summary table", {
  #   df_sum = summarize_data(data)
  #   df_pairs = get_pair_stats(data, verbose = FALSE, by_method = FALSE)
  #   .qc_summary_table(df_sum, df_pairs)
  # })
  summary_df = summary(data)
  .qc_title_page(title, n_samples, summary_df)
  if(is.data.frame(summary_df)) {
    .qc_table_pages(summary_df, fontsize = if(multi) 8 else 9)
  }

  ## split samples into two batches for the read fraction plots if there are many samples
  if(n_samples > 12) {
    half = ceiling(n_samples / 2)
    batches = list(sample_names[1:half], sample_names[(half + 1):n_samples])
  } else {
    batches = list(sample_names)
  }
  read_fraction_plots = list()
  for(chain in c("beta", "alpha")) {
    for(i in seq_along(batches)) {
      local({
        chain = chain
        batch = batches[[i]]
        nm = paste0("Paired by read fraction (", chain, ")")
        title = paste("Fraction of", chain, "chains paired, by read fraction")
        if(length(batches) > 1) {
          half_label = c("first half of samples", "second half of samples")[i]
          nm = paste0(nm, " - ", half_label)
          title = paste0(title, " (", half_label, ")")
        }
        read_fraction_plots[[nm]] <<- function() {
          plot_paired_by_read_fraction_range(data, chain = chain, samples = batch) +
            ggtitle(title)
        }
      })
    }
  }

  ## dataset-level plots
  plots = c(list(
    "Number of reads" = function() plot_n_reads(data) +
      ggtitle("Number of reads per sample"),
    "Pairs called by each algorithm" = function() plot_paired(data, chain = "paired") +
      ggtitle("Number of alpha-beta pairs by pairing algorithm")),
    read_fraction_plots,
    list("Number of partners" = function() plot_num_partners(data) +
      facet_wrap(~Group, ncol = min(3, n_samples)) +
      ggtitle("Number of partners per chain (functional pairs)"),
    "Number of partners (step plot)" = function() {
      plot_n_paired_step(data, chain = "both", max_x = 20) +
        ggtitle("Cumulative fraction of pairs by number of partners") +
        .qc_one_column_legend(n_samples, base_size)
    }
  ))
  .qc_msg("Plotting dataset-level QC plots...")
  if(multi && n_samples > 9) {
    ## with many samples the number of partners plot gets a full page to itself
    i_np = which(names(plots) == "Number of partners")
    .qc_draw_plots(plots[seq_len(i_np - 1)], per_page = 2, base_size = base_size)
    .qc_draw_plots(plots[i_np], per_page = 1, base_size = base_size)
    .qc_draw_plots(plots[-seq_len(i_np)], per_page = 2, base_size = base_size)
  } else {
    .qc_draw_plots(plots, per_page = if(multi) 2 else 1, base_size = base_size)
  }

  ## clonotype indices plots, both on one page
  clonotype_plots = list(
    "Clonotype indices (beta)" = function() plot_clonotype_indices(data, chain = "beta") +
      ggtitle("Fraction of reads from the most frequent clonotypes (beta)"),
    "Clonotype indices (alpha)" = function() plot_clonotype_indices(data, chain = "alpha") +
      ggtitle("Fraction of reads from the most frequent clonotypes (alpha)")
  )
  .qc_draw_plots(clonotype_plots, per_page = 2, base_size = base_size)

  ## sample overlap heatmaps
  if(n_samples > 1) {
    if(all(sapply(c("ComplexHeatmap", "dendsort"), requireNamespace, quietly = TRUE))) {
      .qc_msg("Plotting sample overlap heatmaps...")
      n_top = 200
      heatmaps = lapply(c("paired", "beta", "alpha"), function(chain) {
        function() {
          olap = plot_sample_overlap(data, n_seq = n_top, chain = chain, return_data = TRUE)
          .qc_overlap_heatmap(olap, title = paste("Overlap of top", n_top, chain, "chains between samples"),
                              base_size = base_size)
        }
      }) %>% setNames(paste0("Sample overlap (", c("paired", "beta", "alpha"), ")"))
      .qc_draw_plots(heatmaps, per_page = if(multi) 2 else 1, base_size = base_size)
    } else {
      .qc_msg("Skipping sample overlap heatmaps (requires the 'ComplexHeatmap' and 'dendsort' packages)")
    }
  }

  ## rank plots, both on one page
  .qc_msg("Plotting rank plots...")
  rank_plots = list(
    "Rank plot (beta)" = function() plot_ranks(data, chain = "beta") +
      ggtitle("Clonotype rank plot (beta)") + .qc_one_column_legend(n_samples, base_size),
    "Rank plot (alpha)" = function() plot_ranks(data, chain = "alpha") +
      ggtitle("Clonotype rank plot (alpha)") + .qc_one_column_legend(n_samples, base_size)
  )
  .qc_draw_plots(rank_plots, per_page = 2, base_size = base_size)

  ## per-sample plots
  if(per_sample) {
    for(sample in sample_names) {
      .qc_msg("Plotting per-sample QC: ", sample)
      sample_plots = list(
        "Pairing status vs. rank" = function() plot_paired_vs_rank(data, sample = sample, n_max = n_max),
        "Read fraction vs. pairing status" = function() plot_read_fraction_vs_pair_status(data, sample = sample, n_max = n_max),
        "Pairs with each other" = function() plot_pairs_with_eachother(data, sample = sample, n_max = n_max)
      )
      ## title each plot with the sample and plot name
      sample_plots = lapply(names(sample_plots), function(nm) {
        fn = sample_plots[[nm]]
        function() fn() + ggtitle(paste0(sample, ": ", nm))
      }) %>% setNames(paste0(sample, ": ", names(sample_plots)))
      if(multi) {
        .qc_draw_plots(sample_plots, per_page = 3, base_size = base_size, header = paste("Sample:", sample))
      } else {
        .qc_section_page(paste("Sample:", sample))
        .qc_draw_plots(sample_plots, per_page = 1, base_size = base_size)
      }
    }
  }

  .qc_msg("QC report written to: ", file)
  invisible(file)
}

#' Evaluate an expression, returning a "qc_error" object on failure
#' @noRd
.qc_try = function(name, expr) {
  tryCatch(expr, error = function(e) {
    warning(paste0("QC report: '", name, "' failed: ", conditionMessage(e)), call. = FALSE)
    structure(list(name = name, message = conditionMessage(e)), class = "qc_error")
  })
}

#' Sample overlap heatmap without annotation bars, with row/column names sized to fit
#'
#' Called inside the viewport the heatmap is drawn in, so the font size is chosen from
#' the space available for each row/column name.
#' @noRd
.qc_overlap_heatmap = function(olap, title, base_size = 10) {
  n = nrow(olap)
  ## cluster on Jaccard similarity, as in plot_sample_overlap()
  n_seq = diag(olap)
  jaccard = olap / (outer(n_seq, n_seq, "+") - olap)
  hc = stats::hclust(stats::as.dist(1 - jaccard), method = "average")
  dend = dendsort::dendsort(stats::as.dendrogram(hc), type = "average")

  ## heatmap body takes roughly 60% of the smaller viewport side; fit one name per cell
  vp_h = grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  vp_w = grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  side = 0.6 * min(vp_h, vp_w)
  cell_in = side / n
  fontsize = max(4, min(base_size, 0.8 * cell_in * 72))

  ComplexHeatmap::Heatmap(olap,
                          name = "Overlap",
                          show_row_names = TRUE,
                          show_column_names = TRUE,
                          row_names_gp = grid::gpar(fontsize = fontsize),
                          column_names_gp = grid::gpar(fontsize = fontsize),
                          cluster_rows = dend,
                          cluster_columns = dend,
                          ## square heatmap body
                          width = grid::unit(side, "inches"),
                          height = grid::unit(side, "inches"),
                          column_title = title,
                          column_title_gp = grid::gpar(fontsize = base_size + 2, fontface = "bold"))
}

#' Legend of the "color" aesthetic in one column, with text shrunk if needed to fit
#'
#' Must be called inside the viewport the plot is drawn in, since the text size is
#' chosen from the height of that viewport. The size is set in the guide's own theme so
#' that it isn't overridden by .qc_text_theme().
#' @noRd
.qc_one_column_legend = function(n_items, base_size) {
  vp_h_pt = grid::convertHeight(grid::unit(1, "npc"), "points", valueOnly = TRUE)
  legend_size = max(4, min(base_size - 1, 0.8 * vp_h_pt / (1.4 * (n_items + 2))))
  guides(color = guide_legend(ncol = 1, theme = theme(
    legend.text = element_text(size = legend_size),
    legend.key.height = grid::unit(1.3 * legend_size, "pt"),
    legend.key.spacing.y = grid::unit(0, "pt")
  )))
}

#' Theme adjustments to give QC plots readable text sizes
#' @noRd
.qc_text_theme = function(base_size) {
  theme(
    text = element_text(size = base_size),
    axis.title = element_text(size = base_size),
    axis.text = element_text(size = base_size - 1),
    strip.text = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size - 1),
    plot.title = element_text(size = base_size + 2, face = "bold")
  )
}

#' Draw a list of plots, `per_page` plots stacked vertically on each page
#'
#' Each element of `plot_fns` is a function returning a ggplot or ComplexHeatmap object.
#' Plots always take up 1/per_page of the page height so that they are the same size
#' on every page, including a partially filled last page.
#' @noRd
.qc_draw_plots = function(plot_fns, per_page = 1, base_size = 14, header = NULL) {
  pages = split(seq_along(plot_fns), ceiling(seq_along(plot_fns) / per_page))
  header_h = if(is.null(header)) 0 else 0.5
  for(page in pages) {
    grid::grid.newpage()
    ## outer margins of 0.4 inches
    grid::pushViewport(grid::viewport(width = grid::unit(1, "npc") - grid::unit(0.8, "inches"),
                                      height = grid::unit(1, "npc") - grid::unit(0.8, "inches")))
    if(!is.null(header)) {
      grid::grid.text(header, x = 0, y = grid::unit(1, "npc") - grid::unit(header_h / 2, "inches"),
                      just = "left", gp = grid::gpar(fontsize = base_size + 6, fontface = "bold"))
    }
    grid::pushViewport(grid::viewport(y = 0, just = "bottom",
                                      height = grid::unit(1, "npc") - grid::unit(header_h, "inches"),
                                      layout = grid::grid.layout(nrow = per_page, ncol = 1)))
    for(i in seq_along(page)) {
      nm = names(plot_fns)[page[i]]
      grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      res = .qc_try(nm, {
        obj = plot_fns[[page[i]]]()
        ## ComplexHeatmap objects are drawn with draw(), ggplots with print()
        if(inherits(obj, "Heatmap") || inherits(obj, "HeatmapList")) {
          ComplexHeatmap::draw(obj, newpage = FALSE)
        } else {
          print(obj + .qc_text_theme(base_size), newpage = FALSE)
        }
        TRUE
      })
      if(inherits(res, "qc_error")) .qc_error_text(res, base_size)
      grid::popViewport()
    }
    grid::popViewport(2)
  }
}

#' Error message in place of a failed plot
#' @noRd
.qc_error_text = function(err, base_size = 12) {
  grid::grid.rect(gp = grid::gpar(col = "grey70", fill = NA))
  grid::grid.text(paste0(err$name, "\n\nThis plot could not be created:\n", err$message),
                  gp = grid::gpar(fontsize = base_size, col = "firebrick"))
}

#' Page with a section header
#' @noRd
.qc_section_page = function(text) {
  grid::grid.newpage()
  grid::grid.text(text, gp = grid::gpar(fontsize = 24, fontface = "bold"))
}

#' Title page of the report
#' @noRd
.qc_title_page = function(title, n_samples, summary_df) {
  grid::grid.newpage()
  lines = c(
    paste("TIRTLtools version", as.character(utils::packageVersion("TIRTLtools"))),
    paste("Created:", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste("Number of samples:", n_samples)
  )
  if(is.data.frame(summary_df)) {
    lines = c(lines,
      paste("Total alpha reads:", format(sum(summary_df$`Number of reads (TCR-alpha)`), big.mark = ",")),
      paste("Total beta reads:", format(sum(summary_df$`Number of reads (TCR-beta)`), big.mark = ",")),
      paste("Average pairs per sample:", format(mean(summary_df$`Paired TCRs`), big.mark = ","))
    )
  }
  grid::grid.text(title, y = 0.7, gp = grid::gpar(fontsize = 28, fontface = "bold"))
  grid::grid.text(paste(lines, collapse = "\n"), y = 0.45, gp = grid::gpar(fontsize = 14))
}

#' Draw a data frame as a table over as many pages as needed
#'
#' Works with any columns: column widths are measured from the text, long column
#' names are wrapped onto several lines, numeric columns are right-aligned, the font
#' is shrunk if the table is too wide for the page, and rows are split across pages
#' to fit the page height.
#' @noRd
.qc_table_pages = function(df, fontsize = 9, title = "Summary of samples", min_fontsize = 5) {
  df = as.data.frame(df)
  if(ncol(df) == 0) return(invisible(NULL))
  is_num = vapply(df, is.numeric, logical(1))
  cells = lapply(df, .qc_format_column)
  header = colnames(df)

  ## start the first page before measuring text so that measuring doesn't open a blank page
  grid::grid.newpage()
  margin = 0.5 ## inches
  pad = 0.15 ## inches between columns
  page_w = grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) - 2 * margin
  page_h = grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE) - 2 * margin

  ## text width in inches at a given font size
  str_w = function(x, fs, bold = FALSE) {
    if(length(x) == 0) return(0)
    gp = grid::gpar(fontsize = fs, fontface = if(bold) "bold" else "plain")
    vapply(x, function(s) grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = gp)),
                                             "inches", valueOnly = TRUE), numeric(1))
  }
  ## wrap a header on spaces so that no line is wider than `w`
  wrap_header = function(h, w, fs) {
    words = strsplit(h, " ", fixed = TRUE)[[1]]
    lines = character(0)
    cur = ""
    for(word in words) {
      cand = if(cur == "") word else paste(cur, word)
      if(cur != "" && str_w(cand, fs, bold = TRUE) > w) {
        lines = c(lines, cur)
        cur = word
      } else {
        cur = cand
      }
    }
    paste(c(lines, cur), collapse = "\n")
  }
  ## column widths: fit the full header when there is room; otherwise at least the widest
  ## cell and the widest single header word, with any spare width given to the columns
  ## whose headers need it most
  get_col_w = function(fs) {
    min_w = vapply(seq_along(cells), function(j) {
      words = strsplit(header[j], "[ \n]")[[1]]
      max(str_w(cells[[j]], fs), str_w(words, fs, bold = TRUE))
    }, numeric(1))
    full_w = pmax(min_w, str_w(header, fs, bold = TRUE))
    avail = page_w - pad * (length(cells) - 1)
    if(sum(full_w) <= avail) return(full_w)
    extra = full_w - min_w
    spare = max(0, avail - sum(min_w))
    if(sum(extra) > 0) min_w = min_w + extra * spare / sum(extra)
    min_w
  }

  ## shrink the font until the table fits the page width
  col_w = get_col_w(fontsize)
  while(sum(col_w) + pad * (length(col_w) - 1) > page_w && fontsize > min_fontsize) {
    fontsize = fontsize - 0.5
    col_w = get_col_w(fontsize)
  }
  header_wrapped = vapply(seq_along(header), function(j) wrap_header(header[j], col_w[j], fontsize), character(1))

  x_left = cumsum(c(0, utils::head(col_w, -1) + pad))
  x_right = x_left + col_w
  line_h = 1.2 * fontsize / 72 ## inches
  row_h = 1.6 * fontsize / 72
  n_header_lines = max(lengths(strsplit(header_wrapped, "\n", fixed = TRUE)))
  title_h = 0.5
  header_h = n_header_lines * line_h + 0.1
  rows_per_page = max(1, floor((page_h - title_h - header_h) / row_h))
  n_row = nrow(df)
  pages = if(n_row == 0) list(integer(0)) else split(seq_len(n_row), ceiling(seq_len(n_row) / rows_per_page))

  inch = function(x) grid::unit(x, "inches")
  for(p in seq_along(pages)) {
    rows = pages[[p]]
    if(p > 1) grid::grid.newpage()
    grid::pushViewport(grid::viewport(x = inch(margin), y = inch(margin), width = inch(page_w),
                                      height = inch(page_h), just = c("left", "bottom")))
    page_title = title
    if(length(pages) > 1) page_title = paste0(title, " (", p, "/", length(pages), ")")
    grid::grid.text(page_title, x = 0, y = inch(page_h), just = c("left", "top"),
                    gp = grid::gpar(fontsize = 16, fontface = "bold"))
    header_bottom = page_h - title_h - header_h
    for(j in seq_along(header_wrapped)) {
      grid::grid.text(header_wrapped[j], x = inch(if(is_num[j]) x_right[j] else x_left[j]),
                      y = inch(header_bottom + 0.05), just = c(if(is_num[j]) "right" else "left", "bottom"),
                      gp = grid::gpar(fontsize = fontsize, fontface = "bold", lineheight = 1))
    }
    grid::grid.lines(x = inch(c(0, max(x_right))), y = inch(header_bottom))
    for(i in seq_along(rows)) {
      y_i = header_bottom - (i - 0.5) * row_h
      for(j in seq_along(cells)) {
        grid::grid.text(cells[[j]][rows[i]], x = inch(if(is_num[j]) x_right[j] else x_left[j]), y = inch(y_i),
                        just = if(is_num[j]) "right" else "left", gp = grid::gpar(fontsize = fontsize))
      }
    }
    grid::popViewport()
  }
  invisible(NULL)
}

#' Format one table column as character for printing
#' @noRd
.qc_format_column = function(x) {
  if(is.numeric(x)) {
    whole = all(is.na(x) | (is.finite(x) & x == round(x)))
    out = if(whole) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
          else format(x, big.mark = ",", digits = 3, trim = TRUE)
  } else {
    out = as.character(x)
  }
  out[is.na(x)] = ""
  out
}
