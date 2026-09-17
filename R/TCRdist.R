#' GPU/CPU implementation of TCRdist, a distance/similarity metric for pairs of TCRs
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' An efficient, batched, chunked version of TCRdist that supports NVIDIA GPUs, Apple
#' Silicon GPUs, and a parallel CPU-only C++ backend that does not require python.
#'
#' @details
#' This function calculates pairwise TCRdist (Dash et al., Nature 2017) for a set of TCRs
#' (or between two sets of TCRs) and returns a sparse output with the TCRdist and indices of all pairs
#' that have TCRdist less than or equal to a desired cutoff (default cutoff is 90).
#'
#' TCR encoding, chunking (row-by-row and column-by-column), sparsification, and output
#' assembly all happen in R. For each chunk of \code{chunk_size} TCRs, R calls a single
#' minimal compute kernel for the one step that benefits from parallel hardware --
#' summing substitution-matrix penalties across encoded TCR features -- using either
#' python (`cupy` for NVIDIA GPUs, `mlx` for Apple Silicon GPUs, or `numpy` otherwise) or,
#' with \code{backend = "cpp"}, a parallel C++ implementation via RcppParallel that needs
#' no python at all. R then filters the resulting chunk by \code{tcrdist_cutoff} and
#' appends it to the running edge list.
#'
#' @param tcr1 a data frame with one TCR per row. It must have the columns "va", "vb", "cdr3a", and "cdr3b"
#' @param tcr2 (optional) another data frame of TCRs. If supplied, TCRdist will be calculated
#' for every combination of one TCR from tcr1 and one TCR from tcr2. Otherwise, it will calculate TCRdist
#' for every pair of TCRs in tcr1.
#' @param remove_MAIT whether to remove TCRs from MAIT cells (default is FALSE)
#' @param params (optional) a table of valid parameters for amino acids and va/vb segments.
#' (default is NULL, which uses TIRTLtools::params)
#' @param submat (optional) a substitution matrix with mismatch penalties for each
#' combination of amino acids or va/vb segments (default is NULL, which uses TIRTLtools::submat).
#' @param tcrdist_cutoff (optional) discard all TCRdist values above this cutoff. If not supplied by the user, this will default to 90 for dual-chain TCRdist or 45 for single-chain TCRdist.
#' @param chunk_size (optional) The chunk size to use in calculation of TCRdist (default 1000). If set at n,
#' it will calculate pairwise TCRdist for n x n TCRs at once. This may be as high as allowable by GPU memory.
#' @param write_to_tsv (optional) write the results to a tab-separated file ".tsv" (default is FALSE, does not write .tsv file)
#' @param output_folder (optional) folder to write output ".tsv" files to, if \code{write_to_tsv} is TRUE (default is the current directory).
#' @param backend (optional) the backend to use for the chunk computation (default "auto").
#' One of "auto" (pick a GPU backend if available, otherwise numpy), "cpp" (a parallel C++
#' implementation via RcppParallel -- fast, CPU-only, and does not require python at all),
#' "cpu" (numpy, via python), "cupy" (NVIDIA GPU, via python), or "mlx" (Apple Silicon GPU, via python).
#'
#' @return
#' If write_to_csv is TRUE (default is FALSE), the function will write output .tsv files and return NULL.
#' Otherwise, it will return a list with entries:
#'
#' \code{$edges_df} - a data frame with three columns: "node1_0index", "node2_0index", and "TCRdist".
#' The first two columns contain the indices (0-indexed) of the TCRs for each pair.
#' The last column contains the TCRdist if it is below the cutoff. The output is sparse in that it only contains
#' pairs that have TCRdist <= cutoff.
#'
#' \code{$nodes_df} - a data frame of the TCRs supplied to the function. It contains an additional column
#' "tcr_index" with the (0-indexed) index of each TCR. If \code{tcr2} was supplied, this is
#' \code{bind_rows(tcr1, tcr2)}: tcr1's rows are numbered first, and tcr2's \code{tcr_index}
#' values continue on immediately after the highest \code{tcr_index} in tcr1.
#'
#' @family tcr_similarity
#' @seealso \code{\link{cluster_tcrs}()}, \code{\link{plot_clusters}()}, and \code{\link{identify_non_functional_seqs}()}
#'
#' @export
#' @examples
#' load_example_data(dataset = "SJTRC_minimal")
#' df = get_all_tcrs(SJTRC_minimal, chain="paired", remove_duplicates = TRUE)
#' result = TCRdist(df, tcrdist_cutoff = 90)
#' edge_df = result[['edges_df']] %>%
#'   data.table::as.data.table() ### table of TCRdist values <= cutoff
#' node_df = result[['nodes_df']] %>%
#'   data.table::as.data.table() ### table of input data with indices
#'
#' edge_df ## sparse 3-column output: node1, node2, TCRdist
#' ## note that indices start at 0 and are found in node_df$tcr_index
#'
TCRdist = function(
    tcr1,
    tcr2 = NULL,
    remove_MAIT = FALSE,
    params = NULL,
    submat = NULL,
    tcrdist_cutoff = NULL,
    chunk_size = 1000,
    write_to_tsv = FALSE,
    output_folder = ".",
    backend = c("auto", "cupy", "mlx", "cpp")
    ) {
  ## former parameters, moved to hard-coded
  print_res = TRUE
  return_data = !write_to_tsv
  only_lower_tri = TRUE
  print_chunk_size = 10
  chunk_size_col = NULL

  backend = match.arg(backend)

  use_cpp = identical(backend, "cpp")
  if (!use_cpp) py_require( packages = .get_py_deps_new() )

  has_a = ifelse("cdr3a" %in% colnames(tcr1), TRUE, FALSE)
  has_b = ifelse("cdr3b" %in% colnames(tcr1), TRUE, FALSE)
  if (missing(tcrdist_cutoff) || is.null(tcrdist_cutoff)) {
    if (has_a && has_b) {
      tcrdist_cutoff = 90
      cli::cli_alert_info("Both {.field cdr3a} and {.field cdr3b} found — using {.val tcrdist_cutoff = 90}")
    } else if (has_a || has_b) {
      tcrdist_cutoff = 45
      cli::cli_alert_info("Only one of {.field cdr3a}/{.field cdr3b} found — using {.val tcrdist_cutoff = 45}")
    } else {
      cli::cli_abort("Neither {.field cdr3a} nor {.field cdr3b} found in the data frame")
    }
  }

  compare_to_self = is.null(tcr2)
  tcr1 = prep_for_tcrdist(tcr1, params = params, remove_MAIT = remove_MAIT)
  if (!compare_to_self) tcr2 = prep_for_tcrdist(tcr2, params = params, remove_MAIT = remove_MAIT)

  if (is.null(submat)) submat = TIRTLtools::submat
  if (is.null(params)) params = TIRTLtools::params
  params_vec = stats::setNames(params$value, params$feature)

  tcr1 = as.data.frame(tcr1)
  tcr1$tcr_index = seq_len(nrow(tcr1)) - 1L
  if (!compare_to_self) {
    tcr2 = as.data.frame(tcr2)
    tcr2$tcr_index = max(tcr1$tcr_index) + seq_len(nrow(tcr2))
  }

  tcr1_enc = .encode_tcrs_new(tcr1, params_vec)
  tcr2_enc = if (compare_to_self) tcr1_enc else .encode_tcrs_new(tcr2, params_vec)

  chunk_size = as.integer(chunk_size)
  if (is.null(chunk_size_col)) {
    chunk_size_col = chunk_size
  } else {
    chunk_size_col = as.integer(min(chunk_size_col, nrow(tcr2_enc)))
  }

  if (use_cpp) {
    if (print_res) cli::cli_alert_info("Using {.val cpp} backend (parallel C++ via RcppParallel, no python/GPU)")
  } else {
    backend_selected = .select_tcrdist_backend(backend)
    py_mod = reticulate::import_from_path(
      "TCRdist_core",
      path = system.file("python/TCRdist_new/", package = "TIRTLtools"),
      convert = TRUE,
      delay_load = TRUE
    )
  }

  n1 = nrow(tcr1_enc)
  n2 = nrow(tcr2_enc)
  starts1 = seq.int(0L, n1 - 1L, by = chunk_size)
  starts2 = seq.int(0L, n2 - 1L, by = chunk_size_col)

  chunk_pairs = list()
  for (ch1 in starts1) {
    for (ch2 in starts2) {
      if (compare_to_self && only_lower_tri && ch1 < ch2) next
      chunk_pairs[[length(chunk_pairs) + 1L]] = c(ch1, ch2)
    }
  }
  n_chunks = length(chunk_pairs)
  if (print_res) cli::cli_alert_info("Number of chunks: {n_chunks}")
  milestones = .get_progress_milestones(n_chunks, print_chunk_size)

  if (write_to_tsv) {
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    edge_file = file.path(output_folder, "TCRdist_df.tsv")
    if (file.exists(edge_file)) file.remove(edge_file)
    utils::write.table(tcr1, file.path(output_folder, "tcr1.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    if (!compare_to_self) {
      utils::write.table(tcr2, file.path(output_folder, "tcr2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }

  start_time = Sys.time()
  first_write = TRUE
  edge_list = vector("list", n_chunks)
  for (i in seq_len(n_chunks)) {
    ch1 = chunk_pairs[[i]][1]
    ch2 = chunk_pairs[[i]][2]
    end1 = min(ch1 + chunk_size, n1)
    end2 = min(ch2 + chunk_size_col, n2)
    rows1 = (ch1 + 1L):end1
    rows2 = (ch2 + 1L):end2

    if (use_cpp) {
      chunk_mat = tcrdist_chunk_cpp(
        tcr1_enc = tcr1_enc[rows1, , drop = FALSE],
        tcr2_enc = tcr2_enc[rows2, , drop = FALSE],
        submat = submat
      )
    } else {
      chunk_mat = py_mod$tcrdist_chunk(
        tcr1_enc = tcr1_enc[rows1, , drop = FALSE],
        tcr2_enc = tcr2_enc[rows2, , drop = FALSE],
        submat = submat,
        backend = backend_selected
      )
    }

    edges_tmp = .sparsify_chunk(
      chunk_mat = chunk_mat,
      ch1 = ch1,
      ch2 = ch2,
      tcrdist_cutoff = tcrdist_cutoff,
      compare_to_self = compare_to_self,
      only_lower_tri = only_lower_tri
    )

    if (write_to_tsv) {
      utils::write.table(edges_tmp, edge_file, sep = "\t", row.names = FALSE,
                          col.names = first_write, append = !first_write, quote = FALSE)
      first_write = FALSE
    }
    if (return_data) edge_list[[i]] = edges_tmp

    if (print_res && i %in% milestones) {
      percent = round(100 * i / n_chunks)
      elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cli::cli_alert_info("{percent}% done — time taken so far: {round(elapsed, 2)} seconds")
    }
  }
  if (print_res) {
    total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cli::cli_alert_success("Total time taken: {round(total_time, 2)} seconds")
  }

  if (!return_data) return(invisible(NULL))

  edges_df = data.table::rbindlist(edge_list) |>
    as_tibble() |>
    arrange(node1_0index, node2_0index)

  nodes_df = as_tibble(tcr1) |>
    mutate(source = "tcr1") |>
    select(tcr_index, source, everything())

  attr(edges_df, "pandas.index") <- NULL
  attr(nodes_df, "pandas.index") <- NULL

  if (!compare_to_self) {
    edges_df$node2_0index = edges_df$node2_0index + nrow(nodes_df)
  }

  if (!compare_to_self) {
    tcr2 = tcr2 |> mutate(source = "tcr2")
    nodes_df = nodes_df |> bind_rows(as_tibble(tcr2))
  }

  out = list(edges_df = edges_df, nodes_df = nodes_df)
  return(out)
}

### Encode a prepped TCR data frame into an integer matrix of features
### (same feature layout as the python process_TCRs() in inst/python/TCRdist/TCRdist_gpu.py:
### trimmed+padded cdr3a, va, trimmed+padded cdr3b, vb), using a named vector
### mapping feature (amino acid or V-segment string) -> integer code.
.encode_tcrs_new = function(df, params_vec) {
  use_alpha = all(c("cdr3a", "va") %in% colnames(df))
  use_beta = all(c("cdr3b", "vb") %in% colnames(df))
  if (!use_alpha) cli::cli_alert_warning("cdr3a and/or va not found in input data frame — alpha chain will not be used")
  if (!use_beta) cli::cli_alert_warning("cdr3b and/or vb not found in input data frame — beta chain will not be used")
  if (!use_alpha && !use_beta) cli::cli_abort("Neither alpha chain (cdr3a, va) nor beta chain (cdr3b, vb) found in the data frame")

  pad_len = 29L
  trim_cols = 4:27 ## drop first 3 and last 2 amino acids of the 29-character padded cdr3 (matches python slice(3,-2))

  mats = list()
  if (use_alpha) {
    cdr3a_padded = .pad_center_vec(df$cdr3a, pad_len)
    cdr3a_mat = do.call(rbind, strsplit(cdr3a_padded, ""))[, trim_cols, drop = FALSE]
    cdr3a_int = matrix(params_vec[as.vector(cdr3a_mat)], nrow = nrow(cdr3a_mat))
    va_int = matrix(params_vec[df$va], ncol = 1)
    mats = c(mats, list(cdr3a_int, va_int))
  }
  if (use_beta) {
    cdr3b_padded = .pad_center_vec(df$cdr3b, pad_len)
    cdr3b_mat = do.call(rbind, strsplit(cdr3b_padded, ""))[, trim_cols, drop = FALSE]
    cdr3b_int = matrix(params_vec[as.vector(cdr3b_mat)], nrow = nrow(cdr3b_mat))
    vb_int = matrix(params_vec[df$vb], ncol = 1)
    mats = c(mats, list(cdr3b_int, vb_int))
  }
  encoded = do.call(cbind, mats)
  storage.mode(encoded) = "integer"
  return(encoded)
}

### Turn a dense chunk_mat (n1 x n2 of TCRdist values) into a sparse edge data frame,
### applying the cutoff and (for self-comparisons) the lower-triangle/diagonal filter.
.sparsify_chunk = function(chunk_mat, ch1, ch2, tcrdist_cutoff, compare_to_self, only_lower_tri) {
  idx = which(chunk_mat <= tcrdist_cutoff, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(node1_0index = integer(0), node2_0index = integer(0), TCRdist = integer(0)))
  }
  edges_tmp = data.frame(
    node1_0index = idx[, 1] - 1L + ch1,
    node2_0index = idx[, 2] - 1L + ch2,
    TCRdist = chunk_mat[idx]
  )
  if (compare_to_self) {
    if (only_lower_tri) {
      edges_tmp = edges_tmp[edges_tmp$node1_0index > edges_tmp$node2_0index, , drop = FALSE]
    } else {
      edges_tmp = edges_tmp[edges_tmp$node1_0index != edges_tmp$node2_0index, , drop = FALSE]
    }
  }
  return(edges_tmp)
}

### Progress milestones (chunk indices at which to print progress), analogous
### to the milestone logic in inst/python/TCRdist/TCRdist_gpu.py::TCRdist_batch()
.get_progress_milestones = function(n_chunks, print_chunk_size) {
  if (n_chunks == 0) return(integer(0))
  print_chunk_size = as.integer(print_chunk_size)
  if (is.na(print_chunk_size) || print_chunk_size < 1 || print_chunk_size > 99) print_chunk_size = 10L
  pcts = seq(print_chunk_size, 100, by = print_chunk_size)
  milestones = unique(pmax(1L, as.integer(floor(n_chunks * pcts / 100))))
  return(milestones)
}

### Pick which array backend the python kernel should use.
.select_tcrdist_backend = function(backend = c("auto", "cpp", "cupy", "mlx")) {
  backend = match.arg(backend)
  #if (backend == "cpu") return("numpy")
  #if (backend %in% c("cupy", "mlx")) return(backend)
  if (backend == "cpp") return("cpp")
  if (backend %in% c("auto", "cupy", "mlx")) {
    if (.has_nvidia_gpu()) return("cupy")
    if (.is_apple_silicon()) return("mlx")
  }
}

### Minimal python dependencies for TCRdist_new() -- just numpy plus a GPU array
### library when available. (TCRdist_new() does not need pandas/scipy since all
### chunking, sparsification, and output assembly happen in R.)
.get_py_deps_new = function(with_versions = FALSE) {
  base_deps = c("numpy==2.3.5")
  if (.has_nvidia_gpu()) base_deps = c(base_deps, "cupy==13.6.0")
  if (.is_apple_silicon()) base_deps = c(base_deps, "mlx==0.29.4")
  if (with_versions) return(base_deps)
  return(gsub("==.*", "", base_deps))
}

.fix_py_to_r_df_list = function(res) {
  res = lapply(res, .fix_py_to_r_df)
  return(res)
}

.fix_py_to_r_df = function(df) {
  mat = df$values
  ncols = length(df$columns$values)
  cols = sapply(0:(ncols-1), function(i) df$columns$values[i])
  ll = lapply(1:ncols, function(i) sapply(mat[,i], function(x) x[[1]]))
  df_tmp = as.data.frame(ll)
  colnames(df_tmp) = cols
  return(df_tmp)
}
