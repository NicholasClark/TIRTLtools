
## plot_tshell_old ---------------------------------------------------

# plot_tshell_old = function(
#     ll_alpha,
#     ll_beta,
#     nuc_seq,
#     chain,
#     n_plot = 10,
#     #plot_overlap_only = TRUE,
#     interactive = FALSE
# ) {
#   if(!chain %in% c("alpha", "beta")) stop("'chain' needs to be 'alpha' or 'beta'")
#
#   ll_alpha$row_data$wa = Matrix::colSums(ll_alpha$matrix!= 0)
#   ll_beta$row_data$wb = Matrix::colSums(ll_beta$matrix!= 0)
#
#   alpha_mat = ll_alpha$matrix
#   beta_mat = ll_beta$matrix
#
#   alpha_mat_rf = normalize_rows(alpha_mat)
#   beta_mat_rf = normalize_rows(beta_mat)
#
#   #if("data.frame" %in% class(test$rows_alpha))
#   if(chain == "alpha") idx = which(ll_alpha$row_data$targetSequences == nuc_seq)[1]
#   if(chain == "beta") idx = which(ll_beta$row_data$targetSequences == nuc_seq)[1]
#
#   cor_df = ll_beta$row_data
#
#   idx = which(ll_alpha$row_data$targetSequences == nuc_seq)[1]
#
#   wa1 = ll_alpha$row_data$wa[idx]
#   alpha_vec = extract_col_dense(alpha_mat_rf, idx)
#
#   # cor_df = well_data[[paste("rows", other, sep = "_")]]
#   # mat_dense_t = t(as.matrix(mat_other))
#   # well_sums = Matrix::rowSums(mat_dense_t)
#   # mat_dense_t = mat_dense_t/well_sums
#   #
#   # cor_df$w_other = Matrix::rowSums(mat_other != 0)
#   # cor_df$w_chain = sum(!is.na(df_single$readFraction))
#
#   #mat_dense_t[mat_dense_t <= 0] = NA
#   #mat_log10 = log10(dense_mat)
#
#   res = sparse_col_cor(beta_mat_rf, alpha_vec)
#   cor_df = bind_cols(cor_df, res)
#
#   #cor_df$r_log_log = cor(mat_log10, log10(df_single$readFraction), method = "pearson", use = "pairwise.complete.obs")
#   cor_df$wij = sparse_overlap(beta_mat_rf, alpha_vec)
#   cor_df$wa = wa1
#   cor_df$p_adj = p.adjust(cor_df$p, method = "BH")
#
#   cor_df = cor_df %>% select(chain, row_id, aaSeqCDR3, r, t, p, p_adj, wa, wb, wij, n, targetSequences, v,j, everything())
#
#   cor_df2 = cor_df %>% filter(wij > 0.5*max(cor_df$wij)) %>% arrange(desc(r))
#
#   cor_df3 = cor_df2[1:n_plot,]
#
#   meta_gg = ll_beta$row_data[match(cor_df3$row_id, ll_beta$row_data$row_id),] %>%
#     mutate(chain_name = paste("beta", row_id))
#   mat_gg = as.matrix(beta_mat_rf[,cor_df3$row_id])
#
#   colnames(mat_gg) = meta_gg$targetSequences
#   rownames(mat_gg) = ll_beta$col_data$well
#   gg_df = lapply(1:n_plot, function(i) {
#     df_tmp = tibble(readFraction1 = mat_gg[,i],
#                     readFraction2 = alpha_vec,
#                     well = ll_beta$col_data$well,
#                     chain_name = meta_gg$chain_name[i]
#     )
#   }) %>% bind_rows()
#   gg_df$chain_name = factor(gg_df$chain_name, levels = meta_gg$chain_name)
#   gg_df$readFraction1[gg_df$readFraction1 == 0] = NA
#   #print(cor_df3)
#   gg = ggplot(gg_df, aes(x=readFraction1, y=readFraction2, label = well)) +
#     geom_point(alpha = 0.3) +
#     facet_wrap(~chain_name, scales = "free") +
#     scale_x_log10() +
#     scale_y_log10() +
#     theme_bw()
#   #print(gg)
#   gg2 = ggplot(cor_df2, aes(x=r, y=-log10(p), aa = aaSeqCDR3, v = v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
#   gg3 = ggplot(cor_df2, aes(x=r, y=-log10(p_adj), aa = aaSeqCDR3, v= v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
#   if(interactive) {
#     gg = plotly::ggplotly(gg)
#     gg2 = plotly::ggplotly(gg2)
#     gg3 = plotly::ggplotly(gg3)
#   }
#
#   return(list(cor_df = cor_df, cor_df2 = cor_df2, cor_df3 = cor_df3, gg = gg, gg2 = gg2, gg3 = gg3))
# }

## sparse_col_cor_old ----------------------

# sparse_col_cor_old <- function(X, y) {
#   stopifnot(is(X, "sparseMatrix"), length(y) == nrow(X))
#
#   n <- nrow(X)
#
#   # --- y stats ---
#   y_c  <- y - mean(y)          # centered y; sum(y_c) = 0
#   sd_y <- sqrt(sum(y_c^2) / (n - 1))
#
#   # --- Numerator: t(X_centered) %*% y_c ---
#   # Since sum(y_c) = 0:
#   # sum_i (X_ij - mean_j)(y_i - mean_y) = sum_i X_ij * y_c_i  (cross term vanishes)
#   num <- as.numeric(Matrix::crossprod(X, y_c))   # avoids forming t(X) densely
#
#   # --- Column std of X (sparse-aware, O(nnz)) ---
#   col_means  <- Matrix::colMeans(X)
#   col_sq_means <- Matrix::colMeans(X ^ 2)        # X^2 stays sparse for sparse X
#   sd_x <- sqrt(pmax(col_sq_means - col_means^2, 0) * n / (n - 1))
#
#   # --- Pearson r ---
#   r <- num / ((n - 1) * sd_x * sd_y)
#   r[sd_x == 0] <- NA            # constant columns → NA
#
#   setNames(r, colnames(X))
#   return(r)
# }

## extract_rows_by_index_parquet_old -------------------------------------

# extract_rows_by_index_parquet_old <- function(parquet_path, row_indices) {
#   # DuckDB uses 1-based row numbers via rowid (actually 0-based internally)
#   # We'll use a generated row number approach
#
#   con <- DBI::dbConnect(duckdb::duckdb())
#   on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
#
#   # Convert indices to a comma-separated string for SQL IN clause
#   indices_str <- paste(row_indices - 1L, collapse = ", ")  # convert to 0-based
#
#   query <- sprintf("
#     SELECT * EXCLUDE (rn)
#     FROM (
#       SELECT *, row_number() OVER () - 1 AS rn
#       FROM read_parquet('%s')
#     )
#     WHERE rn IN (%s)
#   ", parquet_path, indices_str)
#
#   DBI::dbGetQuery(con, query)
# }

## load_well_data_to_sparse_old --------------------------------------------

# load_well_data_to_sparse_old = function(folder_path, wells = get_well_subset(1:16,1:24),
#                           well_pos=3, chain = c("both", "alpha", "beta"), nproc = 1L,
#                           columns = NULL,
#                           max_files = Inf) {
#   chain = chain[1]
#
#   files = list.files(path = folder_path, full.names = F)
#   files_full = file.path(folder_path, files)
#   files_no_ending = gsub("\\.tsv", "", files)
#   #files_no_ending = gsub("\\.", "_", files_no_ending) ## replace periods with underscores if necessary
#
#   wells_all = sapply(strsplit(files_no_ending, "_"), function(x) x[[well_pos]])
#   ## make sure wells in the string are 2-3 characters long and only alphanumeric -- necessary for some data
#   wells_all = sapply(wells_all, function(x) substr(x, 1, 3))
#   wells_all = gsub("[^[:alnum:] ]", "", wells_all)
#
#   chains_all = sapply(strsplit(files_no_ending, "_"), function(x) rev(x)[[1]])
#   chains_all[chains_all=="TRA"] = "alpha"
#   chains_all[chains_all=="TRB"] = "beta"
#
#   wells_keep = wells_all %in% wells
#   if(chain %in% c("alpha","beta")) {
#     chains_keep = chains_all == chain
#   } else {
#     chains_keep = TRUE
#   }
#   keep = which(wells_keep & chains_keep)
#
#   meta_tmp = tibble(well = wells_all[keep], chain = chains_all[keep], file_short = files[keep], file_full = files_full[keep])
#   print(paste(dim(meta_tmp)[1], "files total"))
#   n_files = min(max_files, dim(meta_tmp)[1])
#   print(paste("Loading", n_files, "files"))
#   ll = parallel::mclapply(1:n_files, function(i) {
#     print(i)
#     file_full_tmp = meta_tmp$file_full[i]
#     file_short_tmp = meta_tmp$file_short[i]
#     well_tmp = meta_tmp$well[i]
#     chain_tmp = meta_tmp$chain[i]
#     #print(file_full_tmp)
#     df_tmp = data.table::fread(file_full_tmp, select = columns) %>%
#       mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
#       mutate(well_chain = paste(well, chain))
#     return(df_tmp)
#   }, mc.cores = as.integer(nproc))
#   #names(ll) = files_no_ending[keep]
#   return(ll)
# }

# Python/Basilisk helpers ----------------------

## install_python_env ----------------

#' Install python environment
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function install a python environment via the `basilisk` package that is
#' needed for some functions: \code{\link{TCRdist}()}, \code{\link{cluster_tcrs}()}, and \code{\link{run_pairing}()}.
#'
#' It will run automatically when needed, but is available to the user if needed.
#'
#' @param force whether to force installation (default if FALSE)
#'
#' @return returns NULL
#'
#' @family python-deps
#' @noRd

# install_python_env <- function(force = FALSE) {
#   message("Installing/updating the Python environment for TIRTLtools...")
#   install_python_env_silently(force = force)
#   message("Done.")
#   return(invisible(NULL))
# }

## install_python_env_silently -----------------
# Internal: trigger basilisk env creation without doing real work
# install_python_env_silently <- function(force = FALSE) {
#   if (env_is_installed() && !force) {
#     return(invisible(TRUE))
#   }
#   basilisk::basiliskRun(
#     env = TIRTLtools_py_env,
#     fun = function() NULL
#   )
#   #mark_env_installed()
#   invisible(TRUE)
# }

## ensure_python_env -------------
# ensure_python_env <- function(ask = TRUE) {
#
#   if (env_is_installed()) {
#     return(invisible(TRUE))
#   }
#
#   # Non-interactive or explicit "don't ask": auto-install
#   if (!interactive() || !isTRUE(ask)) {
#     install_python_env_silently()
#     return(invisible(TRUE))
#   }
#
#   # Interactive + ask = TRUE: prompt the user
#   version = packageVersion("TIRTLtools")
#   message("This function requires installation of a Python environment.")
#   message(paste("This is the first time you are using the Python backend for TIRTLtools v", version, sep = ""))
#
#
#   ans <- utils::menu(
#     choices = c("Yes, install now", "No, cancel"),
#     title = "Allow TIRTLtools to install its Python environment now?"
#   )
#
#   if (ans != 1) {
#     stop("Installation of the Python environment was cancelled by the user.")
#   }
#
#   install_python_env_silently()
#   invisible(TRUE)
# }

