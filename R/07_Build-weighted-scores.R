#LASSO_CWC_final_coeffs = tar_read(LASSO_CWC_final_coeffs)
#betas_df = LASSO_CWC_final_coeffs$coefficients
#abund_df = tar_read(Final_metabs_long)
#LASSO_res = LASSO_CWC_final_coeffs$coefficients
#mapping_file = tar_read(Mapping_file)
#id_col = "idno"
#metabolite_col = "Metabolite"
#beta_col = "Estimate"
#time_col = "exam"
#score_name = "RepeatedMeasures_score"
#normalize_names  = FALSE
#verbose = TRUE
#na_rm = TRUE              
#min_non_missing = 1


build_weighted_score_with_LASSO <- function(
    betas_df,
    abund_df,
    mapping_file,
    LASSO_res,
    id_col,
    metabolite_col,
    beta_col,
    time_col,
    score_name,
    normalize_names  = FALSE,
    verbose = TRUE,
    na_rm = TRUE,              # ignore NA abundances (treat as zero contribution)
    min_non_missing       # require at least this many non-missing metabolites per row
) {
  
  #----------------Helper -----------------------------#
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  
  #----------------Select betas and Ps -----------------------------#
  df_b <- as.data.frame(betas_df)|>
    tibble::rownames_to_column(var = "Metabolite") |>
  dplyr::select(metabolite_col, dplyr::all_of(beta_col))
                
  names(df_b) <- c("target", "beta")
  
  
  #----------------Select hits -----------------------------#
  
  LASSO_sig <- as.data.frame(LASSO_res)|>
    tibble::rownames_to_column(var = "Metabolite") |>
    dplyr::filter(Estimate !=0) |>
    dplyr::filter(Metabolite %in% mapping_file$Metabolite)
  
  if (nrow(df_b) == 0L) stop("No metabolites pass the p-value threshold.")
  # QC info
  Sig_targets_n <- dim(df_b)[1]
  Sig_targets_names <- as.data.frame(df_b$target)
  
  #----------------Check for duplicates in raw data -----------------------------#
  dupe_counts <- as.data.frame(table(df_b$target), stringsAsFactors = FALSE)
  names(dupe_counts) <- c("target", "n_rows"); dupe_counts$n_rows <- as.integer(dupe_counts$n_rows)
  n_with_duplicates <- sum(dupe_counts$n_rows > 1L)
  n_total_dupe_rows <- sum(pmax(dupe_counts$n_rows - 1L, 0L))
  if (n_with_duplicates > 0L) {
    dup_preview <- dupe_counts[dupe_counts$n_rows > 1L, ]
    dup_preview <- dup_preview[order(-dup_preview$n_rows, dup_preview$target), ]
    preview_str <- paste(utils::head(sprintf("%s (n=%d)", dup_preview$target, dup_preview$n_rows), 10L), collapse = "; ")
    msg("Found %d metabolites with duplicate rows among P<%g; collapsing %d extra rows (mean beta).",
        n_with_duplicates, n_total_dupe_rows)
    msg("Duplicates (up to 10): %s%s", preview_str, if (nrow(dup_preview) > 10L) " ..." else "")
  }
  
  # collapse dupes
  df_b <- aggregate(beta ~ target, data = df_b, FUN = function(z) mean(z, na.rm = TRUE))
  
  #----------------Make scores -----------------------------#
  target_cols_in_abund <- setdiff(names(abund_df), c(id_col, time_col))
  common_targets <- intersect(df_b$target, target_cols_in_abund)
  
  
  missing_targets <- setdiff(df_b$target, target_cols_in_abund)
  nonsig_targets   <- setdiff(target_cols_in_abund, df_b$target)
  
  
  # align to common_targets 
  df_b_common <- df_b[match(common_targets, df_b$target), , drop = FALSE]
  X <- as.matrix(abund_df[, common_targets, drop = FALSE]) #is abundances
  w <- as.numeric(df_b_common$beta) #is betas
  X2 <-scale(X) #Center and scale before use
  
  if (na_rm) {
    # elementwise multiply scaled vars, then rowMeans with na.rm=TRUE
    WX <- t( t(X2) * w )
    score <- rowMeans(WX, na.rm = TRUE)
    
    
    # optional: guard against rows that are almost entirely missing (below missing threshold set in function)
    n_nonmiss <- rowSums(!is.na(X))
    too_sparse <- n_nonmiss < min_non_missing
    if (any(too_sparse)) {
      msg("Warning: %d rows had < %d non-missing metabolites; setting score to NA for those rows.",
          sum(too_sparse), min_non_missing)
      score[too_sparse] <- NA_real_
      too_sparse_excl <- abund_df[too_sparse,]
      too_sparse_excl <- too_sparse_excl[,names(too_sparse)==id_col | names(too_sparse)==time_col]
    }
    
  } else {
    # any NA -> NA
    #score <- as.numeric(X %*% w)
    score <- NA
  }
  
  
  # build output: include id and (optionally) time
  if (requireNamespace("rlang", quietly = TRUE)) {
    out <- tibble::tibble(
      !!rlang::sym(id_col) := abund_df[[id_col]],
      !!rlang::sym(score_name) := score
    )
    if (!is.null(time_col)) {
      out[[time_col]] <- abund_df[[time_col]]
      out <- out[, c(id_col, time_col, score_name)]
    }
  } else {
    out <- tibble::tibble(tmp_id = abund_df[[id_col]], tmp_score = score)
    names(out) <- c(id_col, score_name)
    if (!is.null(time_col)) {
      out[[time_col]] <- abund_df[[time_col]]
      out <- out[, c(id_col, time_col, score_name)]
    }
  }
  
  
  #----------------Score info -----------------------------#
  Score_info = list(Sig_targets_n = Sig_targets_n,
                    Sig_targets_names = Sig_targets_names,
                    Duplicate_targets_n = n_with_duplicates,
                    Included_targets = common_targets,
                    Missing_targets = missing_targets,
                    too_sparse_excl = too_sparse_excl,
                    min_non_missing = min_non_missing)
  
  #----------------All Outputs -----------------------------#
  list(score_info = Score_info, 
       scores = out)
}