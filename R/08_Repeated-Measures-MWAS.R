#rm(list=ls())
#library(targets)
#library(tidyverse)
#library(lmerTest)
#cleaned_metabs = tar_read(Metabs_long_clean)
#metab_mapping = tar_read(Mapping_file)
#traits_db = tar_read(Traits_long)
#numeric_covariates = c("age", "PA", "egfr", "energy", "income")
#factor_covariates = c("gender", "race", "DM",  "smoking", "site")
#predictor = "redmeat"
#lasso_res = targets::tar_read(LASSO_CWC_final_coeffs_db)
#metab <- "Amide_X1.5.AG...1.deoxyglucose"
#lmer(scale(Amide_X1.5.AG...1.deoxyglucose) ~ scale(totalfruit) + scale(age) +  scale(PA) + scale(eGFR) + scale(energy) + sex + race + edu + diabetes + smoking + (1 | idno), data=MWAS_file)


repeated_measures_MWAS_function <- function(cleaned_metabs, metab_mapping,
                                       traits_db, numeric_covariates,
                                       factor_covariates, predictor, lasso_res) {
  
  
  #---------------------------------------------------------------------------------------#
  #--------------------------------1. Build data table------------------------------------#
  #---------------------------------------------------------------------------------------#
  common_keys <- intersect(names(cleaned_metabs), names(traits_db))
  
  MWAS_file <- traits_db |>
    dplyr::left_join(cleaned_metabs, by = common_keys)
  
  inc_metabs <- lasso_res |>
    tibble::rownames_to_column(var = "Metabolite") |>
    dplyr::filter(Metabolite %in% names(cleaned_metabs) & Estimate !=0) |>
    dplyr::pull(Metabolite)
  
  
  #metabs <- metab_mapping |>
   # dplyr::filter(Metabolite %in% names(MWAS_file)) |>
    #dplyr::pull(Metabolite)
  
  stopifnot(all(c(predictor, factor_covariates, numeric_covariates) %in% names(MWAS_file)))
  
  # keep subjects with all exams 

  
  dat <- MWAS_file |>
    dplyr::group_by(idno) |>
    dplyr::filter(exam == 1 | exam ==5) |>
    dplyr::select(idno, exam, all_of(inc_metabs), all_of(predictor), all_of(numeric_covariates), 
                  all_of(factor_covariates)) |>
    dplyr::filter(n_distinct(exam) >= 2) |>
    dplyr::ungroup() |>
    dplyr::mutate(idno = as.factor(idno),
                  exam = as.factor(exam))
  

 
  #---------------------------------------------------------------------------------------#
  #--------------------------------2. Run MWAS--------------------------------------------#
  #---------------------------------------------------------------------------------------#
  out <- data.frame(
    Metab = inc_metabs,
    Nobs = NA_integer_,
    Beta = NA_real_,
    SE = NA_real_,
    P = NA_real_
    )
  
  for (k in seq_along(inc_metabs)) {
    
    metab <- inc_metabs[k]
    
    # RHS: predictor * time + other covariates
    fixed_terms <- c(
      paste0("scale(", predictor, ")"),
      paste0("scale(", numeric_covariates, ")"),
      factor_covariates
    )
    
    rhs_fixed <- paste(fixed_terms, collapse = " + ")
    fml_fixed <- as.formula(paste0("scale(rcompanion::blom(", metab, ")) ~ ", rhs_fixed) )
    fml_mixed <- as.formula(paste0("scale(rcompanion::blom(", metab, ")) ~ ", rhs_fixed, " + (1 | idno)"))
    
    res <- fit_mixed_or_fixed(fml_mixed, fml_fixed, dat)
    fit <- res$fit
    
    out$Nobs[k]  <- stats::nobs(fit)
    out$Model[k] <- res$model_type
    
    # term names to extract
    pred_term   <- paste0("scale(", predictor, ")")
    pred_est   <- extract_term(fit, pred_term)
    
    out[k, c("Beta", "SE", "P")] <- pred_est

  }
  
  out$P_fdr <- p.adjust(out$P, method = "fdr")
  
  #---------------------------------------------------------------------------------------#
  #--------------------------------3. output  --------------------------------------------#
  #---------------------------------------------------------------------------------------#
  
  out
}