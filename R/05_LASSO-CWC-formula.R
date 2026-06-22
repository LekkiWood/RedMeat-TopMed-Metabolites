#Metabs_table <- tar_read(Metabs_long_clean)
#Metabs_table <- Metabs_table[,1:150]
#Traits_table <- tar_read(Traits_long)
#missing_threshold = 0.005

LASSO_CWC_formula_function <- function(Metabs_table, Traits_table, missing_threshold) {

 
  ##############################
  ## 1. Create data
  ##############################
  
  ##############################
  ## 1. Create data
  ##############################
  
  common_keys <- intersect(names(Metabs_table), names(Traits_table))
  Analysis_df <- dplyr::left_join(Metabs_table, Traits_table, by = common_keys)
  
  
  #####Select vars and deal with missing data
  
  Penalized <- Analysis_df |>
    dplyr::filter(!is.na(redmeat_cwc_z) & (exam==1 | exam==5)) |> 
    dplyr::group_by(idno) |>
    dplyr::mutate(
      redmeat_pm  = ifelse(exam == 1 | exam==5, mean(redmeat, na.rm=TRUE), NA_real_), 
      redmeat_cwc = ifelse(exam == 1 | exam==5, redmeat - redmeat_pm, NA_real_) ) |>
    dplyr::select(-redmeat_pm) |>
    dplyr::ungroup()
  
  
  # Calculate missing value proportions and identify columns to remove
  cols_to_keep1 <- Metabs_table |>
    #dplyr::select(-subject_id, -TOM_ID, -sidno, -redmeat, -age, -egfr, -DM, -income, -PA, -ldl, -insulin, -education, -energy, -gender, -exam, -smoking, -race, -gender) |>
    dplyr::select(-sidno, -exam, -subject_id, -TOM_ID, -idno) |>
                    dplyr::summarise(across(everything(), ~ mean(is.na(.)))) |>
    tidyr::pivot_longer(everything(), names_to = "variable", values_to = "missing_rate") |>
    dplyr::filter(missing_rate < missing_threshold) |>
    dplyr::pull(variable)
  
  cols_to_keep2 <- c("sidno", "age", "egfr", "DM", "income", "PA", "education", "redmeat_cwc_z", "energy", "gender", "exam", "race", "smoking")
  
  cols_to_keep <- c(cols_to_keep2, cols_to_keep1)
  
  Penalized_na <- Penalized[cols_to_keep]

  
  Penalized_na <- Penalized_na |>
    dplyr::mutate(across(c(sidno, gender, exam, race), as.factor)) |>
    # Create dummy variables for categorical variables
    dplyr::mutate(
      gender_num  = ifelse(gender == "Male", 1, 0),    # adjust to match your coding
      exam_num    = as.numeric(exam),                  # if exam is ordinal (1,5), this is fine
      is_Chinese  = ifelse(race == "Chinese-American", 1, 0),
      is_African  = ifelse(race == "African-American", 1, 0),
      is_Hispanic = ifelse(race == "Hispanic", 1, 0),
      is_nonsmoker  = ifelse(smoking == "Never", 1, 0),
      is_formersmoker  = ifelse(smoking == "Former", 1, 0)
      
    ) |>
    # Now scale *only* numeric variables that are not dummy-coded and not predictor (already scaled)
    dplyr::mutate(across(
      where(is.numeric) &
        !matches("gender_num|exam_num|is_Chinese|is_African|is_Hispanic|is_nonsmoker|is_formersmoker|redmeat_cwc_z"),
      ~ as.numeric(scale(.x))
    ))
  
  gc()
  
  
  lasso_names_all <- setdiff(names(Penalized_na),
                             c("idno","TOM_ID","subject_id",
                               "hdl","TG","glucose","processedredmeat",
                               "redmeat_z","redmeat_pm_z","redmeat",
                               "redmeat_gm_z","redmeat_pm_cgm_z","redmeat_cgm_z",
                               "educ", "race", "gender", "exam", "site", "smoking", "DM", "race", "sidno", "redmeat", "redmeat_cwc", "redmeat_cwc_z"))
  


  f_all <- (paste(lasso_names_all, collapse = " + "))
  redmeat_formula_all <- as.formula(paste("redmeat_cwc ~", f_all))
 
  

list(LASSO_formula  = redmeat_formula_all,
     LASSO_names = lasso_names_all) 

  
}
