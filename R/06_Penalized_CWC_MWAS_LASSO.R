#rm(list=ls())
#Metabs_table   = tar_read(Metabs_long_clean)  
#Traits_table   = tar_read(Traits_long)
#lambda         = 1 # full grid directly
#kk             = 5
#formula        = as.formula("redmeat_cwc ~ age + egfr + income + PA + education + energy + Amide_X1.5.AG...1.deoxyglucose")
#log_file       = "progress/test.log"
#lasso_names_all = c("redmeat_cwc", "age", "egfr", "income", "PA", "education", "energy", "Amide_X1.5.AG...1.deoxyglucose")
#missing_threshold = 0.005

parallel_CWC_LASSO <- function(Metabs_table, Traits_table, lambda, kk, formula, log_file, lasso_names_all, missing_threshold) {
  
  ## --- setup ---
  options(future.globals.maxSize = 10 * 1024^3)
  dir.create("progress", showWarnings = FALSE, recursive = TRUE)
  if (is.null(log_file) || !nzchar(log_file)) log_file <- "progress/lasso_progress.log"
  suppressMessages(requireNamespace("pryr", quietly = TRUE))
  set.seed(11042012)
  
  ## 1) Join & prepare ---------------------------------------------------------
  ##############################
  ## 1. Create data
  ##############################
  
  common_keys <- intersect(names(Metabs_table), names(Traits_table))
  Analysis_df <- dplyr::left_join(Metabs_table, Traits_table, by = common_keys)
  
  
  #####Select vars and deal with missing data
  
  Penalized <- Analysis_df |>
    dplyr::group_by(idno) |>
    dplyr::mutate(
      redmeat_pm  = ifelse(exam == 1 | exam==5, mean(redmeat, na.rm=TRUE), NA_real_), 
      redmeat_cwc = ifelse(exam == 1 | exam==5, redmeat - redmeat_pm, NA_real_) ) |>
    dplyr::select(-redmeat_pm) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(redmeat_cwc) & (exam==1 | exam==5))
  
  
  # Calculate missing value proportions and identify columns to remove
  cols_to_keep1 <- Metabs_table |>
    #dplyr::select(-subject_id, -TOM_ID, -sidno, -redmeat, -age, -egfr, -DM, -income, -PA, -ldl, -insulin, -education, -energy, -gender, -exam, -smoking, -race, -gender) |>
    dplyr::select(-sidno, -exam, -subject_id, -TOM_ID, -idno) |>
    dplyr::summarise(across(everything(), ~ mean(is.na(.)))) |>
    tidyr::pivot_longer(everything(), names_to = "variable", values_to = "missing_rate") |>
    dplyr::filter(missing_rate < missing_threshold) |>
    dplyr::pull(variable)
  
  cols_to_keep2 <- c("sidno", "age", "egfr", "DM", "income", "PA", "education", "redmeat_cwc", "energy", "gender", "exam", "race", "smoking")
  
  cols_to_keep <- c(cols_to_keep2, cols_to_keep1)
  
  Penalized_na <- Penalized[cols_to_keep]
  
  Penalized_na <- Penalized_na[complete.cases(Penalized_na),]

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
    # Now scale *only* numeric variables that are not dummy-coded
    dplyr::mutate(across(
      where(is.numeric) &
        !matches("gender_num|exam_num|is_Chinese|is_African|is_Hispanic|is_nonsmoker|is_formersmoker"),
      ~ as.numeric(scale(.x))
    ))
  
  gc()
  
  
  lasso_names_all <- setdiff(names(Penalized_na),
                             c("idno","TOM_ID","subject_id",
                               "hdl","TG","glucose","processedredmeat",
                               "redmeat_z","redmeat_pm_z",
                               "redmeat_gm_z","redmeat_pm_cgm_z","redmeat_cgm_z",
                               "educ", "race", "gender", "exam", "site", "smoking", "DM", "race"))

  
  Penalized_na <- Penalized_na |>
    dplyr::select(all_of(lasso_names_all))


  
  ## 2) Parallel setup ---------------------------------------------------------
  # Choose your parallel plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multicore, workers = max(1, 5))
  
  
  N<-dim(Penalized_na)[1]
  ind<-sample(N,N)
  family = gaussian(link="identity")
  
  ## set number of folds
  nk <- floor(N/kk)
  
  Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))
  
  ## first fit good starting model
  PQL <- nlme::lme(redmeat_cwc ~ 1, random=~1 | sidno, data=Penalized_na)
  
  Delta.start<-c(as.numeric(PQL$coef$fixed),rep(0,length(lasso_names_all)),as.numeric(t(PQL$coef$random$sidno)))
  
  Q.start<-as.numeric(lme4::VarCorr(PQL)[1,1])
  
  
  #How long to run
  start_lasso_full_tuning <- Sys.time()
  
  for(j in 1:length(lambda))
  {
    print(paste("Iteration ", j,sep=""))
    
    for (i in 1:kk)
    {
      if (i < kk)
      {
        indi <- ind[(i-1)*nk+(1:nk)]
      }else{
        indi <- ind[((i-1)*nk+1):N]
      }
      
      penalized.train <- Penalized_na[-indi,]
      penalized.test <- Penalized_na[indi,]
      
   
      # ---- Start message ----
     
      msg_start <- paste0(format(Sys.time(), "%a %b %d %X %Y"), " lambda=", lambda[j], " fold=", kk, " START"
      )
      cat(msg_start, "\n", file = log_file, append = TRUE)
      message(msg_start)
      # -----------------------
      
      
      glm2 <- try(glmmLasso::glmmLasso(formula, rnd = list(sidno=~1),  
                            family = family, data = penalized.train, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
                            control = list(print.iter=TRUE))
                  ,silent=FALSE) 
      
      
      if(!inherits(glm2, "try-error"))
      {  
        penalized.test2 <- as.data.frame(penalized.test)
        penalized.test2$sidno <- factor(
          as.character(penalized.test2$sidno),
          levels = levels(penalized.train$sidno))
          
        y.hat<-predict(glm2,penalized.test2)    
        
        Devianz_ma[j,i]<-sum(family$dev.resids(penalized.test$redmeat_cwc,y.hat,wt=rep(1,length(y.hat))))
      }
    }
    print(sum(Devianz_ma[j,]))
    
    # ---- End message  ----
    
    msg_done <- paste0(format(Sys.time(), "%a %b %d %X %Y"), " lambda=", lambda[j], " fold=", i, " DONE"
    )
    cat(msg_done, "\n", file = log_file, append = TRUE)
    message(msg_done)
    # ---------------------------------------
  }
  
  
  end_lasso_full_tuning <- Sys.time()
  

  
  time_lasso_full_tuning <- end_lasso_full_tuning - start_lasso_full_tuning ####Time difference o
  

  Devianz_vec<-apply(Devianz_ma,1,sum)
  
  opt2<-which.min(Devianz_vec)
  
  lasso_final <- glmmLasso(formula, rnd = list(sidno=~1),  
                                     family = family, data = Penalized_na, lambda=lambda[opt2],switch.NR=FALSE,final.re=FALSE,
                                     control = list(print.iter=TRUE))
  
  
  
  
  

  
  list(opt_lambda = lambda[opt2],
       Deviance = Devianz_ma,
       Runtime = time_lasso_full_tuning,
       final_model = lasso_final,
       final_coeffs = summary(lasso_final))
}
