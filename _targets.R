############################################
#                                          #
#            Overview                      #
#                                          #
###########################################
# 
# 
# 
#                   ┌─────────────────────────────────────┐
#                   │ Raw Metabolomics Input Files        │
#                   │  (path_amide, sample_infoo, etc.)   │
#                   └──────────────┬──────────────────────┘
#                                  │
#                                  ▼
#                        ┌──────────────────┐
#                      │ build_metabs_out │
#                      │ (build_metabs_function) │
#                      └──────────────┬──────────┘
#                                     │
#                                     ▼
#   ┌──────────────────────────────┬────────────────────────────────────┐
#   │                              │                                    │
#   ▼                              ▼                                    ▼
#Metabs_long_clean      MetabBuilding_QC, Mapping_file, etc.      CV_and_missingness_file
#   │
#   │
#   │
#   ▼
#┌──────────────────────────────┐
#│ build_traits_out             │
#│ (build_traits)               │
#└──────────────┬───────────────┘
#               │
#               ▼
#         Traits_long
#               │
#               │
#               ▼
#┌──────────────────────────────────────────┐
#│ LASSO_formula (LASSO_formula_function)   │
#│ → Builds redmeat_formula_all             │
#└──────────────┬───────────────────────────┘
#               │
#               ▼
#   ┌──────────────────────────────────────────┐
#   │ lambda_chunks (split from λ = 1:600)     │
#   │ kk = 5 (cross-validation folds)          │
#   └───────────────────┬──────────────────────┘
#                       │
#                       ▼
#          ┌────────────────────────────────────┐
#          │ lasso_chunk_results                │
#          │ (parallel_LASSO)                   │
#          │ pattern = map(lambda_chunk = λ_chunk) │
#          │ Runs each λ subset in parallel,     │
#          │ logs → progress/lasso_chunk_X.log   │
#          └───────────────────┬────────────────┘
#                              │
#                              ▼
#                ┌───────────────────────────────┐
#                │ lasso_combined                │
#                │ (combine_results)             │
#                │ Merges Deviance matrices,     │
#                │ finds global best λ, refits   │
#                │ final model                   │
#                └──────────────┬────────────────┘
#                               │
#                               ▼
#      ┌─────────────────────────────────────────────────────┐
#      │                    LASSO Outputs                   │
#      │────────────────────────────────────────────────────│
#      │ LASSO_optimal_lambda ← lasso_combined$global_opt_λ  │
#      │ LASSO_deviance       ← lasso_combined$Devianz_all   │
#      │ LASSO_optimal_model_output ← lasso_combined$final_model │
#      │ LASSO_coeffs         ← lasso_combined$final_coeffs   │
#      │ lasso_diagnostics    ← diagnostic plot of Deviance   │
#      └─────────────────────────────────────────────────────┘


############################################
#                                         #
#            Script.                      #
#                                         #
###########################################

library('targets')
library('tarchetypes')
library('crew')

tar_option_set(
  controller = crew_controller_local(workers = 3),                # or fewer if RAM is tight
  memory = "transient",       # don't keep targets in memory
  garbage_collection = TRUE   # free memory aggressively
)



Sys.setenv(VROOM_CONNECTION_SIZE = as.character(10 * 1024 * 1024))

#set.seed(11042012)
tar_option_set(packages = c("dplyr", "tidyr", "tibble", "readr", "data.table", "bit64", "foreign", "quarto", "misty", "lmerTest", 
                            "rlang", "purrr", "EnvStats", "mnormt", "MASS", "nlme", "glmmLasso", "dplyr", "future.apply", "progressr", 
                            "furrr", "purrr", "digest", "glmnet", "doMC", "future", "nlme"))

tar_source("/media/Analyses/RedMeat-TopMed-Metabolites/R")

tar_option_set(seed = 11042012)


list(
  
  
  #--------------------------------------------------------------#
  #---------- Make metabolite info.  ----------------------------#
  #--------------------------------------------------------------#
  
  
  #Abundance tables
  tar_target(path_amide, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/25_0107_TOPMed_MESA_Amide-neg_rev031325.csv" , format = "file"),
  tar_target(path_C8, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C8-pos_checksums_rev031325.csv", format = "file"),
  tar_target(path_C18, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C18-neg_checksums_rev031325.csv", format = "file"),
  tar_target(path_HILIC,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_HILIC-pos_checksums_rev031325.csv", format = "file"),
  
  #Info tables
  tar_target(path_amide_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_AmideNeg_SampleInfo_20250329.txt", format = "file"),
  tar_target(path_C8_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C8Pos_SampleInfo_20250329.txt", format = "file"),
  tar_target(path_C18_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C18Neg_SampleInfo_20250329.txt", format = "file"),
  tar_target(path_HILIC_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_HILIC-Pos_SampleInfo_20250329.txt", format = "file"),
  
  #Bridging file
  tar_target(path_bridge,"/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv", format = "file"),
  
  
  # Compute metabolite information
  tar_target(build_metabs_out,   # <- rename
             build_metabs_function(
               path_amide = path_amide,
               path_C8 = path_C8,
               path_C18 = path_C18,
               path_HILIC = path_HILIC,
               path_amide_info = path_amide_info,
               path_C8_info = path_C8_info,
               path_C18_info = path_C18_info,
               path_HILIC_info = path_HILIC_info,
               path_bridge = path_bridge,
               QC_label = "QC-pooled_ref")
  ),
  
  #QC output for metab table
  tar_target(MetabBuilding_QC, build_metabs_out$MetabBuilding_QC),
  #Metabs file - uncleaned
  tar_target(Metabs_long, build_metabs_out$Metabs_long),
  #QC output for mapping file
  tar_target(MappingQC_info, build_metabs_out$MappingQC_info),
  #Mapping file
  tar_target(Mapping_file, build_metabs_out$Mapping_final),
  #CV and missingness
  tar_target(CV_and_missingness_file, build_metabs_out$CV_missing),
  #Duplicate info
  tar_target(Duplicates_file, build_metabs_out$Duplicates_info_file),
  #Final_metabs_file
  tar_target(Metabs_long_clean, build_metabs_out$Metabs_long_clean),
  #--------------------------------------------------------------#
  #---------- Make traits ---------------------------------------#
  #--------------------------------------------------------------#
  
  #Inputs
  tar_target(path_e1_ffq_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1_FFQ_20160520.dta", format = "file"),
  tar_target(path_e1_nutr_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1_Nutrients_20151208.dta", format = "file"),
  tar_target(path_e5_ffq_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_FFQ_20140130.dta", format = "file"),
  tar_target(path_e5_nutr_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_Nutrient_20140416.dta", format = "file"),
  tar_target(path_exam1_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta", format = "file"),
  tar_target(path_exam5_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_FinalLabel_20140613.dta", format = "file"),
  tar_target(path_exam6_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe6_FinalLabel_20220513.dta", format = "file"),
  tar_target(path_cognition_file, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESA_MA_E7_worse_cogdx_02212025.csv", format = "file"),
  
  # Compute table of traits
  tar_target(build_traits_out,   # <- rename
             build_traits(
               path_bridge,
               path_e1_ffq_file,
               path_e1_nutr_file,
               path_e5_ffq_file,
               path_e5_nutr_file, 
               path_exam1_file, 
               path_exam5_file,
               path_exam6_file, 
               path_cognition_file 
             )
  ),
  
  #QC output
  tar_target(Build_traits_QCinfo, build_traits_out$QC_info_out),
  #Traits table
  tar_target(Traits_long, build_traits_out$Covs),
  
  #--------------------------------------------------------------------------#
  # -----------------------Run ML LASSO  ------------------------------------#
  # -----------------------Repeated Measures  -------------------------------#
  #--------------------------------------------------------------------------#
  
  #---------Set up LASSO  ---------#
  
  #---------Create formula & info ---------#
  
  tar_target(LASSO_info, LASSO_formula_function(Metabs_table = Metabs_long_clean, 
                                                Traits_table = Traits_long,
                                                missing_threshold = 0.005)),
  
  tar_target(LASSO_formula, LASSO_info$LASSO_formula),
  tar_target(LASSO_vars, LASSO_info$lasso_names),
  
  #---------Run 0-90 in one chunk ---------# 
  tar_target(
    lasso_results_0,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 0, to = 90, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_0_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_0_optimal_lambda, lasso_results_0$opt_lambda),
  
  #---------Run 100-190 in one chunk ---------# 
  tar_target(
    lasso_results_1,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 100, to = 190, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_1_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_1_optimal_lambda, lasso_results_1$opt_lambda),
  
  #---------Run 200-290 in one chunk ---------# 
  tar_target(
    lasso_results_2,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 200, to = 290, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_2_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_2_optimal_lambda, lasso_results_2$opt_lambda),
  
  #---------Run 300-390 in one chunk ---------# 
  tar_target(
    lasso_results_3,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 300, to = 390, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_3_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_3_optimal_lambda, lasso_results_3$opt_lambda),
  
  #---------Run 400-490 in one chunk ---------# 
  tar_target(
    lasso_results_4,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 400, to = 490, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_4_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_4_optimal_lambda, lasso_results_4$opt_lambda),
  
  
  #---------Run 500-590 in one chunk ---------# 
  tar_target(
    lasso_results_5,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 500, to = 590, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_5_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_5_optimal_lambda, lasso_results_5$opt_lambda),
  
  #---------Run 600-690 in one chunk ---------# 
  tar_target(
    lasso_results_6,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 600, to = 690, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_6_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_6_optimal_lambda, lasso_results_6$opt_lambda),
  
  #---------Run 700-790 in one chunk ---------# 
  tar_target(
    lasso_results_7,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 700, to = 790, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_7_full.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_7_optimal_lambda, lasso_results_7$opt_lambda),
  
  #-----------------------------------#
  #---------Run final model ----------# 
  #-----------------------------------#
  
  tar_target(
    lasso_final_results,
    parallel_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = c(LASSO_0_optimal_lambda, LASSO_1_optimal_lambda, LASSO_2_optimal_lambda, LASSO_3_optimal_lambda, LASSO_4_optimal_lambda, 
                         LASSO_5_optimal_lambda, LASSO_6_optimal_lambda, LASSO_7_optimal_lambda),
      kk             = 5,
      formula        = LASSO_formula,
      log_file       = "progress/lasso_final.log",
      lasso_names_all = LASSO_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_final_optimal_lambda, lasso_final_results$opt_lambda),
  tar_target(LASSO_final_deviance,       lasso_final_results$Deviance),
  tar_target(LASSO_final_optimal_model,  lasso_final_results$final_model),
  tar_target(LASSO_final_coeffs,         lasso_final_results$final_coeffs),
  tar_target(LASSO_final_coeffs_db,     as.data.frame(LASSO_final_coeffs$coefficients)),
  
  #---------Save coeffs ---------#
  
  #save CV file for use
  tar_target(LASSO_filename, paste0("Redmeat_LASSO_coeffs_", Sys.Date(), ".csv")),
  
  
  #export CV as csv
  tar_target(LASSO_csv,
             {
               out_dir  <- "outputs"
               dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
               out_path <- file.path(out_dir, LASSO_filename)
               write.csv(LASSO_final_coeffs_db, out_path)
               out_path
             },
             format = "file"
  ),

  #--------------------------------------------------------------------------#
  # --------------Run CWC LASSO  (within effects) ---------------------------#
  #--------------------------------------------------------------------------#
  
  #---------Set up LASSO  ---------#
  
  #---------Create formula & info ---------#
  
  tar_target(LASSO_CWC_info, LASSO_CWC_formula_function(Metabs_table = Metabs_long_clean, 
                                                        Traits_table = Traits_long,
                                                        missing_threshold = 0.005)),
  
  tar_target(LASSO_CWC_formula, LASSO_CWC_info$LASSO_formula),
  tar_target(LASSO_CWC_vars, LASSO_CWC_info$LASSO_names),
  
  #---------Run 0-90 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_0,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 0, to = 90, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_0_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_0_optimal_lambda, LASSO_CWC_results_0$opt_lambda),
  
  #---------Run 100-190 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_1,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 100, to = 190, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_1_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_1_optimal_lambda, LASSO_CWC_results_1$opt_lambda),
  
  #---------Run 200-290 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_2,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 200, to = 290, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_2_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_2_optimal_lambda, LASSO_CWC_results_2$opt_lambda),
  
  #---------Run 300-390 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_3,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 300, to = 390, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_3_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_3_optimal_lambda, LASSO_CWC_results_3$opt_lambda),
  
  #---------Run 400-490 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_4,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 400, to = 490, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_4_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_4_optimal_lambda, LASSO_CWC_results_4$opt_lambda),
  
  
  #---------Run 500-590 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_5,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 500, to = 590, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_5_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_5_optimal_lambda, LASSO_CWC_results_5$opt_lambda),
  
  #---------Run 600-690 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_6,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 600, to = 690, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_6_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_6_optimal_lambda, LASSO_CWC_results_6$opt_lambda),
  
  #---------Run 700-790 in one chunk ---------# 
  tar_target(
    LASSO_CWC_results_7,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = seq(from = 700, to = 790, by = 10),  # full grid directly
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_7_full.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_7_optimal_lambda, LASSO_CWC_results_7$opt_lambda),
  
  #-----------------------------------#
  #---------Run final model ----------# 
  #-----------------------------------#
  
  tar_target(
    LASSO_CWC_final_results,
    parallel_CWC_LASSO(
      Metabs_table   = Metabs_long_clean,   
      Traits_table   = Traits_long,
      lambda         = c(LASSO_CWC_0_optimal_lambda, LASSO_CWC_1_optimal_lambda, LASSO_CWC_2_optimal_lambda, LASSO_CWC_3_optimal_lambda, LASSO_CWC_4_optimal_lambda, 
                         LASSO_CWC_5_optimal_lambda, LASSO_CWC_6_optimal_lambda, LASSO_CWC_7_optimal_lambda),
      kk             = 5,
      formula        = LASSO_CWC_formula,
      log_file       = "progress/LASSO_CWC_final.log",
      lasso_names_all = LASSO_CWC_vars,
      missing_threshold = 0.005
    )
  ),
  
  
  #---------Save output ---------#
  
  
  tar_target(LASSO_CWC_final_optimal_lambda, LASSO_CWC_final_results$opt_lambda),
  tar_target(LASSO_CWC_final_deviance,       LASSO_CWC_final_results$Deviance),
  tar_target(LASSO_CWC_final_optimal_model,  LASSO_CWC_final_results$final_model),
  tar_target(LASSO_CWC_final_coeffs,         LASSO_CWC_final_results$final_coeffs),
  tar_target(LASSO_CWC_final_coeffs_db,     as.data.frame(LASSO_CWC_final_coeffs$coefficients)),
  
  #---------Save coeffs ---------#
  
  #save CV file for use
  tar_target(LASSO_CWC_LASSO_filename, paste0("Redmeat_CWC_LASSO_coeffs_", Sys.Date(), ".csv")),
  
  
  #export CV as csv
  tar_target(LASSO_CWC_LASSO_csv,
             {
               out_dir  <- "outputs"
               dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
               out_path <- file.path(out_dir, LASSO_CWC_LASSO_filename)
               write.csv(LASSO_CWC_final_coeffs_db, out_path)
               out_path
             },
             format = "file"
  ),
  
  #--------------------------------------------------------------------------#
  # --------------Build scores  ---------------------------------------------#
  #--------------------------------------------------------------------------#
  
  
  #---------Repeated Measures score--------#
  
  tar_target(build_RM_score, build_weighted_score_with_LASSO(betas_df = LASSO_final_coeffs$coefficients,
                                                              abund_df = Metabs_long_clean,
                                                              LASSO_res = LASSO_final_coeffs$coefficients,
                                                              mapping_file = Mapping_file,
                                                              id_col = "idno",
                                                              metabolite_col = "Metabolite",
                                                              beta_col = "Estimate",
                                                              time_col = "exam",
                                                              score_name = "RepeatedMeasures_score",
                                                              normalize_names  = FALSE,
                                                              verbose = TRUE,
                                                              na_rm = TRUE,              
                                                              min_non_missing = 1)      
  ),
  
  tar_target(RM_score_info, build_RM_score$score_info),
  tar_target(RM_scores, build_RM_score$scores),
  
  #---------Within-Person score--------#
  
  tar_target(build_WP_score, build_weighted_score_with_LASSO(betas_df = LASSO_CWC_final_coeffs$coefficients,
                                                             abund_df = Metabs_long_clean,
                                                             LASSO_res = LASSO_CWC_final_coeffs$coefficients,
                                                             mapping_file = Mapping_file,
                                                             id_col = "idno",
                                                             metabolite_col = "Metabolite",
                                                             beta_col = "Estimate",
                                                             time_col = "exam",
                                                             score_name = "WithinPerson_score",
                                                             normalize_names  = FALSE,
                                                             verbose = TRUE,
                                                             na_rm = TRUE,              
                                                             min_non_missing = 1)      
  ),
  
  tar_target(WP_score_info, build_WP_score$score_info),
  tar_target(WP_scores, build_WP_score$scores),

#----------------------------------------------------------------------------#
#--------------------------------Quarto file---------------------------------#
#----------------------------------------------------------------------------#

tarchetypes::tar_quarto(
  quarto_file,
  path = "/media/Analyses/RedMeat-TopMed-Metabolites/RedMeat-TopMed-Metabolites_temp4ASN.qmd",
  quiet = FALSE
)

)

