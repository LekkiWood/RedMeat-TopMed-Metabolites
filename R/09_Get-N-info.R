 
#' 
#' 
#'
#----------------------------------------------------------#
#---------------Practice ---- -----------------------------#
#----------------------------------------------------------#
#library(targets)
#library(tidyverse)
#metabs_db = targets::tar_read(Metabs_long_clean)
#traits_db = targets::tar_read(Traits_long) = targets::tar_read(Traits_long)

#----------------------------------------------------------#
#---------------Function -----------------------------#
#----------------------------------------------------------#

get_N_function <- function(traits_db,
                           metabs_db)
{
  
  #----------------Build QC info------------------------#
  
  ####Those with MESA metabolites data at E1, E5, or E6
  
  ids_E1 <- traits_db |>
    tidyr::drop_na(age) |>
    dplyr::filter(exam==1) |>
    dplyr::select(idno, exam, sidno)
  
  ids_E5 <- traits_db |>
    tidyr::drop_na(age) |>
    dplyr::filter(exam==5) |>
    dplyr::select(idno, exam, sidno)
  
  ids_E6 <- traits_db |>
    tidyr::drop_na(age) |>
    dplyr::filter(exam==6) |>
    dplyr::select(idno, exam, sidno)
    
  
  ids_with_metabs <- metabs_db |>
    dplyr::filter(exam == 1 | exam ==5 | exam == 6) |>
    dplyr::select(sidno, exam, subject_id)
  
  ids_with_E1metabs <- metabs_db |>
    dplyr::filter(exam == 1) |>
    dplyr::select(sidno, exam, subject_id)
  
  ids_with_E5metabs <- metabs_db |>
    dplyr::filter(exam == 5) |>
    dplyr::select(sidno, exam, subject_id)
  
  ids_with_E6metabs <- metabs_db |>
    dplyr::filter(exam == 6) |>
    dplyr::select(sidno, exam, subject_id)
  
  
  ids_with_E1E5metabs <- metabs_db |>
    dplyr::group_by(sidno) |>
    dplyr::filter(exam == 1 | exam ==5) |>
    dplyr::filter(n_distinct(exam) >= 2) |>
    dplyr::ungroup() 
  
  
  ids_with_metabs_and_diet <- traits_db |>
    dplyr::filter(sidno %in%  ids_with_E1E5metabs$sidno) |>
    tidyr::drop_na(redmeat_cwc_z) |>
    dplyr::select(sidno, exam, sidno)
  
  ids_with_2diet <- traits_db |>
    dplyr::filter(sidno %in%  ids_with_E1E5metabs$sidno) |>
    tidyr::drop_na(redmeat_cwc_z) |>
    dplyr::group_by(sidno) |>
    dplyr::filter(exam == 1 | exam ==5) |>
    dplyr::filter(n_distinct(exam) >= 2) |>
    dplyr::ungroup() |>
    dplyr::pull(sidno)
  
  ids_with_metabs_and_2diet <- metabs_db |>
    dplyr::filter(sidno %in%  ids_with_2diet) |>
    dplyr::select(sidno, exam, subject_id)
  
  ids_with_metabs_and_E6metabs <- metabs_db |>
    dplyr::filter(exam==6) |>
    dplyr::select(sidno, exam, subject_id)
  
  ####Those with MESA diet data
  
 
  N_info <- list(ids_E1 =   ids_E1,
                 ids_E5 =   ids_E5,
                 ids_E6 =   ids_E6,
                 
                 ids_with_metabs = ids_with_metabs,
                 ids_with_E1metabs = ids_with_E1metabs,
                 ids_with_E5metabs = ids_with_E5metabs,
                 ids_with_E6metabs = ids_with_E6metabs,
                 
                 ids_with_E1E5metabs =  ids_with_E1E5metabs,
                 ids_with_metabs_and_diet = ids_with_metabs_and_diet,
                 ids_with_metabs_and_2diet = ids_with_metabs_and_2diet,
                 ids_with_metabs_and_E6metabs = ids_with_metabs_and_E6metabs
                     
    )
    
  
  
  
  
  #----------------Outputs -----------------------------#

  
}