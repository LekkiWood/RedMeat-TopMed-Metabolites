#' @param: bridge = bridge file 
#' 
#' 
#'
#----------------------------------------------------------#
#---------------Practice ---- -----------------------------#
#----------------------------------------------------------#

#E1_ffq_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1_FFQ_20160520.dta")
#E1_nutr_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1_Nutrients_20151208.dta")
#E5_ffq_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_FFQ_20140130.dta")
#E5_nutr_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_Nutrient_20140416.dta")

#Standard exam data
#E1_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta")
#E5_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_FinalLabel_20140613.dta")
#E6_raw <- foreign::read.dta("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe6_FinalLabel_20220513.dta")

#Cognition data
#cog_raw <- data.table::fread("/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESA_MA_E7_worse_cogdx_02212025.csv")

#Bridging file
#bridge <- data.table::fread("/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv") |>
#  dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
#  dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`)

#----------------------------------------------------------#
#---------------Function -----------------------------#
#----------------------------------------------------------#

build_traits <- function(
    path_bridge = path_bridge,
    path_e1_ffq_file = path_e1_ffq_file, 
    path_e1_nutr_file = path_e1_nutr_file,
    path_e5_ffq_file = path_e5_ffq_file, 
    path_e5_nutr_file = path_e5_nutr_file,
    path_exam1_file = path_exam1_file,
    path_exam5_file = path_exam5_file,
    path_exam6_file = path_exam6_file,
    path_cognition_file = path_cognition_file,
    min_kcal = 600, max_kcal = 6000,
    verbose = TRUE,
    return_diagnostics = TRUE
) 
{
  
  #----------------------------------------------------------#
  #---------------Read in files -----------------------------#
  #----------------------------------------------------------#
  
  #Diet-data
  E1_ffq_raw <- foreign::read.dta(path_e1_ffq_file)
  E1_nutr_raw <- foreign::read.dta(path_e1_nutr_file)
  E5_ffq_raw <- foreign::read.dta(path_e5_ffq_file)
  E5_nutr_raw <- foreign::read.dta(path_e5_nutr_file)
  
  #Standard exam data
  E1_raw <- foreign::read.dta(path_exam1_file)
  E5_raw <- foreign::read.dta(path_exam5_file)
  E6_raw <- foreign::read.dta(path_exam6_file)
  
  #Cognition data
  cog_raw <- read.csv(path_cognition_file)
  
  #Bridging file
  bridge <- data.table::fread(path_bridge) |>
    dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
    dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`)

  #----------------------------------------------------------#
  #---------------Covariates.   -----------------------------#
  #----------------------------------------------------------#
  
##E1
  
  E1 <- E1_raw |>
    dplyr::select(idno, age1c, gender1, race1c, site1c, egfr1c, dm031c, cig1c, income1, pavcm1c, ldl1, hdl1, trig1, glucos1c, inslnr1t) |>
    dplyr::rename(age = age1c, gender = gender1, race = race1c, site = site1c, egfr = egfr1c, DM = dm031c, smoking = cig1c, income = income1, PA = pavcm1c, ldl = ldl1, hdl = hdl1, TG = trig1, glucose = glucos1c, insulin = inslnr1t) |>
    dplyr::mutate(DM = dplyr::case_when(DM=="NORMAL" ~ 0,
                                        DM=="IFG" ~ 0,
                                        DM=="Untreated DIABETES" ~ 1,
                                        DM=="Treated DIABETES" ~ 1,
                                        TRUE ~ NA_real_)) |>
      dplyr::mutate(gender = factor(gender, labels=c("Female", "Male")),
                  race = factor(race, labels=c("White", "Chinese-American", "African-American", "Hispanic")),
                  site = factor(site, labels=c("WFU", "COL", "JHU", "UMN", "NWU", "UCLA")),
                  smoking = factor(smoking, labels=c("Never", "Former", "Current"), ordered=TRUE))
  
  E1 <- E1 |>
    dplyr::mutate(income = dplyr::case_when(income =="1: < $5,000" ~ 0,
                              income =="2: $5,000-$7,999" ~ 1,
                              income =="3: $8,000-$11,999" ~ 2,
                              income =="4: $12,000-$15,99" ~ 3,
                              income =="5: $16,000-$19,999" ~ 4,
                              income =="6: $20,000-$24,999" ~ 5,
                              income =="7: $25,000-$29,999" ~ 6,
                              income =="8: $30,000-$34,999" ~ 7,
                              income =="9: $35,000-$39,999" ~ 8,
                              income =="10: $40,000-$49,999" ~ 9,
                              income =="11: $50,000-$74,999" ~ 10,
                              income =="12: $75,000-$99,999" ~ 11,
                              income =="13: $100,000+" ~ 12,
                              TRUE ~ NA_real_,)) |>
    dplyr::mutate(exam = 1) 
  
##E5

E5 <- E5_raw |>
  dplyr::select(idno, age5c, gender1, race1c, site5c, egfr5c, dm035c, cig5c, income5, pavcm5c, ldl5, hdl5, trig5, glucose5, insulin5) |>
  dplyr::rename(age = age5c, gender = gender1, race = race1c, site = site5c, egfr = egfr5c, DM = dm035c, smoking = cig5c, income = income5, PA = pavcm5c, ldl = ldl5, hdl = hdl5, TG = trig5, glucose = glucose5, insulin = insulin5) |>
  dplyr::mutate(DM = dplyr::case_when(DM=="NORMAL" ~ 0,
                                      DM=="IMPAIRED FASTING GLUCOSE" ~ 0,
                                      DM=="UNTREATED DIABETES" ~ 1,
                                      DM=="TREATED DIABETES" ~ 1,
                                      TRUE ~ NA_real_)) |>
  dplyr::mutate(gender = factor(gender, labels=c("Female", "Male")),
                race = factor(race, labels=c("White", "Chinese-American", "African-American", "Hispanic")),
                site = factor(site, labels=c("WFU", "COL", "JHU", "UMN", "NWU", "UCLA")),
                smoking = factor(smoking, labels=c("Never", "Former", "Current"), ordered=TRUE))
              

E5 <- E5 |>
  dplyr::mutate(income = dplyr::case_when(income =="1: < $5000" ~ 0,
                                          income =="2: $5000 - $7999" ~ 1,
                                          income =="3: $8000 - $11999" ~ 2,
                                          income =="4: $12000 - $15999" ~ 3,
                                          income =="5: $16000 - $19999" ~ 4,
                                          income =="6: $20000 - $24999" ~ 5,
                                          income =="7: $25000 - $29999" ~ 6,
                                          income =="8: $30000 - $34999" ~ 7,
                                          income =="9: $35000 - $39999" ~ 8,
                                          income =="10: $40000 - $49999" ~ 9,
                                          income =="11: $50000 - $74999" ~ 10,
                                          income =="12: $75000 - $99999" ~ 11,
                                          income =="13: $100,000 - $124,999" ~ 12,
                                          income =="14: $125,000 - $149,999" ~ 12,
                                          income =="15: $150,000 or more9" ~ 12,
                                          TRUE ~ NA_real_,)) |>
  dplyr::mutate(exam = 5) 


##E6

E6 <- E6_raw |>
  dplyr::select(idno, age6c, gender1, race1c, site6c, egfr6c, dm036c, cig6c, income6, pavcm6c, ldl6, hdl6, trig6, glucose6, insulin6) |>
  dplyr::rename(age = age6c, gender = gender1, race = race1c, site = site6c, egfr = egfr6c, DM = dm036c, smoking = cig6c, income = income6, PA = pavcm6c, ldl = ldl6, hdl = hdl6, TG = trig6, glucose = glucose6, insulin = insulin6) |>
  dplyr::mutate(DM = dplyr::case_when(DM=="NORMAL" ~ 0,
                                      DM=="IFG" ~ 0,
                                      DM=="Untreated DIABETES" ~ 1,
                                      DM=="Treated DIABETES" ~ 1,
                                      TRUE ~ NA_real_)) |>
  dplyr::mutate(gender = factor(gender, labels=c("Female", "Male")),
                race = factor(race, labels=c("White", "Chinese-American", "African-American", "Hispanic")),
                site = factor(site, labels=c("WFU", "COL", "JHU", "UMN", "NWU", "UCLA")),
                smoking = factor(smoking, labels=c("Never", "Former", "Current"), ordered=TRUE))

E6 <- E6 |>
  dplyr::mutate(income = dplyr::case_when(income =="< $5000" ~ 0,
                                          income =="$5000 - $7999" ~ 1,
                                          income =="$8000 - $11999" ~ 2,
                                          income =="$12000 - $15999" ~ 3,
                                          income =="$16000 - $19999" ~ 4,
                                          income =="$20000 - $24999" ~ 5,
                                          income =="$25000 - $29999" ~ 6,
                                          income =="$30000 - $34999" ~ 7,
                                          income =="$35000 - $39999" ~ 8,
                                          income =="$40000 - $49999" ~ 9,
                                          income =="$50000 - $74999" ~ 10,
                                          income =="$75000 - $99999" ~ 11,
                                          income =="100,000 - $124,999" ~ 12,
                                          income =="$125,000 - $149,999" ~ 12,
                                          income =="$150,000 OR MORE" ~ 12,
                                          TRUE ~ NA_real_,)) |>
  dplyr::mutate(exam = 6) 



##Bind

Covs <- rbind(E1, E5)
Covs <- rbind(Covs, E6)

##Add in educ

Educ_info <- E1_raw |>
  dplyr::select(idno, educ1) |>
  dplyr::rename(education= educ1) |>
  dplyr::mutate(education = dplyr::case_when(education=="0: NO SCHOOLING" ~ 0,
                             education=="1: GRADES 1-8" ~ 1,
                             education=="2: GRADES 9-11" ~ 2,
                             education=="3: COMPLETED HIGH SCHOOL/GED" ~ 3,
                             education=="4: SOME COLLEGE BUT NO DEGREE" ~ 4,
                             education=="5: TECHNICAL SCHOOL CERTIFICATE" ~ 5,
                             education=="6: ASSOCIATE DEGREE" ~ 6,
                             education=="7: BACHELOR'S DEGREE" ~ 7,
                             education=="8: GRADUATE OR PROFESSIONAL SCHOOL" ~ 8,
                             TRUE ~ NA_real_)) 

Covs <- dplyr::full_join(Covs, Educ_info, dplyr::join_by(idno)) 

#----------------------------------------------------------#
#---------------Cog-status   -----------------------------#
#----------------------------------------------------------#

cog_status <- cog_raw |>
  dplyr::select(idno, cog_dxma6) |>
  dplyr::rename(cogstatus = cog_dxma6)|>
  dplyr::mutate(exam = 6,
               cogstatus = dplyr::case_when(cogstatus=="NI" ~ 0,
                               cogstatus=="MCI" ~ 1,
                               cogstatus=="PD" ~ 2,
                               TRUE ~ NA_real_)) 

cog_status <- cog_status |>
  dplyr::mutate(cogstatus_binary = dplyr::case_when(cogstatus== 0 ~ 0,
                                            cogstatus== 1 ~ 1,
                                            cogstatus==2 ~ 1,
                                            TRUE ~ NA_real_)) 

cog_status <- cog_status |>
  dplyr::mutate(cogstatus = factor(cogstatus, labels=c("Not impaired", "Mild impairment", "Probably dementia"), ordered=TRUE))  |>
  dplyr::filter(!is.na(cogstatus))



Covs <- dplyr::full_join(Covs, cog_status, dplyr::join_by(idno, exam)) 

#----------------------------------------------------------#
#--------------- Diet-vars   -----------------------------#
#----------------------------------------------------------#

Diet1A <- E1_ffq_raw |>
  dplyr::select(idno, fgredmeat1c, fghfprocmeat1c) |>
  dplyr::rename(redmeat = fgredmeat1c, 
                processedredmeat = fghfprocmeat1c) |>
  dplyr::mutate(exam=1)


Diet1B <- E1_nutr_raw |>
  dplyr::select(idno, enrgyn1c) |>
  dplyr::rename(energy = enrgyn1c) |>
  dplyr::mutate(exam=1)

Diet1 <- dplyr::full_join(Diet1A, Diet1B, dplyr::join_by(idno, exam)) 

# Exam 5
Diet5A <- E5_ffq_raw |>
  dplyr::select(idno, fgredmeat5c, fghfprocmeat5c) |>
  dplyr::rename(redmeat = fgredmeat5c, 
                processedredmeat = fghfprocmeat5c) |>
  dplyr::mutate(exam=5)


Diet5B <- E5_nutr_raw |>
  dplyr::select(idno, enrgyn5c) |>
  dplyr::rename(energy = enrgyn5c) |>
  dplyr::mutate(exam=5)

Diet5 <- dplyr::full_join(Diet5A, Diet5B, dplyr::join_by(idno, exam)) 

Diet <- rbind(Diet1, Diet5) 

Diet <- Diet |>
  mutate(redmeat_z = ifelse(exam == 1 | exam==5, (redmeat - mean(redmeat, na.rm=TRUE)) / sd(redmeat, na.rm=TRUE), NA_real_)) |> 
  group_by(idno) |>
  mutate(
    redmeat_pm_z  = ifelse(exam == 1 | exam==5, mean(redmeat_z, na.rm=TRUE), NA_real_), 
    redmeat_cwc_z = ifelse(exam == 1 | exam==5, redmeat_z - redmeat_pm_z, NA_real_) 
  ) |>
  ungroup() |>
  mutate(
    redmeat_gm_z     = ifelse(exam == 1 | exam==5, mean(redmeat_z, na.rm=TRUE), NA_real_), # should be ~0
    redmeat_pm_cgm_z = ifelse(exam == 1 | exam==5, redmeat_pm_z - redmeat_gm_z, NA_real_),
    redmeat_cgm_z    = ifelse(exam == 1 | exam==5, redmeat_z - redmeat_gm_z, NA_real_)
  )


Covs <- dplyr::full_join(Covs, Diet, dplyr::join_by(idno, exam))
#----------------Bridge------------------------#

Covs <- dplyr::left_join(Covs, bridge, dplyr::join_by(idno))

#----------------Build QC info------------------------#

QC_info <- list(
  filenames = list(E1_ffq_raw = path_e1_ffq_file,
                   E1_nutr_raw = path_e1_nutr_file, 
                   E5_ffq_raw = path_e5_ffq_file, 
                   E5_nutr_raw = path_e5_nutr_file,
                   E1_raw = path_exam1_file,
                   E5_raw = path_exam5_file,
                   E6_raw = path_exam6_file,
                   cog_raw = path_cognition_file,
                   bridge = path_bridge
                   
  )
 
)



#----------------Outputs -----------------------------#

list(
  
  QC_info_out = QC_info,
  Covs = Covs
  
)

}
