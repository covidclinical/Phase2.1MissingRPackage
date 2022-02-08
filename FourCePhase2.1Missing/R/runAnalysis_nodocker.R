#' Run analysis
#'
#' @return NULL. Result files are written out to the `getProjectOutputDirectory()` directory.
#'
#' @keywords 4CE
#' @export(getPrivateSummaryRepositoryUrl)
#' @export(getProjectOutputDirectory)
#' @export(getPublicSummaryRepositoryUrl)
#' @export(runAnalysis)
#' @export(smallCountObfuscation)
#' @export(submitAnalysis)
#' @export(validateAnalysis)
#' @import(dplyr)
#' @import(tidyr)
#' @import(tidyverse)
#' @import(tibble)
#' @import(ggplot2)
#' @import(furrr)
#' @import(forcats)
#' @import(tidytext)
#' @import(readr)
#' @import(Rcpp)
#' @import(DT)
#' @import(naniar)
#' @import(cowplot)
#' @import(ggridges)
#' @import(RColorBrewer)
#' @import(mSTEM)
#' @import(Rcpp)
#' @import(stringr)
#' @import(matrixStats)
#' @import(lubridate)
#' @import(beepr)
#' @import(RcppRoll)
#' @import(mSTEM)
#' @import(MCMCpack)
#' @import(kableExtra)

 
runAnalysis_nodocker <- function(data_dir = "~/4ceData/Input",dateFormat="%d-%b-%y",time="all",siteid = "upenn") {
  for (r_file in list.files('R', full.names = TRUE, pattern = '.R$')) source(r_file)
  #install.packages('stm')
  #library(stm)
 


 
  my_dir = data_dir
 
 if (siteid != "VA1" && siteid != "VA2" && siteid != "VA3" && siteid != "VA4" && siteid != "VA5"){
  
  demo_raw <-
  readr::read_csv(
    file.path(my_dir, "LocalPatientSummary.csv"),
    col_types = list(patient_num = readr::col_character()),
    na = c("1900-01-01", "1/1/1900")
  ) %>%
  mutate(
    across(ends_with("_date") & where(is.character), lubridate::mdy),     
    last_discharge_date = if_else(
      !is.na(death_date) & death_date < last_discharge_date,
      death_date,
      last_discharge_date
    )
  )
  }else{
  demo_raw = demographics
  }
  
  obs_raw <-
    readr::read_csv(
      file.path(data_dir, "LocalPatientObservations.csv"),
      col_types = list(patient_num = readr::col_character())
    )
  
  #Read in thrombotic icd codes. Add column "truncated code" for first three characters of each code. Extract unique truncated codes. 
  thrombo_codes <- c("I74", "I75", "I76", "I21", "I22", "I23", "I26", "I27", "Z86", "I63", "I67", "I81", "I82","444","445",
"410","415","V12","434","437","452","453","D65","P60","286","776")
  #Extract information for patients with codes for thrombotic events

patient_obs <- readr::read_csv(file.path(data_dir, "LocalPatientObservations.csv"),
                                 col_types = list(
                                   patient_num = readr::col_character(),
                                   value = readr::col_double()
                                 ),
                                 na = "-999")
 
patient_obs$patient_num <- as.character(patient_obs$patient_num)
 
  te_patients <- patient_obs %>%
    filter(
      concept_type == "DIAG-ICD10",
      concept_code %in% thrombo_codes,
      days_since_admission >= 0
    )
  
  clin_raw <-
    readr::read_csv(
      file.path(data_dir, "LocalPatientClinicalCourse.csv"),
      col_types = list(patient_num = readr::col_character())
    )

lab_bounds1 <- c("1742-6","1751-7","1920-8","1975-2","1988-5","2160-0","2276-4","2532-0","3255-7","33959-8","48065-7","49563-0","6690-2","731-0","751-8","48066-5","5902-2","6598-7")
lab_bounds <- clin_raw[1:length(lab_bounds1),][1]

lab_bounds$LOINC<- lab_bounds1
lab_bounds$short_name = c("ALT","Albumin","AST","Bilirubin","CRP","Creatinine","Ferritin","LDH","Fibrinogen","Procalcitonin","FEU","Troponin_normal","Leukocytes","Lymphocyte","Neutrophil","DDU","PT","Troponin_high")



 
 
  if (time == "phase_1") {
    demo_raw = demo_raw[demo_raw$admission_date <= "2020-07-30",]
    demo_raw = demo_raw[demo_raw$admission_date >= "2020-01-01",]
    demo_raw = demo_raw[!is.na(demo_raw$admission_date),]
    patient_obs = patient_obs[patient_obs$patient_num %in% unique(demo_raw$patient_num),]
    obs_raw = obs_raw[obs_raw$patient_num %in% unique(demo_raw$patient_num),]
    clin_raw = clin_raw[clin_raw$patient_num %in% unique(demo_raw$patient_num),]
  }
  
  if (time == "phase_2") {
    demo_raw = demo_raw[demo_raw$admission_date > "2020-07-30",]
    demo_raw = demo_raw[demo_raw$admission_date <= "2021-09-30",]
    demo_raw = demo_raw[!is.na(demo_raw$admission_date),]
    patient_obs = patient_obs[patient_obs$patient_num %in% unique(demo_raw$patient_num),]
    obs_raw = obs_raw[obs_raw$patient_num %in% unique(demo_raw$patient_num),]
    clin_raw = clin_raw[clin_raw$patient_num %in% unique(demo_raw$patient_num),]
  }
  
  if (time == "all") {
    demo_raw = demo_raw[demo_raw$admission_date <= "2021-09-30",]
    demo_raw = demo_raw[demo_raw$admission_date >= "2020-01-01",]
    demo_raw = demo_raw[!is.na(demo_raw$admission_date),]
    patient_obs = patient_obs[patient_obs$patient_num %in% unique(demo_raw$patient_num),]
    obs_raw = obs_raw[obs_raw$patient_num %in% unique(demo_raw$patient_num),]
    clin_raw = clin_raw[clin_raw$patient_num %in% unique(demo_raw$patient_num),]
  }
  
  valid_days = c()
  ct = 0
  ct2 = 0
  for (i in 1:length(unique(clin_raw$patient_num))){
    test_patient = clin_raw[clin_raw$patient_num==unique(clin_raw$patient_num)[i],]
    j = 1
    while(test_patient$in_hospital[j] == 1 && j != nrow(test_patient)) {
      ct = ct + 1
      j=j+1
    }
    if (j == nrow(test_patient)) {
      valid_days[i] = test_patient$days_since_admission[j]
    }else{
     if (test_patient$in_hospital[1] == 0){
      valid_days[i] = 0
     }else{
    valid_days[i] = test_patient$days_since_admission[(j-1)]
    }
    }
  }
  

  patient_obs_wide <- patient_obs %>%
    left_join(lab_bounds, by = c("concept_code" = "LOINC")) %>%
    filter(days_since_admission >= 0, !is.na(short_name)) %>%
    dplyr::select(-concept_code) %>%
    pivot_wider(
      id_cols = c(patient_num, days_since_admission),
      names_from = short_name,
      values_from = value,
      values_fn = mean
    ) %>%
    left_join(dplyr::select(demo_raw, patient_num, severe), by = "patient_num") %>%
    mutate(severity = factor(severe) %>%
             fct_recode("severe" = "1", "nonsevere" = "0")) %>%
    dplyr::select(-severe)
  ct = 0
  if ("Troponin_high" %in% colnames(patient_obs_wide)){
    if ("Troponin_normal" %in% colnames(patient_obs_wide)) {
      ct = ct + 1
    }else{
      colnames(patient_obs_wide)[which(colnames(patient_obs_wide)=="Troponin_high")] <- "Troponin"
    }
    
  }
  if ("Troponin_normal" %in% colnames(patient_obs_wide)) {
    if ("Troponin_high" %in% colnames(patient_obs_wide)) {
      patient_obs_wide$Troponin = patient_obs_wide$Troponin_high +   patient_obs_wide$Troponin_normal
      for (i in 1:length(patient_obs_wide$Troponin_high)){
        if (!is.na(patient_obs_wide$Troponin_high[i]) | !is.na(patient_obs_wide$Troponin_normal[i])) {
          patient_obs_wide$Troponin[i] = 1
        }
        else {
          patient_obs_wide$Troponin[i] = patient_obs_wide$Troponin_high[i]
        }
      }
    }else{
      colnames(patient_obs_wide)[which(colnames(patient_obs_wide)=="Troponin_normal")] <- "Troponin"
    }
  }
  drops = c("Troponin_high","Troponin_normal")
  patient_obs_wide= patient_obs_wide[,!(names(patient_obs_wide) %in% drops)]
  ct = 0
  if ("DDU" %in% colnames(patient_obs_wide)){
    if ("FEU" %in% colnames(patient_obs_wide)) {
      ct = ct + 1
    }else{
      colnames(patient_obs_wide)[which(colnames(patient_obs_wide)=="DDU")] <- "D_Dimer"
    }
    
  }
  if ("FEU" %in% colnames(patient_obs_wide)) {
    if ("DDU" %in% colnames(patient_obs_wide)) {
      patient_obs_wide$D_Dimer = patient_obs_wide$FEU +   patient_obs_wide$DDU
      for (i in 1:length(patient_obs_wide$DDU)){
        if (!is.na(patient_obs_wide$DDU[i]) | !is.na(patient_obs_wide$FEU[i])) {
          patient_obs_wide$D_Dimer[i] = 1
        }
        else {
          patient_obs_wide$D_Dimer[i] = patient_obs_wide$DDU[i]
        }
      }
    }else{
      colnames(patient_obs_wide)[which(colnames(patient_obs_wide)=="FEU")] <- "D_Dimer"
    }
  }
  
  test_obs_wide = patient_obs_wide[FALSE,]
  for (i in 1:length(valid_days)){
    if ((unique(clin_raw$patient_num))[i] %in% unique(patient_obs_wide$patient_num)) {
      test_patient = patient_obs_wide[patient_obs_wide$patient_num==unique(clin_raw$patient_num)[i],]
      adj_patient = test_patient[test_patient$days_since_admission <= valid_days[i],]
      test_obs_wide = rbind(test_obs_wide, adj_patient)
    }
  }
  
  patient_obs_wide = test_obs_wide
  
  test_clin_raw = clin_raw[FALSE,]
  for (i in 1:length(valid_days)){
    test_patient = clin_raw[clin_raw$patient_num==unique(clin_raw$patient_num)[i],]
    adj_patient = test_patient[test_patient$days_since_admission <= valid_days[i],]
    test_clin_raw = rbind(test_clin_raw, adj_patient)
    }

  
  clin_raw = test_clin_raw
  
  test_obs_raw = obs_raw[FALSE,]
  for (i in 1:length(valid_days)){
    if ((unique(clin_raw$patient_num))[i] %in% unique(obs_raw$patient_num)) {
      test_patient = obs_raw[obs_raw$patient_num==unique(clin_raw$patient_num)[i],]
      adj_patient = test_patient[test_patient$days_since_admission <= valid_days[i],]
      test_obs_raw = rbind(test_obs_raw, adj_patient)
    }
  }
  
  obs_raw = test_obs_raw
  

  
  
  
  

  
  
  drops = c("FEU","DDU")
  patient_obs_wide= patient_obs_wide[,!(names(patient_obs_wide) %in% drops)]
  #Specify lab ordering so that analyses are easy to consolidate across sites
  #If your site has additional labs, please add them at the end of allLabs:
  allLabs <- c("Fibrinogen","Procalcitonin","D_Dimer","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
  labsPresent <- allLabs[ allLabs %in% colnames(patient_obs_wide)]
  #Extract unique lab names
  lab_names <- intersect(lab_bounds$short_name, names(patient_obs_wide))
  patient_obs_wide <- patient_obs_wide[,c("patient_num", "days_since_admission", labsPresent, "severity")]
  Severe_pats <- demo_raw[demo_raw$severe == 1,]
  wide <- patient_obs_wide %>% mutate(severe = patient_num %in% unique(Severe_pats$patient_num)) 
  Severe_patients <- wide[wide$severe==1,]
  Nonsevere_patients <- wide[wide$severe==0,]
  
  
  
  ####### Quantify Missingness
  
  
  avgs = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    for (j in 1:ncol(avgs)) {
      dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
      if (nrow(dat) != 0){
        num_miss = sum(is.na(dat[j+2]))
        avgs[ct,j] = num_miss
        props[ct,j] = num_miss/(nrow(dat[j+2]))
        if ( j == ncol(avgs)){
          ct = ct + 1
        }
      }
    }
  }
  dat_prop = colMeans(props)
  dat_sd_prop = apply(props,2,sd)
  dat_num = colMeans(avgs)
  dat_sd_num = apply(avgs,2,sd)
  site = rep(siteid,length(dat_prop))
  df_prop = data.frame(Labs = colnames(avgs), Avg_Prop_Miss = dat_prop,sd = dat_sd_prop,Site=site)
  df_num = data.frame(Labs = colnames(avgs), Avg_Num_Miss = dat_num,sd = dat_sd_num,Site=site)
  df_prop$Labs <- factor(df_prop$Labs,levels = labsPresent)
  df_num$Labs <- factor(df_prop$Labs,levels = labsPresent )
  prop_plot = ggplot(data = df_prop,
                     aes(x = Labs,
                         y = Avg_Prop_Miss,
                         fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = pmin(Avg_Prop_Miss + sd,1)), width = .2,position=position_dodge(.9)) 
  prop_plot 
  num_plot = ggplot(data = df_num,
                    aes(x = Labs,
                        y = Avg_Num_Miss,
                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(x = "Labs",y = "Average Number Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Num_Miss-sd,0),ymax = Avg_Num_Miss + sd), width = .2,position=position_dodge(.9)) 
  num_plot
  require(gridExtra)
  grid.arrange(prop_plot,num_plot,ncol=1)
  props_quant= colQuantiles(as.matrix(props))
  samplesize = nrow(props)
  
  num_missing = sapply(patient_obs_wide, function(x) sum(is.na(x)))[3:18]
  prop_missing = sapply(patient_obs_wide, function(x) sum(is.na(x))/length(x))[3:18]
  
  
  #first 30 days: avg prop missing and avg num missing
  
  avgs = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  ct = 1
  for (i in 1:30) {
    for (j in 1:ncol(avgs)) {
      dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
      if (nrow(dat) != 0){
        num_miss = sum(is.na(dat[j+2]))
        avgs[ct,j] = num_miss
        props[ct,j] = num_miss/(nrow(dat[j+2]))
        if ( j == ncol(avgs)){
          ct = ct + 1
        }
      }
    }
  }
  dat_prop = colMeans(props)
  dat_sd_prop = apply(props,2,sd)
  dat_num = colMeans(avgs)
  dat_sd_num = apply(avgs,2,sd)
  site = rep(siteid,length(dat_prop))
  df_prop2 = data.frame(Labs = colnames(avgs), Avg_Prop_Miss = dat_prop,sd = dat_sd_prop,Site=site)
  df_num2 = data.frame(Labs = colnames(avgs), Avg_Num_Miss = dat_num,sd = dat_sd_num,Site=site)
  df_prop2$Labs <- factor(df_prop2$Labs,levels = labsPresent)
  df_num2$Labs <- factor(df_prop2$Labs,levels = labsPresent )
  prop_plot2 = ggplot(data = df_prop2,
                      aes(x = Labs,
                          y = Avg_Prop_Miss,
                          fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = pmin(Avg_Prop_Miss + sd,1)), width = .2,position=position_dodge(.9)) 
  prop_plot2 
  num_plot2 = ggplot(data = df_num2,
                     aes(x = Labs,
                         y = Avg_Num_Miss,
                         fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(x = "Labs",y = "Average Number Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Num_Miss-sd,0),ymax = Avg_Num_Miss + sd), width = .2,position=position_dodge(.9)) 
  num_plot2
  require(gridExtra)
  grid.arrange(prop_plot2,num_plot2,ncol=1)
  props_quant2= colQuantiles(as.matrix(props))
  
  
  #First 20 days: proportion missing, num missing
  
  patient_obs_first30 = patient_obs_wide[patient_obs_wide$days_since_admission <= 30,]
  num_missing2 = sapply(patient_obs_first30, function(x) sum(is.na(x)))[3:18]
  prop_missing2 = sapply(patient_obs_first30, function(x) sum(is.na(x))/length(x))[3:18]
  
  
  ## Stratifications by sex and severity
  
  
  props_M = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_F = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_M = (demo_raw[which(demo_raw$sex=="male"),])$patient_num
  indices_F = (demo_raw[which(demo_raw$sex=="female"),])$patient_num
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_M){
      
      for (j in 1:ncol(props_M)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_M[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_M)){
            ct = ct + 1
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_F){
      for (j in 1:ncol(props_F)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_F[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_F)){
            ct = ct + 1
          }
        }
      }
    }
  }

  
  dat_prop_M = colMeans(props_M)
  dat_prop_F = colMeans(props_F)
  dat_sd_M= apply(props_M,2,sd)
  dat_sd_F= apply(props_F,2,sd)
  diff = dat_prop_M-dat_prop_F
  sex_M = rep("Male",length(dat_prop_M))
  sex_F = rep("Female",length(dat_prop_F))
  df_prop_M = data.frame(Labs = colnames(props_M), Avg_Prop_Miss = dat_prop_M,sd = dat_sd_M,Sex=sex_M)
  df_prop_M$Labs <- factor(df_prop_M$Labs,levels =labsPresent)
  df_prop_F = data.frame(Labs = colnames(props_F), Avg_Prop_Miss = dat_prop_F,sd = dat_sd_F,Sex=sex_F)
  df_prop_F$Labs <- factor(df_prop_F$Labs,levels = labsPresent)
  df_prop_diff_sex = data.frame(Labs = colnames(props_F), Avg_Prop_Miss = diff)
  df_prop_diff_sex$Labs <- factor(df_prop_diff_sex$Labs,levels = labsPresent)
  df_prop_Sex = rbind(df_prop_M,df_prop_F)
  site = rep(siteid,nrow(df_prop_Sex))
  df_prop_Sex$Site = site
  prop_plot_Sex = ggplot(data = df_prop_Sex,
                         aes(x = Labs,
                             y = Avg_Prop_Miss,
                             fill = Sex)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Sex", x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Sex")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Sex
  Site = rep(siteid,nrow(df_prop_diff_sex))
  df_prop_diff_sex$site = Site
  Site = rep(siteid,nrow(df_prop_diff_sex))
  df_prop_diff_sex$site = Site
  prop_plot_diff_sex = ggplot(data = df_prop_diff_sex,
                              aes(x = reorder(Labs,Avg_Prop_Miss),
                                  y = Avg_Prop_Miss,
                                  fill = site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Sex",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_sex
  #Stratify by severity!
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_S = (demo_raw[which(demo_raw$severe==1),])$patient_num
  indices_NS = (demo_raw[which(demo_raw$severe==0),])$patient_num
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      for (j in 1:ncol(props_S)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_S[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_S)){
            ct = ct + 1
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      for (j in 1:ncol(props_NS)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_NS)){
            ct = ct + 1
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S$Labs <- factor(df_prop_S$Labs,levels =labsPresent)
  df_prop_NS = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS$Labs <- factor(df_prop_NS$Labs,levels = labsPresent)
  df_prop_diff_severity = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity$Labs <- factor(df_prop_diff_severity$Labs,levels = labsPresent)
  df_prop_Severity = rbind(df_prop_S,df_prop_NS)
  site = rep(siteid,nrow(df_prop_Severity))
  df_prop_Severity$Site = site
  prop_plot_Severity = ggplot(data = df_prop_Severity,
                              aes(x = Labs,
                                  y = Avg_Prop_Miss,
                                  fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity
  Site = rep(siteid,nrow(df_prop_diff_severity))
  df_prop_diff_severity$site = Site
  prop_plot_diff_severity = ggplot(data = df_prop_diff_severity,
                                   aes(x = reorder(Labs,Avg_Prop_Miss),
                                       y = Avg_Prop_Miss,
                                       fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity
  props_M_quant = colQuantiles(as.matrix(props_M))
  props_F_quant = colQuantiles(as.matrix(props_F))
  props_S_quant= colQuantiles(as.matrix(props_S))
  props_NS_quant = colQuantiles(as.matrix(props_NS))
  ssM = nrow(props_M)
  ssF = nrow(props_F)
  ssS = nrow(props_S)
  ssNS = nrow(props_NS)
  
  
  
  
  
  
  
  
  
  
  
  #Stratify by sex-- first 30 days!
  
  props_M = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_F = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_M = (demo_raw[which(demo_raw$sex=="male"),])$patient_num
  indices_F = (demo_raw[which(demo_raw$sex=="female"),])$patient_num
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_M){
      
      for (j in 1:ncol(props_M)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_M[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_M)){
            ct = ct + 1
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_F){
      for (j in 1:ncol(props_F)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_F[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_F)){
            ct = ct + 1
          }
        }
      }
    }
  }
  dat_prop_M = colMeans(props_M)
  dat_prop_F = colMeans(props_F)
  dat_sd_M= apply(props_M,2,sd)
  dat_sd_F= apply(props_F,2,sd)
  diff = dat_prop_M-dat_prop_F
  sex_M = rep("Male",length(dat_prop_M))
  sex_F = rep("Female",length(dat_prop_F))
  df_prop_M2 = data.frame(Labs = colnames(props_M), Avg_Prop_Miss = dat_prop_M,sd = dat_sd_M,Sex=sex_M)
  df_prop_M2$Labs <- factor(df_prop_M2$Labs,levels =labsPresent)
  df_prop_F2 = data.frame(Labs = colnames(props_F), Avg_Prop_Miss = dat_prop_F,sd = dat_sd_F,Sex=sex_F)
  df_prop_F2$Labs <- factor(df_prop_F2$Labs,levels = labsPresent)
  df_prop_diff_sex2 = data.frame(Labs = colnames(props_F), Avg_Prop_Miss = diff)
  df_prop_diff_sex2$Labs <- factor(df_prop_diff_sex2$Labs,levels = labsPresent)
  df_prop_Sex2 = rbind(df_prop_M2,df_prop_F2)
  site = rep(siteid,nrow(df_prop_Sex2))
  df_prop_Sex2$Site = site
  prop_plot_Sex2= ggplot(data = df_prop_Sex2,
                         aes(x = Labs,
                             y = Avg_Prop_Miss,
                             fill = Sex)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Sex", x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Sex")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Sex2
  Site = rep(siteid,nrow(df_prop_diff_sex2))
  df_prop_diff_sex2$site = Site
  Site = rep(siteid,nrow(df_prop_diff_sex2))
  df_prop_diff_sex2$site = Site
  prop_plot_diff_sex2 = ggplot(data = df_prop_diff_sex2,
                               aes(x = reorder(Labs,Avg_Prop_Miss),
                                   y = Avg_Prop_Miss,
                                   fill = site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Sex",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_sex2
  
  
  #Stratify by severity!
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_S = (demo_raw[which(demo_raw$severe==1),])$patient_num
  indices_NS = (demo_raw[which(demo_raw$severe==0),])$patient_num
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      
      for (j in 1:ncol(props_S)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_S[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_S)){
            ct = ct + 1
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      for (j in 1:ncol(props_NS)) {
        dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
        if (nrow(dat) != 0){
          num_miss = sum(is.na(dat[j+2]))
          props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
          if ( j == ncol(props_NS)){
            ct = ct + 1
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S2 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S2$Labs <- factor(df_prop_S2$Labs,levels =labsPresent)
  df_prop_NS2 = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS2$Labs <- factor(df_prop_NS2$Labs,levels = labsPresent)
  df_prop_diff_severity2 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity2$Labs <- factor(df_prop_diff_severity2$Labs,levels = labsPresent)
  df_prop_Severity2 = rbind(df_prop_S2,df_prop_NS2)
  site = rep(siteid,nrow(df_prop_Severity2))
  df_prop_Severity2$Site = site
  prop_plot_Severity2 = ggplot(data = df_prop_Severity2,
                               aes(x = Labs,
                                   y = Avg_Prop_Miss,
                                   fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity2
  Site = rep(siteid,nrow(df_prop_diff_severity2))
  df_prop_diff_severity2$site = Site
  prop_plot_diff_severity2 = ggplot(data = df_prop_diff_severity2,
                                    aes(x = reorder(Labs,Avg_Prop_Miss),
                                        y = Avg_Prop_Miss,
                                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity2
  props_M_quant2 = colQuantiles(as.matrix(props_M))
  props_F_quant2 = colQuantiles(as.matrix(props_F))
  props_S_quant2= colQuantiles(as.matrix(props_S))
  props_NS_quant2 = colQuantiles(as.matrix(props_NS))
  
  #Stratify by patients over 21 (severity)!
  
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_age = (demo_raw[which(demo_raw$age_group != "00to02" & demo_raw$age_group != "03to05" & demo_raw$age_group != "06to11" & demo_raw$age_group != "12to17"),])$patient_num
  
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        
        for (j in 1:ncol(props_S)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_S[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_S)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        for (j in 1:ncol(props_NS)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_NS)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S3 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S3$Labs <- factor(df_prop_S3$Labs,levels =labsPresent)
  df_prop_NS3 = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS3$Labs <- factor(df_prop_NS3$Labs,levels = labsPresent)
  df_prop_diff_severity3 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity3$Labs <- factor(df_prop_diff_severity3$Labs,levels = labsPresent)
  df_prop_Severity3 = rbind(df_prop_S3,df_prop_NS3)
  site = rep(siteid,nrow(df_prop_Severity3))
  df_prop_Severity3$Site = site
  prop_plot_Severity3 = ggplot(data = df_prop_Severity3,
                               aes(x = Labs,
                                   y = Avg_Prop_Miss,
                                   fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity3
  Site = rep(siteid,nrow(df_prop_diff_severity3))
  df_prop_diff_severity3$site = Site
  prop_plot_diff_severity3 = ggplot(data = df_prop_diff_severity3,
                                    aes(x = reorder(Labs,Avg_Prop_Miss),
                                        y = Avg_Prop_Miss,
                                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity3
  props_S_quant3= colQuantiles(as.matrix(props_S))
  props_NS_quant3 = colQuantiles(as.matrix(props_NS))
  
  
  
  
  
  
  #Stratify by patients over 18 (severity) in first 30 days!
  
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))

  indices_age = (demo_raw[which(demo_raw$age_group != "00to02" & demo_raw$age_group != "03to05" & demo_raw$age_group != "06to11" & demo_raw$age_group != "12to17"),])$patient_num
  
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        
        for (j in 1:ncol(props_S)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_S[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_S)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        for (j in 1:ncol(props_NS)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_NS)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S4 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S4$Labs <- factor(df_prop_S4$Labs,levels =labsPresent)
  df_prop_NS4 = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS4$Labs <- factor(df_prop_NS4$Labs,levels = labsPresent)
  df_prop_diff_severity4 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity4$Labs <- factor(df_prop_diff_severity4$Labs,levels = labsPresent)
  df_prop_Severity4 = rbind(df_prop_S4,df_prop_NS4)
  site = rep(siteid,nrow(df_prop_Severity4))
  df_prop_Severity4$Site = site
  prop_plot_Severity4 = ggplot(data = df_prop_Severity4,
                               aes(x = Labs,
                                   y = Avg_Prop_Miss,
                                   fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity4
  Site = rep(siteid,nrow(df_prop_diff_severity4))
  df_prop_diff_severity4$site = Site
  prop_plot_diff_severity4 = ggplot(data = df_prop_diff_severity4,
                                    aes(x = reorder(Labs,Avg_Prop_Miss),
                                        y = Avg_Prop_Miss,
                                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity4
  props_S_quant4= colQuantiles(as.matrix(props_S))
  props_NS_quant4 = colQuantiles(as.matrix(props_NS))
  
  
  
  #Stratify by patients under 18 (severity)!
  
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))

  indices_age = (demo_raw[which(demo_raw$age_group == "00to02" | demo_raw$age_group == "03to05" | demo_raw$age_group == "06to11" | demo_raw$age_group == "12to17"),])$patient_num
  
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        
        for (j in 1:ncol(props_S)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_S[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_S)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:length(unique(patient_obs_wide$patient_num))) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        for (j in 1:ncol(props_NS)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_NS)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S3 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S3$Labs <- factor(df_prop_S3$Labs,levels =labsPresent)
  df_prop_NS3 = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS3$Labs <- factor(df_prop_NS3$Labs,levels = labsPresent)
  df_prop_diff_severity3 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity3$Labs <- factor(df_prop_diff_severity3$Labs,levels = labsPresent)
  df_prop_Severity3 = rbind(df_prop_S3,df_prop_NS3)
  site = rep(siteid,nrow(df_prop_Severity3))
  df_prop_Severity3$Site = site
  prop_plot_Severity3 = ggplot(data = df_prop_Severity3,
                               aes(x = Labs,
                                   y = Avg_Prop_Miss,
                                   fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity3
  Site = rep(siteid,nrow(df_prop_diff_severity3))
  df_prop_diff_severity3$site = Site
  prop_plot_diff_severity3 = ggplot(data = df_prop_diff_severity3,
                                    aes(x = reorder(Labs,Avg_Prop_Miss),
                                        y = Avg_Prop_Miss,
                                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity3
  props_S_quant5= colQuantiles(as.matrix(props_S))
  props_NS_quant5 = colQuantiles(as.matrix(props_NS))
  
  
  
  
  
  
  #Stratify by patients over 18 (severity) in first 30 days!
  
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_age = (demo_raw[which(demo_raw$age_group != "00to02" & demo_raw$age_group != "03to05" & demo_raw$age_group != "06to11" & demo_raw$age_group != "12to17"),])$patient_num
  
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        
        for (j in 1:ncol(props_S)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_S[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_S)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        for (j in 1:ncol(props_NS)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_NS)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S4 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S4$Labs <- factor(df_prop_S4$Labs,levels =labsPresent)
  df_prop_NS4 = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS4$Labs <- factor(df_prop_NS4$Labs,levels = labsPresent)
  df_prop_diff_severity4 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity4$Labs <- factor(df_prop_diff_severity4$Labs,levels = labsPresent)
  df_prop_Severity4 = rbind(df_prop_S4,df_prop_NS4)
  site = rep(siteid,nrow(df_prop_Severity4))
  df_prop_Severity4$Site = site
  prop_plot_Severity4 = ggplot(data = df_prop_Severity4,
                               aes(x = Labs,
                                   y = Avg_Prop_Miss,
                                   fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity4
  Site = rep(siteid,nrow(df_prop_diff_severity4))
  df_prop_diff_severity4$site = Site
  prop_plot_diff_severity4 = ggplot(data = df_prop_diff_severity4,
                                    aes(x = reorder(Labs,Avg_Prop_Miss),
                                        y = Avg_Prop_Miss,
                                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity4
  props_S_quant6= colQuantiles(as.matrix(props_S))
  props_NS_quant6 = colQuantiles(as.matrix(props_NS))
  
  
  #Stratify by patients under 18 (severity) in first 30 days!
  
  props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
  indices_age = (demo_raw[which(demo_raw$age_group == "00to02" | demo_raw$age_group == "03to05" | demo_raw$age_group == "06to11" | demo_raw$age_group == "12to17"),])$patient_num
  
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_S){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        
        for (j in 1:ncol(props_S)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_S[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_S)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  ct = 1
  for (i in 1:30) {
    if (unique(patient_obs_wide$patient_num)[i] %in% indices_NS){
      if (unique(patient_obs_wide$patient_num)[i] %in% indices_age) {
        for (j in 1:ncol(props_NS)) {
          dat <- filter(patient_obs_wide,patient_num==demo_raw$patient_num[i])
          if (nrow(dat) != 0){
            num_miss = sum(is.na(dat[j+2]))
            props_NS[ct,j] = num_miss/(nrow(dat[j+2]))
            if ( j == ncol(props_NS)){
              ct = ct + 1
            }
          }
        }
      }
    }
  }
  dat_prop_S = colMeans(props_S)
  dat_prop_NS = colMeans(props_NS)
  dat_sd_prop_S = apply(props_S,2,sd)
  dat_sd_prop_NS = apply(props_NS,2,sd)
  diff = dat_prop_S-dat_prop_NS
  severe = rep("Severe",length(dat_prop_S))
  non_severe = rep("Non-Severe",length(dat_prop_NS))
  df_prop_S4 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = dat_prop_S,sd = dat_sd_prop_S,severity = severe)
  df_prop_S4$Labs <- factor(df_prop_S4$Labs,levels =labsPresent)
  df_prop_NS4 = data.frame(Labs = colnames(props_NS), Avg_Prop_Miss = dat_prop_NS,sd = dat_sd_prop_NS,severity = non_severe)
  df_prop_NS4$Labs <- factor(df_prop_NS4$Labs,levels = labsPresent)
  df_prop_diff_severity4 = data.frame(Labs = colnames(props_S), Avg_Prop_Miss = diff)
  df_prop_diff_severity4$Labs <- factor(df_prop_diff_severity4$Labs,levels = labsPresent)
  df_prop_Severity4 = rbind(df_prop_S4,df_prop_NS4)
  site = rep(siteid,nrow(df_prop_Severity4))
  df_prop_Severity4$Site = site
  prop_plot_Severity4 = ggplot(data = df_prop_Severity4,
                               aes(x = Labs,
                                   y = Avg_Prop_Miss,
                                   fill = severity)) +
    geom_bar(stat = "identity",width=0.8,position=position_dodge()) + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean +/- Sd Proportion Missing Per Patient",fill="Severity")+
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) +
    geom_errorbar(aes(ymin=pmax(Avg_Prop_Miss-sd,0),ymax = Avg_Prop_Miss + sd), width = .1,position=position_dodge(.9)) 
  prop_plot_Severity4
  Site = rep(siteid,nrow(df_prop_diff_severity4))
  df_prop_diff_severity4$site = Site
  prop_plot_diff_severity4 = ggplot(data = df_prop_diff_severity4,
                                    aes(x = reorder(Labs,Avg_Prop_Miss),
                                        y = Avg_Prop_Miss,
                                        fill = Site)) +
    geom_bar(stat = "identity") + 
    labs(title = "Stratified by Severity",x = "Labs",y = "Mean Difference in Proportion Missing Per Patient",fill="Site") +
    theme(axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1)) 
  prop_plot_diff_severity4
  props_S_quant7= colQuantiles(as.matrix(props_S))
  props_NS_quant7 = colQuantiles(as.matrix(props_NS))
  
  
  
  
  
  
  
  ## temporal trends in pairs of labs
  
  demo_raw$date_diff <- as.Date(as.character(demo_raw$last_discharge_date), format="%Y-%m-%d")-
    as.Date(as.character(demo_raw$admission_date), format="%Y-%m-%d")
  dates_diff = c()
  for (i in 1:nrow(demo_raw)){
    dates_diff = c(dates_diff,rep(demo_raw$date_diff[i],sum(patient_obs_wide$patient_num==demo_raw$patient_num[i])))
  }
  patient_obs_wide = patient_obs_wide[patient_obs_wide$days_since_admission <= dates_diff,]
  dat = patient_obs_wide
  min_corr = 0
  prop_greater = 0
  te_labs2 <- dat %>% 
    filter(days_since_admission < 100,
           days_since_admission >= 0) %>% 
    mutate(te = patient_num %in% unique(te_patients$patient_num))
  #Calculate correlation matrices between lab missing indicators for each time point
  sim_mat = matrix(rep(0),ncol(te_labs2)-4,ncol(te_labs2)-4)
  mat_list = list()
  for (h in 0:max(te_labs2$days_since_admission)){
    data1 = te_labs2[which(te_labs2$days_since_admission==h),]
    na_stats <- data1 %>% 
      dplyr::select(- c(days_since_admission,patient_num,te,severity)) %>% 
      is.na() %>% 
      `!` 
    for (i in 1:ncol(na_stats)){
      na_stats[,i] <- as.integer(as.logical(na_stats[,i]))
    }
    sim_mat = matrix(rep(0),ncol(na_stats),ncol(na_stats))
    for (i in 1:ncol(na_stats)){
      for (j in 1:ncol(na_stats)){
        if (nrow(na_stats) >= 5) {
          sim_mat[i,j] <- cor(na_stats[,i],na_stats[,j],method="spearman")
        }
      }
    }
    rownames(sim_mat) <- colnames(na_stats)
    colnames(sim_mat) <- colnames(na_stats)
    mat_list[[h+1]] = sim_mat
  }
  #Add correlation values to a dataframe if a lab combination passes the threshold
  df <- data.frame(matrix(ncol=0,nrow=length(mat_list)))
  df_ind <- data.frame(matrix(ncol=0,nrow=length(mat_list)))
  names = list()
  final_names = c()
  count = 0
  for (i in 1:ncol(na_stats)){
    for (j in 1:ncol(na_stats)){
      name1 <- c(row.names(mat_list[[1]])[i],colnames(mat_list[[1]])[j])
      name1 <- do.call(paste, c(as.list(name1), sep = ", "))
      name2 <- c(row.names(mat_list[[1]])[j],colnames(mat_list[[1]])[i])
      name2 <- do.call(paste, c(as.list(name2), sep = ", "))
      if (i != j && !(name1 %in% final_names) && !(name2 %in% final_names)) {
        test_vec = c()
        indices = c()
        ct = 0
        for (h in 1:length(mat_list)){
          if (!(is.na(mat_list[[h]][i,j]))) {
            ct = ct+1
            test_vec[ct] = mat_list[[h]][i,j]
            indices[ct] = h
          }else{
            ct = ct + 1
            test_vec[ct] = test_vec[ct-1]
            indices[ct] = h
          }
        }

         count = count+1
         final_names[count] = name1
         df[,count] = test_vec
         df_ind[,count] = indices
      }
    }
  }
  
  colnames(df) <- final_names
  results <- list()
  results[[1]] <- df
  results[[2]] <- mat_list
  temporal = df
  #Get df for severe
  dat = patient_obs_wide
  dat = dat[dat$severity=='severe',]
  min_corr = 0
  prop_greater = 0
  te_labs2 <- dat %>% 
    filter(days_since_admission < 100,
           days_since_admission >= 0) %>% 
    mutate(te = patient_num %in% unique(te_patients$patient_num))
  #Calculate correlation matrices between lab missing indicators for each time point
  sim_mat = matrix(rep(0),ncol(te_labs2)-4,ncol(te_labs2)-4)
  mat_list = list()
  for (h in 0:max(te_labs2$days_since_admission)){
    data1 = te_labs2[which(te_labs2$days_since_admission==h),]
    na_stats <- data1 %>% 
      dplyr::select(- c(days_since_admission,patient_num,te,severity)) %>% 
      is.na() %>% 
      `!` 
    for (i in 1:ncol(na_stats)){
      na_stats[,i] <- as.integer(as.logical(na_stats[,i]))
    }
    sim_mat = matrix(rep(0),ncol(na_stats),ncol(na_stats))
    for (i in 1:ncol(na_stats)){
      for (j in 1:ncol(na_stats)){
        if (nrow(na_stats) >= 5) {
          sim_mat[i,j] <- cor(na_stats[,i],na_stats[,j],method="spearman")
        }
      }
    }
    rownames(sim_mat) <- colnames(na_stats)
    colnames(sim_mat) <- colnames(na_stats)
    mat_list[[h+1]] = sim_mat
  }
  #Add correlation values to a dataframe if a lab combination passes the threshold
  df <- data.frame(matrix(ncol=0,nrow=length(mat_list)))
  df_ind <- data.frame(matrix(ncol=0,nrow=length(mat_list)))
  names = list()
  final_names = c()
  count = 0
  for (i in 1:ncol(na_stats)){
    for (j in 1:ncol(na_stats)){
      name1 <- c(row.names(mat_list[[1]])[i],colnames(mat_list[[1]])[j])
      name1 <- do.call(paste, c(as.list(name1), sep = ", "))
      name2 <- c(row.names(mat_list[[1]])[j],colnames(mat_list[[1]])[i])
      name2 <- do.call(paste, c(as.list(name2), sep = ", "))
      if (i != j && !(name1 %in% final_names) && !(name2 %in% final_names)) {
        test_vec = c()
        indices = c()
        ct = 0
        for (h in 1:length(mat_list)){
          if (!(is.na(mat_list[[h]][i,j]))) {
            ct = ct+1
            test_vec[ct] = mat_list[[h]][i,j]
            indices[ct] = h
          }else{
            ct = ct + 1
            test_vec[ct] = test_vec[ct-1]
            indices[ct] = h
          }
        }
        
         count = count+1
         final_names[count] = name1
         df[,count] = test_vec
         df_ind[,count] = indices
      }
    }
  }
  
  colnames(df) <- final_names
  results <- list()
  results[[1]] <- df
  results[[2]] <- mat_list
  temporal_severe = df
  #get df for non-severe
  dat = patient_obs_wide
  dat = dat[dat$severity=='nonsevere',]
  min_corr = 0
  prop_greater = 0
  te_labs2 <- dat %>% 
    filter(days_since_admission < 100,
           days_since_admission >= 0) %>% 
    mutate(te = patient_num %in% unique(te_patients$patient_num))
  #Calculate correlation matrices between lab missing indicators for each time point
  sim_mat = matrix(rep(0),ncol(te_labs2)-4,ncol(te_labs2)-4)
  mat_list = list()
  for (h in 0:max(te_labs2$days_since_admission)){
    data1 = te_labs2[which(te_labs2$days_since_admission==h),]
    na_stats <- data1 %>% 
      dplyr::select(- c(days_since_admission,patient_num,te,severity)) %>% 
      is.na() %>% 
      `!` 
    for (i in 1:ncol(na_stats)){
      na_stats[,i] <- as.integer(as.logical(na_stats[,i]))
    }
    sim_mat = matrix(rep(0),ncol(na_stats),ncol(na_stats))
    for (i in 1:ncol(na_stats)){
      for (j in 1:ncol(na_stats)){
        if (nrow(na_stats) >= 5) {
          sim_mat[i,j] <- cor(na_stats[,i],na_stats[,j],method="spearman")
        }
      }
    }
    rownames(sim_mat) <- colnames(na_stats)
    colnames(sim_mat) <- colnames(na_stats)
    mat_list[[h+1]] = sim_mat
  }
  #Add correlation values to a dataframe if a lab combination passes the threshold
  df <- data.frame(matrix(ncol=0,nrow=length(mat_list)))
  df_ind <- data.frame(matrix(ncol=0,nrow=length(mat_list)))
  names = list()
  final_names = c()
  count = 0
  for (i in 1:ncol(na_stats)){
    for (j in 1:ncol(na_stats)){
      name1 <- c(row.names(mat_list[[1]])[i],colnames(mat_list[[1]])[j])
      name1 <- do.call(paste, c(as.list(name1), sep = ", "))
      name2 <- c(row.names(mat_list[[1]])[j],colnames(mat_list[[1]])[i])
      name2 <- do.call(paste, c(as.list(name2), sep = ", "))
      if (i != j && !(name1 %in% final_names) && !(name2 %in% final_names)) {
        test_vec = c()
        indices = c()
        ct = 0
        for (h in 1:length(mat_list)){
          if (!(is.na(mat_list[[h]][i,j]))) {
            ct = ct+1
            test_vec[ct] = mat_list[[h]][i,j]
            indices[ct] = h
          }else{
            ct = ct + 1
            test_vec[ct] = test_vec[ct-1]
            indices[ct] = h
          }
        }
        
         count = count+1
         final_names[count] = name1
         df[,count] = test_vec
         df_ind[,count] = indices
      }
    }
  }
  
  colnames(df) <- final_names
  results <- list()
  results[[1]] <- df
  results[[2]] <- mat_list
  temporal_nonsevere = df
  
  
  
  ## Individual temporal analyses
  
  binby = 3
  wide =  patient_obs_wide %>% 
    filter(days_since_admission <= 150,
           days_since_admission >= 0)
  wide = wide[c("patient_num","days_since_admission","severity", labsPresent)]
  total_days = max(wide$days_since_admission)
  n_bins = as.integer(total_days/binby)
  df <- data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat = wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_severe = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dat[which(dat$severity=='severe'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_severe[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_nonsevere = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat = dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct + binby
    dat = dat[which(dat$severity=='nonsevere'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_nonsevere[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df2 <- data.frame(matrix(ncol=0,nrow=nrow(df)*ncol(df)*3))
  count = 0
  for (i in 1:ncol(df)){
    df2[(1+count):(3*nrow(df)+count),1] = rep(colnames(dat)[i],3*nrow(df))
    count = count+3*nrow(df)
  }
  bin1 = "0-3"
  bin2 = "4-6"
  bin3 = "7-9"
  bin4 = "10-12"
  bin5 = "13-15"
  bin6 = "16-18"
  bin7 = "19-21"
  bin8 = "22-24"
  bin9 = "25-27"
  bin10= "28-30"
  bin11 = "31-33"
  bin12 = "34-36"
  bin13 = "37-39"
  bin14 = "40-42"
  bin15 = "43-45"
  bin16 = "46-48"
  bin17 = "49-51"
  bin18 = "52-54"
  bin19 = "55-57"
  bin20 = "58-60"
  bin21 = "61-63"
  bin22 = "64-66"
  bin23 = "67-69"
  bin24 = "70-72"
  bin25 = "73-75"
  bin26 = "76-78"
  bin27 = "79-81"
  bin28 = "82-84"
  bin29 = "85-87"
  bin30 = "88-90"
  bin31 = "91-93"
  bin32 = "94-96"
  bin33 = "97-99"
  bin34 = "100-102"
  bin35 = "103-105"
  bin36 = "106-108"
  bin37 = "109-111"
  bin38 = "112-114"
  bin39 = "115-117"
  bin40 = "118-120"
  bin41 = "121-123"
  bin42 = "124-126"
  bin43 = "127-129"
  bin44 = "130-132"
  bin45 = "133-135"
  bin46 = "136-138"
  bin47 = "139-141"
  bin48 = "142-144"
  bin49 = "145-147"
  bin50 = "148-150"
  
  bins = c(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13,bin14,bin15,bin16,bin17,bin18,bin19
           ,bin20,bin21,bin22,bin23,bin24,bin25,bin26,bin27,bin28,bin29,bin30,bin31,bin32,bin33,bin34,bin35,bin36,bin37,bin38
           ,bin39,bin40,bin41,bin42,bin43,bin44,bin45,bin46,bin47,bin48,bin49,bin50)
  bins = bins[1:n_bins]
  df2[,2] = rep(c(bins,bins,bins),ncol(df))
  c1 = rep("All patients", nrow(df))
  c2 = rep("Severe patients",nrow(df))
  c3 = rep("Non-severe patients",nrow(df))
  ctotal = c(c1,c2,c3)
  df2[,3] = rep(ctotal,ncol(df))
  ct = 0
  for (i in 1:ncol(df)){
    df2[(ct+1):(nrow(df)+ct),4] = df[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_severe[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_nonsevere[,i]
    ct= ct + nrow(df)
  }
  colnames(df2) = c("lab","days_since_admission","name","proportion")
  df2$days_since_admission = factor(df2$days_since_admission,levels = bins)
  level_order = labsPresent
  longhaul_comb_plot = df2 %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab,level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    facet_wrap(~name, nrow = 1) +
    labs(
      title = "Proportion of missing lab values over time",
      x = "Binned days since admission",
      y = NULL, fill = "Prop Missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6),
      axis.ticks.y = element_blank()
    )
  longhaul_comb_plot
  longhaul_comb3 = df2
  df2_all = df2[df2$name=="All patients",] %>%
    dplyr::select(-c(name))
  longhaul_plot = df2_all %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab,level=level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    labs(
      title = "Proportion of missing lab values over time",
      x = "Binned days since admission",
      y = NULL, fill = "Prop missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5,size=9),
      axis.ticks.y = element_blank()
    )
  longhaul_plot
  longhaul3 = df2_all
  
  
  #now bin by 7 days:
  
  
  binby = 7
  wide =  patient_obs_wide %>% 
    filter(days_since_admission <= 150,
           days_since_admission >= 0)
  wide = wide[c("patient_num","days_since_admission","severity", labsPresent)]
  total_days = max(wide$days_since_admission)
  n_bins = as.integer(total_days/binby)
  df <- data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat = wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_severe = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dat[which(dat$severity=='severe'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_severe[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_nonsevere = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat = dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct + binby
    dat = dat[which(dat$severity=='nonsevere'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_nonsevere[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df2 <- data.frame(matrix(ncol=0,nrow=nrow(df)*ncol(df)*3))
  count = 0
  for (i in 1:ncol(df)){
    df2[(1+count):(3*nrow(df)+count),1] = rep(colnames(dat)[i],3*nrow(df))
    count = count+3*nrow(df)
  }
  bin1 = "0-7"
  bin2 = "8-14"
  bin3 = "15-21"
  bin4 = "22-28"
  bin5 = "29-35"
  bin6 = "36-42"
  bin7 = "43-49"
  bin8 = "50-56"
  bin9 = "57-63"
  bin10= "64-70"
  bin11 = "71-77"
  bin12 = "78-84"
  bin13 = "85-91"
  bin14 = "92-98"
  bin15 = "99-105"
  bin16 = "106-112"
  bin17 = "113-119"
  bin18 = "120-126"
  bin19 = "127-133"
  bin20 = "134-140"
  bin21 = "141-147"
  bin22 = "148-154"
  
  bins = c(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13,bin14,bin15,bin16,bin17,bin18,bin19
           ,bin20,bin21,bin22)
  bins = bins[1:n_bins]
  df2[,2] = rep(c(bins,bins,bins),ncol(df))
  c1 = rep("All patients", nrow(df))
  c2 = rep("Severe patients",nrow(df))
  c3 = rep("Non-severe patients",nrow(df))
  ctotal = c(c1,c2,c3)
  df2[,3] = rep(ctotal,ncol(df))
  ct = 0
  for (i in 1:ncol(df)){
    df2[(ct+1):(nrow(df)+ct),4] = df[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_severe[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_nonsevere[,i]
    ct= ct + nrow(df)
  }
  colnames(df2) = c("lab","days_since_admission","name","proportion")
  df2$days_since_admission = factor(df2$days_since_admission,levels = bins)
  level_order = labsPresent
  longhaul_comb_plot = df2 %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab,level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    facet_wrap(~name, nrow = 1) +
    labs(
      title = "Proportion of missing lab values over time",
      x = "Binned days since admission",
      y = NULL, fill = "Prop Missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6),
      axis.ticks.y = element_blank()
    )
  longhaul_comb_plot
  longhaul_comb7 = df2
  df2_all = df2[df2$name=="All patients",] %>%
    dplyr::select(-c(name))
  longhaul_plot = df2_all %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab,level=level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    labs(
      title = "Proportion of missing lab values over time",
      x = "Binned days since admission",
      y = NULL, fill = "Prop missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5,size=9),
      axis.ticks.y = element_blank()
    )
  longhaul_plot
  longhaul7 = df2_all
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #Calculate proportion of existing lab values in days since admission up to 30 days. Stratify by severe vs. non-severe.
  binby = 1
  wide =  patient_obs_wide %>% 
    filter(days_since_admission <= 30,
           days_since_admission >= 0)
  total_days = max(wide$days_since_admission)
  df <- data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days)) {
    dat = wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_severe = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dat[which(dat$severity=='severe'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_severe[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_nonsevere = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat = dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct + binby
    dat = dat[which(dat$severity=='nonsevere'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_nonsevere[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df2 <- data.frame(matrix(ncol=0,nrow=nrow(df)*ncol(df)*3))
  count = 0
  for (i in 1:ncol(df)){
    df2[(1+count):(3*nrow(df)+count),1] = rep(colnames(dat)[i],3*nrow(df))
    count = count+3*nrow(df)
  }
  if (total_days < 30){                        
   days = seq(1,total_days,1)
  }else{
   days = seq(1,30,1)
  }
  df2[,2] = rep(c(days,days,days),ncol(df))
  c1 = rep("All patients", nrow(df))
  c2 = rep("Severe patients",nrow(df))
  c3 = rep("Non-severe patients",nrow(df))
  ctotal = c(c1,c2,c3)
  df2[,3] = rep(ctotal,ncol(df))
  ct = 0
  for (i in 1:ncol(df)){
    df2[(ct+1):(nrow(df)+ct),4] = df[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_severe[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_nonsevere[,i]
    ct= ct + nrow(df)
  }
  colnames(df2) = c("lab","days_since_admission","name","proportion")
  df2$days_since_admission = factor(df2$days_since_admission,levels = bins)
  level_order = labsPresent
  if (total_days < 30){                        
   days = seq(1,total_days,1)
  }else{
   days = seq(1,30,1)
  }
  df2[,2] = rep(c(days,days,days),ncol(df))
  shorthaul_comb= df2
  shorthaul_comb_plot = df2 %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab, level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    facet_wrap(~name, nrow = 1) +
    labs(title = "Proportion of missing lab values over time",
         x = "Binned days since admission",
         y = NULL, fill = "Prop missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6),
      axis.ticks.y = element_blank()
    )
  shorthaul_comb_plot
  df2_all = df2[df2$name=="All patients",] %>%
    dplyr::select(-c(name))
  shorthaul= df2_all
  shorthaul_plot = df2_all %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab, level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    labs(
      title = "Proportion of missing lab values over time",
      x = "Binned days since admission",
      y = NULL, fill = "Missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5,size=9),
      axis.ticks.y = element_blank()
    )
  shorthaul_plot 
  
  
  
  #Calculate proportion of existing lab values in days since admission up to 7 days. Stratify by severe vs. non-severe.
  binby = 1
  wide =  patient_obs_wide %>% 
    filter(days_since_admission <= 7,
           days_since_admission >= 0)
  total_days = max(wide$days_since_admission)
  df <- data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days)) {
    dat = wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_severe = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct+binby
    dat = dat[which(dat$severity=='severe'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_severe[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df_nonsevere = data.frame(matrix(ncol=0,nrow=total_days/binby))
  ct = 0
  for (i in 1:(total_days/binby)) {
    dat = dat =wide[which(wide$days_since_admission >= (ct+1) & wide$days_since_admission <= (binby+1+ct)),]
    ct = ct + binby
    dat = dat[which(dat$severity=='nonsevere'),]
    dat = dplyr::select(dat,-c(days_since_admission,patient_num,severity))
    for (j in 1:(ncol(dat))){
      df_nonsevere[i,j] = (sum(is.na(dat[,j])))/nrow(dat[,j])
    }
  }
  df2 <- data.frame(matrix(ncol=0,nrow=nrow(df)*ncol(df)*3))
  count = 0
  for (i in 1:ncol(df)){
    df2[(1+count):(3*nrow(df)+count),1] = rep(colnames(dat)[i],3*nrow(df))
    count = count+3*nrow(df)
  }
  if (total_days < 7){                        
   days = seq(1,total_days,1)
  }else{
   days = seq(1,7,1)
  }
  df2[,2] = rep(c(days,days,days),ncol(df))
  c1 = rep("All patients", nrow(df))
  c2 = rep("Severe patients",nrow(df))
  c3 = rep("Non-severe patients",nrow(df))
  ctotal = c(c1,c2,c3)
  df2[,3] = rep(ctotal,ncol(df))
  ct = 0
  for (i in 1:ncol(df)){
    df2[(ct+1):(nrow(df)+ct),4] = df[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_severe[,i]
    ct = ct + nrow(df)
    df2[(ct+1):(nrow(df)+ct),4] = df_nonsevere[,i]
    ct= ct + nrow(df)
  }
  colnames(df2) = c("lab","days_since_admission","name","proportion")
  df2$days_since_admission = factor(df2$days_since_admission,levels = bins)
  level_order = labsPresent
  if (total_days < 7){                        
   days = seq(1,total_days,1)
  }else{
   days = seq(1,7,1)
  }
  df2[,2] = rep(c(days,days,days),ncol(df))
  shorthaul_comb2= df2
  shorthaul_comb_plot = df2 %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab, level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    facet_wrap(~name, nrow = 1) +
    labs(title = "Proportion of missing lab values over time",
         x = "Binned days since admission",
         y = NULL, fill = "Prop missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6),
      axis.ticks.y = element_blank()
    )
  shorthaul_comb_plot
  df2_all = df2[df2$name=="All patients",] %>%
    dplyr::select(-c(name))
  shorthaul= df2_all
  shorthaul_plot = df2_all %>%
    ggplot(aes(x = days_since_admission, fill = proportion, y = factor(lab, level_order))) +
    geom_tile(colour = "white", size = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "darkblue", high = "lightgrey") +
    labs(
      title = "Proportion of missing lab values over time",
      x = "Binned days since admission",
      y = NULL, fill = "Missing"
    ) +
    
    theme(
      plot.margin=unit(c(1.02,1.6,1,0),"cm"),
      panel.grid.major = element_blank(),
      legend.position = c(1.078, 0.5),
      legend.key.height = unit(1,'cm'),
      legend.key.width=unit(.4,'cm'),
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5,size=9),
      axis.ticks.y = element_blank()
    )
  shorthaul_plot 
  
  
  
  
  
  
  
  
  
  
  
  
  ####### LDA Topic Modeling
  print("begin topic modeling")
  
  #Read in thrombotic icd codes. Add column "truncated code" for first three characters of each code. Extract unique truncated codes.
  thrombo_codes <- c("I74", "I75", "I76", "I21", "I22", "I23", "I26", "I27", "Z86", "I63", "I67", "I81", "I82","444","445",
"410","415","V12","434","437","452","453","D65","P60","286","776")
  #Generate list of neuro codes (obtained from neuro group)
  neuro_codes <- c("R41","R27","R42","G44","G03","G04","G72","M60","G61","G65","R43","G93","F29","G40","G45","G46","I60","I61","I62","I67","H54","40","298","307","320","321","322","323","330","331","339","345","348","357","359","369","430","431","432")
  ARDs_codes <- c("J80","518")
  
  #Isolate patients with thrombotic events and remove the days_since_admission column. 
  thrombo_patients <- patient_obs %>%
    filter(
      concept_type == "DIAG-ICD10",
      concept_code %in% thrombo_codes,
      days_since_admission >= 0
    ) %>%
    dplyr::select(-days_since_admission) %>%
    distinct()
  neuro_patients <- patient_obs %>%
    filter(
      concept_type == "DIAG-ICD10",
      concept_code %in% neuro_codes,
      days_since_admission >= 0
    ) %>%
    dplyr::select(-days_since_admission) %>%
    distinct()
  ARDs_patients <- patient_obs %>%
    filter(
      concept_type == "DIAG-ICD10",
      concept_code %in% ARDs_codes,
      days_since_admission >= 0
    ) %>%
    dplyr::select(-days_since_admission) %>%
    distinct()
  #Isolate information for unique patients. Sum individual patients with thrombotic events. 
  demo_df <- demo_raw %>%
    mutate(
      TE = patient_num %in% unique(thrombo_patients$patient_num),
      deceased = deceased == 1,
      severe = severe == 1
    )
  demo_df$TE %>% sum()
  demo_df <- demo_df %>%
    mutate(
      Neuro = patient_num %in% unique(neuro_patients$patient_num),
      deceased = deceased == 1,
      severe = severe == 1
    )
  demo_df$Neuro %>% sum()
  demo_df <- demo_df %>%
    mutate(
      ARDs = patient_num %in% unique(ARDs_patients$patient_num),
      deceased = deceased == 1,
      severe = severe == 1
    )
  demo_df$ARDs %>% sum()
  #Generate a table that contains information about each lab for each individual patient. Each observation is identified by their patient number and days since admission, and each lab name and value is included. If there is more than one lab value for that identifier, they are averaged.
  
  
  upper_day1 <- 10
  err_xgbs <- vector("numeric")
  err_svms <- vector("numeric")
  all_evals <- list()
  #filter out observations with days since admission >= a threshold (upper_day, in this case 9 days) and days since admission  <= 0. Take out the days_since_admission column. Obtain number of lab values for each TE patient.
  missing_by_patient <- patient_obs_wide %>%
    filter(
      days_since_admission <= upper_day1, # 10-day window
      days_since_admission >= 0
    ) %>%
    group_by(patient_num) %>%
    summarise(across(-days_since_admission, function(x) {
      sum(!is.na(x))
    })) %>%
    left_join(demo_df, by = "patient_num")
  ct = 0
  if ("Troponin_high" %in% colnames(missing_by_patient)){
    if ("Troponin_normal" %in% colnames(missing_by_patient)){
      ct = ct + 1
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="Troponin_high")] <- "Troponin"
    }
  }
  if ("Troponin_normal" %in% colnames(missing_by_patient)) {
    if ("Troponin_high" %in% colnames(missing_by_patient)) {
      missing_by_patient$Troponin = missing_by_patient$Troponin_normal +   missing_by_patient$Troponin_high
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="Troponin_normal")] <- "Troponin"
    }
  }
  if ("DDU" %in% colnames(missing_by_patient)){
    if ("FEU" %in% colnames(missing_by_patient)) {
      ct = ct + 1
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="DDU")] <- "D_Dimer"
    }
  }
  if ("FEU" %in% colnames(missing_by_patient)) {
    if ("D-Dimer" %in% colnames(missing_by_patient)) {
      missing_by_patient$D_Dimer = missing_by_patient$FEU +   missing_by_patient$D_Dimer
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="FEU")] <- "D_Dimer"
    }
  }
  drops = c("Troponin_high","Troponin_normal","FEU","DDU")
  missing_by_patient= missing_by_patient[,!(names(missing_by_patient) %in% drops)]
  allLabs <- c("Fibrinogen","Procalcitonin","D_Dimer","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
  labsPresent <- allLabs[ allLabs %in% colnames(missing_by_patient)]
  lab_names = labsPresent
  #Obtain total number of observations per lab for TE patients
  uniq_vals <-
    apply(missing_by_patient, 2, function(x) {
      length(unique(x))
    })
  #Remove patients that have no lab values
  missing_by_patient <- missing_by_patient[, uniq_vals > 1] %>%
    mutate(sum_labs = rowSums(across(all_of(lab_names)))) %>% 
    filter(sum_labs > 0)
  #Reduce dataframe to matrix with labs and values
  x_mat <- missing_by_patient[, labsPresent]
  
  print("data processed. Calculate exclusivity")
  
  exclusivity <- function(mod.out, M = 10, frexw = .7) {
    w <- frexw
    if (length(mod.out$beta$logbeta) != 1) stop("Exclusivity calculation only designed for models without content covariates")
    tbeta <- t(exp(mod.out$beta$logbeta[[1]]))
    s <- rowSums(tbeta)
    mat <- tbeta / s # normed by columns of beta now.
    ex <- apply(mat, 2, rank) / nrow(mat)
    fr <- apply(tbeta, 2, rank) / nrow(mat)
    frex <- 1 / (w / ex + (1 - w) / fr)
    index <- apply(tbeta, 2, order, decreasing = TRUE)[1:M, ]
    out <- vector(length = ncol(tbeta))
    for (i in 1:ncol(frex)) {
      out[i] <- sum(frex[index[, i], i])
    }
    out
  }
  x_dfm <- x_mat %>%
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "lab", values_to = "n") %>%
    tidytext::cast_dfm(id, lab, n)
  
  topic_diagnostics = c()
  beta_mat = c()
  beta_mat2 = c()
  TE_results <- c()
  Severity_results <- c()
  Neuro_results <- c()
  ARDs_results <- c()
  
  TE_results2 <- c()
  Severity_results2 <- c()
  Neuro_results2 <- c()
  ARDs_results2 <- c()
  
                        
  print("find K part 1")                      
  
  try(system.time(many_models <- data.frame(K = seq(2, 8, 1)) %>%
                mutate(topic_model = furrr::future_map(K, ~ stm(
                  x_dfm,
                  K = .,
                  seed = TRUE,
                  verbose = FALSE
                )))))
  heldout <- try(make.heldout(x_dfm))
  k_result <- try(many_models %>%
    mutate(
      exclusivity = map(topic_model, exclusivity),
      semantic_coherence = map(topic_model, semanticCoherence, x_dfm),
      eval_heldout = map(topic_model, eval.heldout, heldout$missing),
      residual = map(topic_model, checkResiduals, x_dfm),
      bound = map_dbl(topic_model, function(x) max(x$convergence$bound)),
      lfact = map_dbl(topic_model, function(x) lfactorial(x$settings$dim$K)),
      lbound = bound + lfact,
      iterations = map_dbl(topic_model, function(x) length(x$convergence$bound))
    ))
  topic_diagnostics = try(k_result %>%
    transmute(K,
              `Lower bound` = lbound,
              Residuals = map_dbl(residual, "dispersion"),
              `Semantic coherence` = map_dbl(semantic_coherence, mean),
              `Held-out likelihood` = map_dbl(eval_heldout, "expected.heldout")
    ) %>%
    gather(Metric, Value, -K)) 
  topic_diagnostics_wide = try(data.frame(matrix(ncol=0,nrow=7)))
  topic_diagnostics_wide$K= try(topic_diagnostics$K[1:7])
  topic_diagnostics_wide$LB= try(topic_diagnostics$Value[1:7])
  topic_diagnostics_wide$Res = try(topic_diagnostics$Value[8:14])
  topic_diagnostics_wide$SC = try(topic_diagnostics$Value[15:21])
  topic_diagnostics_wide$HL = try(topic_diagnostics$Value[22:28])
  ind = try(c(which.max(topic_diagnostics_wide$LB)+1,which.max(topic_diagnostics_wide$SC)+1,which.max(topic_diagnostics_wide$HL)+1,which.min(topic_diagnostics_wide$Res)+1))
  Modes <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
  }
  try(Modes(ind))
                      
  print("plot diagnostics") 
  
  topic_diagnostics_plot = try(k_result %>%
    transmute(K,
              `Lower bound` = lbound,
              Residuals = map_dbl(residual, "dispersion"),
              `Semantic coherence` = map_dbl(semantic_coherence, mean),
              `Held-out likelihood` = map_dbl(eval_heldout, "expected.heldout")
    ) %>%
    gather(Metric, Value, -K) %>%
    ggplot(aes(K, Value, color = Metric)) +
    geom_line() +
    geom_point() +
    guides(color = FALSE) +
    scale_x_continuous(breaks = seq(2, 12, 2)) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(
      x = "K (number of topics)",
      y = NULL,
      title = "Model diagnostics by number of topics"
    ))
  topic_diagnostics_plot
  
  print("fit stm 1") 
  
  #Beta matrix contains log probabilities of labs in topics. Generate heat map of beta values for each lab.
  
  K <- try(Modes(ind)[1])
  K2  = 5
  stmfit <- try(stm(x_dfm, K = K, verbose = FALSE, init.type = "Spectral", seed = TRUE))
  stmfit2 <- try(stm(x_dfm, K = K2, verbose = FALSE, init.type = "Spectral", seed = TRUE))
  stmfit_beta <- try(stmfit$beta)
  stmfit_beta2 <- try(stmfit2$beta)
  
  K=try(Modes(ind)[1])
  beta_mat <- try(exp(stmfit$beta$logbeta[[1]]))
  colnames(beta_mat) <- try(stmfit$vocab)
  rownames(beta_mat) <- try(paste("Topic", 1:K))
  try(heatmap((beta_mat), main = "Ordered",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25)))
  try(legend(x="topleft", legend=c("min", "ave", "max"), 
         fill=colorRampPalette(brewer.pal(8, "Blues"))(3)))
  labs = try(stmfit$vocab)
  
  beta_mat2 <- try(exp(stmfit2$beta$logbeta[[1]]))
  colnames(beta_mat2) <- try(stmfit2$vocab)
  rownames(beta_mat2) <- try(paste("Topic", 1:K2))
  try(heatmap((beta_mat2), main = "Ordered",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25)))
  try(legend(x="topleft", legend=c("min", "ave", "max"), 
         fill=colorRampPalette(brewer.pal(8, "Blues"))(3)))
  labs = try(stmfit2$vocab)
  
  
  theta <- try(stmfit$theta)
  theta2 <- try(stmfit2$theta)
  corrFunc <- function(var1, var2) {
    result = cor.test(var1,var2)
    data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
               stringsAsFactors=FALSE)
  }
  ## Pairs of variables for which we want correlations
  theta <- try(stmfit$theta)
  K = try(Modes(ind)[1])
  TE_results <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$TE))[c("estimate","p.value","statistic","method")]))
  Severity_results <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$severe))[c("estimate","p.value","statistic","method")]))
  Neuro_results <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$Neuro))[c("estimate","p.value","statistic","method")]))
  ARDs_results <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$ARDs))[c("estimate","p.value","statistic","method")]))
  
  TE_results_wilcox <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$TE == TRUE],(theta[,x])[missing_by_patient$TE == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Severity_results_wilcox <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$severe == TRUE],(theta[,x])[missing_by_patient$severe == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Neuro_results_wilcox <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$Neuro == TRUE],(theta[,x])[missing_by_patient$Neuro == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  ARDs_results_wilcox <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$ARDs == TRUE],(theta[,x])[missing_by_patient$ARDs == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  
                            
  
  TE_results2 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$TE))[c("estimate","p.value","statistic","method")]))
  Severity_results2 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$severe))[c("estimate","p.value","statistic","method")]))
  Neuro_results2 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$Neuro))[c("estimate","p.value","statistic","method")]))
  ARDs_results2 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$ARDs))[c("estimate","p.value","statistic","method")]))
  
  TE_results_wilcox2 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$TE == TRUE],(theta2[,x])[missing_by_patient$TE == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Severity_results_wilcox2 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$severe == TRUE],(theta2[,x])[missing_by_patient$severe == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Neuro_results_wilcox2 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$Neuro == TRUE],(theta2[,x])[missing_by_patient$Neuro == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  ARDs_results_wilcox2 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$ARDs == TRUE],(theta2[,x])[missing_by_patient$ARDs == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
                             
  
  theta <- try(stmfit$theta)
  K = try(Modes(ind))
  #TE 
  topic_df <- try(data.frame(theta, te = missing_by_patient$TE) %>% 
    pivot_longer(- te, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  
  TE_boxplot <- try(topic_df %>% 
    ggplot(aes(x = te, y = value)) +
    labs(
      x = "Thrombotic Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  TE_boxplot
  #Severity
  topic_df <- try(data.frame(theta, Severity = missing_by_patient$severe) %>% 
    pivot_longer(- Severity, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Severity_boxplot <- try(topic_df %>% 
    ggplot(aes(x = Severity, y = value)) +
    labs(
      x = "Severity",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Severity_boxplot
  #Neuro
  topic_df <- try(data.frame(theta, Neuro = missing_by_patient$Neuro) %>% 
    pivot_longer(- Neuro, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Neuro_boxplot <- try(topic_df %>% 
    ggplot(aes(x = Neuro, y = value)) +
    labs(
      x = "Neurological Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Neuro_boxplot
  #ARDs
  topic_df <- try(data.frame(theta, ARDs = missing_by_patient$ARDs) %>% 
    pivot_longer(- ARDs, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  ARDs_boxplot <- try(topic_df %>% 
    ggplot(aes(x = ARDs, y = value)) +
    labs(
      x = "Acute Respiratory Distress",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  ARDs_boxplot
  
  
  
  
  
  
  
  
  
  
  topic_df2 <- try(data.frame(theta2, te = missing_by_patient$TE) %>% 
    pivot_longer(- te, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  
  TE_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = te, y = value)) +
    labs(
      x = "Thrombotic Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  TE_boxplot2
  #Severity
  topic_df2 <- try(data.frame(theta2, Severity = missing_by_patient$severe) %>% 
    pivot_longer(- Severity, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Severity_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = Severity, y = value)) +
    labs(
      x = "Severity",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Severity_boxplot2
  #Neuro
  topic_df2 <- try(data.frame(theta2, Neuro = missing_by_patient$Neuro) %>% 
    pivot_longer(- Neuro, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Neuro_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = Neuro, y = value)) +
    labs(
      x = "Neurological Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Neuro_boxplot2
  #ARDs
  topic_df2 <- try(data.frame(theta2, ARDs = missing_by_patient$ARDs) %>% 
    pivot_longer(- ARDs, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  ARDs_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = ARDs, y = value)) +
    labs(
      x = "Acute Respiratory Distress",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  ARDs_boxplot2
  
  
  
  
  
  
  
  
  
  
  
  
  print("calculate exclusivity 2") 
  
  
  ## Repeat with 30 days
  
  upper_day2 <- 30
  err_xgbs <- vector("numeric")
  err_svms <- vector("numeric")
  all_evals <- list()
  #filter out observations with days since admission >= a threshold (upper_day, in this case 9 days) and days since admission  <= 0. Take out the days_since_admission column. Obtain number of lab values for each TE patient.
  if (siteid != "penn"){
  missing_by_patient <- patient_obs_wide %>%
    filter(
      days_since_admission <= upper_day2, # 30-day window
      days_since_admission >= 0
    ) %>%
    group_by(patient_num) %>%
    summarise(across(-days_since_admission, function(x) {
      sum(!is.na(x))
    })) %>%
    left_join(demo_df, by = "patient_num")
  ct = 0
  if ("Troponin_high" %in% colnames(missing_by_patient)){
    if ("Troponin_normal" %in% colnames(missing_by_patient)){
      ct = ct + 1
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="Troponin_high")] <- "Troponin"
    }
  }
  if ("Troponin_normal" %in% colnames(missing_by_patient)) {
    if ("Troponin_high" %in% colnames(missing_by_patient)) {
      missing_by_patient$Troponin = missing_by_patient$Troponin_normal +   missing_by_patient$Troponin_high
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="Troponin_normal")] <- "Troponin"
    }
  }
  if ("DDU" %in% colnames(missing_by_patient)){
    if ("FEU" %in% colnames(missing_by_patient)) {
      ct = ct + 1
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="DDU")] <- "D_Dimer"
    }
  }
  if ("FEU" %in% colnames(missing_by_patient)) {
    if ("D-Dimer" %in% colnames(missing_by_patient)) {
      missing_by_patient$D_Dimer = missing_by_patient$FEU +   missing_by_patient$D_Dimer
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="FEU")] <- "D_Dimer"
    }
  }
  drops = c("Troponin_high","Troponin_normal","FEU","DDU")
  missing_by_patient= missing_by_patient[,!(names(missing_by_patient) %in% drops)]
  allLabs <- c("Fibrinogen","Procalcitonin","D_Dimer","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
  labsPresent <- allLabs[ allLabs %in% colnames(missing_by_patient)]
  lab_names = labsPresent
  #Obtain total number of observations per lab for TE patients
  uniq_vals <-
    apply(missing_by_patient, 2, function(x) {
      length(unique(x))
    })
  #Remove patients that have no lab values
  missing_by_patient <- missing_by_patient[, uniq_vals > 1] %>%
    mutate(sum_labs = rowSums(across(all_of(lab_names)))) %>% 
    filter(sum_labs > 0)
  #Reduce dataframe to matrix with labs and values
  x_mat <- missing_by_patient[, labsPresent]
  
  }else{
   missing_by_patient <- patient_obs_wide %>%
    filter(
      days_since_admission <= upper_day2, # 30-day window
      days_since_admission >= 0
    ) %>%
    group_by(patient_num) %>%
    summarise(across(-days_since_admission, function(x) {
      sum(!is.na(x))
    })) %>%
    left_join(demo_df, by = "patient_num")
  ct = 0
  if ("Troponin_high" %in% colnames(missing_by_patient)){
    if ("Troponin_normal" %in% colnames(missing_by_patient)){
      ct = ct + 1
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="Troponin_high")] <- "Troponin"
    }
  }
  if ("Troponin_normal" %in% colnames(missing_by_patient)) {
    if ("Troponin_high" %in% colnames(missing_by_patient)) {
      missing_by_patient$Troponin = missing_by_patient$Troponin_normal +   missing_by_patient$Troponin_high
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="Troponin_normal")] <- "Troponin"
    }
  }
  if ("DDU" %in% colnames(missing_by_patient)){
    if ("FEU" %in% colnames(missing_by_patient)) {
      ct = ct + 1
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="DDU")] <- "D_Dimer"
    }
  }
  if ("FEU" %in% colnames(missing_by_patient)) {
    if ("D-Dimer" %in% colnames(missing_by_patient)) {
      missing_by_patient$D_Dimer = missing_by_patient$FEU +   missing_by_patient$D_Dimer
    }else{
      colnames(missing_by_patient)[which(colnames(missing_by_patient)=="FEU")] <- "D_Dimer"
    }
  }
  drops = c("Troponin_high","Troponin_normal","FEU","DDU")
  missing_by_patient= missing_by_patient[,!(names(missing_by_patient) %in% drops)]
  allLabs <- c("Fibrinogen","Procalcitonin","D_Dimer","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
  labsPresent <- allLabs[ allLabs %in% colnames(missing_by_patient)]
  lab_names = labsPresent
  #Obtain total number of observations per lab for TE patients
  uniq_vals <-
    apply(missing_by_patient, 2, function(x) {
      length(unique(x))
    })
  #Remove patients that have no lab values
  missing_by_patient <- missing_by_patient[, uniq_vals > 1] %>%
    mutate(sum_labs = rowSums(across(all_of(lab_names)))) %>% 
    filter(sum_labs > 0)
  #Reduce dataframe to matrix with labs and values
  x_mat <- missing_by_patient[, labsPresent]
   }
  
  exclusivity <- function(mod.out, M = 10, frexw = .7) {
    w <- frexw
    if (length(mod.out$beta$logbeta) != 1) stop("Exclusivity calculation only designed for models without content covariates")
    tbeta <- t(exp(mod.out$beta$logbeta[[1]]))
    s <- rowSums(tbeta)
    mat <- tbeta / s # normed by columns of beta now.
    ex <- apply(mat, 2, rank) / nrow(mat)
    fr <- apply(tbeta, 2, rank) / nrow(mat)
    frex <- 1 / (w / ex + (1 - w) / fr)
    index <- apply(tbeta, 2, order, decreasing = TRUE)[1:M, ]
    out <- vector(length = ncol(tbeta))
    for (i in 1:ncol(frex)) {
      out[i] <- sum(frex[index[, i], i])
    }
    out
  }
  x_dfm <- x_mat %>%
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "lab", values_to = "n") %>%
    tidytext::cast_dfm(id, lab, n)
  
  topic_diagnostics2 = c()
  beta_mat3 = c()
  beta_mat4 = c()
  TE_results3 <- c()
  Severity_results3 <- c()
  Neuro_results3 <- c()
  ARDs_results3 <- c()
  
  TE_results4 <- c()
  Severity_results4 <- c()
  Neuro_results4 <- c()
  ARDs_results4 <- c()
  
                        
  print("get optimal k 2")                       
  
  try(system.time(many_models <- data.frame(K = seq(2, 8, 1)) %>%
                mutate(topic_model = furrr::future_map(K, ~ stm(
                  x_dfm,
                  K = .,
                  seed = TRUE,
                  verbose = FALSE
                )))))
  heldout <- try(make.heldout(x_dfm))
  k_result <- try(many_models %>%
    mutate(
      exclusivity = map(topic_model, exclusivity),
      semantic_coherence = map(topic_model, semanticCoherence, x_dfm),
      eval_heldout = map(topic_model, eval.heldout, heldout$missing),
      residual = map(topic_model, checkResiduals, x_dfm),
      bound = map_dbl(topic_model, function(x) max(x$convergence$bound)),
      lfact = map_dbl(topic_model, function(x) lfactorial(x$settings$dim$K)),
      lbound = bound + lfact,
      iterations = map_dbl(topic_model, function(x) length(x$convergence$bound))
    ))
  topic_diagnostics = try(k_result %>%
    transmute(K,
              `Lower bound` = lbound,
              Residuals = map_dbl(residual, "dispersion"),
              `Semantic coherence` = map_dbl(semantic_coherence, mean),
              `Held-out likelihood` = map_dbl(eval_heldout, "expected.heldout")
    ) %>%
    gather(Metric, Value, -K)) 
  topic_diagnostics_wide = try(data.frame(matrix(ncol=0,nrow=7)))
  topic_diagnostics_wide$K= try(topic_diagnostics$K[1:7])
  topic_diagnostics_wide$LB= try(topic_diagnostics$Value[1:7])
  topic_diagnostics_wide$Res = try(topic_diagnostics$Value[8:14])
  topic_diagnostics_wide$SC = try(topic_diagnostics$Value[15:21])
  topic_diagnostics_wide$HL = try(topic_diagnostics$Value[22:28])
  ind = try(c(which.max(topic_diagnostics_wide$LB)+1,which.max(topic_diagnostics_wide$SC)+1,which.max(topic_diagnostics_wide$HL)+1,which.min(topic_diagnostics_wide$Res)+1))
  Modes <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
  }
  try(Modes(ind))
                      
  print("plot diagnostics")  
  
  topic_diagnostics_plot = try(k_result %>%
    transmute(K,
              `Lower bound` = lbound,
              Residuals = map_dbl(residual, "dispersion"),
              `Semantic coherence` = map_dbl(semantic_coherence, mean),
              `Held-out likelihood` = map_dbl(eval_heldout, "expected.heldout")
    ) %>%
    gather(Metric, Value, -K) %>%
    ggplot(aes(K, Value, color = Metric)) +
    geom_line() +
    geom_point() +
    guides(color = FALSE) +
    scale_x_continuous(breaks = seq(2, 12, 2)) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(
      x = "K (number of topics)",
      y = NULL,
      title = "Model diagnostics by number of topics"
    ))
  topic_diagnostics_plot
  
  print("stm fit 2")  
  
  #Beta matrix contains log probabilities of labs in topics. Generate heat map of beta values for each lab.
  
  K <- try(Modes(ind)[1])
  K2  = 5
  stmfit <- try(stm(x_dfm, K = K, verbose = FALSE, init.type = "Spectral", seed = TRUE))
  stmfit2 <- try(stm(x_dfm, K = K2, verbose = FALSE, init.type = "Spectral", seed = TRUE))
  stmfit_beta <- try(stmfit$beta)
  stmfit_beta2 <- try(stmfit2$beta)
  
  K=try(Modes(ind)[1])
  beta_mat3 <- try(exp(stmfit$beta$logbeta[[1]]))
  colnames(beta_mat3) <- try(stmfit$vocab)
  rownames(beta_mat3) <- try(paste("Topic", 1:K))
  try(heatmap((beta_mat3), main = "Ordered",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25)))
  try(legend(x="topleft", legend=c("min", "ave", "max"), 
         fill=colorRampPalette(brewer.pal(8, "Blues"))(3)))
  labs = try(stmfit$vocab)
  
  beta_mat4 <- try(exp(stmfit2$beta$logbeta[[1]]))
  colnames(beta_mat4) <- try(stmfit2$vocab)
  rownames(beta_mat4) <- try(paste("Topic", 1:K2))
  try(heatmap((beta_mat4), main = "Ordered",
          col= colorRampPalette(brewer.pal(8, "Blues"))(25)))
  try(legend(x="topleft", legend=c("min", "ave", "max"), 
         fill=colorRampPalette(brewer.pal(8, "Blues"))(3)))
  labs = try(stmfit2$vocab)
  
  
  theta <- try(stmfit$theta)
  theta2 <- try(stmfit2$theta)
  corrFunc <- function(var1, var2) {
    result = cor.test(var1,var2)
    data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
               stringsAsFactors=FALSE)
  }
  ## Pairs of variables for which we want correlations
  theta <- try(stmfit$theta)
  K = try(Modes(ind)[1])
  TE_results3 <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$TE))[c("estimate","p.value","statistic","method")]))
  Severity_results3 <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$severe))[c("estimate","p.value","statistic","method")]))
  Neuro_results3 <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$Neuro))[c("estimate","p.value","statistic","method")]))
  ARDs_results3 <- try(lapply(1:K, function(x) cor.test(theta[, x], as.numeric(missing_by_patient$ARDs))[c("estimate","p.value","statistic","method")]))
  
  TE_results_wilcox3 <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$TE == TRUE],(theta[,x])[missing_by_patient$TE == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Severity_results_wilcox3 <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$severe == TRUE],(theta[,x])[missing_by_patient$severe == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Neuro_results_wilcox3 <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$Neuro == TRUE],(theta[,x])[missing_by_patient$Neuro == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  ARDs_results_wilcox3 <- try(lapply(1:K, function(x) wilcox.test((theta[,x])[missing_by_patient$ARDs == TRUE],(theta[,x])[missing_by_patient$ARDs == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  
  
  TE_results4 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$TE))[c("estimate","p.value","statistic","method")]))
  Severity_results4 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$severe))[c("estimate","p.value","statistic","method")]))
  Neuro_results4 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$Neuro))[c("estimate","p.value","statistic","method")]))
  ARDs_results4 <- try(lapply(1:K2, function(x) cor.test(theta2[, x], as.numeric(missing_by_patient$ARDs))[c("estimate","p.value","statistic","method")]))
  
  TE_results_wilcox4 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$TE == TRUE],(theta2[,x])[missing_by_patient$TE == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Severity_results_wilcox4 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$severe == TRUE],(theta2[,x])[missing_by_patient$severe == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  Neuro_results_wilcox4 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$Neuro == TRUE],(theta2[,x])[missing_by_patient$Neuro == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  ARDs_results_wilcox4 <- try(lapply(1:K2, function(x) wilcox.test((theta2[,x])[missing_by_patient$ARDs == TRUE],(theta2[,x])[missing_by_patient$ARDs == FALSE],conf.int=TRUE,exact=TRUE)[c("p.value","estimate")]))
  
  
  theta <- try(stmfit$theta)
  K = try(Modes(ind))
  #TE 
  topic_df <- try(data.frame(theta, te = missing_by_patient$TE) %>% 
    pivot_longer(- te, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  
  TE_boxplot <- try(topic_df %>% 
    ggplot(aes(x = te, y = value)) +
    labs(
      x = "Thrombotic Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  TE_boxplot
  #Severity
  topic_df <- try(data.frame(theta, Severity = missing_by_patient$severe) %>% 
    pivot_longer(- Severity, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Severity_boxplot <- try(topic_df %>% 
    ggplot(aes(x = Severity, y = value)) +
    labs(
      x = "Severity",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Severity_boxplot
  #Neuro
  topic_df <- try(data.frame(theta, Neuro = missing_by_patient$Neuro) %>% 
    pivot_longer(- Neuro, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Neuro_boxplot <- try(topic_df %>% 
    ggplot(aes(x = Neuro, y = value)) +
    labs(
      x = "Neurological Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Neuro_boxplot
  #ARDs
  topic_df <- try(data.frame(theta, ARDs = missing_by_patient$ARDs) %>% 
    pivot_longer(- ARDs, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  ARDs_boxplot <- try(topic_df %>% 
    ggplot(aes(x = ARDs, y = value)) +
    labs(
      x = "Acute Respiratory Distress",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  ARDs_boxplot
  
  
  
  
  
  
  
  
  
  
  topic_df2 <- try(data.frame(theta2, te = missing_by_patient$TE) %>% 
    pivot_longer(- te, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  
  TE_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = te, y = value)) +
    labs(
      x = "Thrombotic Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  TE_boxplot2
  #Severity
  topic_df2 <- try(data.frame(theta2, Severity = missing_by_patient$severe) %>% 
    pivot_longer(- Severity, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Severity_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = Severity, y = value)) +
    labs(
      x = "Severity",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Severity_boxplot2
  #Neuro
  topic_df2 <- try(data.frame(theta2, Neuro = missing_by_patient$Neuro) %>% 
    pivot_longer(- Neuro, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  Neuro_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = Neuro, y = value)) +
    labs(
      x = "Neurological Event",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  Neuro_boxplot2
  #ARDs
  topic_df2 <- try(data.frame(theta2, ARDs = missing_by_patient$ARDs) %>% 
    pivot_longer(- ARDs, names_to = 'topic') %>% 
    mutate(topic = gsub('X', 'Topic ', topic)))
  ARDs_boxplot2 <- try(topic_df2 %>% 
    ggplot(aes(x = ARDs, y = value)) +
    labs(
      x = "Acute Respiratory Distress",
      y = "Topic Value"
    ) +
    geom_boxplot(alpha = 0.2) +
    facet_wrap(~ topic, scales = 'free') +
    NULL)
  ARDs_boxplot2
  
  
  
  
  
  
  
  
  
  
  
  ####### Generate Table 1
                                     
  #helper functions
                                     
   count_stats <- function(df, count_var, neg_var, group_var, ...) {
   # summary statistics for survival/severity status
   # count values are obfuscated

   count_var <- sym(count_var)
   neg_var <- sym(neg_var)
   Count_var <- sym(stringr::str_to_title(count_var))
   group_var <- sym(group_var)

   df %>%
     select(group_var, count_var) %>%
     group_by(!!group_var) %>%
     summarise(!!neg_var := sum(!!count_var == 0),
               Total = n(),
               .groups = 'drop') %>%
     blur_it(c(neg_var, 'Total'), ...) %>%
     mutate(!!count_var := Total - !!neg_var) %>%
     # mutate_at(vars(all_of(count_var), all_of(neg_var)),
     #           ~ replace(., is.nan(.)|is.na(.), 0)) %>%
     transmute(
       !!group_var,
       !!neg_var :=
         concat(!!neg_var, if_else(Total == 0, 0, !!neg_var/Total)),
       !!Count_var :=
         concat(!!count_var, if_else(Total == 0, 0, !!count_var/Total))) %>%
     pivot_longer(- !!group_var) %>%
     pivot_wider(names_from = !!group_var,
                 values_from = value)
  }

  continuous_stats <- function(df, cont_var, name, group_var,...) {
   # summary statistics for length of stay
   # count values are obfuscated
   med <- sym(paste('Median', name, '[Min, Max]'))
   mea <- sym(paste('Mean', name, '(SD)'))
   cont_var <- sym(cont_var)
   group_var <- sym(group_var)

   df %>%
     dplyr::select(group_var, cont_var) %>%
     group_by(!!group_var) %>%
     summarise(median_var = median(!!cont_var, na.rm = TRUE),
               min_var = min(!!cont_var, na.rm = TRUE),
               max_var = max(!!cont_var, na.rm = TRUE),
               mean_var = mean(!!cont_var, na.rm = TRUE),
               sd_var = sd(!!cont_var, na.rm = TRUE),
               .groups = 'drop') %>%
     mutate_at(vars(sd_var, median_var, mean_var),
               ~ replace(., is.nan(.)|is.na(.), 0)) %>%
     transmute(
       !!group_var,
       !!med := concat_median(median_var, min_var, max_var),
       !!mea := concat_mean(mean_var, sd_var)) %>%
     pivot_longer(- !!group_var) %>%
     pivot_wider(names_from = !!group_var, values_from = value)
  }

  demo_stats <- function(var, df, group_var, ...){
   group_var <- sym(group_var)
   svar <- sym(var)
   df %>%
     group_by(!!svar) %>%
     dplyr::count(!!group_var, name = 'n_var') %>%
     as.data.frame() %>%
     blur_it('n_var', ...) %>%
     group_by(!!group_var) %>%
     mutate(both_neuro = sum(n_var)) %>%
     ungroup() %>%
     mutate(prop = n_var/both_neuro,
            pres = concat(n_var, n_var/both_neuro)) %>%
     pivot_wider(- c(n_var, both_neuro), names_from = !!group_var,
                 values_from = c(n_var, prop, pres),
                 values_fill = list(n_var = 0, prop = 0, pres = '0 (0%)')) %>%
     mutate(variable = paste(var, !!svar, sep = '.')) %>%
     dplyr::select(- !!svar) %>%
     dplyr::select(variable, everything())
  }

  blur_it <- function(df, vars, blur_abs, mask_thres){
   # Obfuscate count values.
   # If blurring range is +/-3, or blur_abs = 3,
   # the count receive a small addition of a random number from -3 to 3.
   # If a count is less than mask_thres, set that count to 0.

   for (var in vars){
     var <- sym(var)
     blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
     df <- df %>%
       mutate(!!var := !!var + blur_vec,
              !!var := ifelse(abs(!!var) < mask_thres, 0, !!var))
   }
   df
  }

  get_tables <- function(
   # neuro_types,
                        demo_df,
                        # scores_unique,
                        # comorb_names_elix,
                        blur_abs,
                        mask_thres,
                        group_var = 'neuro_post',
                        vars_to_obfs = c('sex',
                                         'age_group',
                                         'race',
                                         'Severity',
                                         'Survival')) {

   # total_patients <- length(unique(demo_df$patient_num))
   demo_obfus_table <- lapply(
     vars_to_obfs,
     demo_stats,
     df = demo_df,
     blur_abs = blur_abs,
     mask_thres = mask_thres,
     group_var = group_var
   ) %>%
     bind_rows()

   # suppressWarnings(
   #   other_obfus_table <-
   #     bind_rows(
   #       continuous_stats(demo_df, 'n_stay', 'length of stay', group_var),
   #       count_stats(demo_df, 'severe', 'Nonsevere', group_var, blur_abs, mask_thres),
   #       continuous_stats(demo_df, 'time_to_severe', 'time to severe', group_var),
   #       count_stats(demo_df, 'deceased', 'Alive', group_var, blur_abs, mask_thres),
   #       continuous_stats(demo_df, 'time_to_death', 'time to death', group_var),
   #       continuous_stats(demo_df, 'n_readmissions', 'number of readmissions', group_var),
   #       continuous_stats(
   #         demo_df,
   #         'time_to_first_readmission',
   #         'time to first readmission',
   #         group_var
   #       ),
   #       continuous_stats(
   #         scores_unique,
   #         'elixhauser_score',
   #         'Elixhauser score',
   #         group_var
   #       )
   #     )
   # )
   # elix_obfus_table1 <-
   #   Reduce(
   #     function(...)
   #       left_join(..., by = c("Comorbidity", "Abbreviation")),
   #     lapply(
   #       neuro_types,
   #       list_table1,
   #       df = scores_unique,
   #       num_pats = total_patients,
   #       comorb_names = comorb_names_elix,
   #       group_var = group_var,
   #       blur_abs = blur_abs,
   #       mask_thres = mask_thres
   #     )
   #   ) %>%
   #   mutate(
   #     n_Total = rowSums(select(., starts_with('n_'))),
   #     prop_Total = n_Total / total_patients
   #   ) %>%
   #   arrange(desc(n_Total))

   list(
     demo_table = demo_obfus_table
     # other_obfus_table = other_obfus_table,
     # elix_obfus_table1 = elix_obfus_table1
   )
  }

  get_table1 <- function(
   df, num_pats, comorbidities,
   pat_col = 'patients', ...
  )
  {
   col_char <- paste('n', pat_col, sep = '_')
   npat_col <- sym(col_char)
   proppat_col = sym(paste('prop', pat_col, sep = '_'))
   comorbidities_map = comorbidities$Abbreviation

   df %>%
     dplyr::select(all_of(comorbidities_map)) %>%
     colSums() %>%
     data.frame(n_patients = .) %>%
     tibble::rownames_to_column("Abbreviation") %>%
     blur_it('n_patients', ...) %>%
     mutate(prop_patients = n_patients/num_pats) %>%
     rename(!!npat_col := n_patients,
            !!proppat_col := prop_patients) %>%
     right_join(comorbidities, ., by = "Abbreviation")
  }

  list_table1 <- function(x, df, num_pats, comorb_names, group_var, ...){
   group_var <- sym(group_var)
   get_table1(
     df %>% filter(!!group_var == x),
     num_pats,
     comorbidities = comorb_names,
     pat_col = x, ...)
  }

  tidywhere <- function (fn){
   # https://github.com/r-lib/tidyselect/issues/201#issuecomment-650547846
   # where is not exported in tidyselect
   predicate <- purrr::as_mapper(fn)
   function(x, ...) {
     out <- predicate(x, ...)
     out
   }
  }
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     
  
  
  thrombo_codes <- c("I74", "I75", "I76", "I21", "I22", "I23", "I26", "I27", "Z86", "I63", "I67", "I81", "I82","444","445",
"410","415","V12","434","437","452","453","D65","P60","286","776")
  
  #Generate list of neuro codes (obtained from neuro group)
  
  neuro_codes <- c("R41","R27","R42","G44","G03","G04","G72","M60","G61","G65","R43","G93","F29","G40","G45","G46","I60","I61","I62","I67","H54","40","298","307","320","321","322","323","330","331","339","345","348","357","359","369","430","431","432")
  
  ARDs_codes <- c("J80","518")
  
  
  thrombo_patient <- obs_raw %>% 
    filter(concept_code %in% thrombo_codes,
           days_since_admission >= 0)
  
  thrombo_patient %>% 
    dplyr::count(patient_num) %>% 
    pull(n) %>%
    hist(breaks = 36)
  
  thrombo_patient_vec <- unique(thrombo_patient$patient_num)
  
  neuro_patient <- obs_raw %>% 
    filter(concept_code %in% neuro_codes,
           days_since_admission >= 0)
  
  neuro_patient %>% 
    dplyr::count(patient_num) %>% 
    pull(n) %>%
    hist(breaks = 36)
  
  neuro_patient_vec <- unique(neuro_patient$patient_num)
  
  ARDs_patient <- obs_raw %>% 
    filter(concept_code %in% ARDs_codes,
           days_since_admission >= 0)
  
  ARDs_patient %>% 
    dplyr::count(patient_num) %>% 
    pull(n) %>%
    hist(breaks = 36)
  
  ARDs_patient_vec <- unique(ARDs_patient$patient_num)
  
  
  
  Severe_patient = patient_obs_wide[patient_obs_wide$severity=="severe",]
  
  Severe_patient_vec <- unique(Severe_patient$patient_num)
  
  demo_raw_male <- demo_raw[demo_raw$sex == "male",]
  Male_patient = patient_obs_wide[patient_obs_wide$patient_num %in% demo_raw_male$patient_num,]
  Male_patient_vec <- unique(Male_patient$patient_num)
  
  #####
  # start analysis


  demo_processed <- demo_raw  %>%
    mutate(
      across(
        ends_with("_date") & tidywhere(is.character),
        ~ lubridate::parse_date_time(.x, orders = c("mdy", "ymd"))
      ),
      last_discharge_date = pmin(death_date, last_discharge_date, na.rm = TRUE),
      total_stay = lubridate::interval(admission_date, last_discharge_date) %/% lubridate::days(1)
    ) %>%
    mutate(
      time_to_severe = severe_date - admission_date,
      time_to_severe = ifelse(time_to_severe < 0, NA, time_to_severe),
      time_to_death = death_date - admission_date,
      time_to_death = ifelse(time_to_death < 0, NA, time_to_death),
      sex = as.factor(sex),
      race = as.factor(race),
      age_group = as.factor(age_group),
      Severity = as.factor(severe) %>%
        fct_recode(Severe = "1", `Non-severe` = "0"),
      Survival = as.factor(deceased) %>%
        fct_recode(Alive = "0", Deceased = "1"),
      n_stay = total_stay,
      thrombo = patient_num %in% thrombo_patient_vec,
      neuro = patient_num %in% neuro_patient_vec,
      ARDs = patient_num %in% ARDs_patient_vec,
      Severe = patient_num %in% Severe_patient_vec,
      Sex2 = patient_num %in% Male_patient_vec
    ) 
  
  
  CurrSiteId = siteid
  
  obfus_tables_thrombo <- get_tables(
    demo_processed,
    blur_abs = 0,
    mask_thres = 10,
    group_var = 'thrombo'
  ) %>%
    lapply(function(x) mutate(x, site = CurrSiteId))
  
  obfus_tables_neuro <- get_tables(
    demo_processed,
    blur_abs = 0,
    mask_thres = 10,
    group_var = 'neuro'
  ) %>%
    lapply(function(x) mutate(x, site = CurrSiteId))
  
  obfus_tables_ARDs <- get_tables(
    demo_processed,
    blur_abs = 0,
    mask_thres = 10,
    group_var = 'ARDs'
  ) %>%
    lapply(function(x) mutate(x, site = CurrSiteId))
  
  obfus_tables_Severe <- get_tables(
    demo_processed,
    blur_abs = 0,
    mask_thres = 10,
    group_var = 'Severe'
  ) %>%
    lapply(function(x) mutate(x, site = CurrSiteId))
  
  obfus_tables_Sex <- get_tables(
    demo_processed,
    blur_abs = 0,
    mask_thres = 10,
    group_var = 'Sex2'
  ) %>%
    lapply(function(x) mutate(x, site = CurrSiteId))
  
  
  combined_table_wide <- 
    obfus_tables_thrombo$demo_table %>% 
    mutate(
      site = toupper(site),
      variable = tolower(variable),
      Demo_var = sub("\\..*", "", variable),
      Demo_var_i = sub(".*\\.", "", variable) %>%
        gsub("_", " ", .) %>%
        str_to_title() %>%
        recode(
          `00to02` = "0-2",
          `03to05` = "3-5",
          `06to11` = "6-11",
          `12to17` = "12-17",
          `18to25` = "18-25",
          `26to49` = "26-49",
          `50to69` = "50-69",
          `70to79` = "70-79",
          `80plus` = "80+",
          race.white = "White",
          race.american_indian = "American Indian",
          `Hawaiian Pacific Islander` = "Hawaiian/Pacific Islander",
          `Hispanic Latino` = "Hispanic/Latino"
        )
    )
  
  
  combined_table_wide_thrombo <- 
    obfus_tables_thrombo$demo_table %>% 
    mutate(
      site = toupper(site),
      variable = tolower(variable),
      Demo_var = sub("\\..*", "", variable),
      Demo_var_i = sub(".*\\.", "", variable) %>%
        gsub("_", " ", .) %>%
        str_to_title() %>%
        recode(
          `00to02` = "0-2",
          `03to05` = "3-5",
          `06to11` = "6-11",
          `12to17` = "12-17",
          `18to25` = "18-25",
          `26to49` = "26-49",
          `50to69` = "50-69",
          `70to79` = "70-79",
          `80plus` = "80+",
          race.white = "White",
          race.american_indian = "American Indian",
          `Hawaiian Pacific Islander` = "Hawaiian/Pacific Islander",
          `Hispanic Latino` = "Hispanic/Latino"
        )
    )
  
  combined_table_wide_neuro <- 
    obfus_tables_neuro$demo_table %>% 
    mutate(
      site = toupper(site),
      variable = tolower(variable),
      Demo_var = sub("\\..*", "", variable),
      Demo_var_i = sub(".*\\.", "", variable) %>%
        gsub("_", " ", .) %>%
        str_to_title() %>%
        recode(
          `00to02` = "0-2",
          `03to05` = "3-5",
          `06to11` = "6-11",
          `12to17` = "12-17",
          `18to25` = "18-25",
          `26to49` = "26-49",
          `50to69` = "50-69",
          `70to79` = "70-79",
          `80plus` = "80+",
          race.white = "White",
          race.american_indian = "American Indian",
          `Hawaiian Pacific Islander` = "Hawaiian/Pacific Islander",
          `Hispanic Latino` = "Hispanic/Latino"
        )
    )
  
  combined_table_wide_ARDs <- 
    obfus_tables_ARDs$demo_table %>% 
    mutate(
      site = toupper(site),
      variable = tolower(variable),
      Demo_var = sub("\\..*", "", variable),
      Demo_var_i = sub(".*\\.", "", variable) %>%
        gsub("_", " ", .) %>%
        str_to_title() %>%
        recode(
          `00to02` = "0-2",
          `03to05` = "3-5",
          `06to11` = "6-11",
          `12to17` = "12-17",
          `18to25` = "18-25",
          `26to49` = "26-49",
          `50to69` = "50-69",
          `70to79` = "70-79",
          `80plus` = "80+",
          race.white = "White",
          race.american_indian = "American Indian",
          `Hawaiian Pacific Islander` = "Hawaiian/Pacific Islander",
          `Hispanic Latino` = "Hispanic/Latino"
        )
    )
  
  
  combined_table_wide_Severe <- 
    obfus_tables_Severe$demo_table %>% 
    mutate(
      site = toupper(site),
      variable = tolower(variable),
      Demo_var = sub("\\..*", "", variable),
      Demo_var_i = sub(".*\\.", "", variable) %>%
        gsub("_", " ", .) %>%
        str_to_title() %>%
        recode(
          `00to02` = "0-2",
          `03to05` = "3-5",
          `06to11` = "6-11",
          `12to17` = "12-17",
          `18to25` = "18-25",
          `26to49` = "26-49",
          `50to69` = "50-69",
          `70to79` = "70-79",
          `80plus` = "80+",
          race.white = "White",
          race.american_indian = "American Indian",
          `Hawaiian Pacific Islander` = "Hawaiian/Pacific Islander",
          `Hispanic Latino` = "Hispanic/Latino"
        )
    )
  
  combined_table_wide_Sex <- 
    obfus_tables_Sex$demo_table %>% 
    mutate(
      site = toupper(site),
      variable = tolower(variable),
      Demo_var = sub("\\..*", "", variable),
      Demo_var_i = sub(".*\\.", "", variable) %>%
        gsub("_", " ", .) %>%
        str_to_title() %>%
        recode(
          `00to02` = "0-2",
          `03to05` = "3-5",
          `06to11` = "6-11",
          `12to17` = "12-17",
          `18to25` = "18-25",
          `26to49` = "26-49",
          `50to69` = "50-69",
          `70to79` = "70-79",
          `80plus` = "80+",
          race.white = "White",
          race.american_indian = "American Indian",
          `Hawaiian Pacific Islander` = "Hawaiian/Pacific Islander",
          `Hispanic Latino` = "Hispanic/Latino"
        )
    )
  
  
  
  
  ordered_vars <- c(
    "all", "sex", "age_group", "race",
    "readmitted", "survival"
  )
  
  tableOne_inter<- combined_table_wide %>%
    group_by(variable, Demo_var, Demo_var_i) %>%
    summarise(across(
      starts_with("n_var"),
      function(x) sum(x, na.rm = TRUE)
    ), .groups = "drop")
  
  tableOne_inter_thrombo <- combined_table_wide_thrombo %>%
    group_by(variable, Demo_var, Demo_var_i) %>%
    summarise(across(
      starts_with("n_var"),
      function(x) sum(x, na.rm = TRUE)
    ), .groups = "drop")
  
  tableOne_inter_neuro <- combined_table_wide_neuro %>%
    group_by(variable, Demo_var, Demo_var_i) %>%
    summarise(across(
      starts_with("n_var"),
      function(x) sum(x, na.rm = TRUE)
    ), .groups = "drop")
  
  tableOne_inter_ARDs <- combined_table_wide_ARDs %>%
    group_by(variable, Demo_var, Demo_var_i) %>%
    summarise(across(
      starts_with("n_var"),
      function(x) sum(x, na.rm = TRUE)
    ), .groups = "drop")
  
  tableOne_inter_Severe <- combined_table_wide_Severe %>%
    group_by(variable, Demo_var, Demo_var_i) %>%
    summarise(across(
      starts_with("n_var"),
      function(x) sum(x, na.rm = TRUE)
    ), .groups = "drop")
  
  tableOne_inter_Sex <- combined_table_wide_Sex %>%
    group_by(variable, Demo_var, Demo_var_i) %>%
    summarise(across(
      starts_with("n_var"),
      function(x) sum(x, na.rm = TRUE)
    ), .groups = "drop")
  
  row_order <- c(
    "All Patients", "Female", "Male",
    "0-2", "3-5", "6-11", "12-17", "18-25",
    "26-49", "50-69", "70-79", "80+",
    "American Indian", "Asian", "Black",
    "Hawaiian/Pacific Islander",
    "Hispanic/Latino", "Other", "White",
    "Median Elixhauser score [Min, Max]",
    "Median length of stay [Min, Max]",
    "Alive", "Deceased", "Median time to death [Min, Max]"
  )
  
  get_table_n <- function(var) {
    tableOne_raw %>%
      filter(`Table 1` == var) %>%
      pull(N)
  }
  
  tableOne_raw <- tableOne_inter %>%
    filter(Demo_var == "sex") %>%
    summarise(across(where(is.numeric), sum)) %>%
    data.frame(
      variable = "all",
      Demo_var = "all",
      Demo_var_i = "All Patients",
      .
    ) %>%
    bind_rows(tableOne_inter) %>%
    mutate(
      N = rowSums(across(where(is.numeric))),
      Demo_var = factor(Demo_var, levels = ordered_vars)
    ) %>%
    group_by(Demo_var) %>%
    mutate(across(
      starts_with("n_var"),
      function(x) {
        paste0(x, " (", round(x / sum(x, na.rm = TRUE) * 100, 1), "%", ")")
      }
    )) %>%
    ungroup() %>%
    arrange(Demo_var) %>%
    dplyr::select(Demo_var_i, N, starts_with("n_var")) %>%
    `colnames<-`(gsub(
      x = names(.),
      pattern = "n_var_",
      replacement = ""
    )) %>%
    rename("Table 1" = Demo_var_i) 
  
  
  tableOne_raw_thrombo <- tableOne_inter_thrombo %>%
    filter(Demo_var == "sex") %>%
    summarise(across(where(is.numeric), sum)) %>%
    data.frame(
      variable = "all",
      Demo_var = "all",
      Demo_var_i = "All Patients",
      .
    ) %>%
    bind_rows(tableOne_inter_thrombo) %>%
    mutate(
      N = rowSums(across(where(is.numeric))),
      Demo_var = factor(Demo_var, levels = ordered_vars)
    ) %>%
    group_by(Demo_var) %>%
    mutate(across(
      starts_with("n_var"),
      function(x) {
        paste0(x, " (", round(x / sum(x, na.rm = TRUE) * 100, 1), "%", ")")
      }
    )) %>%
    ungroup() %>%
    arrange(Demo_var) %>%
    dplyr::select(Demo_var_i, N, starts_with("n_var")) %>%
    `colnames<-`(gsub(
      x = names(.),
      pattern = "n_var_",
      replacement = ""
    )) %>%
    rename("Table 1" = Demo_var_i) 
  
  tableOne_raw_neuro <- tableOne_inter_neuro %>%
    filter(Demo_var == "sex") %>%
    summarise(across(where(is.numeric), sum)) %>%
    data.frame(
      variable = "all",
      Demo_var = "all",
      Demo_var_i = "All Patients",
      .
    ) %>%
    bind_rows(tableOne_inter_neuro) %>%
    mutate(
      N = rowSums(across(where(is.numeric))),
      Demo_var = factor(Demo_var, levels = ordered_vars)
    ) %>%
    group_by(Demo_var) %>%
    mutate(across(
      starts_with("n_var"),
      function(x) {
        paste0(x, " (", round(x / sum(x, na.rm = TRUE) * 100, 1), "%", ")")
      }
    )) %>%
    ungroup() %>%
    arrange(Demo_var) %>%
    dplyr::select(Demo_var_i, N, starts_with("n_var")) %>%
    `colnames<-`(gsub(
      x = names(.),
      pattern = "n_var_",
      replacement = ""
    )) %>%
    rename("Table 1" = Demo_var_i) 
  
  tableOne_raw_ARDs <- tableOne_inter_ARDs %>%
    filter(Demo_var == "sex") %>%
    summarise(across(where(is.numeric), sum)) %>%
    data.frame(
      variable = "all",
      Demo_var = "all",
      Demo_var_i = "All Patients",
      .
    ) %>%
    bind_rows(tableOne_inter_ARDs) %>%
    mutate(
      N = rowSums(across(where(is.numeric))),
      Demo_var = factor(Demo_var, levels = ordered_vars)
    ) %>%
    group_by(Demo_var) %>%
    mutate(across(
      starts_with("n_var"),
      function(x) {
        paste0(x, " (", round(x / sum(x, na.rm = TRUE) * 100, 1), "%", ")")
      }
    )) %>%
    ungroup() %>%
    arrange(Demo_var) %>%
    dplyr::select(Demo_var_i, N, starts_with("n_var")) %>%
    `colnames<-`(gsub(
      x = names(.),
      pattern = "n_var_",
      replacement = ""
    )) %>%
    rename("Table 1" = Demo_var_i) 
  
  tableOne_raw_Severe <- tableOne_inter_Severe %>%
    filter(Demo_var == "sex") %>%
    summarise(across(where(is.numeric), sum)) %>%
    data.frame(
      variable = "all",
      Demo_var = "all",
      Demo_var_i = "All Patients",
      .
    ) %>%
    bind_rows(tableOne_inter_Severe) %>%
    mutate(
      N = rowSums(across(where(is.numeric))),
      Demo_var = factor(Demo_var, levels = ordered_vars)
    ) %>%
    group_by(Demo_var) %>%
    mutate(across(
      starts_with("n_var"),
      function(x) {
        paste0(x, " (", round(x / sum(x, na.rm = TRUE) * 100, 1), "%", ")")
      }
    )) %>%
    ungroup() %>%
    arrange(Demo_var) %>%
    dplyr::select(Demo_var_i, N, starts_with("n_var")) %>%
    `colnames<-`(gsub(
      x = names(.),
      pattern = "n_var_",
      replacement = ""
    )) %>%
    rename("Table 1" = Demo_var_i) 
  
  
  tableOne_raw_Sex <- tableOne_inter_Sex %>%
    filter(Demo_var == "sex") %>%
    summarise(across(where(is.numeric), sum)) %>%
    data.frame(
      variable = "all",
      Demo_var = "all",
      Demo_var_i = "All Patients",
      .
    ) %>%
    bind_rows(tableOne_inter_Sex) %>%
    mutate(
      N = rowSums(across(where(is.numeric))),
      Demo_var = factor(Demo_var, levels = ordered_vars)
    ) %>%
    group_by(Demo_var) %>%
    mutate(across(
      starts_with("n_var"),
      function(x) {
        paste0(x, " (", round(x / sum(x, na.rm = TRUE) * 100, 1), "%", ")")
      }
    )) %>%
    ungroup() %>%
    arrange(Demo_var) %>%
    dplyr::select(Demo_var_i, N, starts_with("n_var")) %>%
    `colnames<-`(gsub(
      x = names(.),
      pattern = "n_var_",
      replacement = ""
    )) %>%
    rename("Table 1" = Demo_var_i) 
  
  
  
  tableOne_compiled_thrombo <- tableOne_raw_thrombo %>%
    mutate(
      N = case_when(
        grepl("Median length of stay", `Table 1`) ~ get_table_n("All Patients"),
        grepl("death", `Table 1`) ~ get_table_n("Deceased"),
        grepl("Median Elixhauser score", `Table 1`) ~ get_table_n("All Patients"),
        TRUE ~ N
      )
    ) %>% 
    rename('No thrombotic events' = `FALSE`,
           'Has thrombotic events' = `TRUE`)
  
  tableOne_compiled_neuro<- tableOne_raw_neuro %>%
    mutate(
      N = case_when(
        grepl("Median length of stay", `Table 1`) ~ get_table_n("All Patients"),
        grepl("death", `Table 1`) ~ get_table_n("Deceased"),
        grepl("Median Elixhauser score", `Table 1`) ~ get_table_n("All Patients"),
        TRUE ~ N
      )
    ) %>% 
    rename('No neurological events' = `FALSE`,
           'Has neurological events' = `TRUE`)
  
  tableOne_compiled_ARDs <- tableOne_raw_ARDs %>%
    mutate(
      N = case_when(
        grepl("Median length of stay", `Table 1`) ~ get_table_n("All Patients"),
        grepl("death", `Table 1`) ~ get_table_n("Deceased"),
        grepl("Median Elixhauser score", `Table 1`) ~ get_table_n("All Patients"),
        TRUE ~ N
      )
    ) %>% 
    rename('No ARDs events' = `FALSE`,
           'Has ARDs events' = `TRUE`)
  
  tableOne_compiled_Severe <- tableOne_raw_Severe %>%
    mutate(
      N = case_when(
        grepl("Median length of stay", `Table 1`) ~ get_table_n("All Patients"),
        grepl("death", `Table 1`) ~ get_table_n("Deceased"),
        grepl("Median Elixhauser score", `Table 1`) ~ get_table_n("All Patients"),
        TRUE ~ N
      )
    ) %>% 
    rename('Non-Severe' = `FALSE`,
           'Severe' = `TRUE`)
  
  tableOne_compiled_Sex <- tableOne_raw_Sex %>%
    mutate(
      N = case_when(
        grepl("Median length of stay", `Table 1`) ~ get_table_n("All Patients"),
        grepl("death", `Table 1`) ~ get_table_n("Deceased"),
        grepl("Median Elixhauser score", `Table 1`) ~ get_table_n("All Patients"),
        TRUE ~ N
      )
    ) %>% 
    rename('Female' = `FALSE`,
           'Male' = `TRUE`)
  
  tableOne_compiled_two = cbind(tableOne_compiled_thrombo,tableOne_compiled_neuro[,3:4])
  
  tableOne_compiled_all = cbind(tableOne_compiled_two,tableOne_compiled_ARDs[,3:4],tableOne_compiled_Severe[,3:4],tableOne_compiled_Sex[,3:4])
  
  tableOne_compiled_all= tableOne_compiled_all[1:(nrow(tableOne_compiled_all)-2),]
  
  
  
  if (time == "phase_1") {
    save(df_prop,df_num,df_prop2,df_num2,df_prop_Sex,df_prop_Sex2,df_prop_Severity,df_prop_Severity2,num_missing,prop_missing,num_missing2,prop_missing2,temporal,temporal_severe,temporal_nonsevere,longhaul_comb3,longhaul3,longhaul_comb7,longhaul7,shorthaul_comb,shorthaul_comb2,shorthaul,props_quant,props_quant2,props_M_quant,props_M_quant2,props_F_quant,props_F_quant2,props_S_quant,props_S_quant2,props_NS_quant,props_NS_quant2,props_S_quant3,props_S_quant4,props_NS_quant3,props_NS_quant4,props_S_quant5,props_NS_quant5,props_S_quant6,props_NS_quant6,props_S_quant7,props_NS_quant7,
         topic_diagnostics,topic_diagnostics2,TE_results,Severity_results,Neuro_results,ARDs_results,TE_results2,Severity_results2,Neuro_results2,ARDs_results2,TE_results3,Severity_results3,Neuro_results3,ARDs_results3,TE_results4,Severity_results4,Neuro_results4,ARDs_results4,TE_results_wilcox,Neuro_results_wilcox,Severity_results_wilcox,ARDs_results_wilcox,TE_results_wilcox2,Neuro_results_wilcox2,Severity_results_wilcox2,ARDs_results_wilcox2,TE_results_wilcox3,Neuro_results_wilcox3,Severity_results_wilcox3,ARDs_results_wilcox3,TE_results_wilcox4,Neuro_results_wilcox4,Severity_results_wilcox4,ARDs_results_wilcox4,beta_mat,beta_mat2,beta_mat3,beta_mat4,tableOne_compiled_all,
         file = paste(paste(getProjectOutputDirectory(),"/",sep=""),paste(siteid,"results_phase1.RData",sep="_"),sep=""))
  }
  
  if (time == "phase_2") {
    save(df_prop,df_num,df_prop2,df_num2,df_prop_Sex,df_prop_Sex2,df_prop_Severity,df_prop_Severity2,num_missing,prop_missing,num_missing2,prop_missing2,temporal,temporal_severe,temporal_nonsevere,longhaul_comb3,longhaul3,longhaul_comb7,longhaul7,shorthaul_comb,shorthaul_comb2,shorthaul,props_quant,props_quant2,props_M_quant,props_M_quant2,props_F_quant,props_F_quant2,props_S_quant,props_S_quant2,props_NS_quant,props_NS_quant2,props_S_quant3,props_S_quant4,props_NS_quant3,props_NS_quant4,props_S_quant5,props_NS_quant5,props_S_quant6,props_NS_quant6,props_S_quant7,props_NS_quant7,
         topic_diagnostics,topic_diagnostics2,TE_results,Severity_results,Neuro_results,ARDs_results,TE_results2,Severity_results2,Neuro_results2,ARDs_results2,TE_results3,Severity_results3,Neuro_results3,ARDs_results3,TE_results4,Severity_results4,Neuro_results4,ARDs_results4,TE_results_wilcox,Neuro_results_wilcox,Severity_results_wilcox,ARDs_results_wilcox,TE_results_wilcox2,Neuro_results_wilcox2,Severity_results_wilcox2,ARDs_results_wilcox2,TE_results_wilcox3,Neuro_results_wilcox3,Severity_results_wilcox3,ARDs_results_wilcox3,TE_results_wilcox4,Neuro_results_wilcox4,Severity_results_wilcox4,ARDs_results_wilcox4,beta_mat,beta_mat2,beta_mat3,beta_mat4,tableOne_compiled_all,
         file = paste(paste(getProjectOutputDirectory(),"/",sep=""),paste(siteid,"results_phase2.RData",sep="_"),sep=""))
  }
  
  if (time == "all") {
    save(df_prop,df_num,df_prop2,df_num2,df_prop_Sex,df_prop_Sex2,df_prop_Severity,df_prop_Severity2,num_missing,prop_missing,num_missing2,prop_missing2,temporal,temporal_severe,temporal_nonsevere,longhaul_comb3,longhaul3,longhaul_comb7,longhaul7,shorthaul_comb,shorthaul_comb2,shorthaul,props_quant,props_quant2,props_M_quant,props_M_quant2,props_F_quant,props_F_quant2,props_S_quant,props_S_quant2,props_NS_quant,props_NS_quant2,props_S_quant3,props_S_quant4,props_NS_quant3,props_NS_quant4,props_S_quant5,props_NS_quant5,props_S_quant6,props_NS_quant6,props_S_quant7,props_NS_quant7,
         topic_diagnostics,topic_diagnostics2,TE_results,Severity_results,Neuro_results,ARDs_results,TE_results2,Severity_results2,Neuro_results2,ARDs_results2,TE_results3,Severity_results3,Neuro_results3,ARDs_results3,TE_results4,Severity_results4,Neuro_results4,ARDs_results4,TE_results_wilcox,Neuro_results_wilcox,Severity_results_wilcox,ARDs_results_wilcox,TE_results_wilcox2,Neuro_results_wilcox2,Severity_results_wilcox2,ARDs_results_wilcox2,TE_results_wilcox3,Neuro_results_wilcox3,Severity_results_wilcox3,ARDs_results_wilcox3,TE_results_wilcox4,Neuro_results_wilcox4,Severity_results_wilcox4,ARDs_results_wilcox4,beta_mat,beta_mat2,beta_mat3,beta_mat4,tableOne_compiled_all,
         file = paste(paste(getProjectOutputDirectory(),"/",sep=""),paste(siteid,"results_all.RData",sep="_"),sep=""))
   }
  
}
     
    
    
    



