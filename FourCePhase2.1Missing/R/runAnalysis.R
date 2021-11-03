
#' Runs the analytic workflow for the Missing project
#'
#' @keywords 4CE
#' @export

runAnalysis <- function(data_dir = "/4ceData/Input",dateFormat="%d-%b-%y") {
    
    for (r_file in list.files('/covidclinical/Phase2.1MissingRPackage/helpers', full.names = TRUE, pattern = '.R$')) source(r_file)

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

    ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that 
    ## is mounted to the container
    currSiteId = FourCePhase2.1Data::getSiteId()

    ## run the quality control
    FourCePhase2.1Data::runQC(currSiteId)

    ## DO NOT CHANGE ANYTHING ABOVE THIS LINE
    
    data_dir <- FourCePhase2.1Data::getInputDataDirectoryName()
    currSiteId <- FourCePhase2.1Data::getSiteId()

    #### Quantify missingness across sites
    
    #Read in local patient summary: 
    demo_raw <-
      readr::read_csv(
        file.path(data_dir, "LocalPatientSummary.csv"),
        col_types = list(patient_num = readr::col_character()),
        na = "1900-01-01"
      ) %>%
      mutate(
        across(ends_with("_date") & where(is.character), lubridate::mdy),     
        last_discharge_date = if_else(
          !is.na(death_date) & death_date < last_discharge_date,
          death_date,
          last_discharge_date
        )
      )
    #Read in thrombotic icd codes. Add column "truncated code" for first three characters of each code. Extract unique truncated codes. 
    thrombo_codes <- read_csv("https://raw.githubusercontent.com/covidclinical/Phase2.1AKIRPackage/ac19716a4586f45c398728fcd821ca9d5baffe45/FourCePhase2.1AKI/data-raw/thromb_icd_code.csv") %>%
      mutate(truncated_code = substr(icd_code, 1, 3)) %>%
      pull(truncated_code) %>%
      unique()
    #Extract information for patients with codes for thrombotic events
    patient_obs <- read_csv(
      file.path(data_dir, "LocalPatientObservations.csv"),
      col_types = list(patient_num = readr::col_character()),
      na = "-999"
    )
    te_patients <- patient_obs %>%
      filter(
        concept_type == "DIAG-ICD10",
        concept_code %in% thrombo_codes,
        days_since_admission >= 0
      )
    
    # 4CE long format Labs Data
    # Code, Descriptions and Ranges
    lab_mapping <- read_csv("public-data/loinc-map.csv")
    lab_bounds <- read_csv("public-data/lab_bounds.csv")
    lab_names <- lab_bounds$short_name
    # load('public-data/code.dict.rda')
    
    # Create consistent names for D_Dimer and Troponin
    
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
    
    ### Calculate proportions and numbers missing across all sites
    
    avgs = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    props = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    ct = 1
    for (i in 1:length(unique(patient_obs_wide$patient_num))) {
      for (j in 1:ncol(avgs)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,length(dat_prop))
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


    #first 20 days: avg prop missing and avg num missing

    avgs = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    props = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    ct = 1
    for (i in 1:20) {
      for (j in 1:ncol(avgs)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,length(dat_prop))
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

    patient_obs_first20 = patient_obs_wide[patient_obs_wide$days_since_admission <= 20,]
    num_missing2 = sapply(patient_obs_first20, function(x) sum(is.na(x)))[3:18]
    prop_missing2 = sapply(patient_obs_first20, function(x) sum(is.na(x))/length(x))[3:18]
                           
    ### Do this for sex and severity
                           
    #Stratify by sex!
    props_M = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    props_F = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    indices_M = which(demo_raw$sex=="male")
    indices_F = which(demo_raw$sex=="female")
    ct = 1
    for (i in 1:length(as.integer(unique(patient_obs_wide$patient_num)))) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_M){

      for (j in 1:ncol(props_M)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    for (i in 1:length(as.integer(unique(patient_obs_wide$patient_num)))) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_F){
      for (j in 1:ncol(props_F)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,nrow(df_prop_Sex))
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
    Site = rep(params$site,nrow(df_prop_diff_sex))
    df_prop_diff_sex$site = Site
    Site = rep(params$site,nrow(df_prop_diff_sex))
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
    indices_S = which(demo_raw$severe==1)
    indices_NS = which(demo_raw$severe==0)
    ct = 1
    for (i in 1:length(as.integer(unique(patient_obs_wide$patient_num)))) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_S){

      for (j in 1:ncol(props_S)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    for (i in 1:length(as.integer(unique(patient_obs_wide$patient_num)))) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_NS){
      for (j in 1:ncol(props_NS)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,nrow(df_prop_Severity))
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
    Site = rep(params$site,nrow(df_prop_diff_severity))
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











    #Stratify by sex-- first 20 days!

    props_M = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    props_F = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    indices_M = which(demo_raw$sex=="male")
    indices_F = which(demo_raw$sex=="female")
    ct = 1
    for (i in 1:20) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_M){

      for (j in 1:ncol(props_M)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    for (i in 1:20) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_F){
      for (j in 1:ncol(props_F)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,nrow(df_prop_Sex2))
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
    Site = rep(params$site,nrow(df_prop_diff_sex2))
    df_prop_diff_sex2$site = Site
    Site = rep(params$site,nrow(df_prop_diff_sex2))
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
    indices_S = which(demo_raw$severe==1)
    indices_NS = which(demo_raw$severe==0)
    ct = 1
    for (i in 1:20) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_S){

      for (j in 1:ncol(props_S)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    for (i in 1:20) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_NS){
      for (j in 1:ncol(props_NS)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,nrow(df_prop_Severity2))
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
    Site = rep(params$site,nrow(df_prop_diff_severity2))
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
    indices_S = which(demo_raw$severe==1)
    indices_NS = which(demo_raw$severe==0)
    indices_age = which(demo_raw$age_group != "00to02" & demo_raw$age_group != "03to05" & demo_raw$age_group != "06to11" & demo_raw$age_group != "12to17")

    ct = 1
    for (i in 1:length(as.integer(unique(patient_obs_wide$patient_num)))) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_S){
        if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_age) {

      for (j in 1:ncol(props_S)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    for (i in 1:length(as.integer(unique(patient_obs_wide$patient_num)))) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_NS){
        if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_age) {
      for (j in 1:ncol(props_NS)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,nrow(df_prop_Severity3))
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
    Site = rep(params$site,nrow(df_prop_diff_severity3))
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






    #Stratify by patients over 21 (severity) in first 20 days!

    props_S = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    props_NS = patient_obs_wide[FALSE,] %>% dplyr::select(-c(patient_num,days_since_admission,severity))
    indices_S = which(demo_raw$severe==1)
    indices_NS = which(demo_raw$severe==0)
    indices_age = which(demo_raw$age_group != "00to02" & demo_raw$age_group != "03to05" & demo_raw$age_group != "06to11" & demo_raw$age_group != "12to17")

    ct = 1
    for (i in 1:20) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_S){
        if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_age) {

      for (j in 1:ncol(props_S)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    for (i in 1:20) {
      if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_NS){
        if (as.integer(unique(patient_obs_wide$patient_num))[i] %in% indices_age) {
      for (j in 1:ncol(props_NS)) {
        dat <- filter(patient_obs_wide,patient_num==as.integer(demo_raw$patient_num[i]))
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
    site = rep(params$site,nrow(df_prop_Severity4))
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
    Site = rep(params$site,nrow(df_prop_diff_severity4))
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




    ### Temporal pair analysis
                           
    demo_raw$date_diff <- as.Date(as.character(demo_raw$last_discharge_date), format=dateFormat)-
                  as.Date(as.character(demo_raw$admission_date), format=dateFormat)
    dates_diff = c()
    for (i in 1:nrow(demo_raw)){
      dates_diff = c(dates_diff,rep(demo_raw$date_diff[i],sum(patient_obs_wide$patient_num==as.integer(demo_raw$patient_num[i]))))
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

          if (mean(abs(test_vec) > min_corr) > prop_greater) {
            count = count+1
            final_names[count] = name1
            df[,count] = test_vec
            df_ind[,count] = indices
          }
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

          if (mean(abs(test_vec) > min_corr) > prop_greater) {
            count = count+1
            final_names[count] = name1
            df[,count] = test_vec
            df_ind[,count] = indices
          }
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

          if (mean(abs(test_vec) > min_corr) > prop_greater) {
            count = count+1
            final_names[count] = name1
            df[,count] = test_vec
            df_ind[,count] = indices
          }
        }
        }
    }

    colnames(df) <- final_names
    results <- list()
    results[[1]] <- df
    results[[2]] <- mat_list
    temporal_nonsevere = df
                           
                           
    
                           
    ### Temporal analysis of individual labs
                           
                           
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
    level_order = c("Fibrinogen","Procalcitonin","D_Dimer","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
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
    level_order = c("Fibrinogen","Procalcitonin","D_Dimer","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
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



















    #Calculate proportion of existing lab values in days since admission up to 20 days. Stratify by severe vs. non-severe.
    binby = 1
    wide =  patient_obs_wide %>% 
        filter(days_since_admission <= 20,
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
    days = seq(1,20,1)
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
    level_order = c("Fibrinogen","Procalcitonin","DDU","Troponin","CRP","Ferritin","LDH","Lymphocyte","PT","Albumin","Neutrophil","AST","ALT","Bilirubin","Leukocytes","Creatinine")
    days = seq(1,20,1)
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
                           
    site_results <- paste0(CurrSiteId, "_results")
      assign(site_results, results)
      save(df_prop,df_num,df_prop2,df_num2,df_prop_Sex,df_prop_Sex2,df_prop_Severity,df_prop_Severity2,num_missing,prop_missing,num_missing2,prop_missing2,temporal,temporal_severe,temporal_nonsevere,longhaul_comb3,longhaul3,longhaul_comb7,longhaul7,shorthaul_comb,shorthaul,props_quant,props_quant2,props_M_quant,props_M_quant2,props_F_quant,props_F_quant2,props_S_quant,props_S_quant2,props_NS_quant,props_NS_quant2,props_S_quant3,props_S_quant4,props_NS_quant3,props_NS_quant4,
        file = file.path(
          getProjectOutputDirectory(),
          paste0(CurrSiteId, "_results.rda")
        )
      )
                           
                           
    


































    
    
    
    

    ## To Do: implement analytic workflow, saving results to a site-specific 
    ## file to be sent to the coordinating site later via submitAnalysis()

    ## Save results to appropriately named files for submitAnalysis(), e.g.:
    #write.csv(
    #    matrix(rnorm(100), ncol=5), 
    #    file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_ResultTable.csv"))
    #)

    #write.table(
    #    matrix(rnorm(12), ncol=3), 
    #    file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_ModelParameters.txt"))
    #)
    
}

