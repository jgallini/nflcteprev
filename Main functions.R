
#This program contains all custom functions used in the NFL CTE prevalence analysis

#Functions 
#sratio: calculates selection ratios
#run_analysis: main results run to update prevalences based on selection ratios
#sratio_nonstrat: used for some non-stratified sensitivity analyses


# This "sratio" function calculates the selection ratio estimates using age-stratified
# MCI to Dementia ratios taken from the NFL-LONG study
# It has the option to use the original ratios or a smoothed version of the ratios

# Inputs: data -> a data frame containing individual-level NDI data required to fit the Poisson model
# ratio_type -> should be set to either "ratio" or "smooth_ratio" depending on which MCI to dementia
  # ratios are desired (ratio is raw ratios calculated from NFL-LONG study, smooth_ratio is 
  # these ratios smoothed out using a natural cubic spline, this is the option used in the main analysis)
# its -> scalar number of iterations, 1000 used for main analyses
# lwrbnd_sens -> a logical value (FALSE or TRUE) for whether or not this analysis
  # is the lower bound sensitivity analysis
# sui -> a logical value (FALSE or TRUE) for whether or not this analysis includes suicide
  # in the selection ratio calculation

# Outputs:
# A 3 item list with items srdraw, srmod, and srconf
# 1. srdraw contains <its> selection ratios drawn from random normal distributions defined
# by the mean and variance of the selection ratio model estimate for that iteration (table 3)
# this basically just adds another element of randomness to the selection ratio distribution
# to ensure robust results
# 2. srmod contains the actual selection ratio model estimates from each iteration 
# (no extra layer of randomness here) (table 2)
# 3. srconf contains model confidence interval estimates for the selection ratio term

sratio<- function (data,ratio_type,its,lwrbnd_sens=FALSE,sui=FALSE){
  
  library(stringr)
  library(splines)
  set.seed(123) 
  
  data <- data %>% mutate(age_group = cut(
        Age,breaks = c(-Inf, 55, 60, 65, 70, 75, 80, 85, Inf),
        right  = FALSE,
        labels = c("55-","55-59","60-64","65-69","70-74","75-79","80-84","85+")))

  #mci to dementia ratios taken from NFL-LONG study
  ratio_tbl <- tibble::tribble(
  ~age_group, ~mci_prev, ~dem_prev,
  "55-", 21.9,  7.0,
  "55-59", 24.0,  4.5,
  "60-64", 23.5,  7.8,
  "65-69", 21.9, 10.1,
  "70-74", 18.0,  8.2,
  "75-79", 24.3, 12.1,
  "80-84", 48.7, 20.5,
  "85+",   33.3, 22.2) %>% 
    mutate(age_group = factor(age_group, levels = levels(data$age_group)))

  #calculating age-stratified ratios for MCI to dementia
  ratio_tbl$ratio<-ratio_tbl$mci_prev/ratio_tbl$dem_prev

  #calculating smoothed ratios
  #ratio_tbl$Agenum<-c(50,55,60,65,70,75,80,85)
  ratio_tbl$Agenum_mid<-c(52.5,57,62,67,72,77,82,87.5)
  smoothmod<-lm(ratio ~ ns(Agenum_mid, df = 3), data = ratio_tbl)
  ratio_tbl$smooth_ratio <- predict(smoothmod)

  #tabulating age-stratified dementia cause of death by donor status
  tab_df <- as.data.frame(
    table(Donor = data$Donor,
          DementiaCOD = data$DementiaCOD,
          age_group = data$age_group))
  #attaching mci to dementia ratio to each stratum
  tab_df <- tab_df %>%
    left_join(ratio_tbl %>% select(age_group, all_of(ratio_type)),
              by = "age_group")
  #calculating how many people need to move into symptom category because we're
  #assigning them MCI (within donor strata)
  tab_df <- tab_df %>%
    mutate(mci_to_move = round(ifelse(DementiaCOD == 1,
                                      Freq * .data[[ratio_type]],0)))
  
  #patch here since no dementia cause of death in donor group for 55-59
  #randomly selecting value from binomial distribution (20,5/70), changing this each iteration
  #will index into this vector during each iteration below to patch this issue
  if (min(data$Age,na.rm=TRUE)<60){
  mci5559_vals<-rbinom(its,with(tab_df,Freq[age_group=="55-59"&Donor=="Donor"&DementiaCOD==0]),
                       (with(tab_df,mci_to_move[age_group=="55-59"&Donor=="Non-Donor"&DementiaCOD==1])/
                          with(tab_df,Freq[age_group=="55-59"&Donor=="Non-Donor"&DementiaCOD==0])))}
  
  #people to move, write-around for 55-59 donor group to fix 0 cell issue
  if (min(data$Age,na.rm=TRUE)<60){
  targets <- tab_df %>%
    filter(DementiaCOD == 1, (mci_to_move > 0 | (age_group=="55-59"&Donor=="Donor"))) %>%
    select(age_group, Donor, mci_to_move)
  } else if (min(data$Age,na.rm=TRUE)>=60){
    targets <- tab_df %>%
      filter(DementiaCOD == 1, (mci_to_move > 0)) %>%
      select(age_group, Donor, mci_to_move)
  }
  
  #creating dementia or mci indicator that starts as just dementiaCOD status
  #for everyone, updates using loop below
  data$dem_mci<-data$DementiaCOD
  
  #loop that is reassigning people to symptom group randomly within strata
  #defined by age and donor status
  srmod<-list()
  srdraw<-vector()
  srconf<-list()
  
  for (j in 1:its){
    
    #keeping track of where the loop is
    if (j %% 10 == 0) {print(j)}
    
    #starting over with original data each time
    ndi_new <- data
    
    for (i in seq_len(nrow(targets))) {
      
      age_i   <- targets$age_group[i]
      donor_i <- targets$Donor[i]
      n_move  <- targets$mci_to_move[i]
      if (age_i=="55-59" & donor_i=="Donor"){
        n_move <- mci5559_vals[j]} 
      #overwrite with patch code (random binomial code above) to fix 0 cell issue
      
      # Eligible people: NOT dementia, NOT already MCI
      eligible <- which( ndi_new$age_group == age_i &
                           ndi_new$Donor == donor_i &
                           ndi_new$DementiaCOD == 0 & ndi_new$dem_mci == 0)
      
      if (length(eligible) < n_move) {
        if (j==1) {warning(paste(
          "Not enough eligible subjects in:",
          age_i,
          "Donor =", donor_i,
          "| needed:", n_move,
          "available:", length(eligible))
        )} #only printing these warnings once
        n_move <- length(eligible)
      }
      
      if (n_move > 0) {chosen <- sample(eligible, n_move)
      ndi_new$dem_mci[chosen] <- 1} 
    }
    #using 0 1 coding for poisson model
    donor_mod<-ifelse(data$Donor=="Donor",1,0)
    
    #overwriting dem_mci variable with true values in UNITE data for donors
    #as a lower-bound type of sensitivity analysis (DementiaHX and physdxmci)
    if (lwrbnd_sens==TRUE){
    ndi_new$dem_mci<-ifelse(ndi_new$Donor=="Non-Donor",ndi_new$dem_mci,
                        ifelse(ndi_new$Donor=="Donor"&(ndi_new$DementiaHx==1|ndi_new$phydxmci==1),1,
                            ifelse(ndi_new$Donor=="Donor"&(ndi_new$DementiaHx==0&ndi_new$phydxmci==0),0,NA)))
    }
    
    #overwriting dem_mci with including suicide in "symptoms" group for sensitivity analysis
    if (sui==TRUE){
      ndi_new$dem_mci<-ifelse(ndi_new$cod==1,1,ndi_new$dem_mci)
    }
      
    #now modeling to find selection ratio, poisson regression approximation
    srmod[[j]] <- glm(donor_mod ~ as.factor(dem_mci) + Age + as.factor(poscomb) + 
                        Yrs + as.factor(Race_bin) + PB + as.factor(hof), 
                      family = poisson(link = "log"), data = ndi_new)
    
    # #faster version of calculating confidence intervals, only really using
    #these for poisson table in manuscript of SR model
    b <- coef(srmod[[j]])
    se <- sqrt(diag(vcov(srmod[[j]])))
    srconf[[j]] <- exp(cbind(lower = b - 1.96 * se,
      upper = b + 1.96 * se ))
    
    #getting robust standard errors and effect estimates
    rob_vcov <- sandwich::vcovHC(srmod[[j]], type = "HC0")
    beta_hat <- coef(srmod[[j]])["as.factor(dem_mci)1"]
    se_hat   <- sqrt(diag(rob_vcov))["as.factor(dem_mci)1"]
    
    #drawing randomly from this distribution to allow for more randomness
    srdraw[j] <- exp(rnorm(n = 1, mean = beta_hat, sd = se_hat))
    
  }
  
  return(list(srdraw = srdraw,
       srmod  = srmod, srconf=srconf))
  
}

##############################################################################


#This "run_analysis" function runs main analyses (table 3) by calling the selection ratio function above,
#calculating CTE and dementia counts, and adjusting prevalence estimates

#Inputs: data_ndi-> individual-level NDI data used to fit the Poisson model
# ratio_type -> should be set to either "ratio" or "smooth_ratio" depending on which MCI to dementia
  # ratios are desired (ratio is raw ratios calculated from NFL-LONG study, smooth_ratio is 
  # these ratios smoothed out using a natural cubic spline, this is the option used in the main analysis)
# a_cte -> number of CTE positive brains in brain bank
# c_cte -> number of CTE negative brains in brain bank
# a_dem -> number of dementia cases in brain bank
# c_dem -> number of dementia negative cases in brain bank
# mciunder50 -> a logical value (FALSE or TRUE) for whether or not to include MCI cases
# under 50 years of age, used in sensitivity analyses
# lwrbnd_sens -> a logical value (FALSE or TRUE) for whether or not this analysis
  # is the lower bound sensitivity analysis
# sui -> a logical value (FALSE or TRUE) for whether or not this analysis includes suicide
  # in the selection ratio calculation

#Outputs: vector of randomly drawn selection ratios from distribution (#=its), appropriately adjusted
#CTE prevalence values (#=its), and appropriately adjusted dementia prevalence values (#=its)

run_analysis <- function(data_ndi,#data_unite,
                         ratio_type,
                         a_cte, c_cte, a_dem, c_dem, mciunder50=TRUE,lwrbnd_sens=FALSE,sui=FALSE) {
  
  #deciding which function to run based on whether we're using age-stratified
  #MCI to dementia ratios
  if (ratio_type %in% c("ratio","smooth_ratio")){
  sr_obj <- sratio(data = data_ndi, ratio_type = ratio_type,its = its,lwrbnd_sens=lwrbnd_sens,sui=sui)
  sr_draw <- sr_obj$srdraw
  } else if (ratio_type=="non_strat"){
    sr_draw <- sratio_nonstrat(data=data_ndi,mciunder50 = mciunder50,its=its)
  }
  
  cte_draw <- (a_cte / sr_draw) / ((a_cte / sr_draw) + c_cte)
  dem_draw <- (a_dem / sr_draw) / ((a_dem / sr_draw) + c_dem)
  
  c(`Selection Ratio` = sr_draw, 
    `CTE` = cte_draw,
    `Dementia` = dem_draw)
}

#############################################################################


#This "sratio_nonstart" function that calculates selection ratios using 1 MCI to 
# dementia ratio instead of age stratified ratios like in sratio. This is for sensitivity
#analyses and does not output the full model results like sratio does,
#just the selection ratio estimates

# Inputs: data -> individual level NDI data to fit Poisson model
# mciunder50 -> a logical value (FALSE or TRUE) for whether or not to include MCI cases
  # under 50 years of age, used in sensitivity analyses
#its -> number of iterations, usually 1000

# Outputs: srdraw: a vector of randomly drawn selection ratios from the appropriate 
# normal distributions defined by the model estimates for the selection ratio terms

sratio_nonstrat<-function(data,mciunder50=TRUE,its){
  
  #changing formatting of donor variable just within this function
  data$Donor<-ifelse(data$Donor=="Donor",1,
                     ifelse(data$Donor=="Non-Donor",0,NA))
  
  #overall mci to dementia ratio 
  mrat<-0.238/0.089
  #need to shift people from no dem group over to cognitive problems group
  #who have MCI, calculating total number
  mci_don<-round((mrat)*sum(data$Donor==1&data$DementiaCOD==1,na.rm=TRUE))
  mci_nondon<-round((mrat)*sum(data$Donor==0&data$DementiaCOD==1,na.rm=TRUE))

  #creating weights for each stratum based on literature for how many people 
  #in each age group have MCI, also taking from above NFL-LONG paper
  if (mciunder50==TRUE){
  stratum_weights <- data.frame(
    Age_cat = c("<50","50-54","55-59","60–64","65-69","70–74","75-79","80+"),
    weight  = c(0.219, 0.219, 0.24, 0.235, 0.219, 0.18, 0.243, 0.471))
  # MCI prevalence was 21.9% for ages 50-54 (n=128), 24% for 55-59(n=179), 
  #23.5% for ages 60–64 (n=166), 21.9% for 65–69 (n=169), 18% for 70–74 (n=122), 
  #24.3% for 75–79 (n=107), and 47.1% for 80+ (n=51). 
  #just using same as 50-54 for <50

  #number we need to sample for donors and non-donors
  alloc_don <- stratum_weights %>%
    mutate(n_stratum = round(weight / sum(weight) * mci_don))
  table(data$Age_cat[data$Donor==1 & data$DementiaCOD==0])
  #need to shift 2 people from 55-59 group to 60-64 group, patching this
  alloc_don$n_stratum[alloc_don$Age_cat=="55-59"]<-20
  alloc_don$n_stratum[alloc_don$Age_cat=="60–64"]<-24

  alloc_nondon <- stratum_weights %>%
    mutate(n_stratum = round(weight / sum(weight) * mci_nondon))
  table(data$Age_cat[data$Donor==0 & data$DementiaCOD==0])
  #need to shift 11 people from 50-54 group into <50 group, patching this
  alloc_nondon$n_stratum[alloc_nondon$Age_cat=="50-54"] <- 53
  alloc_nondon$n_stratum[alloc_nondon$Age_cat=="<50"] <- 75
  
  } else if (mciunder50==FALSE){
    stratum_weights <- data.frame(
      Age_cat = c("<50","50-54","55-59","60–64","65-69","70–74","75-79","80+"),
      weight  = c(0, 0.219, 0.24, 0.235, 0.219, 0.18, 0.243, 0.471))
    # MCI prevalence was 21.9% for ages 50-54 (n=128), 24% for 55-59(n=179), 
    #23.5% for ages 60–64 (n=166), 21.9% for 65–69 (n=169), 18% for 70–74 (n=122), 
    #24.3% for 75–79 (n=107), and 47.1% for 80+ (n=51). just using 0 for <50.
    
    #number we need to sample for donors and non-donors
    alloc_don <- stratum_weights %>%
      mutate(n_stratum = round(weight / sum(weight) * mci_don))
    table(data$Age_cat[data$Donor==1 & data$DementiaCOD==0])
    #need to slightly shift these original numbers because we don't have enough donors
    #in all age categories, moving people into closest age category with available
    #people left to sample from
    alloc_don$n_stratum<-c(0,20,20,27,27,19,30,44)
    alloc_nondon <- stratum_weights %>%
      mutate(n_stratum = round(weight / sum(weight) * mci_nondon))
    table(data$Age_cat[data$Donor==0 & data$DementiaCOD==0])
    alloc_nondon$n_stratum <- c(0,53,70,104,72,59,79,154)
    #less shifting needed for non-donors but still need shifting
  }

  #looping through this simulation 1,000 times to bootstrap estimates
  set.seed(92748)
  model_age<-list();meanest<-vector(); 
  bootest<-vector(); expboot<-vector();
  for (i in 1:its){
    
    #keeping track of where the loop is
    if (i %% 10 == 0) {print(i)}
    
    #defining new DementiaCOD or MCI variable (clear cognitive symptoms present variable)
    data$dem_mci<-data$DementiaCOD
    sampled_don <- data %>%
      mutate(data_index = row_number()) %>%   # original row index
      filter(Donor == 1,DementiaCOD==0) %>%
      inner_join(alloc_don, by = "Age_cat") %>%
      group_by(Age_cat) %>%
      group_modify(~ slice_sample(.x, n = unique(.x$n_stratum))) %>%
      ungroup() %>%
      pull(data_index)
    sampled_nondon <- data %>%
      mutate(data_index = row_number()) %>%   # original row index
      filter(Donor == 0,DementiaCOD==0) %>%
      inner_join(alloc_nondon, by = "Age_cat") %>%
      group_by(Age_cat) %>%
      group_modify(~ slice_sample(.x, n = unique(.x$n_stratum))) %>%
      ungroup() %>%
      pull(data_index)
    data$dem_mci[c(sampled_don,sampled_nondon)]<-1
  
    #now running age-adjusted model to get selection ratio, also adjusting for other
    #variables we have on all donors and non-donors
    model_age[[i]] <- glm(Donor ~ as.factor(dem_mci) + Age + as.factor(poscomb) + Yrs + as.factor(Race_bin) + PB + as.factor(hof), 
                        family = poisson(link = "log"), data = data)
  
    #getting robust standard errors
    rob_vcov <- sandwich::vcovHC(model_age[[i]], type = "HC0")
    se_hat   <- sqrt(diag(rob_vcov))["as.factor(dem_mci)1"]
    #mean estimate for selection ratio
    meanest[i]<-summary(model_age[[i]])$coefficients["as.factor(dem_mci)1", "Estimate"]
  
    #randomly drawing selection ratio estimate from distribution
    bootest[i]<-rnorm(1,mean=meanest[i],sd=se_hat)
    expboot[i]<-exp(bootest[i])
  
  }
  
return(srdraw = expboot)
  
}



