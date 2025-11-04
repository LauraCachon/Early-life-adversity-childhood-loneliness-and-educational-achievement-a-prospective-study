## SENSITIVITY ANALYSIS ##

#1. Sex-stratified analysis
#2. Including child's temperament as a confounder
#3. Alternative cut-off points for (Low) High School grades
#4. Complete case analysis 

#1. Sex stratified analysis
  #Almost identical to main analysis, but conducted across two datasets (boys/girls), and therefore not including the variable sex_child1 as a confounder

#Starting from the imp_10it dataset after imputation and data wrangling, create two separate datasets:
male <- filter(imp_10it, sex_child1 == 1)
female <- filter(imp_10it, sex_child1 == 2)


#Binomial models --> Eg: material hardship as predictor (group: boys)

# Create a list of 20 imputed datasets + add variables: male/female
imp_list_m <- lapply(1:20, function(i) {
  data_i <- complete(male_mids, i)
  data_i$mid_high <- ifelse(data_i$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
  data_i$high_only <- ifelse(data_i$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)
  data_i
})


# Main comparison function -males
compare_binomial_models_m <- function(imp_list_m, predictor) {
  # Define common covariates
  covariates <- c(
    "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Add predictor of interest
  formula_str <- paste("~", paste(c(covariates, predictor), collapse = " + "))
  
  # Model 1: mid_high ~ X
  fit_mid_high <- lapply(imp_list_m, function(data) {
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = data)
  })
  
  # Model 2: high_only ~ X
  fit_high_only <- lapply(imp_list_m, function(data) {
    glm(as.formula(paste("high_only", formula_str)), family = "binomial", data = data)
  })
  
  # Pool and extract results
  res_mid_high <- pool_manual(fit_mid_high) #pool_manual function created before (see main analysis)
  res_high_only <- pool_manual(fit_high_only)
  
  # Compare OR for predictor only 
  compa_OR <- res_mid_high %>%
    dplyr::filter(term == predictor) %>%
    dplyr::select(term, or_mh = Odds_Ratio,
                  or_mh_low = OR_Lower_CI,
                  or_mh_high = OR_Upper_CI) %>%
    left_join(
      res_high_only %>%
        dplyr::filter(term == predictor) %>%
        dplyr::select(term, or_ho = Odds_Ratio,
                      or_ho_low = OR_Lower_CI,
                      or_ho_high = OR_Upper_CI),
      by = "term"
    ) %>%
    dplyr::mutate(diff = or_ho - or_mh)
  
  return(list(
    mid_high = res_mid_high,
    high_only = res_high_only,
    comparison = compa_OR
  ))
}


#Predictor: Material Hardship
M_results_mh <- compare_binomial_models_m(imp_list_m, "mh.scale")

# M_results for predictor only
print(M_results_mh$comparison)

#RQ2
#exposure model (PSs) - calculate the weights 
w.imp_male <- weightthem(lon_3l ~ race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                         + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + ppvt.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt
                         , data = male_mids, approach = 'within', method = "gbm", #data mids object
                         estimand = "ATE")


bal.tab(w.imp_male, binary = "std")

#new.names in the stats. analysis (main) script

love.plot(w.imp_male, 
          drop.distance = TRUE, 
          #var.order = "unadjusted",
          abs = TRUE,
          binary = "std",
          line = TRUE, 
          thresholds = c(m = .1),
          var.names = new.names,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"))


#Outcome model (boys)
#eg: ducational attainment 
fits_m <- lapply(seq_along(w.imp_male$models), function(i) {
  data <- complete(w.imp_male, i)
  W <- w.imp_male$models[[i]]
  
  
  ordinal_weightit(edu_att22 ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                   + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                   data = data, #loop, from before
                   link = "logit", weightit = W, vcov = "HC0") 
})

#Difference - ordinal_weightit()
m.imp <- lapply(fits_m, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
  #type = "response")  # optional: gives results on probability scale
  
})

m.imp[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
m.imp #Useful to see! group=outcome level

pooled.m1 <- mice::pool(m.imp, dfcom = Inf)
summary(pooled.m1, conf.int = TRUE, exponentiate = T) #educational attainment




#2. SA Child's Temperament as a confounder
#Add Child's Temperament variable to the dataset, before imputation (all following steps as in the main analysis, adding this variable as confounder)

#select and explore items (6 items, scores 1-5) --> contained in "mydata" object
reverse_shy = c("m2b17c", "m2b17f", "m2b43c", "m2b43f") #years 5 (4pt scale)
mydata[ , reverse_shy] = 6 - mydata[ , reverse_shy]

#merge with data from non-resident mothers (m2b43 items)

mydata <- mydata %>% mutate(shy1= ifelse(is.na(m2b17a), m2b43a, m2b17a))
mydata <- mydata %>% mutate(shy2= ifelse(is.na(m2b17b), m2b43b, m2b17b))
mydata <- mydata %>% mutate(shy3= ifelse(is.na(m2b17c), m2b43c, m2b17c))
mydata <- mydata %>% mutate(shy4= ifelse(is.na(m2b17d), m2b43d, m2b17d))
mydata <- mydata %>% mutate(shy5= ifelse(is.na(m2b17e), m2b43e, m2b17e))
mydata <- mydata %>% mutate(shy6= ifelse(is.na(m2b17f), m2b43f, m2b17f))


summary(mydata$m2b17c)
summary(mydata$shy3) 

items_shy <- c("shy1", "shy2", "shy3", "shy4", "shy5", "shy6")

#sum items 
mydata <- mydata %>%
  mutate(emot_shy1 = rowSums(select(., all_of(items_shy)), na.rm = F)) 

summary(mydata$emot_shy1)


#3. SA- Alternative cut-off points for high school grades variable

#Dichotomize HS grades --> #1= mostly Bs or lower, 0= about half As/Bs(7), mostly As (8) 
imp_10it$hs_binAB <- imp_10it$hsgrades22 #new variable 
imp_10it$hs_binAB <- recode(imp_10it$hs_binAB, "1" = "1", "2" = "1", "3" = "1", "4" = "1", "5" = "1", "6" = "1", "7" = "0", "8" = "0") #6=mostly B's'
summary(imp_10it$hs_binAB) 

#Dichotomize HS grades --> #1= half As/Bs or lower, 0=  mostly As (8)
imp_10it$hs_binA <- imp_10it$hsgrades22 #new variable 
imp_10it$hs_binA <- recode(imp_10it$hs_binA, "1" = "1", "2" = "1", "3" = "1", "4" = "1", "5" = "1", "6" = "1", "7" = "1", "8" = "0") #6=mostly B's'
summary(imp_10it$hs_binA) 

imp_mids <- as.mids(imp_10it) # mids object

#Weights, so that the w.imp_gbmAB object contains the HS variable(s) coded with the new cut off points
w.imp_gbmAB <- weightthem(lon_3l ~ sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                          + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                          data = imp_mids, approach = 'within', method = "gbm", 
                          estimand = "ATE")


#EG: Low high-school grades - mostly Bs or lower 
fits3AB <- lapply(seq_along(w.imp_gbmAB$models), function(i) { 
  data <- complete(w.imp_gbmAB, i)
  
  W <- w.imp_gbmAB$models[[i]]
  
  
  glm_weightit(hs_binAB ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family= binomial, vcov = "HC0") 
  
})

#contrasts
comp.imp3AB <- lapply(fits3AB, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnoravg"
  ) 
  
})


pooled.try3AB <- mice::pool(comp.imp3AB, dfcom = Inf)
summary(pooled.try3AB, conf.int = TRUE, exponentiate = T)


#4. SA- Complete case analysis

ffcw4 <- ffcw3 %>% dplyr::select(-c(PPVT9raw, edu_aux22, neighbour3pcg, sum_no_violent))
summary(ffcw4)
complete <- ffcw4[complete.cases(ffcw4), ]
summary(complete) #1266 observations 

# DATA wrangling #

#1. Dichotomize HS grades
complete$hs_bin <- complete$hsgrades22 #new variable 
complete$hs_bin <- recode(complete$hs_bin, "1" = "1", "2" = "1", "3" = "1", "4" = "1", "5" = "1", "6" = "0", "7" = "0", "8" = "0") #6=mostly B's'
summary(complete$hs_bin) #1=half B, half C or lower, 0=at least mostly Bs

#2. Loneliness cateogry
#Loneliness scores
#recode so that 0=not lonely and 1=sometimes/often lonely
complete <- complete %>%
  mutate(lon5_cbp = recode(lon5_cbp, `2` = 1, `1`= 1, `0`= 0))

complete$lon5_cbp<-as.factor(complete$lon5_cbp)
summary(complete$lon5_cbp)

complete <- complete %>%
  mutate(lon9_cbp = recode(lon9_cbp, `3` = 1, `2`= 1, `1`= 0))

complete$lon9_cbp<-as.factor(complete$lon9_cbp)
summary(complete$lon9_cbp)

# Create a new variable with 4 levels and 3 levels based on lon5_cbp and lon9_cbp
complete <- complete %>%
  mutate(lon_combined = paste(lon5_cbp, lon9_cbp, sep = ","),
         lon_combined = factor(lon_combined, 
                               levels = c("1,1", "1,0", "0,1", "0,0"),
                               labels = c("chronic", "y5_lonely", "y9_lonely", "never")),
         lon_3l = case_when(
           lon_combined == "chronic" ~ 2,
           lon_combined == "y5_lonely" ~ 1,
           lon_combined == "y9_lonely" ~ 1,
           lon_combined == "never" ~ 0
         ), lon_3l = factor(lon_3l, levels = c(0, 1, 2), ordered = TRUE)
  )
summary(complete)

#3. Drop loneliness items

complete <- complete %>%
  dplyr::select(
    -c(lon5_cbp, lon9_cbp, lon_combined))

summary(complete)      

#4. Scales numerical variables

complete <- complete %>%
  mutate(
    mh.scale     = as.numeric(scale(mathard4m)),
    nce.scale    = as.numeric(scale(nce5)),
    chmalt.scale = as.numeric(scale(sum_conflict4)),
    ppvt.scale   = as.numeric(scale(PPVT5raw)),
    
    cm1age.scale = as.numeric(scale(cm1age)), #trying to improve the PS model for RQ2
    
    # Add interaction terms
    mh_ppvt      = mh.scale * ppvt.scale,
    nce_ppvt     = nce.scale * ppvt.scale,
    chmalt_ppvt  = chmalt.scale * ppvt.scale,
    
    age_mh = cm1age.scale * mh.scale
  ) 

# RQ1: Binomial models
complete <- complete %>% mutate(mid_high = ifelse(lon_3l == 0, 0, 1))  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
complete <- complete %>% mutate(high_only = ifelse(lon_3l == 2, 1, 0)) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)

#Main comparison function
compare_binomial_models <- function(predictor) {
  # Define common covariates
  covariates <- c(
    "sex_child1", "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Add predictor of interest
  formula_str <- paste("~", paste(c(covariates, predictor), collapse = " + "))
  
  # Model 1: mid_high ~ X
  fit_mid_high <- 
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = complete)
  
  
  # Model 2: high_only ~ X
  fit_high_only <- 
    glm(as.formula(paste("high_only", formula_str)), family = "binomial", data = complete)
  
  # No pooling now
  res_mid_high <- fit_mid_high
  res_high_only <- fit_high_only
  
  # Compare OR for predictor only 
  compa_OR <- res_mid_high %>%
    dplyr::filter(term == predictor) %>%
    dplyr::select(term, or_mh = Odds_Ratio,
                  or_mh_low = OR_Lower_CI,
                  or_mh_high = OR_Upper_CI) %>%
    left_join(
      res_high_only %>%
        dplyr::filter(term == predictor) %>%
        dplyr::select(term, or_ho = Odds_Ratio,
                      or_ho_low = OR_Lower_CI,
                      or_ho_high = OR_Upper_CI),
      by = "term"
    ) %>%
    dplyr::mutate(diff = or_ho - or_mh)
  
  return(list(
    mid_high = res_mid_high,
    high_only = res_high_only,
    comparison = compa_OR
  ))
}


compare_binomial_models <- function(predictor) {
  # Define common covariates
  covariates <- c(
    "sex_child1", "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Add predictor of interest
  formula_str <- paste("~", paste(c(covariates, predictor), collapse = " + "))
  
  # Model 1: mid_high ~ X
  fit_mid_high <- 
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = complete)
  
  # Model 2: high_only ~ X
  fit_high_only <- 
    glm(as.formula(paste("high_only", formula_str)), family = "binomial", data = complete)
  
  # Get tidy results with OR and CI
  tidy_with_or <- function(model) {
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::rename(
        Odds_Ratio = estimate,
        OR_Lower_CI = conf.low,
        OR_Upper_CI = conf.high
      )
  }
  
  res_mid_high <- tidy_with_or(fit_mid_high)
  res_high_only <- tidy_with_or(fit_high_only)
  
  # Compare OR for predictor only 
  compa_OR <- res_mid_high %>%
    dplyr::filter(term == predictor) %>%
    dplyr::select(term, or_mh = Odds_Ratio,
                  or_mh_low = OR_Lower_CI,
                  or_mh_high = OR_Upper_CI) %>%
    left_join(
      res_high_only %>%
        dplyr::filter(term == predictor) %>%
        dplyr::select(term, or_ho = Odds_Ratio,
                      or_ho_low = OR_Lower_CI,
                      or_ho_high = OR_Upper_CI),
      by = "term"
    ) %>%
    dplyr::mutate(diff = or_ho - or_mh)
  
  return(list(
    mid_high = res_mid_high,
    high_only = res_high_only,
    comparison = compa_OR
  ))
}


compare_binomial_int <- function(predictors) {
  # Define covariates common to all models
  covariates <- c(
    "sex_child1", "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Combine covariates and predictors
  all_terms <- c(covariates, predictors)
  formula_str <- paste("~", paste(all_terms, collapse = " + "))
  
  # Model 1: mid_high
  fit_mid_high <- glm(
    as.formula(paste("mid_high", formula_str)),
    family = "binomial", data = complete
  )
  
  # Model 2: high_only
  fit_high_only <- glm(
    as.formula(paste("high_only", formula_str)),
    family = "binomial", data = complete
  )
  
  # Tidy results with ORs + CIs
  tidy_with_or <- function(model) {
    broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::rename(
        Odds_Ratio = estimate,
        OR_Lower_CI = conf.low,
        OR_Upper_CI = conf.high
      )
  }
  
  res_mid_high <- tidy_with_or(fit_mid_high)
  res_high_only <- tidy_with_or(fit_high_only)
  
  # Comparison dataframe for all predictors of interest
  compa_OR <- res_mid_high %>%
    dplyr::filter(term %in% predictors) %>%
    dplyr::select(term, or_mh = Odds_Ratio,
                  or_mh_low = OR_Lower_CI,
                  or_mh_high = OR_Upper_CI) %>%
    left_join(
      res_high_only %>%
        dplyr::filter(term %in% predictors) %>%
        dplyr::select(term, or_ho = Odds_Ratio,
                      or_ho_low = OR_Lower_CI,
                      or_ho_high = OR_Upper_CI),
      by = "term"
    ) %>%
    dplyr::mutate(diff = or_ho - or_mh)
  
  return(list(
    mid_high = res_mid_high,
    high_only = res_high_only,
    comparison = compa_OR
  ))
}

#A couple of examples:

#1.Predictor: Material Hardship
result_mh.com <- compare_binomial_models("mh.scale")
print(result_mh.com)


#4. All SEF as predictors 
result_sef <- compare_binomial_int(predictors = c("mh.scale", "nce.scale", "chmalt.scale"))
print(result_sef$comparison)


# RQ2 
#Exposure model (GPSs) - calculate the weights 

W <- weightit(lon_3l ~ sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
              + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
              data = complete, method = "gbm", 
              estimand = "ATE")

bal.tab(W, binary = "std")

new.names <- c(sex_child1 = "Sex (Boy/Girl)",
               race_merged = "Race/ethnicity",
               lowbbweight = "Low weight at birth",
               cm1age = "Mother's age",
               race1_mother = "Mother's race/ethnicity",
               mborn = "Mother borned in the US",
               relst1 = "Parents' relationship",
               mother_edu3 = "Mother's education",
               cognit3_mother = "Mother's cognitive ability",
               cm3md_case_con = "Mother's Major Depression",
               cm3alc_case = "Mother's Alcohol Abuse",
               cm3drug_case = "Mother's Drug Abuse",
               cm3gad_case = "Mother's Generalized Anxiety Disorder",
               m_health3 = "Mother's general health status",
               health3 = "Child's general health status",
               disab3 = "Child's physical disability",
               mh.scale = "Material Hardship",
               nce.scale = "Neighborhood Collective Efficacy",
               chmalt.scale = "Child maltreatment",
               ppvt.scale = "Peabody Picture Vocabulary Test",
               mh_ppvt = "Int: Material Hardship x PPVT",
               nce_ppvt = "Int: Neighborhood CE x PPVT",
               chmalt_ppvt = "Int: Child maltreatment x PPVT")

# Pairwise comparisons for treatment groups 0 vs 1
love.plot(W, drop.distance = TRUE,
          abs = TRUE,
          binary = "std",
          line = TRUE,
          thresholds = c(m = .1),
          var.names = new.names,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"),
          stats = "mean.diffs",, which.treat = c("0", "1"))

# Pairwise comparisons for treatment groups 0 vs 2
love.plot(W, drop.distance = TRUE,
          abs = TRUE,
          binary = "std",
          line = TRUE,
          thresholds = c(m = .1),
          var.names = new.names,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"),
          stats = "mean.diffs",, which.treat = c("0", "2"))

# Pairwise comparisons for treatment groups 1 vs 2
love.plot(W,
          drop.distance = TRUE,
          abs = TRUE,
          binary = "std",
          line = TRUE,
          thresholds = c(m = .1),
          var.names = new.names,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"),
          stats = "mean.diffs",
          which.treat = c("1", "2"))


#Outcome model, eg: educational attainment
att_com <- ordinal_weightit(edu_att22 ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                            + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                            data = complete, 
                            link = "logit", weightit = W) 
#contrast
att_ATE <- 
  avg_comparisons(att_com,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg")

att_ATE 

#Aggregate to get one ATE per exposure contrast and exponentiate to get RR

att_ATE_summary <- att_ATE %>%
  group_by(contrast) %>%
  summarise(
    estimate   = mean(estimate, na.rm = TRUE),
    conf.low   = mean(conf.low, na.rm = TRUE),
    conf.high  = mean(conf.high, na.rm = TRUE),
    p.value    = mean(p.value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    RR       = exp(estimate),
    conf.lowRR = exp(conf.low),
    conf.highRR = exp(conf.high))

att_ATE_summary

## THE END ##

