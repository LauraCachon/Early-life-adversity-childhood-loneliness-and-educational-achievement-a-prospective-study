## R code to reproduce main analysis of the paper "Early-life adversity, childhood loneliness, and educational achievement: a prospective study" by Cachon-Alonso et al 2025

## 1. FFCW Data preparation
## 2. Multiple Imputations by Chained Equations
## 3. Binomial logistic regression models (RQ1)
## 4. Propensity score weights and weighted g-computation (RQ2)


### 1.FFCW Data preparation ##

library(tidyverse)
library(haven)
library(dplyr)
options(scipen=999)

#load FFCW data
load("/Users/lacachon/Desktop/FFCW-Paper3/data/data.RData") #waves 1-6 (year 1 to 15)

wave72024v12<-read_dta("/Users/lacachon/Desktop/FFCW-Paper3/data/FF_wave7_2024v1 2.dta") #wave 7 variables
wave7_edu <- wave72024v12 %>% select (idnum, hsgrades22=k7b13, edu_att22=ck7edu, ck7ethrace, cp7kedu,
                                      susp=k7b35b, susp4=k7b39_4, susp5=k7b39_5, susp6=k7b39_6, susp7=k7b39_7, susp8=k7b39_8, susp9=k7b39_9, susp10=k7b39_10, susp11=k7b39_11, susp12=k7b39_12, k7b39_1, k7b39_2, k7b39_3, k7b39_0)
summary(wave7_edu$susp)
allwaves<-full_join(data, wave7_edu, by="idnum")

#Replace with NA
mydata<-allwaves
mydata[mydata < 0] <- NA #Replace all negative values (not in wave, skipped, not asked, refused, etc) with NA


##Inclsuion criteria: 
#Select participants with loneliness complete data 
complete_lon <- mydata %>%
  filter(!is.na(p4l5) & !is.na(p5q3k)) #n=2466 

#Select those who attended a graded high school

complete_lon$hsgrades22<-as.factor(complete_lon$hsgrades22)
complete_lon <- complete_lon %>% filter(hsgrades22 != 9 | is.na(hsgrades22)) #10 people bye

mydata<- complete_lon

#Data wrangling of several variables:
#merge the two child health status measures to try to reduce NAs
summary(mydata$m3b2)
summary(mydata$m3b27)
mydata<- mydata %>% mutate(childhealth3 = ifelse(is.na(m3b2), m3b27, m3b2))
summary(mydata$childhealth3) #15 NAs less

#merge race/ethnicity variables reported at year 15 and 22
mydata<- mydata %>% mutate(race_merged = ifelse(is.na(ck7ethrace), ck6ethrace, ck7ethrace))
summary(mydata$race_merged)

#Replace NAs in NCE5 for mother's reports with data reported by fathers if available

mydata <- mydata %>% mutate(nce5_1= ifelse(is.na(m4i0m1), f4i0m1, m4i0m1))
mydata <- mydata %>% mutate(nce5_2= ifelse(is.na(m4i0m2), f4i0m2, m4i0m2))
mydata <- mydata %>% mutate(nce5_3= ifelse(is.na(m4i0m3), f4i0m3, m4i0m3))
mydata <- mydata %>% mutate(nce5_4= ifelse(is.na(m4i0m4), f4i0m4, m4i0m4))
mydata <- mydata %>% mutate(nce5_5= ifelse(is.na(m4i0m5), f4i0m5, m4i0m5))

mydata <- mydata %>% mutate(nce5_6= ifelse(is.na(m4i0n1), f4i0n1, m4i0n1))
mydata <- mydata %>% mutate(nce5_7= ifelse(is.na(m4i0n2), f4i0n2, m4i0n2))
mydata <- mydata %>% mutate(nce5_8= ifelse(is.na(m4i0n3), f4i0n3, m4i0n3))
mydata <- mydata %>% mutate(nce5_9= ifelse(is.na(m4i0n4), f4i0n4, m4i0n4))


##Create sum scores for certain scales

#1. MATERIAL HARDSHIP #Year 5: 11 items  
mathard_items <- c("m4i23d", "m4i23e", "m4i23f", "m4i23k", "m4i23a", "m4i23i", "m4i23j", "m4i23h", "m4i23b", "m4i23c", "m4i23g") #Year 5: 11 items 

# Recode the value 2 to 0 (reversed coding)
mydata <- mydata %>%
  mutate_at(vars(one_of(mathard_items)), ~ ifelse(. == 2, 0, .)) #now its 0=no, yes=1 


#2. NCE SCALE #Year 5: 9 items in total, scores from 1 to 4 
reverse_cols1 = c('f4i0m3', 'f4i0m4', 'm4i0m3', 'm4i0m4') #years 5 (4pt scale)
mydata[ , reverse_cols1] = 5 - mydata[ , reverse_cols1]

#year 3 (5pt scale) as auxiliary variable, scores from 1 to 5
reverse_cols2 = c('p3k2d', 'p3k2e') 
mydata[ , reverse_cols2] = 6 - mydata[ , reverse_cols2]


## 3. CHILD MALTTREATMENT (CTS) - subscales: psychological aggression, physical assault, and neglect
items_conflict4 <- c("p4g3", "p4g4", "p4g6", "p4g7", "p4g8", "p4g9",  #psychological aggression (6, 10, 8, 14, 9)
                     "p4g10", "p4g11", "p4g13", "p4g14", "p4g15",     #physical assault (7, 4, 11, 13, 3)
                     "p4g16", "p4g17", "p4g18", "p4g19")              #neglect (15, 16, 17, 18, 19)

recode_conflict <- function(x) {
  recode(x,
         `0` = 0,
         `1` = 1,
         `2` = 2,
         `3` = 4, # a score of 3 means 3-5 times
         `4` = 8,
         `5` = 15,
         `6` = 25, # a score of 6 means more than 20 times
         `7` = 0) #a score of 7 means "yes, but not in the past year"
}

# Recoding all relevant columns (e.g., p3j1 to p3j19, p4g1 to p4g19, p5q1a to p5q1n)
mydata <- mydata %>%
  mutate(across(all_of(items_conflict4), 
                ~ as.numeric(.))) %>% # Convert to numeric
  mutate(across(all_of(items_conflict4), 
                recode_conflict)) 

#sum score
mydata <- mydata %>%
  # Sum for items_conflict4
  mutate(sum_conflict4 = rowSums(select(., all_of(items_conflict4)), na.rm = F)) #if any of the item is NA, then final score is NA

#Non-violent discipline subscale(aux variable)
no_violent<- c("p4g1", "p4g5", "p4g12", "p4g2")

#sum here 
mydata <- mydata %>%
  # Sum for items_conflict4
  mutate(sum_no_violent = rowSums(select(., all_of(no_violent)), na.rm = F)) #if any of the item is NA, then final score is NA

#Parental relationship status: combining married and cohabiting at birth
mydata$cm1marf <- as.factor(mydata$cm1marf) #1=married
mydata$cm1cohf <- as.factor(mydata$cm1cohf) #1=cohabite

summary(mydata$cm1marf)
summary(mydata$cm1cohf)

mytable <- xtabs(~cm1cohf+cm1marf, data=mydata) #2x2 table
ftable(mytable) 

mydata <- mydata %>%
  mutate(relst1 = paste(cm1marf, cm1cohf, sep = ","),
         relst1 = factor(relst1, 
                         levels = c("1,0", "0,1", "0,0"),
                         labels = c("married", "cohabite", "single")))

summary(mydata$relst1)

#Truncate Health Status: (1=fair/poor health, 0=good/very good/excellent)
#child health status
mydata$health3 <- factor(ifelse(mydata$childhealth3 >= 4, 1, 0),
                         levels = c(0, 1))

mydata$childhealth3 <- as.factor(mydata$childhealth3)

#mother's health status
mydata$m_health3 <- factor(ifelse(mydata$m3j1 >= 4, 1, 0),
                           levels = c(0, 1))

mydata$m3j1 <- as.factor(mydata$m3j1)
summary(mydata$m3j1)
summary(mydata$m_health3)


## OUTCOMES ##

#DICHOTOMIZING high school grades
mydata$hs_bin <- mydata$hsgrades22 #new variable 
mydata$hs_bin <- recode(mydata$hs_bin, "1" = "1", "2" = "1", "3" = "1", "4" = "1", "5" = "1", "6" = "0", "7" = "0", "8" = "0") #6=mostly B's'
summary(mydata$hs_bin) #1=half B, half C or lower, 0=at least mostly Bs

#Creating the OUT-OF-SCHOOL SUSPENSION variable
items_susp <- paste0("susp", 4:12) #only data for the n=1100 that where suspended at some point

#sum suspensions from 4th grade to 12th grade
mydata <- mydata %>%
  # Sum for items_conflict4
  mutate(sum_susp = rowSums(select(., all_of(items_susp)), na.rm = F)) 

summary(mydata$sum_susp)

mydata$sum_susp_ <- as.factor(mydata$sum_susp)
summary(mydata$sum_susp_) #34 were suspended only before 4th (score of 0) For the rest, 1-9 susp between 4th grade and end of HS (this will be 0 in or analysis)
#12 people with sum_susp = 1, but NAs for when did it happen (this will be NA in our analysis)

#see:

#12 with NAs about when the suspension happened
problem_cases <- mydata %>% filter(susp == 1 & is.na(sum_susp))
prob<- problem_cases %>% select (idnum,
                                 susp, susp4, susp5, susp6, susp7, susp8, susp9, susp10, susp11, susp12, k7b39_3, k7b39_2, k7b39_1, k7b39_0)
#34 suspended before 4th grade
problem_cases2 <- mydata %>% filter(sum_susp == 0)
prob2<- problem_cases2 %>% select (idnum,
                                   susp, susp4, susp5, susp6, susp7, susp8, susp9, susp10, susp11, susp12, k7b39_3, k7b39_2, k7b39_1, k7b39_0)


### recode susp first so that 2=0 (no suspended)
# Recode the value 2 to 0 (reversed coding)
mydata$susp<-ifelse(mydata$susp == 2, 0, mydata$susp)
summary(mydata$susp) #0=no suspended, 1=yes, suspended

mydata$susp_ <- as.factor(mydata$susp)
summary(mydata$susp_)

#merge
mydata$suspHS <- ifelse(mydata$susp == 1, mydata$sum_susp, mydata$susp) 
summary(mydata$suspHS)

mydata$suspHS_ <- as.factor(mydata$suspHS)
summary(mydata$suspHS_) #done! 0=those that never were suspended, and those 34 that were before 4th grade (maybe we could do SA without those 34 to double checking case of reverse causality)

mydata$susp_bin<-ifelse(mydata$suspHS > 0, 1, mydata$suspHS) %>% as.factor
summary(mydata$susp_bin)

summary(mydata$cm3cogsc)

#Create FFCW2: all variables, with sum scores and individual items        
ffcw2 <- mydata %>% 
  select(
    idnum, 
    
    #confounders 
    sex_child1 = cm1bsex, #1=boy, 2=girl
    race_merged, #reported by YA and merged years 15 and 22
    lowbbweight = cm1lbw, #weight at birth in grams
    cm1age, 
    race1_mother = cm1ethrace, 
    mborn = m1h2, #were you born in the US? (yes/no)
    relst1, #married, cohabiting, single
    
    mother_edu3=cm3edu,
    cognit3_mother=cm3cogsc, #WAIS
    cm3md_case_con, #CIDI parental mental health - depression (binary)
    cm3alc_case, cm3drug_case, #CIDI alcohol and drug dependence 
    cm3gad_case, #CIDI anxiety
    m_health3,  #mother's general health status
    health3,    #child's general health status 
    disab3 = p3a2, #Does child have any physical disabilities?
    
    #socio-environmental indicators
    
    m4i23d, m4i23e, m4i23f, m4i23k, m4i23a, m4i23i, m4i23j, m4i23h, m4i23b, m4i23c, m4i23g,
    
    nce5_1, nce5_2, nce5_3, nce5_4, nce5_6, nce5_7, nce5_8, nce5_9, #replacing with father values if mother scores are missing
    
    p3k1a, p3k1b, p3k1c, p3k1d, p3k1e, p3k2a, p3k2b, p3k2c, p3k2d, p3k2e, #auxiliary NCE Y3
    
    sum_conflict4, sum_no_violent,
    
    #verbal ability
    
    PPVT5raw=ch4ppvtraw, 
    PPVT9raw=ch5ppvtraw, #auxiliary variable Y9
    
    #childhood loneliness
    
    lon5_cbp=p4l5, #Child complains of loneliness. Scale: CBP
    lon9_cbp=p5q3k,
    
    #Outcomes
    
    hsgrades22, edu_att22, susp_bin, edu_aux22=cp7kedu) %>% mutate(
      neighbour3pcg = p3k1a + p3k1b + p3k1c + p3k1d + p3k1e + p3k2a + p3k2b + p3k2c + p3k2d + p3k2e, #informal social control and levels of cohesion and trust (two subscales combined)
      nce5 = nce5_1 + nce5_2 + nce5_3 + nce5_4 + nce5_6 + nce5_7 + nce5_8 + nce5_9,
      mathard4m= m4i23d+m4i23e+m4i23f+m4i23k+m4i23a+m4i23i+m4i23j+m4i23h+m4i23b+m4i23c+m4i23g) #11 items!! (other waves less, but wont be used,so)


summary(ffcw2)

# Create FFCW3: all variables, only sum scores for CTS, MH and NCE (24 key variables, idnum, 4 aux variables, lon x2)
ffcw3 <- ffcw2 %>%
  select(
    -c(
      # Items used in neighbour3pcg
      p3k1a, p3k1b, p3k1c, p3k1d, p3k1e,
      p3k2a, p3k2b, p3k2c, p3k2d, p3k2e,
      
      # Items used in nce5
      nce5_1, nce5_2, nce5_3, nce5_4, 
      nce5_6, nce5_7, nce5_8, nce5_9,
      
      # Items used in mathard4m
      m4i23d, m4i23e, m4i23f, m4i23k,
      m4i23a, m4i23i, m4i23j, m4i23h,
      m4i23b, m4i23c, m4i23g
    )
  )

# Convert to factors 
ffcw3[c(
  "sex_child1", "lowbbweight", "mother_edu3", "race_merged", "race1_mother", "health3", "m_health3", "cm3md_case_con",
  "cm3gad_case", "cm3alc_case", "cm3drug_case", "mborn", "relst1", "disab3", "susp_bin")] <- lapply(ffcw3[c(
    "sex_child1", "lowbbweight", "mother_edu3", "race_merged", "race1_mother", "health3", "m_health3", "cm3md_case_con",
    "cm3gad_case", "cm3alc_case", "cm3drug_case", "mborn", "relst1","disab3", "susp_bin")], as.factor)

# Convert to ordered factors
ffcw3[c("lon5_cbp", "lon9_cbp", "edu_att22", "edu_aux22", "hsgrades22")] <- lapply(ffcw3[c(
  "lon5_cbp", "lon9_cbp", "edu_att22","edu_aux22", "hsgrades22")], ordered) 


ffcw3[c("cm1age", "cognit3_mother", "neighbour3pcg", "nce5", "mathard4m", "sum_conflict4", "sum_no_violent", "PPVT5raw", "PPVT9raw")] <- lapply(ffcw3[c("cm1age", "cognit3_mother", "neighbour3pcg", "nce5", "mathard4m", "sum_conflict4", "sum_no_violent", "PPVT5raw", "PPVT9raw")], as.numeric)

library(skimr)
skim(ffcw3) #--> FFCW3 READY FOR THE IMPUTATION (check that variables "type" makes sense)

## 2. MULTIPLE IMPUTATION BY CHAINED EQUATIONS ##
library(mice)

# 1. Precompute interaction variables where both inputs are observed
ffcw3$mh.va    <- with(ffcw3, (mathard4m - 1.14) * (PPVT5raw - 62.18))
ffcw3$nce.va   <- with(ffcw3, (nce5 - 17) * (PPVT5raw - 62.18))
ffcw3$malt.va  <- with(ffcw3, (sum_conflict4 - 39.6) * (PPVT5raw - 62.18))

# 2. Initialize mice for method and predictorMatrix setup
init <- mice(ffcw3, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

# 3. Exclude idnum and passive variables from being predictors
predM[, "idnum"] <- 0

# Set up passive imputation formulas
meth["mh.va"]    <- "~I((mathard4m - 1.14)*(PPVT5raw - 62.18))"
meth["nce.va"]   <- "~I((nce5 - 17)*(PPVT5raw - 62.18))"
meth["malt.va"]  <- "~I((sum_conflict4 - 39.6)*(PPVT5raw - 62.18))"

# Prevent circular prediction
predM[c("mathard4m", "PPVT5raw"), "mh.va"] <- 0
predM[c("nce5", "PPVT5raw"), "nce.va"]     <- 0
predM[c("sum_conflict4", "PPVT5raw"), "malt.va"] <- 0

# Prevent passive variables from predicting others
predM["mh.va", ]    <- 0
predM["nce.va", ]   <- 0
predM["malt.va", ]  <- 0

#Prevent edu_aux from predicting the variables in which it causes a problem
predM["sum_conflict4", "edu_aux22"] <- 0
predM["neighbour3pcg", "edu_aux22"] <- 0

# 4. Run the imputation
imputed_10 <- mice(ffcw3, maxit = 10,
                   method = meth, 
                   predictorMatrix = predM, 
                   m = 20, 
                   seed = 123)

plot(imputed_10) #convergence plot (looks good)

imp_10it <- complete(imputed_10, action = "long", include = TRUE)


# DATA CODING AFTER IMPUTATON #

#1. Dichotomize HS grades
imp_10it$hs_bin <- imp_10it$hsgrades22 #new variable 
imp_10it$hs_bin <- recode(imp_10it$hs_bin, "1" = "1", "2" = "1", "3" = "1", "4" = "1", "5" = "1", "6" = "0", "7" = "0", "8" = "0") #6=mostly B's'
summary(imp_10it$hs_bin) #1=half B, half C or lower, 0=at least mostly Bs

#2. Loneliness cateogry
#Loneliness scores
#recode so that 0=not lonely and 1=sometimes/often lonely
imp_10it <- imp_10it %>%
  mutate(lon5_cbp = recode(lon5_cbp, `2` = 1, `1`= 1, `0`= 0))

imp_10it$lon5_cbp<-as.factor(imp_10it$lon5_cbp)
summary(imp_10it$lon5_cbp)

imp_10it <- imp_10it %>%
  mutate(lon9_cbp = recode(lon9_cbp, `3` = 1, `2`= 1, `1`= 0))

imp_10it$lon9_cbp<-as.factor(imp_10it$lon9_cbp)
summary(imp_10it$lon9_cbp)

# Create a new variable with 4 levels and 3 levels based on lon5_cbp and lon9_cbp
imp_10it <- imp_10it %>%
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
summary(imp_10it)

#3. Drop auxiliary variables, loneliness items

imp_10it <- imp_10it %>%
  dplyr::select(
    -c(neighbour3pcg, edu_aux22, PPVT9raw, sum_no_violent, lon5_cbp, lon9_cbp, lon_combined))

summary(imp_10it)      

#4. Scales numerical variables - IN EACH IMPUTED DATASET
library(dplyr)

imp_10it <- imp_10it %>%
  group_by(.imp) %>%
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
  ) %>%
  ungroup()

imp_mids <- as.mids(imp_10it) #I need a mids object

## STATISTICAL ANALYSIS ##

#### RQ1 ###

## 3. BINOMIAL MODELS ## 
#following https://peopleanalytics-regression-book.org/ord-reg.html#testing-the-proportional-odds-assumption

# Create a list of 20 imputed datasets + add variables
imp_list <- lapply(1:20, function(i) {
  data_i <- complete(imp_mids, i)
  data_i$mid_high <- ifelse(data_i$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
  data_i$high_only <- ifelse(data_i$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)
  data_i
})

#inspect new variables
data_check <- imp_list[[1]] #new variables make sense 

#previous step done just once

#### macro function to run the binomial models and runs the comparison of ORs
library(broom)
library(dplyr)

# Pooling function 
pool_manual <- function(model_list) {
  estimates_list <- lapply(model_list, tidy)
  estimates_df <- bind_rows(estimates_list, .id = "imp")
  
  terms <- unique(estimates_df$term)
  
  pooled_results <- lapply(terms, function(term) {
    term_data <- estimates_df %>% filter(term == !!term)
    q_bar <- mean(term_data$estimate)
    u_bar <- mean(term_data$std.error^2)
    b <- var(term_data$estimate)
    t_var <- u_bar + (1 + 1/length(model_list)) * b
    se_total <- sqrt(t_var)
    df <- (length(model_list) - 1) * (1 + u_bar / ((1 + 1/length(model_list)) * b))^2
    p_value <- 2 * pt(-abs(q_bar / se_total), df = df)
    
    ci_low <- q_bar - qt(0.975, df = df) * se_total
    ci_high <- q_bar + qt(0.975, df = df) * se_total
    
    or <- exp(q_bar)
    or_low <- exp(ci_low)
    or_high <- exp(ci_high)
    
    data.frame(
      term = term,
      estimate = q_bar,
      std.error = se_total,
      Lower_CI = ci_low,
      Upper_CI = ci_high,
      Odds_Ratio = or,
      OR_Lower_CI = or_low,
      OR_Upper_CI = or_high,
      p.value = p_value
    )
  })
  
  bind_rows(pooled_results)
}

# Main comparison function
compare_binomial_models <- function(imp_list, predictor) {
  # Define common covariates
  covariates <- c(
    "sex_child1", "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Add predictor of interest
  formula_str <- paste("~", paste(c(covariates, predictor), collapse = " + "))
  
  # Model 1: mid_high ~ X
  fit_mid_high <- lapply(imp_list, function(data) {
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = data)
  })
  
  # Model 2: high_only ~ X
  fit_high_only <- lapply(imp_list, function(data) {
    glm(as.formula(paste("high_only", formula_str)), family = "binomial", data = data)
  })
  
  # Pool and extract results
  res_mid_high <- pool_manual(fit_mid_high)
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

#for models with 3 predictor and interaction
compare_binomial_int <- function(imp_list, predictors) {
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
  fit_mid_high <- lapply(imp_list, function(data) {
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = data)
  })
  
  # Model 2: high_only
  fit_high_only <- lapply(imp_list, function(data) {
    glm(as.formula(paste("high_only", formula_str)), family = "binomial", data = data)
  })
  
  # Pool results
  res_mid_high <- pool_manual(fit_mid_high)
  res_high_only <- pool_manual(fit_high_only)
  
  # Create comparison dataframe for all predictors
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

#1.Predictor: Material Hardship
result_mh <- compare_binomial_models(imp_list, "mh.scale")

#Results for predictor only
print(result_mh$comparison)

#For full results
print(result_mh$mid_high)
print(result_mh$high_only)

#2. Predictor: NCE
result_nce <- compare_binomial_models(imp_list, "nce.scale")
print(result_nce$comparison)

#3. Predictor: Child maltreatment
result_chmalt <- compare_binomial_models(imp_list, "chmalt.scale")
print(result_chmalt$comparison)

#4. All SEF as predictors 
result_sef <- compare_binomial_int(imp_list, predictors = c("mh.scale", "nce.scale", "chmalt.scale"))
print(result_sef$comparison)

#5. Predictor: Verbal Ability
result_ppvt <- compare_binomial_models(imp_list, "ppvt.scale")
print(result_ppvt$comparison)

#6. Predictor: MH, VA, MH x VA
result_mh.va <- compare_binomial_int(imp_list, predictors = c("mh.scale", "ppvt.scale", "mh.scale:ppvt.scale")) #misses the interaction in the output
print(result_mh.va$comparison)

#7. Predictor: NCE, VA, NCE x VA
result_nce.va <- compare_binomial_int(imp_list, predictors = c("nce.scale", "ppvt.scale", "nce.scale:ppvt.scale"))
print(result_nce.va$comparison)

#8. Predictor CH MALT, VA, CH MALT x VA
result_chmalt.va <- compare_binomial_int(imp_list, predictors = c("chmalt.scale", "ppvt.scale", "chmalt.scale:ppvt.scale"))
print(result_chmalt.va$comparison)


## ANALYSE INTERACTIONS FURTHER

## 1. Simple slopes for binomial model: persistent vs transient/never loneliness (whole sample)
int_mids <- imp_10it #dataset, not mids object

int_mids$persist_tran <- ifelse(int_mids$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (persistent/transient loneliness vs never)
int_mids$persistent_lon <- ifelse(int_mids$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (persistent loneliness vs transient/never) <<- interaction here

#dataset 1
imputed_1 <- subset(int_mids, .imp == 1)
summary(imputed_1)

library(emmeans)
library(ggplot2)

# Fit glm (first imputation)
fit1 <- glm(persistent_lon ~ sex_child1 + race_merged + lowbbweight + cm1age + race1_mother +
              mborn + relst1 + mother_edu3 + cognit3_mother +
              cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case +
              m_health3 + health3 + disab3 + mh.scale * ppvt.scale,
            data = imputed_1, family = "binomial")

zvals <- c(-1, 0, 1)

emtr <- emtrends(
  fit1,
  specs = "ppvt.scale",
  var = "mh.scale",
  at = list(ppvt.scale = zvals),
  nuisance = c("sex_child1", "race_merged", "lowbbweight", "cm1age",
               "race1_mother", "mborn", "relst1", "mother_edu3", "cognit3_mother",
               "cm3md_case_con", "cm3alc_case", "cm3drug_case", "cm3gad_case",
               "m_health3", "health3", "disab3")
) 

#ORs with CIs
# Get confidence intervals
ci_emtr <- confint(emtr)
ci_df <- as.data.frame(ci_emtr)

# Exponentiate log-odds to ORs
ci_df$OR      <- exp(ci_df$mh.scale.trend)
ci_df$OR_low  <- exp(ci_df$asymp.LCL)
ci_df$OR_high <- exp(ci_df$asymp.UCL)

ci_df


#Plot --> gets predicted probabilities 
library(ggeffects)
ggpredict(fit1, terms = c("mh.scale", "ppvt.scale [-1,0,1]")) |>
  plot()

## RQ2

# 4. Propensity scores weights and weighted g-computation

#weighting the imputed data
library(MatchThem)
library(WeightIt)
library(cobalt)
library(marginaleffects)

#exposure model (PSs) - calculate the weights 

#Generalized Boosted Models (best model after trying different weight estimation methods)

w.imp_gbm <- weightthem(lon_3l ~ sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                        + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                        data = imp_mids, approach = 'within', method = "gbm", 
                        estimand = "ATE")

w.imp_gbm
summary(w.imp_gbm)
summary(get.w(w.imp_gbm))

#Covariate balance table
bal.tab(w.imp_gbm)
bal.tab(w.imp_gbm, m.threshold = 0.1)


#Covariate balance plot
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


love.plot(w.imp_gbm, 
          drop.distance = TRUE, 
          #var.order = "unadjusted",
          abs = TRUE,
          stars = "std",
          line = TRUE, 
          thresholds = c(m = .1),
          var.names = new.names1,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"))



## Outcome model: weighted g-computation 

#Outcome 1: educational attainment
fits <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  
  ordinal_weightit(edu_att22 ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                   + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                   data = data, #loop, from before
                   link = "logit", weightit = W, vcov = "HC0") 
}) 

#Difference - ordinal_weightit()
comp.imp <- lapply(fits, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg" 
  ) 
  
  
})


comp.imp[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
comp.imp #Useful to see! group=outcome level

pooled.try <- mice::pool(comp.imp, dfcom = Inf)
summary(pooled.try, conf.int = TRUE, exponentiate = T) #educational attainment (RR results)



#Outcome 2: out-of-school suspension 
fits2 <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  
  glm_weightit(susp_bin ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0") 
}) 


#Contrasts
comp.imp2 <- lapply(fits2, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  
  ) 
  
})


comp.imp2[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
comp.imp2 #Useful to see! group=outcome level

pooled.try2 <- mice::pool(comp.imp2, dfcom = Inf)
summary(pooled.try2, conf.int = TRUE, exponentiate = T) #suspension


#Outcome 3: low high-school grades

fits3 <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  
  glm_weightit(hs_bin ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family= binomial, vcov = "HC0") 
})

comp.imp3 <- lapply(fits3, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg" 
  ) 
  
})


comp.imp3[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
comp.imp3 #Useful to see! group=outcome level

pooled.try3 <- mice::pool(comp.imp3, dfcom = Inf)
summary(pooled.try3, conf.int = TRUE, exponentiate = T)




