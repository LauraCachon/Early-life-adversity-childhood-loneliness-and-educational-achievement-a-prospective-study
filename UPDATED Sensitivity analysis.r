## SENSITIVITY ANALYSIS ##

#1. Sex-stratified analysis (Tables s4-s6, Figure s5)
#2. Including child's temperament as a confounder (Table s7-s8, Figure s6))
#3. Alternative cut-off points for (Low) High School grades (Table s9, s10)
#4. Complete case analysis (Table s11, s12, Figure s7)
#5. Multiple imputation in the original study sample (Table s13-s14, Figure s8)

#1. Sex stratified analysis
#Almost identical to main analysis, but conducted across two datasets (boys/girls), and therefore not including the variable sex_child1 as a confounder
#starting form the imp_10it dataset after imputation and data wrangling (Created in Summary code-main analysis)

# Convert the full dataset first
imp_mids #make sure it is only selected sample

library(mice)
male_mids   <- filter(imp_mids, sex_child1 == 1)
female_mids <- filter(imp_mids, sex_child1 == 2)

## Binomial models

# Create a list of 20 imputed datasets + add variables: male/female
#male
imp_list_m <- lapply(1:20, function(i) {
  data_i <- complete(male_mids, i)
  data_i$mid_high <- ifelse(data_i$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
  data_i$high_only <- ifelse(data_i$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)
  data_i
})

#female
imp_list_f <- lapply(1:20, function(i) {
  data_i <- complete(female_mids, i)
  data_i$mid_high <- ifelse(data_i$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
  data_i$high_only <- ifelse(data_i$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)
  data_i
})


#### macro function to run the binomial models and runs the comparison of ORs
library(broom)
library(dplyr)

# Pooling function --> same for male and female (used in the later functions)
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

#Binomial models: males

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

#for models with 3 predictor and interaction males
compare_binomial_int_m <- function(imp_list_m, predictors) {
  # Define covariates common to all models
  covariates <- c(
    "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Combine covariates and predictors
  all_terms <- c(covariates, predictors)
  formula_str <- paste("~", paste(all_terms, collapse = " + "))
  
  # Model 1: mid_high
  fit_mid_high <- lapply(imp_list_m, function(data) {
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = data)
  })
  
  # Model 2: high_only
  fit_high_only <- lapply(imp_list_m, function(data) {
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
M_results_mh <- compare_binomial_models_m(imp_list_m, "mh.scale")

#M_resultss for predictor only
print(M_results_mh$comparison)

#For full M_resultss
print(M_results_mh$mid_high)
print(M_results_mh$high_only)

#2. Predictor: NCE
M_results_nce <- compare_binomial_models_m(imp_list_m, "nce.scale")
print(M_results_nce$comparison)

#3. Predictor: Child maltreatment
M_results_chmalt <- compare_binomial_models_m(imp_list_m, "chmalt.scale")
print(M_results_chmalt$comparison)

#4. All SEF as predictors 
M_results_sef <- compare_binomial_int_m(imp_list_m, predictors = c("mh.scale", "nce.scale", "chmalt.scale"))
print(M_results_sef$comparison)

#5. Predictor: Verbal Ability
M_results_ppvt <- compare_binomial_models_m(imp_list_m, "ppvt.scale")
print(M_results_ppvt$comparison)

#6. Predictor: MH, VA, MH x VA
M_results_mh.va <- compare_binomial_int_m(imp_list_m, predictors = c("mh.scale", "ppvt.scale", "mh.scale:ppvt.scale")) #misses the int_meraction in the output
print(M_results_mh.va$comparison)

#7. Predictor: NCE, VA, NCE x VA
M_results_nce.va <- compare_binomial_int_m(imp_list_m, predictors = c("nce.scale", "ppvt.scale", "nce.scale:ppvt.scale"))
print(M_results_nce.va$comparison)

#8. Predictor CH MALT, VA, CH MALT x VA
M_results_chmalt.va <- compare_binomial_int_m(imp_list_m, predictors = c("chmalt.scale", "ppvt.scale", "chmalt.scale:ppvt.scale"))
print(M_results_chmalt.va$comparison)

#interaction 
interaction_m <- imp_list_m[[1]] %>%
  mutate(persist_tran   = ifelse(lon_3l == 0, 0, 1),
         persistent_lon = ifelse(lon_3l == 2, 1, 0))

library(emmeans)
library(ggplot2)

# Interaction (males): persistent loneliness MH x VA
int_m <- glm(persistent_lon ~ race_merged + lowbbweight + cm1age + race1_mother +
               mborn + relst1 + mother_edu3 + cognit3_mother +
               cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case +
               m_health3 + health3 + disab3 + mh.scale * ppvt.scale,
             data = interaction_m, family = "binomial")

zvals <- c(-1, 0, 1)

emtr_m <- emtrends(
  int_m,
  specs = "ppvt.scale",
  var = "mh.scale",
  at = list(ppvt.scale = zvals),
  nuisance = c("race_merged", "lowbbweight", "cm1age",
               "race1_mother", "mborn", "relst1", "mother_edu3", "cognit3_mother",
               "cm3md_case_con", "cm3alc_case", "cm3drug_case", "cm3gad_case",
               "m_health3", "health3", "disab3")
) 

#ORs with CIs
# Get confidence intervals
ci_emtr_m <- confint(emtr_m)
ci_df_m <- as.data.frame(ci_emtr_m)

# Exponentiate log-odds to ORs
ci_df_m$OR      <- exp(ci_df_m$mh.scale.trend)
ci_df_m$OR_low  <- exp(ci_df_m$asymp.LCL)
ci_df_m$OR_high <- exp(ci_df_m$asymp.UCL)

ci_df_m

#Plot --> gets predicted probabilities 
library(ggeffects)
ggpredict(int_m, terms = c("mh.scale", "ppvt.scale [-1,0,1]")) |>
  plot()



#Binomial models: female

# Main comparison function - females
compare_binomial_models_f <- function(imp_list_f, predictor) {
  # Define common covariates
  covariates <- c(
    "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Add predictor of interest
  formula_str <- paste("~", paste(c(covariates, predictor), collapse = " + "))
  
  # Model 1: mid_high ~ X
  fit_mid_high <- lapply(imp_list_f, function(data) {
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = data)
  })
  
  # Model 2: high_only ~ X
  fit_high_only <- lapply(imp_list_f, function(data) {
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

#for models with 3 predictor and interaction males
compare_binomial_int_f <- function(imp_list_f, predictors) {
  # Define covariates common to all models
  covariates <- c(
    "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
    "relst1", "mother_edu3", "cognit3_mother", "cm3md_case_con", "cm3alc_case",
    "cm3drug_case", "cm3gad_case", "m_health3", "health3", "disab3"
  )
  
  # Combine covariates and predictors
  all_terms <- c(covariates, predictors)
  formula_str <- paste("~", paste(all_terms, collapse = " + "))
  
  # Model 1: mid_high
  fit_mid_high <- lapply(imp_list_f, function(data) {
    glm(as.formula(paste("mid_high", formula_str)), family = "binomial", data = data)
  })
  
  # Model 2: high_only
  fit_high_only <- lapply(imp_list_f, function(data) {
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
F_results_mh <- compare_binomial_models_f(imp_list_f, "mh.scale")

#F_resultss for predictor only
print(F_results_mh$comparison)

#For full F_resultss
print(F_results_mh$mid_high)
print(F_results_mh$high_only)

#2. Predictor: NCE
F_results_nce <- compare_binomial_models_f(imp_list_f, "nce.scale")
print(F_results_nce$comparison)

#3. Predictor: Child maltreatment
F_results_chmalt <- compare_binomial_models_f(imp_list_f, "chmalt.scale")
print(F_results_chmalt$comparison)

#4. All SEF as predictors 
F_results_sef <- compare_binomial_int_f(imp_list_f, predictors = c("mh.scale", "nce.scale", "chmalt.scale"))
print(F_results_sef$comparison)

#5. Predictor: Verbal Ability
F_results_ppvt <- compare_binomial_models_f(imp_list_f, "ppvt.scale")
print(F_results_ppvt$comparison)

#6. Predictor: MH, VA, MH x VA
F_results_mh.va <- compare_binomial_int_f(imp_list_f, predictors = c("mh.scale", "ppvt.scale", "mh.scale:ppvt.scale")) #misses the int_meraction in the output
print(F_results_mh.va$comparison)

#7. Predictor: NCE, VA, NCE x VA
F_results_nce.va <- compare_binomial_int_f(imp_list_f, predictors = c("nce.scale", "ppvt.scale", "nce.scale:ppvt.scale"))
print(F_results_nce.va$comparison)

#8. Predictor CH MALT, VA, CH MALT x VA
F_results_chmalt.va <- compare_binomial_int_f(imp_list_f, predictors = c("chmalt.scale", "ppvt.scale", "chmalt.scale:ppvt.scale"))
print(F_results_chmalt.va$comparison)



#Sex stratified tables

library(flextable)
library(officer)
library(dplyr)


# TABLE S5: SEX-STRATIFIED (FEMALE) - Binomial models  #

## ------------------------------------------------------------
## 1) Helper: format "OR (low - high)" rounded to 2 decimals
## ------------------------------------------------------------
fmt_or <- function(or, low, high) {
  sprintf("%.2f (%.2f - %.2f)", or, low, high)
}

## ------------------------------------------------------------
## 2) Helper: turn one compare_binomial_*_f() result data frame
##    into a small tibble of predictor rows with formatted OR strings
##    (col1 = or_mh, col2 = or_ho -- rename below if reversed)
## ------------------------------------------------------------
make_rows <- function(df, labels) {
  # labels = character vector matching the row order in df, giving the
  # display name you want in the "Predictors" column
  data.frame(
    Predictor = labels,
    OR1 = fmt_or(df$or_mh, df$or_mh_low, df$or_mh_high),
    OR2 = fmt_or(df$or_ho, df$or_ho_low, df$or_ho_high),
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------
## 3) Build each model's rows from your existing result objects
## ------------------------------------------------------------

model_list <- list(
  "Model 1" = make_rows(
    F_results_mh$comparison,
    c("Material hardship")
  ),
  
  "Model 2" = make_rows(
    F_results_nce$comparison,
    c("Neighborhood Collective Efficacy")
  ),
  
  "Model 3" = make_rows(
    F_results_chmalt$comparison,
    c("Child maltreatment")
  ),
  
  "Model 4" = make_rows(
    F_results_sef$comparison,
    c("Material Hardship", "Neighborhood Collective Efficacy", "Child maltreatment")
  ),
  
  "Model 5" = make_rows(
    F_results_ppvt$comparison,
    c("Verbal ability (PPVT)")
  ),
  
  "Model 6" = make_rows(
    F_results_mh.va$comparison,
    c("Material hardship", "Verbal ability (PPVT)", "Interaction: material hardship x verbal ability")
  ),
  
  "Model 7" = make_rows(
    F_results_nce.va$comparison,
    c("Neighborhood Collective Efficacy", "Verbal ability (PPVT)", "Interaction: neighbourhood x verbal ability")
  ),
  
  "Model 8" = make_rows(
    F_results_chmalt.va$comparison,
    c("Child maltreatment", "Verbal ability (PPVT)", "Interaction: child maltreatment x verbal ability")
  )
  
)

## ------------------------------------------------------------
## 4) Stack into one long data frame with a Model column
##    (repeated model name will be merged visually in step 5)
## ------------------------------------------------------------
table_df <- do.call(rbind, lapply(names(model_list), function(m) {
  d <- model_list[[m]]
  data.frame(Model = m, d, stringsAsFactors = FALSE)
}))

rownames(table_df) <- NULL

## ------------------------------------------------------------
## 5) Build the flextable with merged Model cells + custom headers
## ------------------------------------------------------------
ft <- flextable(table_df) %>%
  set_header_labels(
    Model = "Model",
    Predictor = "Predictors",
    OR1 = "Odds ratio (Recurrent/Transient vs Not reported)",
    OR2 = "Odds Ratio (Recurrent vs Transient/Not reported)"
  ) %>%
  merge_v(j = "Model") %>%              # merge repeated Model cells vertically
  theme_booktabs() %>%
  align(j = c("OR1", "OR2"), align = "center", part = "all") %>%
  valign(j = "Model", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  autofit()

ft  # preview in RStudio Viewer


# TABLE S4: SEX STRATIFIED (MALE) - Binomial models #

## ============================================================
## Same table, for boys (M_results_* objects)
## ============================================================

model_list_m <- list(
  
  "Model 1" = make_rows(
    M_results_mh$comparison,
    c("Material hardship")
  ),
  
  "Model 2" = make_rows(
    M_results_nce$comparison,
    c("Neighborhood Collective Efficacy")
  ),
  
  "Model 3" = make_rows(
    M_results_chmalt$comparison,
    c("Child maltreatment")
  ),
  
  "Model 4" = make_rows(
    M_results_sef$comparison,
    c("Material Hardship", "Neighborhood Collective Efficacy", "Child maltreatment")
  ),
  
  "Model 5" = make_rows(
    M_results_ppvt$comparison,
    c("Verbal ability (PPVT)")
  ),
  
  "Model 6" = make_rows(
    M_results_mh.va$comparison,
    c("Material hardship", "Verbal ability (PPVT)", "Interaction: material hardship x verbal ability")
  ),
  
  "Model 7" = make_rows(
    M_results_nce.va$comparison,
    c("Neighborhood Collective Efficacy", "Verbal ability (PPVT)", "Interaction: neighbourhood x verbal ability")
  ),
  
  "Model 8" = make_rows(
    M_results_chmalt.va$comparison,
    c("Child maltreatment", "Verbal ability (PPVT)", "Interaction: child maltreatment x verbal ability")
  )
)

table_df_m <- do.call(rbind, lapply(names(model_list_m), function(m) {
  d <- model_list_m[[m]]
  data.frame(Model = m, d, stringsAsFactors = FALSE)
}))

rownames(table_df_m) <- NULL

ft_m <- flextable(table_df_m) %>%
  set_header_labels(
    Model = "Model",
    Predictor = "Predictors",
    OR1 = "Odds ratio (Recurrent/Transient vs Not reported)",
    OR2 = "Odds Ratio (Recurrent vs Transient/Not reported)"
  ) %>%
  merge_v(j = "Model") %>%
  theme_booktabs() %>%
  align(j = c("OR1", "OR2"), align = "center", part = "all") %>%
  valign(j = "Model", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  autofit()

ft_m  # preview



## Analyse significant interactions further

## 1. Simple slopes for binomial model
interaction <- imp_10it %>% filter(.imp ==1)

interaction$persist_tran <- ifelse(interaction$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (persistent/transient loneliness vs never)
interaction$persistent_lon <- ifelse(interaction$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (persistent loneliness vs transient/never) <<- interaction here

interaction_m <- interaction %>% filter(sex_child1 == 1)
interaction_f <- interaction %>% filter(sex_child1 == 2)

library(emmeans)
library(ggplot2)

# Interaction (males): persistent loneliness MH x VA
int_m <- glm(persistent_lon ~ race_merged + lowbbweight + cm1age + race1_mother +
               mborn + relst1 + mother_edu3 + cognit3_mother +
               cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case +
               m_health3 + health3 + disab3 + mh.scale * ppvt.scale,
             data = interaction_m, family = "binomial")

zvals <- c(-1, 0, 1)

emtr_m <- emtrends(
  int_m,
  specs = "ppvt.scale",
  var = "mh.scale",
  at = list(ppvt.scale = zvals),
  nuisance = c("race_merged", "lowbbweight", "cm1age",
               "race1_mother", "mborn", "relst1", "mother_edu3", "cognit3_mother",
               "cm3md_case_con", "cm3alc_case", "cm3drug_case", "cm3gad_case",
               "m_health3", "health3", "disab3")
) 

#ORs with CIs
# Get confidence intervals
ci_emtr_m <- confint(emtr_m)
ci_df_m <- as.data.frame(ci_emtr_m)

# Exponentiate log-odds to ORs
ci_df_m$OR      <- exp(ci_df_m$mh.scale.trend)
ci_df_m$OR_low  <- exp(ci_df_m$asymp.LCL)
ci_df_m$OR_high <- exp(ci_df_m$asymp.UCL)

ci_df_m

#Plot --> gets predicted probabilities 
library(ggeffects)
ggpredict(int_m, terms = c("mh.scale", "ppvt.scale [-1,0,1]")) |>
  plot()

### RQ2 ###

#exposure model (PSs) - calculate the weights 
??weightthem
w.imp_male <- weightthem(lon_3l ~ race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                         + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + ppvt.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt
                         , data = male_mids, approach = 'within', method = "gbm", #data mids object
                         estimand = "ATE")

w.imp_female <- weightthem(lon_3l ~ race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                           + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + ppvt.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt
                           , data = female_mids, approach = 'within', method = "gbm", #data mids object
                           estimand = "ATE")


bal.tab(w.imp_male, binary = "std")
bal.tab(w.imp_female, binary = "std")

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

love.plot(w.imp_female, 
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


## OUTCOME MODELS - BOYS ##
## USING "HC0" sandwich variance matrix ##
#1. low high-school grades

fits3_m <- lapply(seq_along(w.imp_male$models), function(i) {
  data <- complete(w.imp_male, i)
  W <- w.imp_male$models[[i]]
  
  
  glm_weightit(hs_bin ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family= binomial, vcov = "HC0") #Increase R in the final models 
})

m.imp3 <- lapply(fits3_m, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"
  ) 
  
})

pooled.m3 <- mice::pool(m.imp3, dfcom = Inf)
summary(pooled.m3, conf.int = TRUE, exponentiate = T)

#2. out-of-school suspension 
fits2_m <- lapply(seq_along(w.imp_male$models), function(i) {
  data <- complete(w.imp_male, i)
  W <- w.imp_male$models[[i]]
  
  
  glm_weightit(susp_bin ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0") 
}) 

#Contrasts
m.imp2 <- lapply(fits2_m, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"
  ) 
  
})

pooled.m2 <- mice::pool(m.imp2, dfcom = Inf)
summary(pooled.m2, conf.int = TRUE, exponentiate = T) #suspension

#3. educational attainment 
fits_m <- lapply(seq_along(w.imp_male$models), function(i) {
  data <- complete(w.imp_male, i)
  W <- w.imp_male$models[[i]]
  
  
  glm_weightit(edu_nonacad ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0")
})


#Difference 
m.imp <- lapply(fits_m, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
})

pooled.m1 <- mice::pool(m.imp, dfcom = Inf)
summary(pooled.m1, conf.int = TRUE, exponentiate = T) #educational attainment


# OUTCOME MDOELS - GILRS #

#1. low high-school grades
fits3_f <- lapply(seq_along(w.imp_female$models), function(i) {
  data <- complete(w.imp_female, i)
  W <- w.imp_female$models[[i]]
  
  
  glm_weightit(hs_bin ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family= binomial, vcov = "HC0") #Increase R in the final models 
})

f.imp3 <- lapply(fits3_f, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
})

pooled.f3 <- mice::pool(f.imp3, dfcom = Inf)
summary(pooled.f3, conf.int = TRUE, exponentiate = T)


#2. out-of-school suspension 
fits2_f <- lapply(seq_along(w.imp_female$models), function(i) {
  data <- complete(w.imp_female, i)
  W <- w.imp_female$models[[i]]
  
  
  glm_weightit(susp_bin ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0") 
}) 

#Contrasts
f.imp2 <- lapply(fits2_f, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
})

pooled.f2 <- mice::pool(f.imp2, dfcom = Inf)
summary(pooled.f2, conf.int = TRUE, exponentiate = T) #suspension

#3. educational attainment 
fits_f <- lapply(seq_along(w.imp_female$models), function(i) {
  data <- complete(w.imp_female, i)
  W <- w.imp_female$models[[i]]
  
  
  glm_weightit(edu_nonacad ~ lon_3l + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0")
}) 


#Difference - ordinal_weightit()
f.imp <- lapply(fits_f, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
  #type = "response")  # optional: gives results on probability scale
  
})

pooled.f1 <- mice::pool(f.imp, dfcom = Inf)
summary(pooled.f1, conf.int = TRUE, exponentiate = T) #educational attainment


#2. SA Child's Temperament as a confounder
#Add Child's Temperament variable to the dataset, before imputation (all following steps as in the main analysis, adding this variable as confounder)

## first,  add Child emotionality and shyness as covariate - ##
emot_SA <- da31622.0001 %>% select(idnum=IDNUM, m2b17a=M2B17A, m2b43a=M2B43A, m2b17b=M2B17B, m2b43b=M2B43B, m2b17c=M2B17C, m2b43c=M2B43C, 
                                   m2b17d=M2B17D, m2b43d=M2B43D, m2b17e=M2B17E, m2b43e=M2B43E, m2b17f=M2B17F, m2b43f=M2B43F) #select items 

emot_SA <- emot_SA %>%
  mutate(across(everything(), ~ as.numeric(sub("^\\(([-0-9]+)\\).*", "\\1", as.character(.)))))


incl_emot <- left_join(mydata, emot_SA, by="idnum")

#Replace with NA
incl_emot[incl_emot < 0] <- NA #Replace all negative values (not in wave, skipped, not asked, refused, etc) with NA

#select and explore items (6 items, scores 1-5) --> contained in "incl_emot" object
reverse_shy = c("m2b17c", "m2b17f", "m2b43c", "m2b43f") #years 5 (4pt scale)
incl_emot[ , reverse_shy] = 6 - incl_emot[ , reverse_shy]

#merge with data from non-resident mothers (m2b43 items)
incl_emot <- incl_emot %>% mutate(shy1= ifelse(is.na(m2b17a), m2b43a, m2b17a))
incl_emot <- incl_emot %>% mutate(shy2= ifelse(is.na(m2b17b), m2b43b, m2b17b))
incl_emot <- incl_emot %>% mutate(shy3= ifelse(is.na(m2b17c), m2b43c, m2b17c))
incl_emot <- incl_emot %>% mutate(shy4= ifelse(is.na(m2b17d), m2b43d, m2b17d))
incl_emot <- incl_emot %>% mutate(shy5= ifelse(is.na(m2b17e), m2b43e, m2b17e))
incl_emot <- incl_emot %>% mutate(shy6= ifelse(is.na(m2b17f), m2b43f, m2b17f))


summary(incl_emot$m2b17c)
summary(incl_emot$shy3) #reduced NAs?

items_shy <- c("shy1", "shy2", "shy3", "shy4", "shy5", "shy6")

#sum items 
incl_emot <- incl_emot %>%
  mutate(emot_shy1 = rowSums(select(., all_of(items_shy)), na.rm = F)) #if any of the item is NA, then final score is NA

summary(incl_emot$emot_shy1)


## STEP 0: MULTIPLE IMPUTATION ##

#Create FFCW2: all variables, with sum scores and individual items        
ffcw2 <- incl_emot %>% 
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
    emot_shy1, #SA-emotionality and shyness
    
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


## check that all alright with ffcw3

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


ffcw3[c("cm1age", "cognit3_mother", "neighbour3pcg", "nce5", "mathard4m", "emot_shy1", "sum_conflict4", "sum_no_violent", "PPVT5raw", "PPVT9raw")] <- lapply(ffcw3[c("cm1age", "cognit3_mother", "neighbour3pcg", "nce5", "mathard4m", "emot_shy1", "sum_conflict4", "sum_no_violent", "PPVT5raw", "PPVT9raw")], as.numeric)

library(skimr)
skim(ffcw3) #--> FFCW3 READY FOR THE IMPUTATION (check that variables "type" makes sense)

summary(ffcw3$emot_shy1)

#multiple imputation 
library(mice)

# 1. Precompute interaction variables where both inputs are observed
ffcw3$mh.va    <- with(ffcw3, (mathard4m - 1.14) * (PPVT5raw - 62.2))
ffcw3$nce.va   <- with(ffcw3, (nce5 - 15.4) * (PPVT5raw - 62.2)) #updated to reflect the new (corrected) mean value
ffcw3$malt.va  <- with(ffcw3, (sum_conflict4 - 39.6) * (PPVT5raw - 62.2))

# 2. Initialize mice for method and predictorMatrix setup
init <- mice(ffcw3, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

# 3. Exclude idnum and passive variables from being predictors
predM[, "idnum"] <- 0

# Set up passive imputation formulas
meth["mh.va"]    <- "~I((mathard4m - 1.14)*(PPVT5raw - 62.2))"
meth["nce.va"]   <- "~I((nce5 - 15.4)*(PPVT5raw - 62.2))"
meth["malt.va"]  <- "~I((sum_conflict4 - 39.6)*(PPVT5raw - 62.2))"

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

# 4. Run the imputation with 10 iterations
imputed_10 <- mice(ffcw3, maxit = 10,
                   method = meth, 
                   predictorMatrix = predM, 
                   m = 20, 
                   seed = 123)
plot(imputed_10)

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

#3. Drop auxiliary variables, loneliness items, and hsgrades22 (didnt but doesnt matter)

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

#5. ## dichotomous educational attainment outcome
#starting from imp_10it 
summary(imp_10it$edu_att22) #1=less than high school, 2=high school or equivalent, 3=some college or technical education, and 4=completed college or graduate school.
imp_10it$edu_nonacad <- recode(imp_10it$edu_att22, "1" = "1", "2" = "1", "3" = "0", "4" = "0") #non-academic education vs HS max 

summary(imp_10it$edu_nonacad)

imp_mids <- as.mids(imp_10it) 

#### RQ1 ###

##  BINOMIAL MODELS ## 
# Create a list of 20 imputed datasets + add variables
imp_list <- lapply(1:20, function(i) {
  data_i <- complete(imp_mids, i)
  data_i$mid_high <- ifelse(data_i$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
  data_i$high_only <- ifelse(data_i$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)
  data_i
})

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
    "sex_child1", "emot_shy1", "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
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
    "sex_child1", "emot_shy1", "race_merged", "lowbbweight", "cm1age", "race1_mother", "mborn",
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
sa.shy_mh <- compare_binomial_models(imp_list, "mh.scale")
print(sa.shy_mh$comparison)

#2. Predictor: NCE
sa.shy_nce <- compare_binomial_models(imp_list, "nce.scale")
print(sa.shy_nce$comparison)

#3. Predictor: Child maltreatment
sa.shy_chmalt <- compare_binomial_models(imp_list, "chmalt.scale")
print(sa.shy_chmalt$comparison)

#4. All SEF as predictors 
sa.shy_sef <- compare_binomial_int(imp_list, predictors = c("mh.scale", "nce.scale", "chmalt.scale"))
print(sa.shy_sef$comparison)

#5. Predictor: Verbal Ability
sa.shy_ppvt <- compare_binomial_models(imp_list, "ppvt.scale")
print(sa.shy_ppvt$comparison)

#6. Predictor: MH, VA, MH x VA
sa.shy_mh.va <- compare_binomial_int(imp_list, predictors = c("mh.scale", "ppvt.scale", "mh.scale:ppvt.scale")) #misses the interaction in the output
print(sa.shy_mh.va$comparison)

#7. Predictor: NCE, VA, NCE x VA
sa.shy_nce.va <- compare_binomial_int(imp_list, predictors = c("nce.scale", "ppvt.scale", "nce.scale:ppvt.scale"))
print(sa.shy_nce.va$comparison)

#8. Predictor CH MALT, VA, CH MALT x VA
sa.shy_chmalt.va <- compare_binomial_int(imp_list, predictors = c("chmalt.scale", "ppvt.scale", "chmalt.scale:ppvt.scale"))
print(sa.shy_chmalt.va$comparison)

#Table- Supplementary material 

library(flextable)
library(officer)
library(dplyr)

fmt_or <- function(or, low, high) sprintf("%.2f (%.2f - %.2f)", or, low, high)

make_rows <- function(df, labels) {
  data.frame(
    Predictor = labels,
    OR1 = fmt_or(df$or_mh, df$or_mh_low, df$or_mh_high),
    OR2 = fmt_or(df$or_ho, df$or_ho_low, df$or_ho_high),
    stringsAsFactors = FALSE
  )
}

model_list <- list(
  "Model 1" = make_rows(sa.shy_mh$comparison,      c("Material hardship")),
  "Model 2" = make_rows(sa.shy_nce$comparison,     c("Neighborhood Collective Efficacy")),
  "Model 3" = make_rows(sa.shy_chmalt$comparison,  c("Child maltreatment")),
  "Model 4" = make_rows(sa.shy_sef$comparison,     c("Material Hardship", "Neighborhood Collective Efficacy", "Child maltreatment")),
  "Model 5" = make_rows(sa.shy_ppvt$comparison,    c("Verbal ability (PPVT)")),
  "Model 6" = make_rows(sa.shy_mh.va$comparison,   c("Material hardship", "Verbal ability (PPVT)", "Interaction: material hardship x verbal ability")),
  "Model 7" = make_rows(sa.shy_nce.va$comparison,  c("Neighborhood Collective Efficacy", "Verbal ability (PPVT)", "Interaction: neighbourhood x verbal ability")),
  "Model 8" = make_rows(sa.shy_chmalt.va$comparison, c("Child maltreatment", "Verbal ability (PPVT)", "Interaction: child maltreatment x verbal ability"))
)

table_df <- do.call(rbind, lapply(names(model_list), function(m) {
  data.frame(Model = m, model_list[[m]], stringsAsFactors = FALSE)
}))
rownames(table_df) <- NULL

ft <- flextable(table_df) %>%
  set_header_labels(
    Model = "Model",
    Predictor = "Predictors",
    OR1 = "Odds ratio (Recurrent/Transient vs Not reported)",
    OR2 = "Odds Ratio (Recurrent vs Transient/Not reported)"
  ) %>%
  merge_v(j = "Model") %>%
  theme_booktabs() %>%
  align(j = c("OR1", "OR2"), align = "center", part = "all") %>%
  valign(j = "Model", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  autofit()

ft


## RQ2 

#https://iqss.github.io/dss-ps/assessing.html 

#weighting the imputed data
library(MatchThem)
library(WeightIt)
library(cobalt)
library(marginaleffects)

w.imp_gbm <- weightthem(lon_3l ~ sex_child1 + emot_shy1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                        + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                        data = imp_mids, approach = 'within', method = "gbm", 
                        estimand = "ATE")


bal.tab(w.imp_gbm)

?love.plot
love.plot(w.imp_gbm, binary = "std")

new.names1 <- c(sex_child1 = "Sex (Boy/Girl)",
                emot_shy1 = "Child's temperament",
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



#var.names in stats. analysis (main) script

love.plot(w.imp_gbm, 
          drop.distance = TRUE, 
          #var.order = "unadjusted",
          abs = TRUE,
          binary = "std",
          line = TRUE, 
          thresholds = c(m = .1),
          
          var.names = new.names1,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"))


## Outcome models ##

#1. low high-school grades
sa_fits3 <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  
  glm_weightit(hs_bin ~ lon_3l + emot_shy1 + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family= binomial, vcov = "HC0") 
})

comp.shy3 <- lapply(sa_fits3, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
})

comp.shy3[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
comp.shy3 #Useful to see! group=outcome level

pooled.try3 <- mice::pool(comp.shy3, dfcom = Inf)
summary(pooled.try3, conf.int = TRUE, exponentiate = T)


#2. out-of-school suspension 
sa_fits2 <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  
  glm_weightit(susp_bin ~ lon_3l + emot_shy1 + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0") 
}) 


#Contrasts
comp.shy2 <- lapply(sa_fits2, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
})


comp.shy2[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
comp.shy2 #Useful to see! group=outcome level

pooled.try2 <- mice::pool(comp.shy2, dfcom = Inf)
summary(pooled.try2, conf.int = TRUE, exponentiate = T) #suspension


#3. Educational attainment 
sa_fits <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  
  glm_weightit(edu_nonacad ~ lon_3l + emot_shy1 + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, weightit = W, family = binomial, vcov = "HC0") 
}) 

#Difference 
comp.shy <- lapply(sa_fits, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg"  #Note: we need the log risk ratio because Rubin’s pooling rules don’t apply to the risk ratio but do to the log risk ratio. 
                  #We will exponentiate the log risk ratio and its confidence interval after pooling.
  ) 
  
})


comp.shy[[1]] #to check which are the contrast --> 2-1 / 3-1 / 3-2
comp.shy #Useful to see! group=outcome level

pooled.try <- mice::pool(comp.shy, dfcom = Inf)
summary(pooled.try, conf.int = TRUE, exponentiate = T) #educational attainment


# Table - Suppl Material

fmt_rr <- function(estimate, low, high) sprintf("%.2f (%.2f - %.2f)", estimate, low, high)

make_ate_rows <- function(pooled_summary) {
  s <- summary(pooled_summary, conf.int = TRUE, exponentiate = TRUE)
  data.frame(
    Comparison = c("Transient vs. not reported",
                   "Recurrent vs. not reported",
                   "Recurrent vs. transient"),
    RR = fmt_rr(s$estimate, s$conf.low, s$conf.high),
    stringsAsFactors = FALSE
  )
}

outcome_list <- list(
  "High school grades"        = make_ate_rows(pooled.try3),
  "Out-of-school suspension"  = make_ate_rows(pooled.try2),
  "Educational attainment"    = make_ate_rows(pooled.try)
)

table_df <- do.call(rbind, lapply(names(outcome_list), function(outcome) {
  data.frame(Outcome = outcome, outcome_list[[outcome]], stringsAsFactors = FALSE)
}))
rownames(table_df) <- NULL

table_df$Outcome <- factor(table_df$Outcome, levels = names(outcome_list))
table_df <- table_df[order(table_df$Outcome), ]
table_df$Outcome <- as.character(table_df$Outcome)

ft <- flextable(table_df) %>%
  set_header_labels(
    Outcome = "Outcome",
    Comparison = "Pairwise comparison",
    RR = "Estimate and CIs \u2013 Risk Ratio"
  ) %>%
  merge_v(j = "Outcome") %>%
  theme_booktabs() %>%
  align(j = "RR", align = "center", part = "all") %>%
  valign(j = "Outcome", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(j = "Outcome", part = "body") %>%
  bold(part = "header") %>%
  autofit()

ft


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


#1. Low high-school grades = mostly Bs or lower 
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


#2. Low high-school grades = halfs As/Bs or lower 
fits3A <- lapply(seq_along(w.imp_gbmAB$models), function(i) { #shouldnt it be w.imp_bm
  data <- complete(w.imp_gbmAB, i)
  
  W <- w.imp_gbmAB$models[[i]]
  
  
  
  glm_weightit(hs_binA ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family= binomial, vcov = "HC0") 
  
})

#contrasts
comp.imp3A <- lapply(fits3A, function(fit) {
  avg_comparisons(fit,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnoravg"
  ) 
  
})

pooled.try3A <- mice::pool(comp.imp3A, dfcom = Inf)
summary(pooled.try3A, conf.int = TRUE, exponentiate = T)



#4. SA- Complete case analysis

ffcw4 <- ffcw3 %>% dplyr::select(-c(PPVT9raw, edu_aux22, neighbour3pcg, sum_no_violent))
summary(ffcw4)
complete <- ffcw4[complete.cases(ffcw4), ]
summary(complete) #1266 observations 

# DATA wrangling #

#1. Dichotomize HS grades
complete$hs_bin <- complete$hsgrades22 
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


### BINOMIAL MODELS ###

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

library(dplyr)
library(broom)

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

library(dplyr)
library(broom)

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


#1.Predictor: Material Hardship
result_mh.com <- compare_binomial_models("mh.scale")
print(result_mh.com)

#2. Predictor: NCE
result_nce <- compare_binomial_models("nce.scale")
print(result_nce$comparison)

#3. Predictor: Child maltreatment
result_chmalt <- compare_binomial_models("chmalt.scale")
print(result_chmalt$comparison)

#4. All SEF as predictors 
result_sef <- compare_binomial_int(predictors = c("mh.scale", "nce.scale", "chmalt.scale"))
print(result_sef$comparison)

#5. Predictor: Verbal Ability
result_ppvt <- compare_binomial_models("ppvt.scale")
print(result_ppvt$comparison)

#6. Predictor: MH, VA, MH x VA
result_mh.va <- compare_binomial_int(predictors = c("mh.scale", "ppvt.scale", "mh.scale:ppvt.scale")) #misses the interaction in the output
print(result_mh.va$comparison)

#7. Predictor: NCE, VA, NCE x VA
result_nce.va <- compare_binomial_int(predictors = c("nce.scale", "ppvt.scale", "nce.scale:ppvt.scale"))
print(result_nce.va$comparison)

#8. Predictor CH MALT, VA, CH MALT x VA
result_chmalt.va <- compare_binomial_int(predictors = c("chmalt.scale", "ppvt.scale", "chmalt.scale:ppvt.scale"))
print(result_chmalt.va$comparison)

#Table- Supplementary material
fmt_or <- function(or, low, high) sprintf("%.2f (%.2f - %.2f)", or, low, high)

make_rows <- function(df, labels) {
  data.frame(
    Predictor = labels,
    OR1 = fmt_or(df$or_mh, df$or_mh_low, df$or_mh_high),
    OR2 = fmt_or(df$or_ho, df$or_ho_low, df$or_ho_high),
    stringsAsFactors = FALSE
  )
}

model_list <- list(
  "Model 1" = make_rows(result_mh$comparison,      c("Material hardship")),
  "Model 2" = make_rows(result_nce$comparison,     c("Neighborhood Collective Efficacy")),
  "Model 3" = make_rows(result_chmalt$comparison,  c("Child maltreatment")),
  "Model 4" = make_rows(result_sef$comparison,     c("Material Hardship", "Neighborhood Collective Efficacy", "Child maltreatment")),
  "Model 5" = make_rows(result_ppvt$comparison,    c("Verbal ability (PPVT)")),
  "Model 6" = make_rows(result_mh.va$comparison,   c("Material hardship", "Verbal ability (PPVT)", "Interaction: material hardship x verbal ability")),
  "Model 7" = make_rows(result_nce.va$comparison,  c("Neighborhood Collective Efficacy", "Verbal ability (PPVT)", "Interaction: neighbourhood x verbal ability")),
  "Model 8" = make_rows(result_chmalt.va$comparison, c("Child maltreatment", "Verbal ability (PPVT)", "Interaction: child maltreatment x verbal ability"))
)

table_df <- do.call(rbind, lapply(names(model_list), function(m) {
  data.frame(Model = m, model_list[[m]], stringsAsFactors = FALSE)
}))
rownames(table_df) <- NULL

ft <- flextable(table_df) %>%
  set_header_labels(
    Model = "Model",
    Predictor = "Predictors",
    OR1 = "Odds ratio (Recurrent/Transient vs Not reported)",
    OR2 = "Odds Ratio (Recurrent vs Transient/Not reported)"
  ) %>%
  merge_v(j = "Model") %>%
  theme_booktabs() %>%
  align(j = c("OR1", "OR2"), align = "center", part = "all") %>%
  valign(j = "Model", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  autofit()

ft  


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


#Outcome model


#1. low high-school grades
hs_com <- glm_weightit(hs_bin ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                       + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                       data = complete,
                       weightit = W, family = binomial) 

#contrast
hs_ATE <- 
  avg_comparisons(hs_com,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg")

hs_ATE 

#average and exp
hs_ATE_summary <- hs_ATE %>%
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

hs_ATE_summary



#2. out-of-school suspension
susp_com <- glm_weightit(susp_bin ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                         + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                         data = complete,
                         weightit = W, family = binomial) 



#contrast
susp_ATE <- 
  avg_comparisons(susp_com,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg")

susp_ATE 

#avergae and exp
susp_ATE_summary <- susp_ATE %>%
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

susp_ATE_summary

#3. educational attainment
att_com <- glm_weightit(edu_nonacad ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
                        + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + ppvt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
                        data = complete, weightit = W, family = binomial, vcov = "HC0") 



#contrast
att_ATE <- 
  avg_comparisons(att_com,
                  variables = list(lon_3l = "pairwise"),
                  comparison = "lnratioavg")

att_ATE 


#Aggregate to get one ATE per exposure contrast 

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

# Table- Supplementary material

fmt_rr <- function(estimate, low, high) sprintf("%.2f (%.2f - %.2f)", estimate, low, high)

make_ate_rows <- function(ate_summary) {
  data.frame(
    Comparison = c("Transient vs. not reported",
                   "Recurrent vs. not reported",
                   "Recurrent vs. transient"),
    RR = fmt_rr(ate_summary$RR, ate_summary$conf.lowRR, ate_summary$conf.highRR),
    stringsAsFactors = FALSE
  )
}

outcome_list <- list(
  "High school grades"        = make_ate_rows(hs_ATE_summary),
  "Out-of-school suspension"  = make_ate_rows(susp_ATE_summary),
  "Educational attainment"    = make_ate_rows(att_ATE_summary)
)

table_df <- do.call(rbind, lapply(names(outcome_list), function(outcome) {
  data.frame(Outcome = outcome, outcome_list[[outcome]], stringsAsFactors = FALSE)
}))
rownames(table_df) <- NULL

table_df$Outcome <- factor(table_df$Outcome, levels = names(outcome_list))
table_df <- table_df[order(table_df$Outcome), ]
table_df$Outcome <- as.character(table_df$Outcome)

ft <- flextable(table_df) %>%
  set_header_labels(
    Outcome = "Outcome",
    Comparison = "Pairwise comparison",
    RR = "Estimate and CIs \u2013 Risk Ratio"
  ) %>%
  merge_v(j = "Outcome") %>%
  theme_booktabs() %>%
  align(j = "RR", align = "center", part = "all") %>%
  valign(j = "Outcome", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(j = "Outcome", part = "body") %>%
  bold(part = "header") %>%
  autofit()

ft

# 5. SA - Multiple Imputation from the original sample
# --> very similar to main analysis, but selecting the analytical sample AFTER running multiple imputation

## 0. Data download and loading ##

#Data can be downloaded from https://doi.org/10.3886/ICPSR31622.v5 after creating a ICPSR account
#Start the download (by opening the “Download” tab and selecting “R”; this will download a folder called ICPSR_31622 containing all documents (data and documentation files)
#Inside the ICPSR_31622 folder, you'll find folder DS0001, which contains file 31622-0001-Data.rda (this is the R file you need)

#Once you have loaded the R file 31622-0001-Data.rda:
#Select and rename variables from ICPSR_31622 (da31622.0001 dataset) (…)

datav5 <- da31622.0001 %>% select(idnum=IDNUM, p4l5=P4L5, p5q3k=P5Q3K, hsgrades22=K7B13, edu_att22=CK7EDU, ck7ethrace=CK7ETHRACE, cp7kedu=CP7KEDU, ck6ethrace=CK6ETHRACE,
                                  susp=K7B35B, susp4=K7B39_4, susp5=K7B39_5, susp6=K7B39_6, susp7=K7B39_7, susp8=K7B39_8, susp9=K7B39_9, susp10=K7B39_10, susp11=K7B39_11, susp12=K7B39_12, k7b39_3=K7B39_3, k7b39_2=K7B39_2, k7b39_1=K7B39_1, k7b39_0=K7B39_0,
                                  m3b2=M3B2, m3b27=M3B27, m4i0m1=M4I0M1, m4i0m2=M4I0M2, m4i0m3=M4I0M3, m4i0m4=M4I0M4, m4i0m5=M4I0M5, m4i0n1=M4I0N1, m4i0n2=M4I0N2, m4i0n3=M4I0N3, m4i0n4=M4I0N4, f4i0m1=F4I0M1, f4i0m2=F4I0M2, 
                                  f4i0m3=F4I0M3, f4i0m4=F4I0M4, f4i0m5=F4I0M5, f4i0n1=F4I0N1, f4i0n2=F4I0N2, f4i0n3=F4I0N3, f4i0n4=F4I0N4, m4i23d=M4I23D, m4i23e=M4I23E, m4i23f=M4I23F, m4i23k=M4I23K, m4i23a=M4I23A, m4i23i=M4I23I, 
                                  m4i23j=M4I23J, m4i23h=M4I23H, m4i23b=M4I23B, m4i23c=M4I23C, m4i23g=M4I23G, p3k2d=P3K2D, p3k2e=P3K2E, p4g3=P4G3, p4g4=P4G4, p4g6=P4G6, p4g7=P4G7, p4g8=P4G8, p4g9=P4G9, p4g10=P4G10, p4g11=P4G11, 
                                  p4g13=P4G13, p4g14=P4G14, p4g15=P4G15, p4g16=P4G16, p4g17=P4G17, p4g18=P4G18, p4g19=P4G19, p4g1=P4G1, p4g5=P4G5, p4g12=P4G12, p4g2=P4G2, cm1marf=CM1MARF, cm1cohf=CM1COHF, m3j1=M3J1, cm3cogsc=CM3COGSC, 
                                  cm1bsex=CM1BSEX, cm1lbw=CM1LBW, cm1age=CM1AGE, cm1ethrace=CM1ETHRACE, m1h2=M1H2, cm3edu=CM3EDU, cm3md_case_con=CM3MD_CASE_CON, cm3alc_case=CM3ALC_CASE, cm3drug_case=CM3DRUG_CASE, cm3gad_case=CM3GAD_CASE, 
                                  p3a2=P3A2, p3k1a=P3K1A, p3k1b=P3K1B, p3k1c=P3K1C, p3k1d=P3K1D, p3k1e=P3K1E, p3k2a=P3K2A, p3k2b=P3K2B, p3k2c=P3K2C, p3k2d=P3K2D, p3k2e=P3K2E, ch4ppvtraw=CH4PPVTRAW, ch5ppvtraw=CH5PPVTRAW)

mydata <- datav5

mydata <- mydata %>%
  mutate(across(everything(), ~ as.numeric(sub("^\\(([-0-9]+)\\).*", "\\1", as.character(.)))))

### 1.FFCW Data preparation ##

#Replace with NA
mydata[mydata < 0] <- NA #Replace all negative values (not in wave, skipped, not asked, refused, etc) with NA

#Data wrangling of several variables (done for the whole sample, while in the main analysis only for the selected sample):
#merge the two child health status measures to try to reduce NAs
summary(mydata$m3b2)
summary(mydata$m3b27)
mydata<- mydata %>% mutate(childhealth3 = ifelse(is.na(m3b2), m3b27, m3b2))
summary(mydata$childhealth3) #15 NAs less

#merge race/ethnicity variables reported at year 15 and 22
mydata<- mydata %>% mutate(race_merged = ifelse(is.na(ck7ethrace), ck6ethrace, ck7ethrace))
summary(mydata$race_merged)


##Create sum scores for certain scales
#1. MATERIAL HARDSHIP #Year 5: 11 items  
mathard_items <- c("m4i23d", "m4i23e", "m4i23f", "m4i23k", "m4i23a", "m4i23i", "m4i23j", "m4i23h", "m4i23b", "m4i23c", "m4i23g") #Year 5: 11 items 

# Recode the value 2 to 0 (reversed coding)
mydata <- mydata %>%
  mutate_at(vars(one_of(mathard_items)), ~ ifelse(. == 2, 0, .)) #now its 0=no, yes=1 


#2. NCE SCALE #Year 5: 9 items in total, scores from 1 to 4 
reverse_cols1 = c('f4i0n3', 'f4i0n4', 'm4i0n3', 'm4i0n4') #year 5 (4pt scale) #previous error, was m3, m4, instead of n3, n4, now correct!
mydata[ , reverse_cols1] = 5 - mydata[ , reverse_cols1]

#year 3 (5pt scale) as auxiliary variable, scores from 1 to 5
reverse_cols2 = c('p3k2d', 'p3k2e') #correct items
mydata[ , reverse_cols2] = 6 - mydata[ , reverse_cols2]

#Replace NAs in NCE5 for mother's reports with data reported by fathers if available
mydata <- mydata %>% mutate(nce5_1= ifelse(is.na(m4i0m1), f4i0m1, m4i0m1))
mydata <- mydata %>% mutate(nce5_2= ifelse(is.na(m4i0m2), f4i0m2, m4i0m2))
mydata <- mydata %>% mutate(nce5_3= ifelse(is.na(m4i0m3), f4i0m3, m4i0m3))
mydata <- mydata %>% mutate(nce5_4= ifelse(is.na(m4i0m4), f4i0m4, m4i0m4))
mydata <- mydata %>% mutate(nce5_5= ifelse(is.na(m4i0m5), f4i0m5, m4i0m5))

mydata <- mydata %>% mutate(nce5_6= ifelse(is.na(m4i0n1), f4i0n1, m4i0n1))
mydata <- mydata %>% mutate(nce5_7= ifelse(is.na(m4i0n2), f4i0n2, m4i0n2))
mydata <- mydata %>% mutate(nce5_8= ifelse(is.na(m4i0n3), f4i0n3, m4i0n3)) #reverse-coded
mydata <- mydata %>% mutate(nce5_9= ifelse(is.na(m4i0n4), f4i0n4, m4i0n4)) #reverse-coded

#These items are summed when creating ffcw2 object

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

## 2. MULTIPLE IMPUTATION BY CHAINED EQUATIONS ## - INCLUDING ALL SAMPLE N=4898
library(mice)

# 1. Precompute interaction variables where both inputs are observed
ffcw3$mh.va    <- with(ffcw3, (mathard4m - 1.14) * (PPVT5raw - 62.2))
ffcw3$nce.va   <- with(ffcw3, (nce5 - 15.4) * (PPVT5raw - 62.2)) #updated to reflect the new (corrected) mean value
ffcw3$malt.va  <- with(ffcw3, (sum_conflict4 - 39.6) * (PPVT5raw - 62.2))

# 2. Initialize mice for method and predictorMatrix setup
init <- mice(ffcw3, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

# 3. Exclude idnum and passive variables from being predictors
predM[, "idnum"] <- 0

# Set up passive imputation formulas
meth["mh.va"]    <- "~I((mathard4m - 1.14)*(PPVT5raw - 62.2))"
meth["nce.va"]   <- "~I((nce5 - 15.4)*(PPVT5raw - 62.2))"
meth["malt.va"]  <- "~I((sum_conflict4 - 39.6)*(PPVT5raw - 62.2))"

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

#Select the study sample N=2456

##Inclusion criteria: 
#Select participants with loneliness complete data 
complete_lon <- mydata %>%
  filter(!is.na(p4l5) & !is.na(p5q3k)) #n=2466 

#Select those who attended a graded high school
complete_lon$hsgrades22<-as.factor(complete_lon$hsgrades22)
complete_lon <- complete_lon %>% filter(hsgrades22 != 9 | is.na(hsgrades22)) #10 people bye
compout<- complete_lon

imp_wei <- imp_10it #after data wrangling
imp_wei$selected <- ifelse(imp_wei$idnum %in% compout$idnum, 1, 0) #Create selected variable

imp_org <- imp_wei %>% filter(selected ==1)
imp_10it <- imp_org

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

#5. dichotomous educational attainment outcome
#starting from imp_10it 
summary(imp_10it$edu_att22) #1=less than high school, 2=high school or equivalent, 3=some college or technical education, and 4=completed college or graduate school.
imp_10it$edu_nonacad <- recode(imp_10it$edu_att22, "1" = "1", "2" = "1", "3" = "0", "4" = "0") #non-academic education vs HS max 

summary(imp_10it$edu_nonacad)

imp_mids <- as.mids(imp_10it) #I need a mids object (important to do it after all postMICE data wrangling, so that all variables are included and in the correct form)

## STATISTICAL ANALYSIS ##

# RQ1 (imputation from original sample-SA) #

# Create a list of 20 imputed datasets + add variables
imp_list <- lapply(1:20, function(i) {
  data_i <- complete(imp_mids, i)
  data_i$mid_high <- ifelse(data_i$lon_3l == 0, 0, 1)  # 2 or 1 vs 0 (high/mid loneliness vs low loneliness)
  data_i$high_only <- ifelse(data_i$lon_3l == 2, 1, 0) # 2 vs 1 or 0 (high loneliness vs mid/low loneliness)
  data_i
})

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

## TABLE s13

library(flextable)
library(officer)
library(dplyr)

fmt_or <- function(or, low, high) sprintf("%.2f (%.2f - %.2f)", or, low, high)

make_rows <- function(df, labels) {
  data.frame(
    Predictor = labels,
    OR1 = fmt_or(df$or_mh, df$or_mh_low, df$or_mh_high),
    OR2 = fmt_or(df$or_ho, df$or_ho_low, df$or_ho_high),
    stringsAsFactors = FALSE
  )
}

model_list <- list(
  "Model 1" = make_rows(result_mh$comparison,      c("Material hardship")),
  "Model 2" = make_rows(result_nce$comparison,     c("Neighborhood Collective Efficacy")),
  "Model 3" = make_rows(result_chmalt$comparison,  c("Child maltreatment")),
  "Model 4" = make_rows(result_sef$comparison,     c("Material Hardship", "Neighborhood Collective Efficacy", "Child maltreatment")),
  "Model 5" = make_rows(result_ppvt$comparison,    c("Verbal ability (PPVT)")),
  "Model 6" = make_rows(result_mh.va$comparison,   c("Material hardship", "Verbal ability (PPVT)", "Interaction: material hardship x verbal ability")),
  "Model 7" = make_rows(result_nce.va$comparison,  c("Neighborhood Collective Efficacy", "Verbal ability (PPVT)", "Interaction: neighbourhood x verbal ability")),
  "Model 8" = make_rows(result_chmalt.va$comparison, c("Child maltreatment", "Verbal ability (PPVT)", "Interaction: child maltreatment x verbal ability"))
)

table_df <- do.call(rbind, lapply(names(model_list), function(m) {
  data.frame(Model = m, model_list[[m]], stringsAsFactors = FALSE)
}))
rownames(table_df) <- NULL

ft <- flextable(table_df) %>%
  set_header_labels(
    Model = "Model",
    Predictor = "Predictors",
    OR1 = "Odds ratio (Recurrent/Transient vs Not reported)",
    OR2 = "Odds Ratio (Recurrent vs Transient/Not reported)"
  ) %>%
  merge_v(j = "Model") %>%
  theme_booktabs() %>%
  align(j = c("OR1", "OR2"), align = "center", part = "all") %>%
  valign(j = "Model", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  autofit()

ft  #preview & copy to Clipboard from here


## RQ2 

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
          var.names = new.names,
          colors = c("red", "blue"),
          shapes = c("triangle filled", "circle filled"),
          sample.names = c("Unweighted", "PS Weighted (GBM)"))



## Outcome model: weighted g-computation 

#Outcome 1: low high-school grades

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

#Outcome 3: educational attainment
fits <- lapply(seq_along(w.imp_gbm$models), function(i) {
  data <- complete(w.imp_gbm, i)
  W <- w.imp_gbm$models[[i]]
  
  glm_weightit(edu_nonacad ~ lon_3l + sex_child1 + race_merged + lowbbweight + cm1age + race1_mother + mborn + relst1 + mother_edu3 + cognit3_mother
               + cm3md_case_con + cm3alc_case + cm3drug_case + cm3gad_case + m_health3 + health3 + disab3 + mh.scale + nce.scale + chmalt.scale + mh_ppvt + nce_ppvt + chmalt_ppvt,
               data = data, #loop, from before
               weightit = W, family = binomial, vcov = "HC0") 
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

## Table s14

fmt_rr <- function(estimate, low, high) sprintf("%.2f (%.2f - %.2f)", estimate, low, high)

make_ate_rows <- function(pooled_summary) {
  s <- summary(pooled_summary, conf.int = TRUE, exponentiate = TRUE)
  data.frame(
    Comparison = c("Transient vs. not reported",
                   "Recurrent vs. not reported",
                   "Recurrent vs. transient"),
    RR = fmt_rr(s$estimate, s$conf.low, s$conf.high),
    stringsAsFactors = FALSE
  )
}

outcome_list <- list(
  "High school grades"        = make_ate_rows(pooled.try3),
  "Out-of-school suspension"  = make_ate_rows(pooled.try2),
  "Educational attainment"    = make_ate_rows(pooled.try)
)

table_df <- do.call(rbind, lapply(names(outcome_list), function(outcome) {
  data.frame(Outcome = outcome, outcome_list[[outcome]], stringsAsFactors = FALSE)
}))
rownames(table_df) <- NULL

table_df$Outcome <- factor(table_df$Outcome, levels = names(outcome_list))
table_df <- table_df[order(table_df$Outcome), ]
table_df$Outcome <- as.character(table_df$Outcome)

ft <- flextable(table_df) %>%
  set_header_labels(
    Outcome = "Outcome",
    Comparison = "Pairwise comparison",
    RR = "Estimate and CIs \u2013 Risk Ratio"
  ) %>%
  merge_v(j = "Outcome") %>%
  theme_booktabs() %>%
  align(j = "RR", align = "center", part = "all") %>%
  valign(j = "Outcome", valign = "top", part = "body") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(j = "Outcome", part = "body") %>%
  bold(part = "header") %>%
  autofit()

ft  #visualize


## THE END ##