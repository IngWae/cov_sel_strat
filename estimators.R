### SUPPORTING FUNCTIONS -------------------------------------------------------
# odds():
# Calculates the odds from a propprtion 
odds <- function(pt){
  mean(pt)/mean(1-pt)
}

# odds_ratio():
# Calculates the odds ratio from two proportions 
odds_ratio <- function(pt, pref){
  odds_t <- pt/(1-pt)
  odds_ref <- pref/(1-pref)
  return(odds_t/odds_ref)
}


##### PROPENSITY SCORES --------------------------------------------------------
# propensity_scores():
# Function to calculate the propensity score for a sample 
# Arguments are a sample, treatment (dmt/smoking) and type of propensity 
# score model (correct/incorrect/nonparametric)
propensity_scores <- function(samp, tr, type){
  # SCENARIO I
  # Correct model specification 
  if(tr == "dmt" & type == "correct"){
    temp_samp <- samp %>% dplyr::select(dmt, age, cardiovasc_dis, 
                                        pulmonary_dis, diabetes, 
                                        progressive_course, edss) %>%
      mutate(age_sqrt = sqrt(age),
             edss_sq = edss^2, 
             cv_edss = edss*cardiovasc_dis, 
             pd_edss = edss*pulmonary_dis, 
             dia_edss = edss*diabetes)
    ps_model <- multinom(dmt ~., data = temp_samp, trace=FALSE)
    return(ps_model$fitted.values)
  } 
  # Incorrect model specification 
  else if(tr == "dmt" & type == "incorrect"){
    temp_samp <- samp %>% dplyr::select(dmt, age, cardiovasc_dis, 
                                        pulmonary_dis, diabetes, 
                                        progressive_course, edss)
    ps_model <- multinom(dmt ~., data = temp_samp, trace=FALSE)
    return(ps_model$fitted.values)
  } 
  # Non parametric
  else if(tr == "dmt" & type == "nonparametric"){
    temp_samp <- data.frame(samp %>% 
                              dplyr::select(age, cardiovasc_dis, pulmonary_dis, 
                                            diabetes, progressive_course, 
                                            edss, dmt)) 
    
    ps_model <- ranger(dmt~., data =temp_samp, probability = T)
    
    ps <- predict(ps_model, type="response", data = temp_samp)$predictions
    
    return(ps)
  } 
  # SCENARIO II
  # Correct model specification 
  else if(tr=="smoking" & type == "correct"){
    temp_samp <- samp %>% 
      dplyr::select(smoking, age, male) %>%
      mutate(age_sq = age^2, male_age = male*age)
    ps_model <- glm(smoking ~ ., data = temp_samp, family = binomial())
    return(ps_model$fitted.values)
  }
  # Incorrect model specification 
  else if(tr == "smoking" & type == "incorrect"){
    temp_samp <- samp %>% dplyr::select(smoking, age, male)
    ps_model <- glm(smoking ~ age + male, data = temp_samp, family = binomial())
    return(ps_model$fitted.values)
  }
  # Non parametric
  else if(tr == "smoking" & type == "nonparametric"){
    temp_samp <- data.frame(samp %>% 
                              dplyr::select(male, age, smoking))
    
    temp_samp$smoking <- as.factor(temp_samp$smoking)
    
    ps_model <- ranger(smoking ~ ., data = temp_samp, probability = T)
    
    ps <- predict(ps_model, data = temp_samp, type="response")$predictions[,"1"]
    
    return(ps)
  }
}


#########    UNADJUSTED OR  ----------------------------------------------------
# unadjusted_OR():
# Function that calculates the unadjusted odds ratio in a given sample and  for 
# a given treatment variable (dmt/smoking)
unadjusted_OR <- function(samp, tr){
  
  # Function to calculate Pr(Y=1) for each treatment group
  p_tr <- function(filtered_temp){
    return(filtered_temp %>% dplyr::select(severity) %>% 
             summarize(mean=mean(severity))%>% pull())
  }
  
 
  # SCENARIO I 
  if(tr == "dmt"){
    temp <- samp %>% 
      dplyr::select(severity, No_treatment, No_risk, 
                    Low_risk, Moderate_to_high_risk)
    # Calculate Pr(Y=1) for each treatment group
    p_NT <- p_tr(temp %>% filter(No_treatment==1)) 
    p_NR <- p_tr(temp %>% filter(No_risk==1)) 
    p_LR <- p_tr(temp %>% filter(Low_risk==1))
    p_MHR <- p_tr(temp %>% filter(Moderate_to_high_risk==1))
    est <- c(odds_ratio(pt = p_NR, pref = p_NT), 
             odds_ratio(pt = p_LR, pref = p_NT),
             odds_ratio(pt = p_MHR, pref = p_NT))
    
    return(est)
    
  } 
  # SCENARIO II
  else if(tr == "smoking"){
    temp <- samp %>% dplyr::select(severity, smoking)
    # Calculate Pr(Y=1) for each treatment group
    p_ntsmok <- p_tr(temp %>% filter(smoking == 0)) 
    p_smok <- p_tr(temp %>% filter(smoking == 1))
    
    est <- c(odds_ratio(pt = p_smok, pref = p_ntsmok))
    return(est)
  }
}

########## ADJUST FOR ALL ----------------------------------------------------
# adjusted_for_all(): 
# Function that calculates the COR conditional on all available covariates in 
# a given sample and treatment variable (dmt/smoking)
adjust_for_all <- function(samp, tr){
  # SCENARIO I
  if(tr == "dmt"){
    # Fit a logistic regression including all covariates 
    covars <- c("age", "male", "cardiovasc_dis", "pulmonary_dis", "diabetes",
                "smoking", "bmi", "progressive_course", "edss")
    
    temp_formula <- as.formula(paste("severity ~ No_risk + Low_risk + Moderate_to_high_risk +", 
                                  paste(covars, collapse = "+")))
  
    temp_model <- do.call("glm", args = list(temp_formula, data = samp, family = binomial()))
    
    # Return the exponential function of the estimated coefficients 
    return(c(exp(coef(summary(temp_model))[2,1]), 
             exp(coef(summary(temp_model))[3,1]), 
             exp(coef(summary(temp_model))[4,1])))
    
  } 
  # SCENARIO II
  else if(tr == "smoking"){
    # Fit a logistic regression including all covariates 
    covars <- c("age", "male", "cardiovasc_dis", "pulmonary_dis", "diabetes",
                "bmi", "progressive_course", "edss",
                "No_risk", "Low_risk", "Moderate_to_high_risk")
    
    temp_formula <- as.formula(paste("severity ~ smoking + ", 
                                  paste(covars, collapse = "+")))
    temp_model <- do.call("glm", 
                          args = list(temp_formula, data = samp, family = binomial()))
    
    # Return the exponential function of the estimated coefficient
    return(exp(coef(summary(temp_model))[2,1]))
  }
}



##### UNIVARIABLE PREFILTERING -------------------------------------------------
# univariable_prefiltering(): 
# Function to calculate the COR through univariate prefiltering in a given sample 
# and for a given treatment variable (dmt/smoking)
univariable_prefiltering <- function(samp, tr){
  # SCENARIO I
  if(tr == "dmt"){
    # Covariates available for selection 
    covars <- c("age", "male",  "cardiovasc_dis", "pulmonary_dis", "diabetes",
                "smoking", "bmi", "progressive_course", "edss")
    # Select covariates that are significant in univariate logistic regression 
    selected_covars <- c()
    for(i in 1:length(covars)){
      temp_formula <- as.formula(paste("severity ~ ", paste(covars[i])))
      temp_model <- do.call("glm", args = list(temp_formula, data = samp, family = binomial()))
      if(coef(summary(temp_model))[2,4] < 0.05){
        selected_covars <- c(selected_covars, covars[i])
      }
    }
    
    # Fit logistic regression including selected covariates 
    temp_formula <- as.formula(paste("severity ~ No_risk + Low_risk + Moderate_to_high_risk +", 
                                  paste(selected_covars, collapse = "+")))
    
    temp_model <- do.call("glm", args = list(temp_formula, data = samp, family = binomial()))
    # Return the exponential function of the estimated coefficients
    return(c(exp(coef(summary(temp_model))[2,1]),  
             exp(coef(summary(temp_model))[3,1]), 
             exp(coef(summary(temp_model))[4,1])))
    
  } 
  # SCENARIO II 
  else if(tr == "smoking"){
    # Covariates available for selection 
    covars <- c("age", "male", "cardiovasc_dis", "pulmonary_dis", "diabetes",
                "bmi", "progressive_course","edss",
                "No_risk", "Low_risk", "Moderate_to_high_risk")
    # Select covariates that are significant in univariate logistic regression 
    selected_covars <- c()
    for(i in 1:length(covars)){
      temp_formula <- as.formula(paste("severity ~ ", paste(covars[i])))
      temp_model <- do.call("glm", args = list(temp_formula, data = samp, 
                                         family = binomial()))
      
      if(coef(summary(temp_model))[2,4] < 0.05){
        selected_covars <- c(selected_covars, covars[i])
      }
    }
    # Fit logistic regression including selected covariates 
    temp_formula <- as.formula(paste("severity ~ smoking + ", 
                                  paste(selected_covars, collapse = "+")))
    temp_model <- do.call("glm", args = list(temp_formula, data = samp, 
                                       family = binomial()))
    # Return the exponential function of the estimated coefficient
    return(exp(coef(summary(temp_model))[2,1]))
  }
}


###### STEPWISE SELECTION -------------------------------------------------
# stepwise_selection(): 
# Function to estimate the COR using stepwise regression in a given sample and 
# for a given treatment variable (dmt/smoking)
stepwise_selection <- function(samp, tr){
  # Sample with outcome, treatment and covariates avaiable for selection 
  temp_samp <- samp %>% dplyr::select(severity, No_risk, Low_risk, Moderate_to_high_risk,
                                      age, male, cardiovasc_dis, pulmonary_dis, diabetes,
                                      smoking, bmi, progressive_course, edss)
  # SCENARIO I 
  if(tr == "dmt"){
    # Base model includes only the treatment variable 
    base <- do.call("glm", args = list(severity ~ No_risk + Low_risk + Moderate_to_high_risk, 
                                       data = temp_samp, family = binomial()))
    # Full model includes all covariates 
    full <- do.call("glm", args = list(severity ~ ., 
                                       data = temp_samp, family = binomial()))
    # Stepwise selection in both directions 
    temp_model <- step(base, direction = "both", 
                       scope = list(lower=base, upper=full), trace = 0)
    # Return the exponential function of the estimated coefficients
    return(c(exp(coef(summary(temp_model))[2,1]), 
             exp(coef(summary(temp_model))[3,1]), 
             exp(coef(summary(temp_model))[4,1])))
    
  } 
  # SCENARIO II 
  else if(tr == "smoking"){
    # Base model includes only the treatment variable 
    base <- do.call("glm", args = list(severity ~ smoking, 
                                       data = temp_samp, family = binomial()))
    # Full model includes all the covariates 
    full <- do.call("glm", args = list(severity ~ ., data = temp_samp, family = binomial()))
    # Stepwise selection in both directions 
    temp_model <- step(base, direction = "both", scope = list(lower=base, upper=full), trace = 0)
    # Return the exponential function of the estimated coefficient
    return(exp(coef(summary(temp_model))[2,1]))
  }
}



#########    REGERESSION IMPUTATION --------------------------------------------
# parametric_reg_imputation(): 
# Function to calculate the adjusted MOR with regression imputation in a given sample and 
# for a given treatment variable (dmt/smoking)
parametric_reg_imputation <- function(samp, tr){
  
  # SCENARIO I 
  if(tr == "dmt"){
    
    x <- levels(samp[[tr]])
    
    # Sample including outcome, treatment and all covariates 
    temp_samp <- data.frame(samp %>% 
                              dplyr::select(severity, smoking, age, male, bmi, diabetes, 
                                            cardiovasc_dis, pulmonary_dis, progressive_course, 
                                            edss, No_treatment, No_risk, Low_risk, Moderate_to_high_risk))
    
    # Fit logistic regressions 
    po_model0 <- glm(severity ~ male + age + smoking + bmi + 
                       diabetes + cardiovasc_dis + pulmonary_dis + 
                       progressive_course + edss, 
                     data = temp_samp %>% filter(No_treatment==1), 
                     family = binomial())
    po_model1 <- glm(severity ~ male + age + smoking + bmi + 
                       diabetes + cardiovasc_dis + pulmonary_dis + 
                       progressive_course + edss, 
                     data = temp_samp %>% filter(No_risk==1), 
                     family = binomial())
    po_model2 <- glm(severity ~ male + age + smoking + bmi + 
                       diabetes + cardiovasc_dis + pulmonary_dis + 
                       progressive_course + edss, 
                     data = temp_samp %>% filter(Low_risk==1), 
                     family = binomial())
    po_model3 <- glm(severity ~ male + age + smoking + bmi + 
                       diabetes + cardiovasc_dis + pulmonary_dis + 
                       progressive_course + edss, 
                     data = temp_samp %>% filter(Moderate_to_high_risk==1), 
                     family = binomial())
    
    
    
    # Matrix for potential outcome predictions
    po <- matrix(nrow=nrow(temp_samp), ncol = length(x))
    colnames(po) <- x
    # Get predictions at each treatment level 
    
    po[,"No_treatment"] <- predict(object = po_model0,  newdata  = temp_samp, type = "response")
    po[,"No_risk"] <- predict(object = po_model1,  newdata  = temp_samp, type = "response")
    po[,"Low_risk"] <- predict(object = po_model2,  newdata  = temp_samp, type = "response")
    po[,"Moderate_to_high_risk"] <- predict(object = po_model3,  newdata  = temp_samp, type = "response")
    
    # Estimates of the marginal probabilities are the mean of the predictions 
    marg_prob <- colMeans(po, na.rm=TRUE)
    
    # Regression imputation estimates 
    or10 <- odds_ratio(marg_prob["No_risk"], marg_prob["No_treatment"])
    or20 <- odds_ratio(marg_prob["Low_risk"], marg_prob["No_treatment"])
    or30 <- odds_ratio(marg_prob["Moderate_to_high_risk"], marg_prob["No_treatment"])
    
    # Return odd ratio estimates 
    return(c(or10, or20, or30))
    
  } 
  # SCENARIO II 
  else if(tr == "smoking"){
    
    temp_samp <- data.frame(samp %>% 
                              dplyr::select(severity, smoking, age, male))
    
    # Fit a logistic regression 
    po_model0 <- glm(severity ~ age + male, data = temp_samp %>% filter(smoking==0), family = binomial())
    po_model1 <-  glm(severity ~ age + male, data = temp_samp %>% filter(smoking==1), family = binomial())
    
    po <- matrix(nrow=nrow(temp_samp), ncol = 2)
    colnames(po) <- c("0","1")
    
    po[,"0"] <- predict(object = po_model0,  newdata  = temp_samp, type = "response")
    po[,"1"] <- predict(object = po_model1,  newdata  = temp_samp, type = "response")
    
    # Estimates of marginal probabilities are the average predictions 
    marg_prob <- colMeans(po, na.rm=TRUE)
    
    # Return odd ratio estimate
    return(odds_ratio(marg_prob["1"], marg_prob["0"]))
  }
  
}


#########    AIPW -------------------------------------------------
# parametric_aipw(): 
# Function to calculate the adjusted MOR in a given sample and 
# for a given treatment variable (dmt/smoking). Requires estimated propensity scores. 
parametric_aipw <- function(ps, samp, tr){
  # SCENARIO I
  if(tr == "dmt"){
    # Sample including outcome, treatment and all covariates 
    temp_samp <- data.frame(samp %>% 
                              dplyr::select(severity, age, male, smoking, bmi, 
                                            cardiovasc_dis, pulmonary_dis, diabetes, 
                                            progressive_course, edss, dmt))
    
    yform <- as.formula(paste("severity ~ age + male + smoking + bmi + 
                              cardiovasc_dis + pulmonary_dis + diabetes + 
                              progressive_course + edss"))
    
    temp_model <- PSweight(ps.estimate = ps, weight = "IPW", 
                           yname = "severity", zname = "dmt", 
                           data = temp_samp, augmentation = T,
                           out.formula = yform, out.method = "glm", 
                           delta = 0.001,
                           family = "binomial")
    
    contrasts_mult <- rbind(c(0,0, 1, -1), c(1,0, 0, -1), c(0,1, 0, -1) )
    # Return odd ratio estimates 
    return(c(exp(summary(temp_model, type = "OR", contrast = contrasts_mult)$estimates[1,1]),
             exp(summary(temp_model, type = "OR", contrast = contrasts_mult)$estimates[2,1]),
             exp(summary(temp_model, type = "OR", contrast = contrasts_mult)$estimates[3,1])))
    
  } 
  # SCENARIO II
  else if(tr == "smoking"){
    # Add column for Pr(T=0|X)
    ps <- cbind(ps, 1-ps)
    colnames(ps) <- c("1", "0")
    
    # Sample including outcome, treatment and all covariates 
    temp_samp <- data.frame(samp %>% dplyr::select(severity, age, male, smoking, bmi, 
                                                   cardiovasc_dis, pulmonary_dis, diabetes, 
                                                   progressive_course, edss,
                                                   No_risk, Low_risk, Moderate_to_high_risk)) 
    # Outcome regression formula (excluding treatment)
    temp_formula <- as.formula(paste("severity ~ age + male"))
    # Fit model 
    temp_model <- PSweight(ps.estimate = ps, weight = "IPW",
                           yname = "severity", zname = "smoking",
                           data = temp_samp,  augmentation = T,
                           out.formula = temp_formula, out.method = "glm",
                           delta = 0.001,
                           family = "binomial")
    # Return odd ratio estimate
    return(c(exp(summary(temp_model, type = "OR", contrast = c(-1,1))$estimate[1])))
  }
}

#########    NONPARAMETRIC ESTIMATION ------------------------------------------
# aipw(): 
# Function to calculate the AIPW odds ratio in a given sample and 
# for a given treatment variable (dmt/smoking). Requires estimated propensity scores 
# and predictions of potential outcomes. 
aipw <- function(samp, ps, pot_out, tr){
  # Ratio formula 
  aipw_calc <- function(y, tr_bin, ps,  pot_out){
    aipw <- mean(y*tr_bin/ps - (tr_bin-ps)/ps*pot_out)
    return(aipw/(1-aipw))
    
  }
  # SCENARIO I 
  if(tr=="dmt"){
    # Propensity score trimming 
    trim <- which(!(rowSums(ps<0.00001)|rowSums(ps>1-0.00001)| rowSums(is.na(ps))))
    samp <- samp[trim,]
    pot_out <- pot_out[trim,]
    ps <- ps[trim,]
    
    # Outcome vector 
    y <- samp$severity
    # Reference and treatment levels 
    ref_level = "No_treatment"
    treatment_level = c("No_risk", "Low_risk", "Moderate_to_high_risk")
    # Binary variable of reference level 
    tr_bin <- samp[[ref_level]]
    # Propensity scores at reference level 
    ps_0 <- ps[,"No_treatment"]
    po_0 <- pot_out[,"No_treatment"]
    # Reference level constitutes the denominator in odds ratio 
    denom <- aipw_calc(y,tr_bin,ps_0,po_0)
    
    # Vector of estimates at each treatment level
    est <- c()
    k = 1 
    for(j in treatment_level){
      # Extract treatment variable, propensity scores and predicted potential outcomes
      # at treatment level j 
      tr_bin <- samp[[j]]
      ps_j <- ps[,j]
      po_j <- pot_out[,j]
      # Calculate odds ratio 
      est[k] <- aipw_calc(y,tr_bin,ps_j,po_j)/denom
      k = k + 1
    }
    
  }  
  # SCENARIO II 
  else if(tr == "smoking"){
    # Propensity score trimming
    trim <- which(!(ps<0.00001|ps>1-0.00001|is.na(ps)))
    samp <- samp[trim,]
    pot_out <- pot_out[trim,]
    ps <- ps[trim]
    
    # Outcome vector 
    y <- samp$severity
    # Treatment vector 
    tr_bin <- samp$smoking
    # Propensity score at reference level 
    ps_0 <- 1-ps
    po_0 <- pot_out[,"0"]
    # Reference level constitutes denominator 
    denom <- aipw_calc(y,tr_bin,ps_0,po_0)
    
    # Propensity scores at treatment level 
    ps_1 <- ps
    po_1 <- pot_out[,"1"]
    # Odds ratio estimate s
    est <- aipw_calc(y, tr_bin, ps_1, po_1)/denom
    
  }
  return(c(est))
}

# nonparametric_est():
# Function to calculate the nonparametric odds ratio estimates
nonparametric_est <- function(samp, ps, tr){
  # SCENARIO I 
  if(tr == "dmt"){
    # Treatment levels 
    x <- levels(samp[[tr]])
    # Remove treatment factor (keep only binary vectors)
    temp_samp <- samp %>% dplyr::select(-dmt )%>% 
      mutate(severity = as.factor(severity))
    
    temp_samp0 <- temp_samp %>% filter(No_treatment==1)
    temp_samp1 <- temp_samp %>% filter(No_risk==1)
    temp_samp2 <- temp_samp %>% filter(Low_risk==1)
    temp_samp3 <- temp_samp %>% filter(Moderate_to_high_risk==1)
    
    # Fit random forest 
    po_model0 <- ranger(severity ~ male + age + smoking + bmi + 
                          diabetes + cardiovasc_dis + pulmonary_dis + 
                          progressive_course + edss, data = temp_samp0, probability = T)
    
    po_model1 <- ranger(severity ~ male + age + smoking + bmi + 
                          diabetes + cardiovasc_dis + pulmonary_dis + 
                          progressive_course + edss, data = temp_samp1, probability = T)
    
    po_model2 <- ranger(severity ~ male + age + smoking + bmi + 
                          diabetes + cardiovasc_dis + pulmonary_dis + 
                          progressive_course + edss, data = temp_samp2, probability = T)
    
    po_model3 <- ranger(severity ~ male + age + smoking + bmi + 
                          diabetes + cardiovasc_dis + pulmonary_dis + 
                          progressive_course + edss, data = temp_samp3, probability = T)
    
    
    
    # Matrix for potential outcome predictions
    po <- matrix(nrow=nrow(temp_samp), ncol = length(x))
    colnames(po) <- x
    # Get predictions at each treatment level 
    
    po[,"No_treatment"] <- predict(object = po_model0,  data  = temp_samp, type = "response")$predictions[,"1"]
    po[,"No_risk"] <- predict(object = po_model1,  data  = temp_samp, type = "response")$predictions[,"1"]
    po[,"Low_risk"] <- predict(object = po_model2,  data  = temp_samp, type = "response")$predictions[,"1"]
    po[,"Moderate_to_high_risk"] <- predict(object = po_model3,  data  = temp_samp, type = "response")$predictions[,"1"]
    
    
    # Estimates of the marginal probabilities are the mean of the predictions 
    marg_prob <- colMeans(po, na.rm=TRUE)
    
    # Regression imputation estimates 
    or10 <- odds_ratio(marg_prob["No_risk"], marg_prob["No_treatment"])
    or20 <- odds_ratio(marg_prob["Low_risk"], marg_prob["No_treatment"])
    or30 <- odds_ratio(marg_prob["Moderate_to_high_risk"], marg_prob["No_treatment"])
    
    # Return regression imputation and AIPW estimates 
    return(c(or10, or20, or30,
             aipw(samp=samp%>% dplyr::select(-dmt), ps=ps, pot_out=po, tr = "dmt")))
    
  } 
  # SCENARIO II 
  else if(tr == "smoking"){
    # Remove DMT factor (keep only dummies) 
    temp_samp <- samp %>% dplyr::select(-dmt) %>% 
      mutate(severity = as.factor(severity))
    
    temp_samp1 <- temp_samp %>% filter(smoking==1)
    temp_samp0 <- temp_samp %>% filter(smoking==0)
    
    # Fit random forest 
    po_model1 <- ranger(severity ~  age + male, data = temp_samp1, probability = T)
    po_model0 <- ranger(severity ~  age + male, data = temp_samp0, probability = T)
    
    # Matrix for potential outcome predictions 
    po <- matrix(nrow=nrow(temp_samp), ncol = 2)
    colnames(po) <- c("0","1")
    
    po[,"0"] <- predict(object = po_model0,  data  = temp_samp, type = "response")$predictions[,"1"]
    po[,"1"] <- predict(object = po_model1,  data  = temp_samp, type = "response")$predictions[,"1"]
    
    # Estimates of marginal probabilities are the average predictions 
    marg_prob <- colMeans(po, na.rm=TRUE)
    
    return(c(odds_ratio(marg_prob["1"], marg_prob["0"]), # Regression imputation estimate
             aipw(samp=samp%>% dplyr::select(-dmt), ps=ps, pot_out=po, tr = "smoking"))) # AIPW estimate
  }
}

######################################
#########     ANALYSIS      ##########
######################################

# estimation(): 
# Function that applies all estimators to a list of samples 
# Takes list of samples, list of correct propensity scores, 
# list of misspecified propensity score, treatment variable of interest
# and whether large sample (ls) properties are of interest
estimation <- function(m, sample_list, ps_corr, ps_err, ps_np=NULL, tr, ls = F){
  if(m == 1){
    print("#### Processing samples under outcome model A.")
  } else if(m == 2){
    print("#### Processing samples under outcome model B.")
  }else if(m == 3){
    print("#### Processing samples under outcome model C.")
  }
  
  # Unadjusted odds ratio 
  print("---------- UA OR")
  UA <- do.call(rbind, lapply(sample_list, unadjusted_OR, tr = tr))
  
  # Regression adjusting for all covariates 
  print("---------- ALL COR")
  ALL <- do.call(rbind, lapply(sample_list, adjust_for_all, tr = tr))
  
  # Regression imputation 
  print("---------- RI")
  RI <- do.call(rbind, lapply(sample_list, parametric_reg_imputation, tr=tr))
  
  
  # AIPW with correctly specified propensity score model 
  print("---------- AIPW")
  invisible(capture.output(AIPW <- do.call(rbind, lapply(seq_along(sample_list), 
                                                            function(x){parametric_aipw(ps_corr[[x]], 
                                                                                        sample_list[[x]], tr = tr)}))))
  
  # AIPW with misspecified propensity score model 
  print("---------- AIPW*")
  invisible(capture.output(AIPW_err <- do.call(rbind, lapply(seq_along(sample_list), 
                                                                function(x){parametric_aipw(ps_err[[x]], 
                                                                                            sample_list[[x]], tr = tr)}))))
  
  # If main simulation, also estimate UPF and SS, as well as nonparametric estimates
  if(ls == F){
    # Regression with Stepwise selection 
    print("---------- SS COR")
    SS <-  do.call(rbind, mclapply(sample_list, stepwise_selection, tr=tr))
    
    # Regression with univariable pre-filtering 
    print("---------- UPF COR")
    UPF <- do.call(rbind, mclapply(sample_list, univariable_prefiltering, tr = tr))
    
    # Non parametric regression imputation and AIPW
    print("---------- NPRI and NPAIPW")
    NP <- do.call(rbind, lapply(seq_along(sample_list), 
                                function(x){nonparametric_est(
                                  samp = sample_list[[x]], 
                                  ps = ps_np[[x]],
                                  tr = tr)}))
    # Separate results 
    if(tr=="dmt"){
      # Non parametric Regression Imputation
      npRI <- NP[,1:3, drop = F]
      # Non parametric AIPW 
      npAIPW <- NP[,4:6, drop = F]
      
    } else if(tr=="smoking"){
      # Non parametric Regression Imputation
      npRI <- NP[,1, drop = F]
      # Non parametric AIPW 
      npAIPW <- NP[,2, drop=F]
    }
    
    
    # Return a list with all estimates 
    return(list(UA = UA, ALL = ALL, SS = SS, UPF = UPF, 
                RI = RI,
                AIPW = AIPW, AIPW_err = AIPW_err,
                npRI = npRI, npAIPW = npAIPW))
    
  }
  # Return a list with all estimates 
  return(list(UA = UA, ALL = ALL, 
              RI = RI,
              AIPW = AIPW, AIPW_err = AIPW_err))
}


# performance_measures():
# Function that calculates bias and MSE 
# Takes a list of estimates, treatment of interest, and true odds ratio 
# as arguments 
performance_measures <- function(est_list, tr, true_or){
  
  if(tr == "dmt"){
    # Bias is the difference between estimate and true value
    bias <- t(apply(est_list, 1, function(x) x-true_or))
    
    return(list(bias = round(colMeans(bias), 2), 
                mse = round(colMeans(bias^2), 2), 
                bias_mcse = round(apply(bias, 2, sd), 2),
                mse_mcse = round(apply(bias^2, 2, sd), 2)))
  }
  else if(tr == "smoking"){
    # Bias 
    bias <- est_list - true_or
    # Return bias and coverage, as well as their respective monte carlo standard error
    return(list(bias = round(mean(bias),2),
                mse = round(mean(bias^2),2), 
                bias_mcse = round(sd(bias),2),
                mse_mcse = round(sd(bias^2),2))) 
  }
}




