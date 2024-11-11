### SUPPORTING FUNCTIONS --------------------------------------------

# Expit function 
expit <- function(x){
  exp(x)/(1+exp(x))
}

# Logit function 
logit <- function(x){
  log(x/(1-x))
}

# Binomial draw for each element in probability vector 
binomial_draw <- function(prob_vec){
  apply(prob_vec, 1, rbinom, n = 1, size =1)
}

# Multinomial draw for each row in weight matrix
multinomial_draw <- function(weight_mat){
  t(apply(weight_mat, 1, rmultinom, n=1, size=1))
}

# Function to get the label of choice in multinomial draw
get_choice <- function(x){
  return(names(x)[which(x==1)])
}

### CONSTANT OPTIMIZATION 
# Optimize the constant shift for the BMI distribution 
constant_bmi <- function(chisq_vec, pY){
  
  cnst_bmi <- function(cnst, chisq_vec,  pY){
    p_hat_underw <- length(which(chisq_vec+cnst<18.5))/length(chisq_vec)
    p_hat_obese <- length(which(chisq_vec+cnst>30))/length(chisq_vec)
    value <- abs(p_hat_obese - pY[1])^2 + abs(p_hat_underw - pY[2])^2
    return(value)
  }
  
  value <- optimize(cnst_bmi,
                    interval = c(13, 20),
                    chisq_vec = chisq_vec,
                    pY=pY,
                    lower =13,
                    upper = 20,
                    tol=.Machine$double.eps^0.5)$minimum
  return(value)
}



### INTERCEPT OPTIMISATION 
# Optimizing intercept for logit model
beta0_logit <- function(betaX, X, pY){
  
  beta0l <- function(beta0, betaX, X,  pY){
    
    # Linear function of covariates
    etaY  <- as.numeric(X %*% c(beta0, betaX))
    
    # Minimize the distance between the average prevalence 
    # and the wanted prevalence 
    value <- abs(mean(expit(etaY))  - pY)
    
    return(value)
  }
  
  value <- optimize(beta0l,
                    interval = logit(pY) + c(-40,40),
                    betaX = betaX,
                    X = X,
                    pY=pY,
                    tol=.Machine$double.eps^0.5)$minimum
  return(value)
}

# Optimizing intercepts for multinomial logit model
beta0_mlogit <- function(betaX, X, pY, var){
  
  beta0ml <- function(beta0, betaX, X, pY){
    # Linear function of covariates 
    etaY <- matrix(0, nrow = nrow(X), ncol = nrow(betaX))
    for(j in 1:nrow(betaX)){
      etaY[,j]  <- as.numeric(X %*% c(beta0[j], betaX[j,]))
    }
    # Denominator for prevalence expression 
    denom <- 1 + rowSums(exp(etaY))
    # Estimated prevalences 
    ratio <- apply(exp(etaY), 2, function(x){ x/denom })
    # Minimize the square root of the sum of squared distances 
    value <- sqrt(sum((colMeans(ratio) - pY)^2) + ((1- sum(colMeans(ratio))-(1-sum(pY))))^2)
    return(value)
  }
  
  if(var == "dc"){
    value <- optim(par = rep(0, length(pY)),
                   fn = beta0ml,
                   betaX = betaX,
                   X = X,
                   pY = pY)$par
  } else if(var =="dmt"){
    value <- optim(par = rep(0, length(pY)),
                   fn = beta0ml,
                   betaX = betaX,
                   X = X,
                   pY = pY, 
                   method = "BFGS")$par
  }
  return(value)
}

### QUANTILE OPTIMIZATION 
# Function to optimize quantile in distribution for probit draws 
quantile_probit <- function(etaY, pY){
  q_probit <- function(qntl, etaY,  pY){
    # Minimize the distance between the average prevalence 
    # and the wanted prevalence 
    value <- abs(mean(pnorm(etaY, mean = quantile(etaY, qntl), sd = sd(etaY)))  - pY)
    return(value)
  }
  
  value <- optimize(q_probit,
                    interval = c(0, 1),
                    etaY = etaY,
                    pY=pY,
                    tol=.Machine$double.eps^0.5)$minimum
  return(value)
}


### COVARIATES ------------------------------------------------------------
# dgp_covars(): 
# Data generating function of covariates and treatments
# - n = sample size
# - parameters = list of target parameters, 
#                if NULL, sample means from Louapre et al are used
dgp_covars <- function(i,n,parameters=NULL){
  
  # Averages and proportions from Louapre et al 
  if(is.null(parameters)){
    parameters = list(p_male = 98/347,
                      mean_age = 44.6, 
                      sd_age = 12.8,
                      p_smok = 33/347,
                      p_bmi = c(24/347, 0.045), # Obesity, underweight
                      p_diabetes = 16/347,
                      p_cardiov = 23/347,
                      p_pulm = 15/347,
                      p_dc = c(276/347,48/347, 17/347), # RRMS, SPMS, PPMS
                      md_edss = 2, 
                      p_dmt = c(63/347,         # None
                                20/347, 33/347, #Interferon beta,Glatiramer
                                33/347, 35/347, #Teriflunomide,Dimethylfumarate
                                57/347, 5/347,  #Natalizumab,Other
                                42/347, 38/347, #Fingolimod,Ocrelizumab
                                17/347, 3/347)  #Rituximab,Cladribine  
    )
  }
  
  # ID (used for joins and other operations)
  id <- 1:n 
  
  #################
  ##    MALE     ##
  #################
  male <- rbinom(n, 1, parameters$p_male)
  
  ###################
  ##      AGE      ##
  ###################
  # Parameters of lognormal distribution 
  varlog <- log(1+(parameters$sd_age/parameters$mean_age)^2)
  meanlog <- log(parameters$mean_age+0.55)-0.5*varlog
  # Draws
  age <- floor(rlnorm(n, meanlog = meanlog, sdlog = sqrt(varlog)))
  
  # Initiate tibble and truncate age 
  samp <- tibble(id=id, male=male, age=age) %>% 
    mutate(age = case_when(
      age >85 ~ 85, 
      age < 18 ~ 18, 
      TRUE ~ age))
  
  ###############
  #### SMOKING ##
  ###############
  # Create square of age and male-age interaction
  covars_smok <- as.matrix(samp %>% 
                             dplyr::select(male, age) %>%
                             mutate(age_sq = age^2, 
                                    male_age = male*age))
  # Slope coefficients
  coef_smok <- c(0.5, 0.001, 0.0002, 0.001)
  # Find intercept
  intrcpt_smok <- beta0_logit(betaX = coef_smok, 
                              X = cbind(1, covars_smok), 
                              pY = parameters$p_smok)
  # Generate probabilities
  prob_smok <- as.matrix(1/(1 + exp( - intrcpt_smok - covars_smok %*% coef_smok)))
  # Draw smoking status and add to tibble
  samp <- samp %>% mutate(smoking = binomial_draw(prob_smok[, 1, drop = FALSE]))
  
  #################
  ####   BMI    ###
  #################
  # Degrees of freedom for chi square distribution 
  dfs <- matrix(0.1*samp$age + 2*samp$male + 6*samp$smoking - 
                  0.15*samp$age*samp$smoking - 3*samp$male*samp$smoking, 
                ncol=1)
  # Draws of chi square deviates
  chisq_dev <- apply(dfs, 1, rchisq, n = 1)
  # Determine constant shift
  constant <-constant_bmi(chisq_dev,parameters$p_bmi)
  # Add BMI and obesity indicator to sample 
  samp <- samp %>% mutate(bmi = chisq_dev+constant) %>% 
    mutate(obesity = ifelse(bmi>30, 1, 0))
  
  #################
  #    DIABETES   #
  #################
  # Covariates   
  covars_di <- as.matrix(samp %>% dplyr::select(age, male, smoking, bmi))
  # Slope coefficients
  coef_di <- c(0.5, 1, 1, 1)/10
  # Determine intercept 
  intrcpt_di <- beta0_logit(betaX = coef_di, X = cbind(1, covars_di), pY = parameters$p_diabetes)
  # Fitted probabilities 
  prob_di <- matrix(1/(1 + exp( - intrcpt_di - covars_di %*% coef_di)), ncol=1)
  # Add diabetes dummy to tibble
  samp <- samp %>% 
    mutate(diabetes =  binomial_draw(prob_di[, 1, drop = FALSE]))
  
  ###############################
  #    CARDIOVASCULAR DISEASE   #
  ###############################
  # Covariates 
  covars_cd <- as.matrix(samp %>% 
                           dplyr::select(age, male, smoking, bmi, 
                                         diabetes))
  # Slope coefficients
  coef_cd <- c(0.5, 2, 4, 0.5, 3)/10 
  intrcpt_cd <- beta0_logit(betaX = coef_cd, X = cbind(1, covars_cd), 
                            pY = parameters$p_cardiov)
  # Fitted probabilities
  prob_cd <- matrix(1/(1 + exp( - intrcpt_cd - covars_cd %*% coef_cd)), ncol = 1)
  # Add dummy for cardiovascular disease to tibble
  samp <- samp %>% 
    mutate(cardiovasc_dis = binomial_draw(prob_cd[, 1, drop = FALSE]))
  
  ###############################
  #      PULMONARY DISEASE      #
  ###############################
  # Covariates
  covars_pd <- as.matrix(samp %>% dplyr::select(age, male, smoking, bmi))
  # Slope coefficients
  coef_pd <- c(0.5, 0.1, 5, 0.1)/10
  # Determine intercept
  intrcpt_pd <- beta0_logit(betaX = coef_pd, X = cbind(1, covars_pd), 
                            pY = parameters$p_pulm)
  
  # Fitted probabilities 
  prob_pd <- matrix(1/(1 + exp( - intrcpt_pd - covars_pd %*% coef_pd)), ncol = 1)
  
  # Add dummy for pulmonary disease to tibble
  samp <- samp %>% 
    mutate(pulmonary_dis = binomial_draw(prob_pd[, 1, drop = FALSE]))
  
  #####################
  ## DISEASE COURSE  ##
  #####################
  # Proportion in sample 
  p_dc <- parameters$p_dc
  # Covariates
  covars_dc <- as.matrix(samp %>% dplyr::select(age, male, smoking, 
                                                bmi, cardiovasc_dis, 
                                                pulmonary_dis, diabetes))
  # Coefficient matrix 
  coef_dc <- rbind(c(3, -1, 5, 3, 3, 3, 2),  # RRMS
                   c(3, -1, 5, 5, 5, 3, 2), # SPMS 
                   c(2, 2, 5, 2, 2, 2, 2))/5e1  # PPMS 
  
  # Intercept 
  intrcpt <- beta0_mlogit(betaX = coef_dc, X = cbind(1, covars_dc), 
                          pY = p_dc, var = "dc")
  
  # Fitted probabilities 
  weights_dc <- exp(cbind(1,covars_dc) %*% t(cbind(intrcpt, coef_dc)))
  denom <- 1 + rowSums(weights_dc)
  prob_dc <- apply(weights_dc, 2, function(x){ x/denom })
  prob_dc <- cbind(1-rowSums(prob_dc), prob_dc)
  
  # Draws of disease course 
  disease_course <- multinomial_draw(prob_dc)
  colnames(disease_course) <- c("CIS", "RRMS", "SPMS", "PPMS")
  disease_course <- data.frame(id= 1:n, disease_course,
                               disease_course=apply(disease_course, 
                                                    1, get_choice))
  # Add disease course and dummy for progressive MS to tibble
  samp <- samp %>% full_join(disease_course, by="id") %>% 
    mutate(progressive_course = ifelse(SPMS==1 | PPMS==1, 1, 0))
  
  #####################
  ####### EDSS  ######
  #####################
  # EDSS scale labels
  edss_scale <- seq(0, 9.5, 0.5)
  # Continuous edss scale
  edss_cont <-  exp(- 2*samp$CIS  - 1*samp$RRMS + 1*samp$SPMS + 2*samp$PPMS) +
    rchisq(n, df=7) 
  # Indicator if observation is above or below median of edss_cont
  edss_bin <- ifelse(edss_cont < median(edss_cont),1,0)
  # Discretize the continuous EDSS scale
  edss_df <- data.frame(id=1:n, edss_cont, edss_bin)
  # Below median
  lower_md <- edss_df %>% 
    filter(edss_bin == 1) %>% 
    mutate(edss_fac = discretize(data.frame(edss_cont), method = "interval", 
                                 breaks = length(which(edss_scale<=parameters$md_edss)))$edss_cont)
  levels(lower_md$edss_fac) <- seq(0,parameters$md_edss, 0.5)
  lower_md$edss_fac <- as.character(lower_md$edss_fac)
  # Above median 
  higher_md <- edss_df %>% 
    filter(edss_bin == 0) %>% 
    mutate(edss_fac = discretize(data.frame(edss_cont), method = "interval", 
                                 breaks = length(which(edss_scale>=parameters$md_edss)))$edss_cont)
  levels(higher_md$edss_fac) <- seq(parameters$md_edss,9.5, 0.5)
  higher_md$edss_fac <- as.character(higher_md$edss_fac)
  
  # Join together
  edss <- rbind(lower_md, higher_md)
  # Remake into factor
  edss$edss_fac <- factor(edss$edss_fac, levels = seq(0, 9.5, 0.5))
  # Make numeric and add to tibble
  edss$edss <- as.numeric(levels(edss$edss_fac))[edss$edss_fac]
  samp <- full_join(samp, edss %>% dplyr::select(id, edss), by = "id")
  
  
  ########################
  ####       DMT      ####
  ########################
  # Names of treatments
  dmts <- c("None",
            "Interferon beta", "Glatiramer",
            "Teriflunomide", "Dimethylfumarate", "Natalizumab", "Other",
            "Fingolimod",  "Ocrelizumab", "Rituximab", "Cladribine" ) #,  "Alemtuzumab")
  
  # Proportions in sample 
  p_dmt <- parameters$p_dmt
  
  # Covariates
  covars_dmt <- as.matrix(samp %>%
                            dplyr::select(cardiovasc_dis, pulmonary_dis, diabetes,
                                          CIS, RRMS, SPMS, PPMS, edss) %>% 
                            mutate(edss_sq = edss^2, 
                                   cv_edss = edss*cardiovasc_dis, 
                                   pd_edss = edss*pulmonary_dis, 
                                   dia_edss = edss*diabetes)) 
  
  # Coefficient matrix
  coef_dmt <- matrix(rep(1, 12*11),
                     ncol = 12, nrow=11, byrow=F)
  rownames(coef_dmt) <- dmts
  colnames(coef_dmt) <- colnames(covars_dmt)
  # CIS are more likely to have no treatment or first line treatments 
  coef_dmt[, c("CIS")] <- c(5, 7, rep(4,2), 1, 3,  rep(-2, 5))
  # RRMS are more likely to have first line treatments or no treatment
  coef_dmt[, c("RRMS")] <- c(5, 7, rep(4,2), 1, 3,  rep(-2, 5))
  # SPMS more likely to have second line treatments 
  coef_dmt[, c("SPMS")] <- c(rep(-2,4), rep(5,7))
  # Few treatments for PPMS, Ocrelizumab is one 
  coef_dmt[c("Ocrelizumab", "Other"), c("PPMS")] <- c(10,5)
  
  # Those with comorbidities are more likely to delay initiation of treatment 
  # Some treatments are used for patients intolerant to other treatments 
  coef_dmt[, c("cardiovasc_dis", "diabetes", "pulmonary_dis")] <- c(7, rep(-5,4), rep(-1,6))
  coef_dmt[c("Cladribine", "Other"), 
           c("cardiovasc_dis", "pulmonary_dis", "diabetes")] <- 3
  
  coef_dmt[c("Teriflunomide", "Dimethylfumarate"),1:3] <- 5
  coef_dmt[c("Cladribine"), 1:3] <- 5
  
  # Those with more disability are more likely to be on second line treatment
  coef_dmt[,c("edss")] <- c(rep(-1, 2), rep(3,4), rep(5, 3), rep(1,2))
  coef_dmt[,c("edss_sq")] <- c(rep(0.1, 2), rep(-0.01,4), rep(-0.5, 3), rep(0.05,2))
  coef_dmt[,c("dia_edss")] <- seq(0.001, 0.05, length.out = 11)
  coef_dmt[,c("pd_edss")] <- seq(-0.001, 0.02, length.out = 11)
  coef_dmt[,c("cv_edss")] <- seq(0.01, -0.01, length.out = 11)
  
  coef_dmt[,c(1:8)] <- coef_dmt[,c(1:8)]/10
  
  # Find intercepts  
  intrcpt <- beta0_mlogit(coef_dmt, cbind(1,covars_dmt), p_dmt, var = "dmt")
  
  # Matrix of fitted probabilities 
  weights_dmt <- exp(cbind(1,covars_dmt) %*% t(cbind(intrcpt, coef_dmt))) # no negative values
  denom <- 1 + rowSums(weights_dmt)
  prob_dmt <- apply(weights_dmt, 2, function(x){ x / denom })
  prob_dmt <- cbind(prob_dmt, 1-rowSums(prob_dmt))
  
  # Draws of disease course, add id to join
  dmt <-  as_tibble(cbind(id=1:n, multinomial_draw(prob_dmt)), 
                    .name_repair = "minimal")
  colnames(dmt) <- c("id",  dmts, "Alemtuzumab")
  
  # Add to tibble and create DMT variable with four categories
  samp <- full_join(samp, dmt, by="id") %>%
    mutate(No_treatment = case_when( None == 1 ~ 1, TRUE ~ 0),
           No_risk = case_when(`Interferon beta` == 1 | 
                                 Glatiramer == 1  ~ 1, 
                               TRUE ~ 0),
           Low_risk = case_when(Teriflunomide == 1 | 
                                  Dimethylfumarate== 1 |
                                  Natalizumab ==1 | 
                                  Other ==1 ~ 1, 
                                TRUE ~ 0),
           Moderate_to_high_risk = case_when(Fingolimod == 1 | 
                                               Ocrelizumab == 1 |
                                               Rituximab == 1 |
                                               Cladribine ==1 | 
                                               Alemtuzumab == 1 ~ 1, 
                                             TRUE ~ 0)) %>%
    mutate(dmt = case_when(No_treatment==1 ~ "No_treatment",
                           No_risk==1 ~ "No_risk",
                           Low_risk ==1 ~ "Low_risk",
                           Moderate_to_high_risk ==1 ~ "Moderate_to_high_risk")) %>%
    mutate(dmt = as.factor(dmt),
           No_treatment = as.integer(No_treatment),
           No_risk = as.integer(No_risk), 
           Low_risk = as.integer(Low_risk), 
           Moderate_to_high_risk = as.integer(Moderate_to_high_risk))
  
  return(samp) 
}

### OUTCOMES ------------------------------------------------
# dgp_outc(): 
# Function generating samples with outcomes 
# - samp: a tibble  created with dgp_covars 
# - marg_prob: vector of marginal probabilities of potential outcomes
# - tr: a string indicating treatment variable, takes values 'dmt' or 'smoking'
# - mod: a string taking values 'AB' or 'C', indicating which outcome model
# - intrctn: a boolean separating outcome model A and B. If set to false,
# data are generated according to outcome model A, otherwise it is generated 
# according to outcome model B
dgp_outc <- function(samp, marg_prob, tr, mod,  intrctn = T){
  
  n <- nrow(samp)
  # SCENARIO I
  if(tr == "dmt"){
    # Covariates for generating the outcome under DMT
    covars <- samp %>% 
      dplyr::select(age, male, smoking, bmi, 
                    cardiovasc_dis, pulmonary_dis, diabetes, 
                    CIS, RRMS, SPMS, PPMS,  
                    edss)
    if(mod=="AB"){
      covars <- as.matrix(covars)
      # Model B
      if(intrctn==T){
        # Coefficient vector for No treatment 
        coef_po_nt <- c(0.1, 2, 0.2, 0.015, 
                        2, 2, 1.5, 
                        1, 5, 0.5, 1.5, 
                        1.75)/7
        # Coefficient vector for No risk 
        coef_po_nr <- c(0.4, 3, 1, 0.01, 
                        1.5, 1, 1.75, 
                        0.25, 5, 1, 5,  
                        1)/5
        # Coefficient vector for Low risk 
        coef_po_lr <- c(0.75, 1.5, 0.5, 0.02, 
                        1.6, 0.5, 0.01,
                        1.2, 1.6, 5, 1.7,
                        0.75)/3
        # Coefficient vector for Moderate to high risk  
        coef_po_mhr <- c(0.2, 0.01, 1.6, 0.025, 
                         1.5, 3.5, 2, 
                         5, 0.5, 0.2, 1, 
                         0.85)/2
      } 
      # Model A 
      else if(intrctn==F){
        coef_po_nt <- c(0.2, 1, 1.6, 0.05, 
                        1.5, 3.5, 0.001, 
                        0.25, 5, 0.2, 1, 
                        1.5)/5
        coef_po_nr <- coef_po_nt
        coef_po_lr <- coef_po_nt
        coef_po_mhr <- coef_po_nt
      }
      
      # Intercepts 
      intrcpt_nt <- beta0_logit(betaX = coef_po_nt, 
                                X = cbind(1, covars), 
                                pY = marg_prob[1])
      intrcpt_nr <- beta0_logit(betaX = coef_po_nr, 
                                X = cbind(1, covars), 
                                pY = marg_prob[2])
      intrcpt_lr <- beta0_logit(betaX = coef_po_lr, 
                                X = cbind(1, covars), 
                                pY = marg_prob[3])
      intrcpt_mhr <- beta0_logit(betaX = coef_po_mhr, 
                                 X = cbind(1, covars), 
                                 pY = marg_prob[4])
      
      # Matrix of fitted probabilities 
      prob_po <- cbind(1/(1 + exp(- intrcpt_nt - covars %*% coef_po_nt)),
                       1/(1 + exp(- intrcpt_nr - covars %*% coef_po_nr)), 
                       1/(1 + exp(- intrcpt_lr - covars %*% coef_po_lr)),
                       1/(1 + exp(- intrcpt_mhr - covars %*% coef_po_mhr)))
      
      # Draws of potential outcomes 
      samp <- samp %>% 
        mutate(Y_nt = binomial_draw(prob_po[, 1, drop = FALSE]), #No treatment
               Y_nr = binomial_draw(prob_po[, 2, drop = FALSE]), #No risk
               Y_lr =  binomial_draw(prob_po[, 3, drop = FALSE]), #Low risk
               Y_mhr = binomial_draw(prob_po[, 4, drop = FALSE])) #Moderate/High risk
      
    }
    else if(mod == "C"){
      covars_expanded <- model.matrix(~ -1 + .^2 +I(age^2) + I(age^3) + 
                                        I(bmi^2) + I(age^0.5) + I(age^4) + 
                                        I(edss^2), data = covars)
      
      # Coefficient vector for No treatment 
      coef_po_nt <- c(c(0.1, 2, 2, 0.015, 
                        2, 2, 1.5, 
                        1, 5, 2, 1.5, 
                        1.5)/7, rep(0.5, 72))
      # Coefficient vector for No risk 
      coef_po_nr <- c(c(0.4, 1, 1, 0.01, 
                        1.5, 1.8, 1.75, 
                        0.25, 0.5, 5, 1,  
                        1.1)/5, rep(0.75, 72))
      # Coefficient vector for Low risk 
      coef_po_lr <- c(c(0.5, 1.5, 0.5, 0.015, 
                        1.6, 2, 1.7,
                        1.2, 1.6, 4, 1.7,
                        1.75)/3, rep(1, 72))
      # Coefficient vector for Moderate to high risk  
      coef_po_mhr <- c(c(0.025, 1, 1, 0.01, 
                         1.5, 2, 1.5, 
                         1, 5, 0.2, 1, 
                         1.25)/2, rep(1.25, 72))
      
      eta_Y_nt <- (covars_expanded%*%coef_po_nt)
      eta_Y_nr <- (covars_expanded%*%coef_po_nr) 
      eta_Y_lr <- (covars_expanded%*%coef_po_lr)
      eta_Y_mhr <- (covars_expanded%*%coef_po_mhr)
      
      q_Y_nt <- quantile_probit(eta = eta_Y_nt,
                                pY = marg_prob[1])
      q_Y_nr <- quantile_probit(eta = eta_Y_nr,
                                pY = marg_prob[2])
      q_Y_lr <- quantile_probit(eta = eta_Y_lr,
                                pY = marg_prob[3])
      q_Y_mhr <- quantile_probit(eta = eta_Y_mhr,
                                 pY = marg_prob[4])
      
      # Matrix of fitted probabilities 
      prob_po <- cbind(pnorm(eta_Y_nt, 
                             mean = quantile(eta_Y_nt, q_Y_nt), 
                             sd = sd(eta_Y_nt)),
                       pnorm(eta_Y_nr, 
                             mean = quantile(eta_Y_nr, q_Y_nr), 
                             sd = sd(eta_Y_nr)),
                       pnorm(eta_Y_lr, 
                             mean = quantile(eta_Y_lr, q_Y_lr), 
                             sd = sd(eta_Y_lr)),
                       pnorm(eta_Y_mhr, 
                             mean = quantile(eta_Y_mhr, q_Y_mhr), 
                             sd = sd(eta_Y_mhr)))
      
      # Draws of potential outcomes 
      samp <- samp %>% 
        mutate(Y_nt = binomial_draw(prob_po[, 1, drop = FALSE]), 
               Y_nr = binomial_draw(prob_po[, 2, drop = FALSE]), 
               Y_lr =  binomial_draw(prob_po[, 3, drop = FALSE]), 
               Y_mhr = binomial_draw(prob_po[, 4, drop = FALSE]))
      
      
    }
    # Observed outcome
    severity <- c()
    for(i in 1:n){
      if(samp$No_treatment[i] == 1){
        severity[i] <- samp$Y_nt[i]
      } else if(samp$No_risk[i] == 1){
        severity[i] <- samp$Y_nr[i]
      } else if(samp$Low_risk[i] == 1){
        severity[i] <- samp$Y_lr[i]
      } else if(samp$Moderate_to_high_risk[i] == 1){
        severity[i] <- samp$Y_mhr[i]
      }
    }
    # Create data frame including all variables necessary for analyses 
    samp <- samp %>% mutate(severity = severity) %>% 
      dplyr::select(male, age, smoking, bmi, 
                    diabetes, cardiovasc_dis, pulmonary_dis,
                    progressive_course, edss,
                    No_treatment, No_risk, Low_risk, Moderate_to_high_risk, 
                    dmt,
                    severity)
    
  }
  # SCENARIO II
  if(tr == "smoking"){
    # Covariates for outcome generating models 
    covars <- samp %>% 
      dplyr::select(age, male, bmi, 
                    cardiovasc_dis, pulmonary_dis, diabetes,
                    CIS, RRMS, SPMS, PPMS,  
                    edss, 
                    `Interferon beta`, Glatiramer, 
                    Teriflunomide, Dimethylfumarate, Natalizumab, Other, 
                    Fingolimod, Ocrelizumab, Rituximab, Cladribine, Alemtuzumab)
    if(mod == "AB"){
      covars <- as.matrix(covars)
      if(intrctn == T){
        # Coefficient vector for Smoking = 0
        coef_po_ntsmok <- c(0.2, 1, 0.01, 
                            1.5, 2, 1.5, 
                            0.25, 0.5, 1.2, 1,  
                            0.3, 
                            0.1, 0.2, 
                            0.7, 0.5, 0.4, 0.2, 
                            0.5, 0.6, 0.4, 0.6, 0.5)/3
        
        # Coefficient vector for Smoking = 1
        coef_po_smok <- c(0.5, 1.5, 0.05,
                          1.4, 1.45, 1.7,
                          0.5, 1.6, 2.1, 1.25,
                          0.8, 
                          0.2, 0.5, 
                          0.5, 0.2, 0.6, 0.4, 
                          0.65, 0.5, 0.5, 0.3, 0.4)/3
      } else if(intrctn == F){
        # Coefficient vector for Smoking = 0
        coef_po_ntsmok <- c(0.5, 1, 0.01, 
                            1.5, 2, 1.5, 
                            0.25, 0.5, 1, 1,  
                            0.3, 
                            0.1, 0.2, 
                            0.4, 0.3, 0.4, 0.2, 
                            0.5, 0.6, 0.4, 0.6, 0.5)/5
        # Coefficient vector for Smoking = 1
        coef_po_smok <- coef_po_ntsmok
      }
      
      # Intercepts 
      intrcpt_ntsmok <- beta0_logit(betaX = coef_po_ntsmok, 
                                    X = cbind(1, covars), 
                                    pY = marg_prob[1])
      intrcpt_smok <- beta0_logit(betaX = coef_po_smok, 
                                  X = cbind(1, covars), 
                                  pY = marg_prob[2])
      
      # Matrix of fitted probabilities 
      prob_po <- cbind(1/(1 + exp(- intrcpt_ntsmok - covars %*% coef_po_ntsmok)),
                       1/(1 + exp(- intrcpt_smok - covars %*% coef_po_smok)))
      
      # Draws for Potential outcomes 
      samp <- samp %>% 
        mutate(Y_ntsmok = binomial_draw(prob_po[, 1, drop = FALSE]), 
               Y_smok = binomial_draw(prob_po[, 2, drop = FALSE])) 
      
    } 
    if(mod == "C" ){
      # Add non-linear terms and interactions 
      covars_expanded <- model.matrix(~ -1 + .^2 + I(age^2) + I(age^3) + 
                                        I(bmi^2) + I(age^0.5) + I(age^4) + 
                                        I(edss^2), 
                                      data = covars) 
      
      # Coefficient vector for Smoking = 0
      coef_po_ntsmok <- c(c(0.2, 1, 0.01, 
                            1.5, 2, 1.5, 
                            0.2, 0.5, 1, 1,  
                            2, 
                            0.1, 0.2, 
                            0.4, 0.3, 0.4, 0.2, 
                            0.5, 0.6, 0.4, 0.6, 0.5)/5, 
                          rep(1, 237))
      
      # Coefficient vector for Smoking = 1
      coef_po_smok <- c(c(0.5, 1.5, 0.06,
                          1.4, 1.2, 1.7,
                          0.5, 1.6, 2, 1.25,
                          2.7, 
                          0.2, 0.5, 
                          0.5, 0.2, 0.6, 0.2, 
                          0.65, 0.5, 0.5, 0.3, 0.4)/5, 
                        rep(1.5, 237))
      eta_Y_ntsmok <- covars_expanded%*%coef_po_ntsmok 
      eta_Y_smok <- covars_expanded%*%coef_po_smok 
      
      # Quantile of normal distribution 
      q_Y_ntsmok <- quantile_probit(eta = eta_Y_ntsmok, pY = marg_prob[1])
      q_Y_smok <- quantile_probit(eta = eta_Y_smok, pY = marg_prob[2])
      
      # Matrix of probabilities   
      prob_po <- cbind(pnorm(eta_Y_ntsmok,
                             mean = quantile(eta_Y_ntsmok, q_Y_ntsmok),
                             sd = sd(eta_Y_ntsmok)),
                       pnorm(eta_Y_smok,
                             mean = quantile(eta_Y_smok, q_Y_smok),
                             sd = sd(eta_Y_smok)))
      
      # Draws of potential outcomes 
      samp <- samp %>% 
        mutate(Y_ntsmok = binomial_draw(prob_po[, 1, drop = FALSE]), 
               Y_smok = binomial_draw(prob_po[, 2, drop = FALSE]))
      
    }
    # The observed outcome is the potential outcome 
    # corresponding to smoking status 
    severity <- ifelse(samp$smoking == 1, samp$Y_smok, samp$Y_ntsmok)
    # Create data frame including all variables necessary for analyses 
    samp <- samp %>% mutate(severity = severity) %>% 
      dplyr::select(male, age, smoking, bmi, 
                    diabetes, cardiovasc_dis,pulmonary_dis,
                    progressive_course, edss,
                    No_treatment, No_risk, Low_risk, Moderate_to_high_risk, 
                    dmt,
                    severity)
  }
  return(samp)
} 






