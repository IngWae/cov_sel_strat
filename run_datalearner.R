# run_datalearner(): 
# Function that runs the full simulation, from sampling to finished figures and tables
# - seed: Desried seed 
# - Rs: Number of replications
# - n: Sample size 
# - marg_prob_dmt: True marginal probabilities of outcome under each level of DMT
# - marg_prob_smok: True marginal probabilities of outcome under each level of Smoking
run_datalearner <- function(seed, Rs, n, marg_prob_dmt, marg_prob_smok, ls=F){
    # True marginal odds ratios 
    true_OR_dmt <- c(odds_ratio(marg_prob_dmt[2], marg_prob_dmt[1]), 
                   odds_ratio(marg_prob_dmt[3], marg_prob_dmt[1]), 
                   odds_ratio(marg_prob_dmt[4], marg_prob_dmt[1]))
  
    true_OR_smok <- c(odds_ratio(marg_prob_smok[2], marg_prob_smok[1]))
  
    # Set seed
    set.seed(seed)
    
    # Generate covariates 
    print("######## GENERATING DATA ######## ")
    print("######## Sampling covariates")
    samples_covars <-  lapply(1:Rs, dgp_covars, n=n)
    
    # Estimate propensity scores in scenario I 
    print("######## Estimating propensity scores in Scenario I")
    # Propensity scores 
    # - Correctly specified propensity score model 
    ps_dmt_corr <- mclapply(samples_covars, propensity_scores, tr = "dmt", type = "correct")
    # - Misspecified propensity score model 
    ps_dmt_err <- mclapply(samples_covars, propensity_scores, tr = "dmt", type = "incorrect")
    # - Nonparametric propensity score model 
    if(ls == F){
      ps_dmt_np <- mclapply(samples_covars, propensity_scores, tr = "dmt", type = "nonparametric")
    }
    
    # Sample under scenario I 
    print("######## Sampling in Scenario I")
    samples_dmt <- list(A = mclapply(samples_covars, dgp_outc, 
                                   marg_prob = marg_prob_dmt, 
                                   tr = "dmt", mod = "AB", intrctn = F), 
                        B = mclapply(samples_covars, dgp_outc,  
                                   marg_prob = marg_prob_dmt,
                                   tr = "dmt", mod = "AB", intrctn = T), 
                        C = mclapply(samples_covars, dgp_outc, 
                                   marg_prob = marg_prob_dmt, 
                                   tr = "dmt", mod = "C", intrctn = T))
    
    # Estimate propensity scores in scenario II 
    print("######## Estimating propensity scores in Scenario II")
    # Propensity scores 
    # - Correctly specified model 
    ps_smok_corr <- mclapply(samples_covars, propensity_scores, tr = "smoking", type = "correct")
    # - Misspecified model 
    ps_smok_err <- mclapply(samples_covars, propensity_scores, tr = "smoking", type = "incorrect")
    # - Non parametric
    if(ls == F){
    ps_smok_np <- mclapply(samples_covars, propensity_scores, tr = "smoking", type = "nonparametric")
    }
    
    # Sample under scenario II
    print("######## Sampling in Scenario II")
    samples_smok <- list(A = mclapply(samples_covars, dgp_outc,
                                    marg_prob = marg_prob_smok, 
                                    tr = "smoking", mod = "AB", intrctn = F), 
                         B = mclapply(samples_covars, dgp_outc, 
                                    marg_prob = marg_prob_smok, 
                                    tr = "smoking", mod = "AB", intrctn = T), 
                         C = mclapply(samples_covars, dgp_outc, 
                                    marg_prob = marg_prob_smok, 
                                    tr = "smoking", mod = "C", intrctn = T))

  # Large sample properties of estimators   
  if(ls == T){
    print("######## LARGE SAMPLE PROPERTIES ########")
    print("######## SCENARIO I ########")
    # Odds ratio estimates with standard deviation and confidence intervals 
    ests_dmt <- sapply(seq_along(samples_dmt), 
                         function(x){estimation(m=x,
                                                samples_dmt[[x]], 
                                                ps_corr = ps_dmt_corr, 
                                                ps_err = ps_dmt_err, 
                                                tr = "dmt", ls = T)},
                         simplify = FALSE, USE.NAMES = TRUE)
    print("######## SCENARIO II ########")
    # Odds ratio estimates with standard deviation and confidence intervals 
    ests_smok <- mclapply(seq_along(samples_smok), 
                          function(x){estimation(m=x,
                                                 samples_smok[[x]], 
                                                 ps_corr = ps_smok_corr, 
                                                 ps_err = ps_smok_err,
                                                 tr = "smoking", ls = T)})
    
    return(list(ls_dmt = ls_table(ests_dmt, tr = "dmt", true_ors = true_OR_dmt),
                ls_smok = ls_table(ests_smok, tr = "smoking", true_ors = true_OR_smok)))
      
      
  } # Main simulation 
    else if(ls == F){
    
      print("######## MONTE CARLO SIMLATION STUDY ########")
    
      print("######## SCENARIO I ########")
      # Odds ratio estimates 
      ests_dmt <- sapply(seq_along(samples_dmt), 
                         function(x){estimation(m=x,
                                                samples_dmt[[x]], 
                                                ps_corr = ps_dmt_corr, 
                                                ps_err = ps_dmt_err, 
                                                ps_np = ps_dmt_np,
                                                tr = "dmt")},
                         simplify = FALSE, USE.NAMES = TRUE)
      
      # Bias and MSE 
      measures_dmt <- mclapply(ests_dmt, lapply, performance_measures,  
                             true_or = true_OR_dmt, tr = "dmt")
      
      print("######## SCENARIO II ########")
      # Odds ratio estimates 
      ests_smok <- mclapply(seq_along(samples_smok), 
                          function(x){estimation(m=x,
                                                 samples_smok[[x]], 
                                                 ps_corr = ps_smok_corr, 
                                                 ps_err = ps_smok_err, 
                                                 ps_np = ps_smok_np,
                                                 tr = "smoking")})
      
      # Bias and MSE 
      measures_smok <- mclapply(ests_smok, lapply, performance_measures, 
                              true_or = true_OR_smok, tr = "smoking")
      
      # Return figures and tables 
      return(list(ol_dmt_corr = overlap_panel(samp = samples_covars[[1]], ps = ps_dmt_corr[[1]], tr = "dmt"),
                  ol_dmt_err = overlap_panel(samp = samples_covars[[1]], ps = ps_dmt_err[[1]], tr = "dmt"),
                  ol_dmt_np = overlap_panel(samp = samples_covars[[1]], ps = ps_dmt_np[[1]], tr = "dmt"),
                  bp_dmt = boxplot_panel(est_list = ests_dmt, tr = "dmt"),
                  res_dmt = tables(measures_dmt, tr = "dmt", type = "measure"),
                  mcse_dmt =  tables(measures_dmt, tr = "dmt", type = "mcse"),
                  ol_smok_corr = overlap_plot(samp = samples_covars[[1]], ps = ps_smok_corr[[1]], tr = "smoking"),
                  ol_smok_err = overlap_plot(samp = samples_covars[[1]], ps = ps_smok_err[[1]], tr = "smoking"),
                  ol_smok_np = overlap_plot(samp = samples_covars[[1]], ps = ps_smok_np[[1]], tr = "smoking"),
                  bp_smok = boxplot_panel(est_list = ests_smok, tr = "smoking"),
                  res_smok = tables(measures_smok, tr = "smoking", type = "measure"),
                  mcse_smok = tables(measures_smok, tr = "smoking", type = "mcse")))
  }
}
