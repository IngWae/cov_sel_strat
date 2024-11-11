##############
### TABLES ###
##############

# ls_tables():
# Function to create table of large sample bias. 
# Takes a large sample estimates, treatment variable (dmt/smoking) and 
# the true odds ratios as arguments. 
ls_table <- function(ls_ests, tr, true_ors){
  # List of estimators 
  estimators <- c("UA", "ALL", "RI", "AIPW", "AIPW*")
  
  if(tr == "smoking"){
    # Initiate  table (a data frame)
    tab <- as.data.frame(matrix(0, nrow = length(estimators), ncol = 4))
    # Fill in table 
    colnames(tab) <- c("Outcome Model:", c("A", "B", "C"))
    tab[,1] <- estimators
    for(j in 1:3){
      for(i in 1:length(estimators)){
        tab[i,j+1] <- round(c(ls_ests[[j]][[i]][1]-true_ors),2)
      }
    }
    # Format latex table
    return(kable(tab, format = "latex", row.names = F, 
                 align = c('l',rep('c',3)), booktabs=T,
                 escape = FALSE, linesep = "\\addlinespace",
                 caption = "Numerical approximations of the asymptotic bias in 
                 samples with n = 1 million. The true MCOR is $\\theta_{1|0} = 1.225$.", 
                 label = "scen2_lsbias") %>% 
             add_header_above(c(" ","Bias" = 3)) %>% 
             kableExtra::group_rows("$\\theta_{1|0}$", 1,5, escape = FALSE))
    
  } else if(tr == "dmt"){
    # Initiate  table (a data frame)
    tab <- as.data.frame(matrix(0, nrow = 3*length(estimators), ncol = 4))
    # Fill in table 
    colnames(tab) <- c("Outcome Model:", c("A", "B", "C"))
    tab[,1] <- estimators
    for(j in 1:3){
      for(i in 1:length(estimators)){
        tab[i,j+1] <- round(c(ls_ests[[j]][[i]][1]- true_ors[1]),2)
        tab[i+length(estimators),j+1] <- round(c(ls_ests[[j]][[i]][2]-true_ors[2]),2)
        tab[i+2*length(estimators),j+1] <- round(c(ls_ests[[j]][[i]][3]-true_ors[3]),2)
    }
    }
    # Format latex table
    return(kable(tab, format = "latex", row.names = F, align = c('l',rep('c',3)), 
                 booktabs=T, escape = FALSE,
                 linesep = "\\addlinespace",
                 caption = "Numerical approximations of the asymptotic bias in 
                 samples with n = 1 million. The true MCORs are $\\theta_{1|0} = 1$, 
                 $\\theta_{2|0} = 1.07$, and $\\theta_{3|0} = 1.42$.", 
                 label = "scen1_lsbias") %>% 
             add_header_above(c(" ", "Bias" = 3)) %>% 
             kableExtra::group_rows("$\\theta_{1|0}$", 1,5, escape = FALSE) %>% 
             kableExtra::group_rows("$\\theta_{2|0}$", 6,10, escape = FALSE) %>% 
             kableExtra::group_rows("$\\theta_{3|0}$", 11,15, escape = FALSE))
  }
}

# tables(): 
# Function to create tables of bias and mse for a list of performance measures. 
# Requires treatment variable to be specified (dmt/smoking)
tables <- function(measure_list, tr, type){
  # List of estimators 
  estimators <- c("UA", "ALL", "SS", "UPF", "RI",  
                  "AIPW", "AIPW*", "npRI", "npAIPW")
  if(tr == "dmt"){
    # Initiate  table (a data frame)
    tab <- as.data.frame(matrix(0, nrow = 3*length(estimators), ncol = 1+3*2))
    # Fill in table 
    colnames(tab) <- c("Outcome Model:", rep(c("A", "B", "C"), 2))
    tab[,1] <- c(rep(estimators,3))
    
    if(type == "measure"){
      for(j in 0:2){
        for(i in 1:length(estimators)){
          tab[i+j*length(estimators),2:4] <- c(measure_list[[1]][[i]]$bias[j+1],
                                               measure_list[[2]][[i]]$bias[j+1],
                                               measure_list[[3]][[i]]$bias[j+1])
          tab[i+j*length(estimators),5:7] <- c(measure_list[[1]][[i]]$mse[j+1],
                                               measure_list[[2]][[i]]$mse[j+1],
                                               measure_list[[3]][[i]]$mse[j+1])
        }
      }
      # Format latex table
      return(kable(tab, format = "latex", row.names = F, 
                   align = rep('c',7), booktabs=T,
                   escape = FALSE, linesep = "\\addlinespace",
                   caption = "Bias and MSE of estimators in Scenario I.", 
                   label = "scen1_bias") %>% 
               add_header_above(c(" ", "Bias" = 3, "MSE" = 3)) %>% 
               kableExtra::group_rows("$\\theta_{1|0}$", 1,9, escape = FALSE) %>% 
               kableExtra::group_rows("$\\theta_{2|0}$", 10,18, escape = FALSE) %>% 
               kableExtra::group_rows("$\\theta_{3|0}$", 19,27, escape = FALSE))
    } else if(type == "mcse"){
      for(j in 0:2){
        for(i in 1:length(estimators)){
          tab[i+j*length(estimators),2:4] <- c(measure_list[[1]][[i]]$bias_mcse[j+1],
                                               measure_list[[2]][[i]]$bias_mcse[j+1],
                                               measure_list[[3]][[i]]$bias_mcse[j+1])
          tab[i+j*length(estimators),5:7] <- c(measure_list[[1]][[i]]$mse_mcse[j+1],
                                               measure_list[[2]][[i]]$mse_mcse[j+1],
                                               measure_list[[3]][[i]]$mse_mcse[j+1])
        }
      }
      # Format latex table
      return(kable(tab, format = "latex", row.names = F, 
                   align = rep('c',7), booktabs=T,
                   escape = FALSE, linesep = "\\addlinespace",
                   caption = "Monte Carlo standard errors of bias and MSE in Scenario I.", 
                   label = "scen1_mcse") %>% 
               add_header_above(c(" ", "Bias" = 3, "MSE" = 3)) %>% 
               kableExtra::group_rows("$\\theta_{1|0}$", 1,9, escape = FALSE) %>% 
               kableExtra::group_rows("$\\theta_{2|0}$", 10,18, escape = FALSE) %>% 
               kableExtra::group_rows("$\\theta_{3|0}$", 19,27, escape = FALSE))
    }
    
    
    } else if(tr == "smoking"){
    # Initiate table (a data frame)
    tab <- as.data.frame(matrix(0, nrow = length(estimators), ncol = 1+3*2))
    # Fill in table 
    colnames(tab) <- c("Outcome Model:", rep(c("A", "B", "C"), 2))
    tab[,1] <- estimators
    
    if(type == "measure"){
      for(i in 1:length(estimators)){
        tab[i,2:4] <- c(measure_list[[1]][[i]]$bias,
                        measure_list[[2]][[i]]$bias,
                        measure_list[[3]][[i]]$bias)
        tab[i,5:7] <- c(measure_list[[1]][[i]]$mse,
                        measure_list[[2]][[i]]$mse,
                        measure_list[[3]][[i]]$mse)
      }
      # Format latex table 
      return(kable(tab, format = "latex", row.names = F, 
                   align = rep('c',7), booktabs=T, 
                   escape = FALSE, linesep = "\\addlinespace",
                   caption = "Bias and MSE of estimators in Scenario I.", 
                   label = "scen2_bias") %>% 
               add_header_above(c(" ", "Bias" = 3, "MSE" = 3)) %>% 
               kableExtra::group_rows("$\\theta_{1|0}$", 1,9, escape = FALSE))
    } else if(type == "mcse"){
      for(i in 1:length(estimators)){
        tab[i,2:4] <- c(measure_list[[1]][[i]]$bias_mcse,
                        measure_list[[2]][[i]]$bias_mcse,
                        measure_list[[3]][[i]]$bias_mcse)
        tab[i,5:7] <- c(measure_list[[1]][[i]]$mse_mcse,
                        measure_list[[2]][[i]]$mse_mcse,
                        measure_list[[3]][[i]]$mse_mcse)
      }
      # Format latex table 
      return(kable(tab, format = "latex", row.names = F, 
                   align = rep('c',7), booktabs=T, 
                   escape = FALSE, linesep = "\\addlinespace",
                   caption = "Monte Carlo standard errors of bias and MSE in Scenario II.", 
                   label = "scen2_mcse") %>% 
               add_header_above(c(" ", "Bias" = 3, "MSE" = 3)) %>% 
               kableExtra::group_rows("$\\theta_{1|0}$", 1,9, escape = FALSE))
    }
      
  }
}


