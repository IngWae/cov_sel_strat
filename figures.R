# list_to_df():
# Function to transform sample list into a data frame for plotting 
# The list is formatted into a long table 
list_to_df <- function(est_list, tr){
   names(est_list)[7] <- "AIPW*"
   # List of estimators 
    estimators <- c("UA", "ALL", "SS", "UPF",   
                    "RI", "AIPW", "AIPW*", 
                    "npRI", "npAIPW")
    
    # Initiate data matrix with the first columns 
    out_df <- cbind("UA", unname(est_list[["UA"]]))
    
    # Loop through estimators to unlist estimates and add to data frame  
    for(j in 2:length(estimators)){
        out_df <- rbind(out_df, cbind(estimators[[j]], unname(est_list[[estimators[[j]]]])))
    }
    if(tr == "dmt"){
      # Set column names 
      colnames(out_df) <- c("estimator", "est1", "est2", "est3")
      out_df <- data.frame(out_df)%>% 
        mutate(est1 = as.numeric(est1),
               est2 = as.numeric(est2),
               est3 = as.numeric(est3))
    } else if(tr == "smoking"){
      # Set column names
      colnames(out_df) <- c("estimator", "est1")
      out_df <- data.frame(out_df) %>% 
        mutate(est1 = as.numeric(est1))
    }
  # Return data frame
  return(out_df)
}

###################
### PLOT THEMES ###
###################
# boxplot_theme(): 
# Function that specifies the theme for the boxplots 
boxplot_theme <- function(){ 
  theme(text = element_text(size=11),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.text=element_text(size=10))
}

# overlap_theme(): 
# Function that specifies the theme for the overlap (density) plots
overlap_theme <- function(){
  theme(text = element_text(size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key=element_blank(), 
        legend.title = element_blank(), 
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill=NA),
        legend.text=element_text(size=12))
}


################
### BOXPLOTS ###
################

# boxplot_function():
# Function to create boxplots from a list of estimates 
boxplot_function <- function(est_list, tr){
  
  if(tr == "dmt"){
    # Create data frame on correct format 
    temp_df <-list_to_df(est_list, tr) %>% 
      mutate(estimator =  factor(estimator, 
                                 levels = c("UA", "ALL", "SS", "UPF", 
                                            "RI", "AIPW", 
                                            "AIPW*", "npRI", "npAIPW")))
    
    # Create the plot of estimates of theta1
    plot_theta1 <- ggplot(data = temp_df, aes(x=estimator, y = est1)) + 
      geom_boxplot(outlier.shape = 1, lwd=0.3, outlier.size = 1, 
                   outlier.color = "black", outlier.fill = "white")+ 
      boxplot_theme()  + theme(plot.margin=grid::unit(c(0,0,0,0),"cm"),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())+
      ylab(" ") + xlab(" ") + 
      geom_hline(yintercept = 1, lty = "dashed") 
    
    # Create the plot of estimates of theta2
    plot_theta2 <- ggplot(data = temp_df, aes(x=estimator, y = est2)) + 
      geom_boxplot(outlier.shape = 1, lwd=0.3, outlier.size = 1, 
                   outlier.color = "black", outlier.fill = "white")+ 
      boxplot_theme() +  theme(plot.margin=grid::unit(c(0,0,0,0),"cm"),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())+
      ylab( " ") + xlab(" ")+ 
      geom_hline(yintercept = 1.065789, lty = "dashed") 
    
    # Create the plot of estimates of theta3
    plot_theta3 <- ggplot(data = temp_df, aes(x=estimator, y = est3)) + 
      geom_boxplot(outlier.shape = 1, lwd=0.3, outlier.size = 1, 
                   outlier.color = "black", outlier.fill = "white")+ 
      boxplot_theme() +  theme(plot.margin=grid::unit(c(0,0,0,0),"cm"))+
      ylab( " ") + xlab(" ")+
      geom_hline(yintercept = 1.421053, lty = "dashed") 
    
    
    return(list(plot_theta1, plot_theta2, plot_theta3))
    
  } else if(tr == "smoking"){
    # Crete data frame on the correct format 
    temp_df <- list_to_df(est_list, tr)  %>% 
      mutate(estimator =  factor(estimator, 
                                 levels = c("UA", "ALL", "SS", "UPF", 
                                            "RI", "AIPW", 
                                            "AIPW*", "npRI", "npAIPW")))
    
    # Create boxplot of estimates of the odds ratio 
    plot_theta1 <- ggplot(data = temp_df, aes(x=estimator, y = est1)) + 
      geom_boxplot(outlier.shape = 1, lwd=0.3, outlier.size = 1, 
                   outlier.color = "black", outlier.fill = "white")+ 
      boxplot_theme() + 
      ylab(" ") + xlab(" ") + 
      geom_hline(yintercept = 1.224647, lty = "dashed") 
    
    return(plot_theta1)
  }
  
}

# boxplot_panel():
# Function to create panel of boxplots created with boxplot_function. 
boxplot_panel <- function(est_list, tr){
  # Create list of boxplots (one entry per outcome model)
  plot_list <- list()
  for(i in 1:3){
    plot_list[[i]] <- boxplot_function(est_list = est_list[[i]], tr = tr)
  }
  # Create grid of plots under treatment DMT 
  if(tr == "dmt"){
     return(grid.arrange(arrangeGrob(arrangeGrob(plot_list[[1]][[1]] + ylab(TeX(sprintf("$\\hat{\\theta}_{1|0}$")))+ 
                              scale_y_continuous(limits = c(0,3), breaks= seq(0,3,by=0.5)),
                              top = textGrob( "Outcome Model A", 
                              gp = gpar(fontsize = 11))),
                              arrangeGrob(plot_list[[2]][[1]] + 
                                scale_y_continuous(limits = c(0,3), breaks=NULL),
                              top = textGrob( "Outcome Model B", 
                              gp = gpar(fontsize = 11))),
                              arrangeGrob(plot_list[[3]][[1]]+ 
                                scale_y_continuous(limits = c(0,3), breaks=NULL),
                              top = textGrob( "Outcome Model C", 
                              gp = gpar(fontsize = 11))),
                              left = textGrob( "No risk vs. No treatment", rot=90,
                              gp = gpar(fontsize = 10)), nrow=1),
                  arrangeGrob(plot_list[[1]][[2]] + ylab(TeX(sprintf("$\\hat{\\theta}_{2|0}$")))+ 
                                scale_y_continuous(limits = c(0.5,2.5), breaks=seq(0.5,2.5,by=0.5)),
                              plot_list[[2]][[2]]+ 
                                scale_y_continuous(limits = c(0.5,2.5), breaks=NULL),
                              plot_list[[3]][[2]]+ 
                                scale_y_continuous(limits = c(0.5,2.5), breaks=NULL),
                              left = textGrob( "Low risk vs. No treatment", rot=90,
                              gp = gpar(fontsize = 10)), nrow=1),
                  arrangeGrob(plot_list[[1]][[3]] + ylab(TeX(sprintf("$\\hat{\\theta}_{3|0}$"))) + 
                                scale_y_continuous(limits = c(0,4), breaks =seq(0,4,by=0.5)),
                              plot_list[[2]][[3]]+ 
                                scale_y_continuous(limits = c(0,4), breaks=NULL),
                              plot_list[[3]][[3]]+ 
                                scale_y_continuous(limits = c(0,4), breaks=NULL),
                              left = textGrob( "Moderate to high risk vs. No treatment", rot=90,
                              gp = gpar(fontsize = 10)), nrow=1),
                  nrow=3, heights = c(2.1/7,1.9/7,3/7)))
  } # Create grid of plots under treatment Smoking
  else if(tr == "smoking"){
    return(grid.arrange(arrangeGrob(plot_list[[1]] +ylab(TeX(sprintf("$\\hat{\\theta}_{1|0}$"))) +
                                                    ggtitle(paste("Outcome Model A")) + 
                                      scale_y_continuous(limits = c(0.5,3), breaks = seq(0,3,by=0.5)), 
                                           plot_list[[2]] + ggtitle(paste("Outcome Model B"))+ 
                                      scale_y_continuous(limits = c(0.5,3), breaks=NULL), 
                                           plot_list[[3]] + ggtitle(paste("Outcome Model C"))+ 
                                      scale_y_continuous(limits = c(0.5,3), breaks=NULL), 
                             left = textGrob( "Smoker vs. Non-Smoker", rot=90,
                                              gp = gpar(fontsize = 10)),
                             nrow=1)))
  }
}

######################
### OVERLAP PLOTS ####
######################
# overlap_plot(): 
# Function to create density plots for overlap assessment 
# Takes a sample with associated propensity score and treatment of 
# interest as arguments 
overlap_plot <- function(samp, ps, tr){
  
  if(tr == "dmt"){
    # Create data frame with the correct format for plotting 
    temp_df <- data.frame(ps = ps, 
                          No_treatment = samp$No_treatment,
                          No_risk = samp$No_risk,
                          Low_risk = samp$Low_risk,
                          Moderate_to_high_risk = samp$Moderate_to_high_risk) %>% 
      mutate(No_treatment = factor(No_treatment, levels = c(0,1)),
             No_risk = factor(No_risk, levels = c(0,1)),
             Low_risk = factor(Low_risk, levels = c(0,1)),
             Moderate_to_high_risk = factor(Moderate_to_high_risk, levels = c(0,1)))
    
    # Density plot for the first treatment level (No treatment)
    plot_tr0 <- ggplot(data = temp_df, aes(x = ps.No_treatment, linetype=No_treatment)) + 
      geom_density(lwd = 0.3,key_glyph = draw_key_path) + 
      overlap_theme() +
      ylab("Density") + xlab("Propensity score") +
      scale_linetype_discrete(limits = c("0", "1"), labels = c("Other","No treatment")) +
      scale_x_continuous(breaks = c(seq(0,1, by = 0.2)), limits = c(0, 1))+
      scale_y_continuous(breaks=seq(0,20, by=5), limits=c(0,20))
    
    # Density plot for the second treatment level (No risk)
    plot_tr1 <- ggplot(data = temp_df, aes(x = ps.No_risk, linetype=No_risk)) + 
      geom_density(lwd = 0.3,key_glyph = draw_key_path) + 
      overlap_theme() +
      ylab("Density") + xlab("Propensity score") +
      scale_linetype_discrete(limits = c("0", "1"), labels = c("Other","No risk")) +
      scale_x_continuous(breaks = c(seq(0,1, by = 0.2)), limits = c(0, 1)) +
      scale_y_continuous(breaks=seq(0,20, by=5), limits=c(0,20))
    
    # Density plot for the third treatment level (Low risk)
    plot_tr2 <- ggplot(data = temp_df, aes(x = ps.Low_risk, linetype=Low_risk)) + 
      geom_density(lwd = 0.3,key_glyph = draw_key_path) + 
      overlap_theme() +
      ylab("Density") + xlab("Propensity score") +
      scale_linetype_discrete(limits = c("0", "1"), labels = c("Other","Low risk")) +
      scale_x_continuous(breaks = c(seq(0,1, by = 0.2)), limits = c(0, 1))+
      scale_y_continuous(breaks=seq(0,15, by=5), limits=c(0,15))
   
     # Density plot for the fourth treatment level (Moderate to high risk)
    plot_tr3 <- ggplot(data = temp_df, aes(x = ps.Moderate_to_high_risk, linetype=Moderate_to_high_risk)) + 
      geom_density(lwd = 0.3,key_glyph = draw_key_path) + 
      overlap_theme() +
      ylab("Density") + xlab("Propensity score") +
      scale_linetype_discrete(limits = c("0", "1"), labels = c("Other","Moderate to high risk")) +
      scale_x_continuous(breaks = c(seq(0,1, by = 0.2)), limits = c(0, 1)) +
      scale_y_continuous(breaks=seq(0,10, by=1), limits=c(0,10))
    
    
    # Return all plots in a list
    return(list(plot_tr0, plot_tr1, plot_tr2, plot_tr3))
    
  } else if(tr == "smoking"){
    # Create data frame with the correct format for plotting 
    temp_df <- data.frame(ps=ps, smoking = samp$smoking) %>% 
      mutate(smoking = factor(smoking, levels = c(0,1)))
    
    # Create density plot for Smokers versus Non-Smokers 
    plot_tr1 <- ggplot(data = temp_df, aes(x = ps, linetype=smoking)) + 
      geom_density(lwd = 0.3,key_glyph = draw_key_path) + 
      overlap_theme() +
      ylab("Density") + xlab("Propensity score") +
      scale_linetype_discrete(limits = c("0", "1"), labels = c("Non-smoker","Smoker")) +
      scale_x_continuous(breaks = seq(0,1, by=0.1), limits = c(0, 1)) #+
     # scale_y_continuous(breaks=seq(0,20, by=5), limits=c(0,20))
    
    # Return plot 
    return(list(plot_tr1))
  }
}

# overlap_panel(): 
# For treatment DMT, overlap plots are arranged in a 2x2 grid 
overlap_panel <- function(samp, ps, tr){
  # Create plots 
  plot_list <- overlap_plot(samp=samp, ps=ps, tr=tr)
  # Create grid 
  return(grid.arrange(plot_list[[1]], plot_list[[2]], 
                                plot_list[[3]], plot_list[[4]], nrow = 2))

}







