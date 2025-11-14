.libPaths(c("/home/users/galanter/Documents/R/win-library/4.2 cluster"))


library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(viridis)

sim_truth <- readRDS("truth.rds")



sim <- readRDS("pre_specified_sim.rds")


generate_results <- function(sim, sim_truth){
  
  sim_vals <- do.call(rbind,sapply(1:length(sim$results_complex),
                                   function(x) sim$results_complex[[x]][2]))
  
  
  sim_vals <- data.frame(select(sim$results,c("n","prop_adjust",
                                              "outcome_adjust",
                                              "adjust_set",  
                                               "scenario",
                                              "misspec")),
                         sim_vals)
  
  names(sim_vals)[7:8] <- 1:2
  
  names(sim_vals)[6] <- "Misspecification"

  sim_vals <- pivot_longer(sim_vals,7:8,names_to  = "Variable",values_to = "Value")
  
  sim_vals$Variable <- as.integer(sim_vals$Variable)
  
  sim_vals$n <- recode(sim_vals$n,
    "500" = "N = 500",
    "1000" = "N = 1000",
    "10000" = "N = 10000"
  )
  
  sim_vals$n <- factor(sim_vals$n, levels= c("N = 500", "N = 1000","N = 10000"))
  
  sim_vals$Set <- sim_vals$adjust_set
  
  sim_vals$Set <- recode(sim_vals$Set,
    "Confounders" = "C",
    "Precision" = "C+P",
    "Instrument" = "C+I",
    "All" = "A"
  )
  
  sim_vals$Set <- factor(sim_vals$Set,
                         c("A","C+P","C","C+I"))
  
  
  #formatting truth values
  
  truth <- sim_truth$results[,-c(1:6,9)]
  
  names(truth)[3:4] <- 1:2
  
  names(truth)[2] <- "Misspecification"
  
  truth <- pivot_longer(truth,3:4,names_to = "Variable",values_to = "Truth") 
  
  truth$Variable <- as.integer(truth$Variable)
  
  #getting final results
  
  sim_vals <- left_join(sim_vals,truth,by = c("scenario",
                                              "Variable","Misspecification"))
  
  sim_vals$diff <- sim_vals$Value - sim_vals$Truth 
  
  sim_vals$diff_2 <- sim_vals$diff^2
  
  
  sim_vals$prop_adjust <- ifelse(sim_vals$prop_adjust,
                                 "Propensity Yes","Propensity No")
  
  sim_vals$outcome_adjust <- ifelse(sim_vals$outcome_adjust == "Oracle",
                                    "Outcome Model Yes","Outcome Model No")
  
  sim_vals <- sim_vals %>% unite(Adjustment,c(prop_adjust,outcome_adjust),
                                 sep = "\n")
  
  print(names(sim_vals))
  
  sim_summary <- sim_vals %>% group_by(n,Adjustment,Set,
                                       scenario,Variable,Misspecification) %>%
    dplyr::summarise(
      MSE = mean(diff_2),
      Bias = mean(diff),
      SD = sd(Value)
    )
  
  
  
   for(m in c("sqr", "sin", "noncenter","none")){


    for(i in 1:2){

      pdf(paste0("pre_specified_misspec_",m,"_value_boxplots_var",i,".pdf"),
          height = 10, width = 12)
      print(ggplot(data = sim_vals[sim_vals$scenario == "C" &
                                     sim_vals$Variable == i &
                                     sim_vals$Misspecification == m,],
                   aes(y = Value,x = Set,color = Set)) +
              geom_boxplot()+ geom_hline(aes(yintercept = Truth),linetype = "dotted",
                                         color = "black") +
              facet_grid(rows = ggplot2::vars(n),cols = ggplot2::vars(Adjustment) ) +
              theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                    text = element_text(size=20),)) +
        scale_color_viridis(discrete = T)
      dev.off()

    }
    
    sim_summary_temp <- filter(sim_summary,Misspecification == m)
    
    print(names(sim_summary_temp))
    
    sim_tab <- pivot_wider(sim_summary_temp, names_from = Adjustment,values_from = c(MSE,Bias,SD))
    
    sim_tab <- sim_tab[,c(-5)]
    
    sim_tab <- arrange(sim_tab, scenario,Variable,n)
    
    sink(paste0("_misspec_",m,"pre_specified_simulation_tables.txt"))

    
    
    print("Scenario C")
    
    for(i in 1:2){
      
      print(kable(sim_tab[sim_tab$scenario == "C" & sim_tab$Variable == i,-c(1)],format = "latex",digits = 2))
      
    }
    
    print("Scenario A")
    
    for(i in 1:2){
      
      print(kable(sim_tab[sim_tab$scenario == "A" & sim_tab$Variable == i,-c(1)],format = "latex",digits = 2))
      
    }
    
    sink()
    
  }
  
}


propensity_plots <- function(sim, r){

  sim_props <-lapply(1:length(sim$results_complex),
                     function(x) sim$results_complex[[x]][4][[1]][1:r,])
    
  sim_props <- do.call(rbind,sim_props)

  
  dat <- select(sim$results,c("n","prop_adjust",
                                              "adjust_set",  
                                              "scenario",
                                              "misspec"))
  
  dat$n <- recode(dat$n,
                       "500" = "N = 500",
                       "1000" = "N = 1000",
                       "10000" = "N = 10000"
  )
  
  dat$n <- factor(dat$n, levels= c("N = 500", "N = 1000","N = 10000"))
  
  dat <- dat %>% rename(Set = adjust_set)
  
  dat$Set <- recode(dat$Set,
                         "Confounders" = "C",
                         "Precision" = "C+P",
                         "Instrument" = "C+I",
                         "All" = "A"
  )
  
  dat$Set <- factor(dat$Set,
                         c("A","C+P","C","C+I"))
  
  
  temp <- lapply(1:nrow(dat),function(x) dat[rep(x,r),])
  
  
  temp <- do.call(rbind,temp)

  
  sim_props <- cbind(temp,sim_props)
  
  names(sim_props)[6:7] <- 1:2
  
  sim_props <- gather(sim_props, key = "Variable", value = "Value",
                      6,7)
  
  sim_props$Variable <- as.integer(sim_props$Variable)
  
  names(sim_props)[5] <- "Misspecification"
  
  names(sim_props)[4] <- "scenario"
  
  names(sim_props)[1] <- "n"
  
  names(sim_props)[3] <- "Set"
  
  for(m in c("sqr", "sin", "noncenter","none")){
    
    for(i in 1:2){
      
      dat <- filter(sim_props,
                    scenario == "C" , 
                    Variable == i ,
                    Misspecification == m)
      
      jpeg(paste0("pre_specified_prop_boxplots_C_misspec_",m,"_var",i,".jpg"),
           height = 10, width = 12, 
           units = "in", res = 72)
      print(ggplot(data = dat,
                   aes(y = Value,x = Set,color = Set)) + 
              geom_boxplot() +facet_grid(rows = ggplot2::vars(n))+
              theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                    text = element_text(size=20),))
      dev.off()
      
    }
    
    
    sim_summary <- filter(sim_props,Misspecification == m) %>% 
      group_by(Variable, Set,n) %>%
      dplyr::summarise(
        min = min(Value),
        Q1 = quantile(Value, probs = 0.25),
        median = median(Value),
        Q4 = quantile(Value, probs = 0.75),
        max = max(Value)
      )
    
    sink(paste0("pre_specfied_prop_summaries_misspec",m,".txt"))
    
    print(kable(sim_summary,format = "latex"))
    
    sink()
    
  }
  
}






#try(generate_results(sim, sim_truth))

try(propensity_plots(sim,5))






