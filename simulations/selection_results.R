.libPaths(c("/home/users/galanter/Documents/R/win-library/4.2 cluster"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

sim_truth <- readRDS("truth.rds")

sim <- readRDS("selection_sim.rds")


generate_results <- function(sim, sim_truth){

  sim_vals <- do.call(rbind,sapply(1:length(sim$results_complex),
                                   function(x) sim$results_complex[[x]][2]))

  sim_vals <- data.frame(sim$results[,4:9],sim_vals)

  names(sim_vals)[7:8] <- 1:2

  names(sim_vals)[5] <- "Misspecification"

  names(sim_vals)[6] <- "Method"


  sim_vals <- pivot_longer(sim_vals,7:8,names_to  = "Variable",values_to = "Value")

  sim_vals$Variable <- as.integer(sim_vals$Variable)

  #formatting truth values

  truth <- sim_truth$results[,-c(1:6,9)]

  names(truth)[3:4] <- 1:2

  names(truth)[2] <- "Misspecification"

  truth <- pivot_longer(truth,3:4,names_to = "Variable",values_to = "Truth")

  truth$Variable <- as.integer(truth$Variable)

  #getting final results

  sim_vals <- left_join(sim_vals,truth,by = c("scenario","Variable","Misspecification"))

  sim_vals$diff <- sim_vals$Value - sim_vals$Truth

  sim_vals$diff_2 <- sim_vals$diff^2

  sim_vals$prop_adjust <- ifelse(sim_vals$prop_adjust,"Y","N")

  sim_vals$Strategy <- paste0("Outcome: ",sim_vals$outcome_adjust)

  table(sim_vals$Strategy)

  sim_vals$n <- recode(sim_vals$n,
                  "500" = "N = 500",
                  "1000" = "N = 1000",
                  "10000" = "N = 10000"
  )

  table(sim_vals$n)

  sim_vals$n <- factor(sim_vals$n, levels= c("N = 500", "N = 1000","N = 10000"))


  sim_vals$Method <- recode(sim_vals$Method,
                            "no_selection" = "A",
                            "oal" = "OAL",
                            "glider" = "GLiDeR",
                            "dip" = "DiPS",
                            "oracle" = "C+P",
                            "causal_ball" = "CB",
                            "hd_balancing" = "HDCB")

  sim_vals$Method <- factor(sim_vals$Method,
                            levels = c("A", "C+P", "OAL",
                                       "GLiDeR", "DiPS", "HDCB", "CB"))


  sim_summary <- sim_vals %>% group_by(n,Strategy,scenario,Variable, Method,Misspecification) %>%
    dplyr::summarise(
      MSE = mean(diff_2),
      Bias = mean(diff),
      SD = sd(Value)
    )



  for(m in c("noncenter","none")){


    for(i in 1:2){

      pdf(paste0("selection_misspec_",m,"_value_boxplots_C_var",i,".pdf"),height = 10, width = 12)
      print(ggplot(data = sim_vals[sim_vals$scenario == "C" & sim_vals$Variable == i &
                                     sim_vals$Misspecification == m,],
                   aes(y = Value,x = Method,color = Method)) +
              geom_boxplot()+ geom_hline(aes(yintercept = Truth),linetype = "dotted",
                                         color = "black") +
              facet_grid(rows = ggplot2::vars(n),cols = ggplot2::vars(Strategy) ) +
              theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                    text = element_text(size=20),))
      dev.off()

    }

    sim_summary_temp <- filter(sim_summary,Misspecification == m)

    sim_tab <- pivot_wider(sim_summary_temp, names_from = Strategy,values_from = c(MSE,Bias,SD))

    sim_tab <- sim_tab[,c(-5)]

    sim_tab <- arrange(sim_tab, scenario,Variable,n,Method)

    sink(paste0("selection_misspec_",m,"_simulation_v8_tables.txt"))


    print("Scenario C")

    for(i in 1:2){

      print(kable(sim_tab[sim_tab$scenario == "C" & sim_tab$Variable == i,-c(1)],format = "latex",digits = 2))

    }

    sink()

  }

}


propensity_plots <- function(sim,r){
  
  saveRDS(sim$results_complex[[3]][4],"sample_props")
  
  sim_props <- lapply(1:nrow(sim$results),function(x) sim$results_complex[[x]][4][[1]][1:r,])
  
  sim_props <- do.call(rbind,sim_props)
  
  dat <- data.frame(sim$results[,c(4,7,8,9)])
  
  temp <- lapply(1:nrow(dat),function(x) dat[rep(x,r),])
  
  temp <- do.call(rbind,temp)
  
  sim_props <- cbind(temp,sim_props)
  
  
  
  names(sim_props)[5:6] <- 1:2
  
  sim_props <- pivot_longer(sim_props, names_to = "Variable", values_to = "Value",
                            cols = 5:6)
  
  sim_props$Variable <- as.integer(sim_props$Variable)
  
  names(sim_props)[3] <- "Misspecification"
  
  names(sim_props)[4] <- "Method"
  
  names(sim_props)[2] <- "scenario"
  
  names(sim_props)[1] <- "n"
  
  sim_props$n <- recode(sim_props$n,
                       "500" = "N = 500",
                       "1000" = "N = 1000",
                       "10000" = "N = 10000"
  )
  
  sim_props$n <- factor(sim_props$n, levels= c("N = 500", "N = 1000","N = 10000"))
  
  
  sim_props$Method <- recode(sim_props$Method,
                            "no_selection" = "A",
                            "oal" = "OAL",
                            "glider" = "GLiDeR",
                            "dip" = "DiPS",
                            "oracle" = "C+P",
                            "causal_ball" = "CB",
                            "hd_balancing" = "HDCB")
  
  sim_props$Method <- factor(sim_props$Method,
                            levels = c("A", "C+P", "OAL", 
                                       "GLiDeR", "DiPS", "HDCB", "CB"))
  
  
  for(m in c("noncenter","none")){
    
    for(i in 1:2){
      
      dat <- filter(sim_props,scenario == "C" , 
                    Variable == i ,
                    Misspecification == m)
      
      
      jpeg(paste0("selection_prop_boxplots_C_misspec_",m,"_var",i,".jpg"),height = 10, width = 12, 
           units = "in", res = 72)
      print(ggplot(data = dat,
                   aes(y = Value,x = Method,color = Method)) + 
              geom_boxplot() +facet_grid(rows = ggplot2::vars(n))+
              theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                    text = element_text(size=20),))
      dev.off()
      
    }
    
    
    sim_summary <- filter(sim_props,Misspecification == m) %>% group_by(Variable, Method,n) %>%
      dplyr::summarise(
        min = min(Value),
        Q1 = quantile(Value, probs = 0.25),
        median = median(Value),
        Q4 = quantile(Value, probs = 0.75),
        max = max(Value)
      )
    
    sink(paste0("selection_prop_summaries","misspec",m,".txt"))
    
    print(kable(sim_summary,format = "latex"))
    
    sink()
    
  }
  
}


selection_plots <- function(sim){
  
  sim_select <- lapply(1:nrow(sim$results),function(x) sim$results_complex[[x]][6])
  
  helper <- function(x){
    
    x <- c(x[[1]][,1],x[[1]][,2])
  }
  
  types <- rep("E",87)
  types[86] <- "C"
  types[c(2,3,87)] <- "P"
  types[c(4,5)] <- "I"
    
  types <- factor(types,levels = c("C","P","I","E"))
   
  
  sim_select <- lapply(sim_select,helper)

  
  sim_select <- lapply(sim_select, function(x) data.frame(sort(rep(1:2,87)),
                                                          rep(1:87,2),
                                                          rep(types,2),x))
  sim_select <- do.call(rbind,sim_select)
  
  dat <- data.frame(sim$results[,c(4,7,8,9)])
  
  temp <- lapply(1:nrow(dat),function(x) dat[rep(x,2*87),])
  
  temp <- do.call(rbind,temp)
  
  sim_select <- cbind(temp,sim_select)
  
  names(sim_select) <- c("n","Scenario","Misspecification","Method",
                         "Variable","Covariate","Type","Selected")
  
  sim_select$n <- recode(sim_select$n,
                        "500" = "N = 500",
                        "1000" = "N = 1000",
                        "10000" = "N = 10000"
  )
  
  sim_select$n <- factor(sim_select$n, levels= c("N = 500", "N = 1000","N = 10000"))
  
  
  sim_select$Method <- recode(sim_select$Method,
                             "no_selection" = "A",
                             "oal" = "OAL",
                             "glider" = "GLiDeR",
                             "dip" = "DiPS",
                             "oracle" = "C+P",
                             "causal_ball" = "CB",
                             "hd_balancing" = "HDCB")
  
  sim_select$Method <- factor(sim_select$Method,
                             levels = c("A", "C+P", "OAL", 
                                        "GLiDeR", "DiPS", "HDCB", "CB"))
  
  
  sim_select <- filter(sim_select, 
                       Method %in% c("OAL",
                                     "GLiDeR", "DiPS",
                                     "CB"))
  
  sim_summary <- sim_select %>% group_by(Variable, Misspecification,Method, Covariate, Type,n) %>%
    dplyr::summarise(
      Proportion = mean(Selected)
    )
  
  for(m in c("none","noncenter")){
    
    for(i in 1:2){
      
      dat2 <- filter(sim_summary,Variable == i, Misspecification == m)
      
      pdf(paste0("selection_selected_proportions_C_misspec_",m,"_var",i,".pdf"),height = 10, width = 12)
      print(ggplot(data = dat2,
                   aes(y = Proportion,x = Covariate, color = Type)) + 
              geom_point() + ylim(0,1)+
              facet_grid(rows = ggplot2::vars(n),cols = ggplot2::vars(Method)) +
              theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                    text = element_text(size=20),))
      dev.off()
      
    }
    
  }
  
}






generate_results(sim, sim_truth)

propensity_plots(sim, 5)

selection_plots(sim)
