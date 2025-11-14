

.libPaths(new = c("/home/users/galanter/Documents/R/win-library/4.2"))


library(SimEngine)

run_on_cluster(
  first = {
    
    select_sim <- new_sim()
    
    
    
    select_sim %<>% set_levels(
      n = c(500, 1000, 10000),
      prop_adjust = c(T,F),
      outcome_adjust = c("Oracle","None"),
      adjust_set = c("Confounders","Precision","Instrument","All"),
      scenario = c("C"),
      misspec = c("none", "sqr", "sin", "noncenter"),
      prop.selection = c(
        "no_selection"
      )
    )
    
    
    
    
    
    select_sim %<>% set_script(function() {
      
      
      source("simulation_functions.R")
      
      
      
      options(dplyr.summarise.inform = FALSE)
      
      tailor_search_inds <- 79 + 6:7
      
      gamma <- c(-15, 17, 8)
      theta0 <- c(0.84, -0.84)
      
      #getting parameters
      
      beta <- switch(
        L$scenario,
        "A" = c(6, -5, -5, -5, 0, 0, rep(0, 79), rep(0, 2)),
        "C" = c(6, 0, -5, 10, 0, 0, rep(0, 79), 5, 0)
      )
      
      alpha <- switch(
        L$scenario,
        "A" = c(0, 0.5, 0, 0, 0.5, 0.5, rep(0, 79), rep(0, 2)),
        "C" = c(0, rep(0, 3), 0.5, -0.25, rep(0, 79), -0.5, 0)
      )
      
      tau <- matrix(0, nrow = 84, ncol = 2)
      
      if (L$scenario == "A")
      {
        tau[1, 1] <- 0.25
      }
      
      #only for scenario c
      adjust_inds <- switch(L$adjust_set,
                            "Confounders" = 85,
                            "Precision" = c(2,3,85,86),
                            "Instrument" = c(4,5,85),
                            "All" = c(2,3,4,5,85,86))
      
      if (L$prop_adjust) {
        prop_inds <- adjust_inds
      } else {
        prop_inds <- NULL
      }
      
      if (L$outcome_adjust == "Oracle") {
        out_inds <- adjust_inds
      } else {
        out_inds <- NULL
      }
      
      
      
      #number of tailoring variables being considered in analysis
      t_num <- length(tailor_search_inds)
      
      Theta_vec <- seq(-2, 2, by = 0.2)
      #7 knots
      knots_vec <- c(-1.8, -1.2, -.6, 0, .6, 1.2, 1.8)
      
      
      batch({df <- data.gen(
        n = L$n ,
        p_n = 84,
        p_t = 2,
        tau = tau,
        beta = beta,
        alpha = alpha,
        gamma = gamma,
        theta0 = theta0,
        smooth = T,
        sd.error = 0.7,
        c = 5,
        misspec = L$misspec
      )})
      
      
      #temporary storage for the list of info for this rep
      ret <- get_all_estimates(
        dat = df,
        Theta_vec = Theta_vec,
        knots_vec = knots_vec,
        tailor_search_inds = tailor_search_inds,
        prop_adjust = L$prop_adjust,
        prop_adjust_inds = prop_inds,
        outcome_adjust = L$outcome_adjust,
        outcome_adjust_inds = out_inds
      )
      
      ranking <- rank(-ret$value)
      
      options(dplyr.summarise.inform = TRUE)
      
      print("rep")
      
      return(list(
        ".complex" = list(
          "ranking" = ranking,
          "value" = ret$value,
          "coefs" = ret$coefs,
          "propensities" = ret$props,
          "weights" = ret$weights,
          "selected" = ret$selected,
          "info" = sessionInfo()
        )
      ))
      
      
    })
    
    select_sim %<>% set_config(
      num_sim = 1000,
      packages = c(
        "splines",
        "dplyr",
        "tidyr",
        "glmnet",
        "R.utils",
        "Ball",
        "grpreg",
        "glmpath",
        "np"
      ),
      n_cores = 1000,
      batch_levels=c("n","scenario","misspec"), 
      return_batch_id=T
    )
    
  },
  
  main = {
    select_sim %<>% run()
  },
  
  last = {
  },
  
  cluster_config = list(js = "slurm")
)