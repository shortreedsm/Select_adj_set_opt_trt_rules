

.libPaths(c("/home/users/galanter/Documents/R/win-library/4.2 cluster"))

library(SimEngine)

run_on_cluster(
  first = {
    
    select_sim <- new_sim()
    
    
    
    select_sim %<>% set_levels(
      n = c(500,1000, 10000),
      prop_adjust = c(T),
      outcome_adjust = c("Lasso", "None","Oracle"),
      scenario = c("C"),
      misspec = c("none", "noncenter"), #no sqr and sin
      prop.selection = c(
        "no_selection",
        "oal",
        "glider",
        "dip",
        "oracle",
        "causal_ball",
        "hd_balancing"
      )
    )
 
    
    
    
    
    
    select_sim %<>% set_script(function() {
      library(
        CBPSmod,
        lib.loc = c(
          .libPaths(),
          "/home/users/galanter/Documents/R/win-library/4.2"
        )
      )
      
      
      source("GLiDeR_Rfunctions_tailoring.R")
      source("causal_ball_script_tailoring.R")
      source("dips_tailoring.R")
      source("simulation_functions.R")
      
      
      
      options(dplyr.summarise.inform = FALSE)
      
      tailor_search_inds <- 80 + 6:7
      
      gamma <- c(-15, 17, 8)
      theta0 <- c(0.84, -0.84)
      
      #getting parameters
      
      beta <- switch(
        L$scenario,
        "A" = c(6, -5, -5, -5, 0, 0, rep(0, 80), rep(0, 2)),
        "C" = c(6, 0, -5, 10, 0, 0, rep(0, 80), 5, 0)
      )
      
      alpha <- switch(
        L$scenario,
        "A" = c(0, 0.5, 0, 0, 0.5, 0.5, rep(0, 80), rep(0, 2)),
        "C" = c(0, rep(0, 3), 0.5, -0.25, rep(0, 80), -0.5, 0)
      )
      
      tau <- matrix(0, nrow = 85, ncol = 2)
      
      if (L$scenario == "A")
      {
        tau[1, 1] <- 0.25
      }
      
      batch({
        
        df <- data.gen(
          n = L$n ,
          p_n = 85,
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
        )
        
        
      })
      
      
      
      adjust_inds <- 1:(length(alpha) - 1)
      
      if (L$prop_adjust) {
        prop_inds <- adjust_inds
      } else {
        prop_inds <- NULL
      }
      
      if (L$outcome_adjust == "Lasso") {
        out_inds <- adjust_inds
      } else if (L$outcome_adjust == "Oracle") {
        if (L$scenario == "A") {
          out_inds <- 1:3 
        } else {
          out_inds <- c(2:3, 86,87)
        }
      } else {
        out_inds <- NULL
      }
      
      #number of tailoring variables being considered in analysis
      t_num <- length(tailor_search_inds)
      
      Theta_vec <- seq(-2, 2, by = 0.2)
      #7 knots
      knots_vec <- c(-1.8, -1.2, -.6, 0, .6, 1.2, 1.8)
      
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
      
      print(L$rep_id)
      
      return(list(
        ".complex" = list(
          "ranking" = ranking,
          "value" = ret$value,
          "coefs" = ret$coefs,
          "propensities" = ret$props,
          "weights" = ret$weights,
          "selected" = ret$selected,
          "bad_props" = ret$bad_props,
          "info" = sessionInfo()
        )
      ))
      
      
    })
    
    select_sim %<>% set_config(
      num_sim = 500,
      batch_levels=c("n", "scenario","misspec"),
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
      n_cores = 1000
    )
    
  },
  
  main = {
    select_sim %<>% run()
  },
  
  last = {
  },
  
  cluster_config = list(js = "slurm")
)