
library( SimEngine)

# Data generation ------
#simulate the data from p_n non tailoring and p_t tailoring variables
#arguments
#   n - number of observations to generate
#   p_n - number of non tailoring covariates
#   p_t - number of tailoring covariates
#   tau - matrix of coefficient of non-tailoring covs on tailoring covs
#    rows are non tailoring and columns are tailoring (p_n x p_t)
#   beta - coefficients of main effects of covariates in outcome model
#    in order intercept, nontailoring, tailoring
#   alpha - coefficients in treatment model
#    in order intercept, nontailoring, tailoring
#   gamma -  tx effect coefficients (ie contrast model)
#    in order intercept, tailoring
#   theta0 - paramater which controls optimal tailoring levels - length p_t
#   smooth - boolean indicator for shape of tailoring term
#   sd.error - positive standard deviation of outcome around mean
#   c - posisitve steepness coefficient for smooth tailoring model
#returns 
#   dat - data matrix of id, tx, covariates, tx, outcome, and counterfactuals
data.gen = function(n, p_n, p_t, tau, beta, alpha, gamma, theta0, smooth, 
                    sd.error=0.7, c = 1, misspec) {
  
  
  
  #get covariates
  
  X_nt = matrix(rnorm(p_n*n, sd = 1), nrow = n)
  
  X_t_means = X_nt %*% tau
  
  X_t = matrix(rnorm(p_t*n, mean = X_t_means,sd = 1), nrow = n)
  
  X = cbind(X_nt,X_t)
  
  X_temp <- X
  
  #transforming precision variables
  if(L$misspec == "sqr") 
  {X_temp[,2:3] <- c(X[,2]^2,X[,3]^2)
  } else if (L$misspec == "sin"){
    X_temp[,2:3] <- c(sin(X[,2]),sin(X[,3]))
  }else if (L$misspec == "noncenter") {
    X_temp[,2:3] <- c((X[,2]+0.5)^2,(X[,3]+0.5)^2)}
  
  
  #covariate names
  colnames(X) <- paste0("x",1:(p_n+p_t))
  
  #error for outcomes
  noise = rnorm(n, mean=0, sd=sd.error)
  
  #function relating theta, the tailoring variables, and the tx effect
  if (smooth==TRUE) {
    expit = function(x) 1/(1+exp(-c*x))
  }else {
    expit = function(x) as.numeric(x>0)
  }
  
  #get tx probabilies
  p = 1/(1+exp(- cbind(1,X) %*% alpha))
  
  #get observed tx
  A = rbinom(n, size=1, prob=p)
  
  #conditional means
  mu_0 = cbind(1,X_temp) %*% beta
  
  mu_1 = mu_0 + cbind(1, sapply(1:p_t,function(x) 
    expit(theta0[x]-X_t[,x]) ))  %*% gamma
  
  #potential outcomes
  Y_0 = mu_0 + noise
  Y_1 = mu_1 + noise
  Y = A * Y_1+ (1 - A) * Y_0
  dat = data.frame(id=1:n, A, X,
                   Y=Y, Y_0 = Y_0, Y_1 = Y_1)
  
  return(dat)
}


# Simulation functions -----

marginal_sim <- new_sim()

marginal_sim %<>% set_levels(
  n = c(1000000),
  prop_adjust = c(T),
  outcome_adjust = c("None"),
  scenario = c("C"),
  misspec = c("none","sqr","sin","noncenter")
  
)





marginal_sim %<>% set_script(function() {
  
  
  options(dplyr.summarise.inform = FALSE)
  
  tailor_search_inds = 80+6:7
  
  gamma = c(-15, 17, 8)
  theta0 = c(0.84,-0.84)
  
  #getting parameters
  
  beta <- switch(L$scenario,
                 "A" = c(6,-5,-5,-5,0,0,rep(0,80),rep(0,2)),
                 "C" = c(6,0,-5,10,0,0,rep(0,80),5,0))
  
  alpha <- switch(L$scenario,
                  "A" = c(0,0.5,0,0,0.5,0.5,rep(0,80),rep(0,2)),
                  "C" = c(0,rep(0,3),0.5,-0.25,rep(0,80),-0.5,0))
  
  tau <- matrix(0, nrow = 85, ncol = 2)
  
  if(L$scenario == "A" ) 
  {tau[1,1] <- 0.25}
  
  adjust_inds <- 1:(length(alpha) -1)
  
  if (L$prop_adjust){
    prop_inds <- adjust_inds
  } else { prop_inds <- NULL}
  
  if (L$outcome_adjust == "Lasso" ){
    out_inds <- adjust_inds
  } else if(L$outcome_adjust == "Oracle"){
    if(L$scenario == "A") {out_inds <- 1:3} else {out_inds <- c(2:3,86)}
  } else { out_inds <- NULL }
  
  #number of tailoring variables being considered in analysis
  t_num <- length(tailor_search_inds)
  
  Theta_vec = seq(-2, 2, by=0.2)
  #7 knots
  knots_vec = c(-1.8,-1.2,-.6,0,.6,1.2,1.8)
  
  set.seed(L$sim_id)
  
  dat = data.gen(n = L$n , p_n = 85, p_t = 2, tau = tau,
                 beta = beta, alpha = alpha, gamma = gamma,
                 theta0 = theta0, smooth = T, 
                 sd.error = 0.7, c = 5, misspec = L$misspec)
  
  
  
  Theta = seq(-2, 2, by=0.2)
  
  truth = data.frame(Theta,0,0)
  
  ret <- rep(NA,2)
  
  names(truth)<- c("Theta","T1","T2")
  
  #get counterfactual average outcomes for each rule
  
  for(k in 1:2){
    for (i in 1:length(Theta)) {
      Y = ifelse(dat[,(k+87)] < Theta[i], dat$Y_1, dat$Y_0)
      truth[i,(k+1)] = mean(Y)
    }
    ret[k]<- (max(truth[,k+1]))
  }
  
  truth2 <- pivot_longer(truth, cols = 2:3, names_to = "Variable", 
                         values_to = "Value")
  
  
  pdf(paste0(L$misspec,"true_marginal_plot.pdf"),height = 10, width = 12)
  print(ggplot( data = truth2, 
                aes(x = Theta, y = Value, color = Variable, linetype = Variable)) +
          geom_line())
  dev.off()
  
  
  return(list(value1 = ret[1],
              value2 = ret[2]))
  
  
})

marginal_sim %<>% set_config(
  num_sim = 1,
  packages = c("dplyr","tidyr","ggplot2"))



marginal_sim %<>% run()