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
data.gen <-
  function(n,
           p_n,
           p_t,
           tau,
           beta,
           alpha,
           gamma,
           theta0,
           smooth,
           sd.error = 0.7,
           c = 1,
           misspec) {

    #get covariates

    X_nt <- matrix(rnorm(p_n * n, sd = 1), nrow = n)

    X_t_means <- X_nt %*% tau

    X_t <- matrix(rnorm(p_t * n, mean = X_t_means, sd = 1), nrow = n)

    X <- cbind(X_nt, X_t)

    X_temp <- X

    #transforming precision variables
    if (L$misspec == "sqr")
    {
      X_temp[, 2:3] <- c(X[, 2] ^ 2, X[, 3] ^ 2)
    } else if (L$misspec == "sin") {
      X_temp[, 2:3] <- c(sin(X[, 2]), sin(X[, 3]))
    } else if (L$misspec == "noncenter") {
      X_temp[, 2:3] <- c((X[, 2] + 0.5) ^ 2, (X[, 3] + 0.5) ^ 2)
    }


    #covariate names
    colnames(X) <- paste0("x", 1:(p_n + p_t))

    #error for outcomes
    noise <- rnorm(n, mean = 0, sd = sd.error)

    #function relating theta, the tailoring variables, and the tx effect
    if (smooth == TRUE) {
      expit <- function(x)
        1 / (1 + exp(-c * x))
    } else {
      expit <- function(x)
        as.numeric(x > 0)
    }

    #get tx probabilies
    p <- 1 / (1 + exp(-cbind(1, X) %*% alpha))

    #get observed tx
    A <- rbinom(n, size = 1, prob = p)

    #conditional means
    mu_0 <- cbind(1, X_temp) %*% beta

    mu_1 <- mu_0 + cbind(1, sapply(1:p_t, function(x)
      expit(theta0[x] - X_t[, x])))  %*% gamma

    #potential outcomes
    Y_0 <- mu_0 + noise
    Y_1 <- mu_1 + noise
    Y <- A * Y_1 + (1 - A) * Y_0
    dat <- data.frame(
      id = 1:n,
      A,
      X,
      Y = Y,
      Y_0 = Y_0,
      Y_1 = Y_1
    )

    return(dat)
  }


#create expanded dataset
create.pseudo <- function(dat, t_ind, Theta_vec) {
  n <- nrow(dat)
  dat.pseudo <- cbind(dat, 0)
  names(dat.pseudo) <- c(names(dat), "regime")
  dat.pseudo <- dat.pseudo[NULL, ]
  A <- dat$A

  P0 <- dat[, t_ind]

  for (ii in 1:(length(Theta_vec))) {
    #regime of interest: A = 1(P0<theta)
    temp.index <-
      (A == 1 & P0 < Theta_vec[ii]) | (A == 0 & P0 >= Theta_vec[ii])
    if (sum(temp.index) == 0)
      next
    dat.p <-
      cbind(dat[temp.index, ], rep(Theta_vec[ii], sum(temp.index)))
    names(dat.p) <- c(names(dat), "regime")
    dat.pseudo <- rbind(dat.pseudo, dat.p)
  }
  dat.pseudo <- dat.pseudo[order(dat.pseudo$id, dat.pseudo$regime), ]
  return(dat.pseudo)
}


#propensity model with no selection
no_selection <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  inds <- union(adjust_inds, t_ind)

  fit.treatment <- glm(dat$A ~ as.matrix(dat[, (2 + inds)]),
                       family = binomial())

  dat$pred.treatment <- fit.treatment$fitted

  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment > 1)
  )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)

  return(list(dat = dat, selected = rep(NA, var_num), bad_props = bad_props))

}


oracle <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  if (L$scenario == "A") {
    inds <- union(1:3, t_ind)

  } else{
    inds <- union(c(2, 3, 86), t_ind)

  }

  fit.treatment <- glm(dat$A ~ as.matrix(dat[, (2 + inds)]),
                       family = binomial())

  dat$pred.treatment <- fit.treatment$fitted

  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment > 1)
  )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)

  return(list(dat = dat, selected = rep(NA, var_num), bad_props = bad_props))

}

#weighted absolute mean difference function to choose lamdba weight in OAL
#arguments
#    betas - vector of beta values from outcome model
#    X - covariate vector
#    A - treatment vector
#    weights - matrix of fitted weights, a column per lamdba
#returns
# wAMD
wAMD <- function(betas, X, A, weights) {
  d <- length(betas)

  exposed_denom <- A %*% weights

  unexposed_denom <- (1 - A) %*% weights

  exposed_num <- matrix(NA, nrow = d,
                        ncol = ncol (weights))

  unexposed_num <- exposed_num

  for (j in 1:d) {
    exposed_num[j, ] <- (X[, j] * A) %*% weights

    unexposed_num[j, ] <- (X[, j] * (1 - A)) %*% weights

  }

  wAMDs <- abs(betas) %*%
    abs(
      sweep(exposed_num, 2, exposed_denom, "/")
      - sweep(unexposed_num, 2, unexposed_denom, "/")
    )

  return(wAMDs)

}


oal <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  inds <- c(setdiff(adjust_inds, t_ind), t_ind)

  X <- as.matrix(dat[, (2 + inds)])

  d <- length(inds)

  n <- nrow(dat)

  #lambda <- n^c(-10,-5,-1,-0.75,-0.5,-0.25,0.25,0.49)

  orig_lambda <- n ^ c(-20, -15, -10, -5, -1, -0.75, 0, 0.49)

  temp <- as.matrix(cbind(dat$A, X))

  outcome_betas <- lm(dat$Y ~ temp)$coefficients[-(1:2)]

  gamma <- 6 - 2 * log(orig_lambda, base = n)

  b_weights <- matrix(NA, nrow = (d),
                      ncol = length(gamma))

  for (j in 1:(d - 1)) {
    b_weights[j, ] <- abs(outcome_betas[j]) ^ (-gamma)
  }

  #tailoring ind has 0 weight
  b_weights[d, ] <- 0

  scale_factors <- apply(b_weights, 2, sum) / d

  lambda <- orig_lambda * scale_factors

  weights <- matrix(NA, nrow = nrow(dat), ncol = length(lambda))

  for (r in 1:length(lambda)) {
    fit.treatment <- glmnet(
      x = X,
      y = dat$A,
      family = binomial(link = "logit"),
      alpha = 1,
      lambda = lambda[r],
      penalty.factor = b_weights[, r]
    )

    fits <- predict(fit.treatment,
                    X,
                    s = lambda[r],
                    type = "response")

    weights[, r] <- 1 / fits * dat$A +
      (1 - dat$A) * 1 / (1 - fits)
  }


  wAMD_vals <- wAMD(
    betas = outcome_betas,
    X = X,
    A = dat$A,
    weights = weights
  )


  best_lambda_ind <- which.min(wAMD_vals)

  #print(log(orig_lambda[best_lambda_ind],n))

  final_mod <- glmnet(
    x = X,
    y = dat$A,
    family = binomial(link = "logit"),
    alpha = 1,
    lambda = lambda[best_lambda_ind],
    penalty.factor = b_weights[, best_lambda_ind]
  )

  fits <-   predict(final_mod,
                    X,
                    s = lambda[best_lambda_ind],
                    type = "response")

  coefs <-   coef(final_mod,
                  s = lambda[best_lambda_ind])[-1]


  #get selected variables in order of data in the DiPS function
  temp_selected <- as.vector(abs(coefs) > 10^(-5))


  selected <- rep(0, var_num)

  #get the true selected indices
  temp_selected <- inds[temp_selected]

  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1

  dat$pred.treatment <- fits

  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment > 1)
  )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)

  return(list(dat = dat, selected = selected, bad_props = bad_props))

}

glider <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  #put tailoring index last, glider will set last index weight to 0
  inds <- c(setdiff(adjust_inds, t_ind), t_ind)

  mod.treatment <- GLiDeR(
    Xorig = dat[, (2 + inds)],
    Yorig = dat$Y ,
    Aorig = dat$A,
    lambda = NULL
  )

  X <- as.matrix(cbind(1, dat[, (2 + inds)]))

  gamma <- c(mod.treatment$gamma0, mod.treatment$gamma)

  dat$pred.treatment <- 1 / (1 + exp(-(X %*% gamma)))

  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment > 1)
  )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)

  #get selected variables in order of data in the DiPS function
  temp_selected <- abs(mod.treatment$gamma) > 10^(-5)

  selected <- rep(0, var_num)

  #get the true selected indices
  temp_selected <- inds[temp_selected]

  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1

  return(list(dat = dat, selected = selected, bad_props = bad_props))

}

causal_ball <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  #set tailoring index as last
  inds <- c(setdiff(adjust_inds, t_ind), t_ind)

  #tell CBS that tailoring index is last

  temp <- CBS(
    X = dat[, (2 + inds)],
    D = dat$A,
    Y = dat$Y,
    t_ind = length(inds)
  )


  dat$pred.treatment <- temp$fitted.ps

  bad_props <- c(mean(dat$pred.treatment < 0),
                 mean(dat$pred.treatment == 0),
                 mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                 mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                 mean(dat$pred.treatment == 1),
                 mean(dat$pred.treatment > 1)
  )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

  temp_selected <- temp$selected

  selected <- rep(0, var_num)

  #get the true selected indices
  temp_selected <- inds[temp_selected == 1]

  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1

  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)

  return(list(dat = dat, selected = selected, bad_props = bad_props))

}

dip <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  inds <- c(setdiff(adjust_inds, t_ind), t_ind)

  temp <-
    DiPS(
      Yi = dat$Y,
      Ti = dat$A,
      Xi = as.matrix(dat[, (2 + inds)]),
      Wi = rep(1, L$n),
      bwd = "plug",
      q = 4,
      alp = 1 / (2 + 4),
      C = 1,
      constEff = TRUE,
      t_ind = length(inds)
    )

  bad_props <- c(mean(temp[[1]] < 0),
                 mean(temp[[1]] == 0),
                 mean(temp[[1]] > 0 & temp[[1]] < 10 ^ (-6)),
                 mean(temp[[1]] < 1 & temp[[1]] > (1 - 10 ^ (-6))),
                 mean(temp[[1]] == 1),
                 mean(temp[[1]] > 1)
  )

  #put weights in [0,1] range
  probs <- pmin(pmax(temp[[1]], 10 ^ (-6)), 1 - 10 ^ (-6))

  dat$weights.trt <- dat$A * 1 / probs +
    (1 - dat$A) * 1 / probs

  dat$pred.treatment <- probs

  #get selected variables in order of data in the DiPS function
  temp_selected <- temp[[2]]

  selected <- rep(0, var_num)

  #get the true selected indices
  temp_selected <- inds[temp_selected]

  #get a vector of all indices with whether they were selected
  selected[temp_selected] <- 1

  return(list(dat = dat, selected = selected, bad_props = bad_props))

}

hd_balancing <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  X <- as.matrix(dat[, (2 + adjust_inds)])

  A <- dat$A

  dat$pred.treatment <- NA


    dat$pred.treatment <- hdCBPS(
      formula = A ~ X,
      y = dat$Y,
      ATT = 0,
      #default interations is 1000
      iterations = 100
    )$fitted.values

    bad_props <- c(mean(dat$pred.treatment < 0),
                   mean(dat$pred.treatment == 0),
                   mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                   mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                   mean(dat$pred.treatment == 1),
                   mean(dat$pred.treatment > 1)
    )

  dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

  dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
    (1 - dat$A) * 1 / (1 - dat$pred.treatment)

  return(list(dat = dat, selected = rep(NA, var_num), bad_props = bad_props))


}


#treatment weights from logistic regression with baseline variables
#arguments:
#  dat - data frame with tx var A at least
#  adjust_inds - vector of columns of any confounders in dat
create.analysis <- function(dat, adjust_inds, t_ind) {
  var_num <- ncol(dat) - 5

  select_fun <- eval(as.name(L$prop.selection))

  #baseline treatment weights
  if (is.null(adjust_inds)) {
    fit.treatment <-
      glm(dat$A ~ as.matrix(dat[, (2 + t_ind)]), family = binomial())

    dat$pred.treatment <- fit.treatment$fitted

    bad_props <- c(mean(dat$pred.treatment < 0),
                   mean(dat$pred.treatment == 0),
                   mean(dat$pred.treatment > 0 & dat$pred.treatment < 10 ^ (-6)),
                   mean(dat$pred.treatment < 1 & dat$pred.treatment > (1 - 10 ^ (-6))),
                   mean(dat$pred.treatment == 1),
                   mean(dat$pred.treatment > 1)
      )

    dat$pred.treatment <- pmin(pmax(dat$pred.treatment, 10 ^ (-6)), 1 - 10 ^ (-6))

    dat$weights.trt <- dat$A * 1 / dat$pred.treatment +
      (1 - dat$A) * 1 / (1 - dat$pred.treatment)

    dat$weights.trt  <- pmin(dat$weights.trt, 10)

    return(list(dat = dat, selected  = rep(0, var_num), bad_props = bad_props))
  }
  else {

    withTimeout({
      ret <- select_fun(dat, adjust_inds, t_ind)},
      timeout = 7200,
      onTimeout = "error")

    ret$dat$weights.trt <- pmin(ret$dat$weights.trt, 10)

    return(ret)
  }

}

#arguments:
# dat: dataframe of observations with IPW weights
# t_ind - the index of the tailoring variable under consideration
# prop_adjust - logical, whether to adjust for covariates in the propensity
# prop_adjust_inds - inds of vars to adjust for in propensity
# outcome_adjust - how to adjust for covariates in the marginal:
#    via "Lasso", via the "Oracle", or "None"
# outcome_adjust_inds - inds of vars to adjust for besides regime in marginal
# Theta_vec - vector of regime values under consideration
# knots - vector of spline knots
#returns:
# v - estimated value of regime correspondng to each entry of Theta_vec
# beta - coefficients of the model
estimate.regime <- function(dat,
                            t_ind,
                            prop_adjust,
                            prop_adjust_inds,
                            outcome_adjust,
                            outcome_adjust_inds,
                            Theta_vec,
                            knots)
{
  #change from v6, IPW weights are specific to variable
  temp <- create.analysis(dat, prop_adjust_inds, t_ind)

  #make sure tailoring variable isn't in margianl model
  outcome_adjust_inds <- setdiff(outcome_adjust_inds,t_ind)

  selected <- temp$selected

  bad_props <- temp$bad_props

  dat <- temp$dat

  dat.pseudo <- create.pseudo(dat, (t_ind + 2), Theta_vec)

  X_marg <- ns(x = dat.pseudo$regime, knots = knots)


  #add in adjustment variables
  if (outcome_adjust == "Lasso") {
    X <- cbind(X_marg, dat.pseudo[, (outcome_adjust_inds + 2)])

    X <- as.matrix(X)

    penalty.factor <- c(rep(0, ncol(X_marg)),
                        rep(1, ncol(dat.pseudo[, (outcome_adjust_inds +
                                                    2)])))

    mod <- glmnet(
      y = dat.pseudo$Y,
      x = X,
      alpha = 1,
      penalty.factor = penalty.factor,
      weights = dat.pseudo$weights.trt
    )


    pred <- predict(mod, newx = X)

    bic <- rep(NA, ncol(pred))

    for (j in 1:ncol(pred)) {
      #gets the weighted mse and number of knots, calculates criteria val
      wmse <- mean(dat.pseudo$weights.trt * (dat.pseudo$Y - pred[, j]) ^
                     2)

      p <- mod$df[j] - 1


      bic[j] <- log(wmse) * nrow(X) + log(nrow(X)) * p
    }

    best_s <- mod$lambda[which.min(bic)]

    v_selected <- as.vector(coef(mod, s = best_s))[-1] != 0

    res <-
      lm(dat.pseudo$Y ~ X[, v_selected], weights = dat.pseudo$weights.trt)

    beta <- as.vector(coef(mod, s = best_s))

    beta[c(TRUE, v_selected)] <- res$coefficients

    v <- rep(NA, length(Theta_vec))

    X_pred_all <- ns(x = Theta_vec, knots = knots)

    for (i in 1:length(Theta_vec)) {
      X_pred <- X_pred_all[i, ]

      X_pred <-
        data.frame(1, t(X_pred), dat[, (outcome_adjust_inds + 2)])

      v[i] <- mean(as.matrix(X_pred) %*% beta)
    }

  } else if (outcome_adjust == "Oracle") {



    X <- as.matrix(cbind(X_marg, dat.pseudo[, (outcome_adjust_inds + 2)]))

    res <- lm(dat.pseudo$Y ~ X, weights = dat.pseudo$weights.trt)

    beta <- res$coefficients

    v <- rep(NA, length(Theta_vec))

    X_pred_all <- ns(x = Theta_vec, knots = knots)

    for (i in 1:length(Theta_vec)) {
      X_pred <- X_pred_all[i, ]

      X_pred <-
        data.frame(1, t(X_pred), dat[, (outcome_adjust_inds + 2)])

      v[i] <- mean(as.matrix(X_pred) %*% beta)
    }

  } else {
    res <- lm(dat.pseudo$Y ~ X_marg, weights = dat.pseudo$weights.trt)

    beta <- res$coefficients

    v <- rep(NA, length(Theta_vec))

    X_pred <- ns(x = Theta_vec, knots = knots)

    X_pred <- as.matrix(cbind(1, X_pred))

    v <- X_pred %*% beta


  }


  return(
    list(
      beta = beta,
      v = v,
      props = dat$pred.treatment,
      weights = dat$weights.trt,
      selected = selected,
      bad_props = bad_props
    )
  )

}


#helper function to run all lasso types within a given rep
get_all_estimates <- function(dat,
                              Theta_vec,
                              knots_vec,
                              tailor_search_inds,
                              prop_adjust,
                              prop_adjust_inds,
                              outcome_adjust,
                              outcome_adjust_inds) {
  t_num <- length(tailor_search_inds)

  value <- rep(NA, t_num)

  coefs <-
    matrix(NA, length(c(knots_vec, outcome_adjust_inds)) + 3, t_num)

  props <- matrix(NA, nrow = nrow(dat), ncol = t_num)

  weights <- matrix(NA, nrow = nrow(dat), ncol = t_num)

  selected <- matrix(NA, nrow = (ncol(dat) - 5), ncol = t_num)

  bad_props <- matrix(NA, nrow = t_num, ncol = 6)

  for (r in 1:t_num) {
    t_ind <- tailor_search_inds[r]
    #get the msm model
    res <- estimate.regime(
      dat,
      t_ind = t_ind,
      prop_adjust = prop_adjust,
      prop_adjust_inds = prop_adjust_inds,
      outcome_adjust = outcome_adjust,
      outcome_adjust_inds = outcome_adjust_inds,
      Theta_vec = Theta_vec,
      knots = knots_vec
    )



    #record the maximum value function
    value[r] <- max(res$v)

    coefs[1:length(res$beta), r] <- res$beta

    props[, r] <- res$props

    weights[, r] <- res$weights

    selected[, r] <- res$selected

    bad_props[r,] <- res$bad_props

  }

  return(
    list(
      value = value,
      coefs = coefs,
      props = props,
      weights = weights,
      selected = selected,
      bad_props = bad_props
    )
  )

}
