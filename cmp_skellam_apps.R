# Code to fit MPCMP-Skellam models after successful simulations were
# conducted using cmp_skellam_sims.R 
k_range <- -10:10
sum_index2_l <- pmax(0, -k_range)
sum_index2_u <- sum_index2_l + 10
sum_index1_l <- sum_index2_l + k_range
sum_index1_u <- sum_index2_u + k_range
m1 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index1_l), length(k_range), 11, byrow = T)
m2 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index2_l), length(k_range), 11, byrow = T)

#### Below code uses block update for team strength parameters ####
fit_mpcmp_skel_apps <- function(data, niters = 1000){
  # Extract data
  dat <- data
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, -10:10)
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise
  beta01 <- log(mean(dat$hgoal))
  beta02 <- log(mean(dat$vgoal))
  beta_p <- c(beta01, beta02) # Prior mean for beta
  beta0 <- c(beta01, beta02)
  # Sum-to-zero constraint on attack and defence strengths
  team_att_str_h <- log(with(dat, tapply(hgoal, homeID, mean))) - beta01
  team_att_str_a <- log(with(dat, tapply(vgoal,awayID, mean))) - beta02
  if (sum(is.infinite(team_att_str_h) > 0)){
    iiInf <- is.infinite(team_att_str_h)
    team_att_str_h[iiInf] <- mean(team_att_str_h[!iiInf])
  }
  if (sum(is.infinite(team_att_str_a) > 0)){
    iiInf <- is.infinite(team_att_str_a)
    team_att_str_a[iiInf] <- mean(team_att_str_a[!iiInf])
  }
  team_att_str <- 0.5*(team_att_str_h + team_att_str_a)
  alpha0 <- team_att_str - mean(team_att_str)
  print(alpha0)
  team_def_str_h <- -log(with(dat, tapply(hgoal, awayID, mean))) + beta01
  team_def_str_a <- -log(with(dat, tapply(vgoal, homeID, mean))) + beta02
  if (sum(is.infinite(team_def_str_h) > 0)){
    iiInf <- is.infinite(team_def_str_h)
    team_def_str_h[iiInf] <- mean(team_def_str_h[!iiInf])
  }
  if (sum(is.infinite(team_att_str_a) > 0)){
    iiInf <- is.infinite(team_def_str_a)
    team_def_str_a[iiInf] <- mean(team_def_str_a[!iiInf])
  }
  team_def_str <- 0.5*(team_def_str_h + team_def_str_a)
  delta0 <- team_def_str - mean(team_def_str)
  print(delta0)
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams) # Initialise at equidispersion for each team
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta block update
    prop_covar <- 0.0001*matrix(c(0.15, 0.7*0.15, 0.7*0.15, 0.15), 2, 2)
    #beta1 <- rmvnorm(1, beta0, prop_covar)[1, ]
    beta1 <- rnorm(2, beta0, 0.005)
    lmu11 <- with(dat, beta1[1] + alpha0[homeID] - delta0[awayID])
    lmu21 <- with(dat, beta1[2] + alpha0[awayID] - delta0[homeID])
    if (max(lmu11) >= log(100)){lmu11 <- pmin(log(100), lmu11)}
    if (min(lmu11) < log(0.01)){lmu11 <- pmax(log(0.01), lmu11)}
    if (max(lmu21) >= log(100)){lmu21 <- pmin(log(100), lmu21)}
    if (min(lmu21) < log(0.01)){lmu21 <- pmax(log(0.01), lmu21)}
    lmu11_index <- floor((lmu11 - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((lmu21 - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((beta0 - beta_p)^2 - 
                              (beta1 - beta_p)^2)/(2*0.01^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      mu10 <- mu11
      mu20 <- mu21
      current <- prop
    }
    # nu update - on log scale (xi)
    xi1 <- rnorm(nteams, xi0, 0.35)
    if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
    if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
    nu11_lp <- exp(xi1[dat$homeID])
    nu21_lp <- exp(xi1[dat$awayID])
    nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
    nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
    prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (xi0^2 - xi1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    xi0[accept] <- xi1[accept]
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.003)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    lmu11 <- with(dat, beta0[1] + alpha1[homeID] - delta0[awayID])
    lmu21 <- with(dat, beta0[2] + alpha1[awayID] - delta0[homeID])
    if (max(lmu11) >= log(100)){lmu11 <- pmin(log(100), lmu11)}
    if (min(lmu11) < log(0.01)){lmu11 <- pmax(log(0.01), lmu11)}
    if (max(lmu21) >= log(100)){lmu21 <- pmin(log(100), lmu21)}
    if (min(lmu21) < log(0.01)){lmu21 <- pmax(log(0.01), lmu21)}
    lmu11_index <- floor((lmu11 - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((lmu21 - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum((alpha0^2 - alpha1^2)/(2*0.25^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      alpha0 <- alpha1
      current <- prop
    }
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.003)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    lmu11 <- with(dat, beta0[1] + alpha0[homeID] - delta1[awayID])
    lmu21 <- with(dat, beta0[2] + alpha0[awayID] - delta1[homeID])
    if (max(lmu11) >= log(100)){lmu11 <- pmin(log(100), lmu11)}
    if (min(lmu11) < log(0.01)){lmu11 <- pmax(log(0.01), lmu11)}
    if (max(lmu21) >= log(100)){lmu21 <- pmin(log(100), lmu21)}
    if (min(lmu21) < log(0.01)){lmu21 <- pmax(log(0.01), lmu21)}
    lmu11_index <- floor((lmu11 - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((lmu21 - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum((delta0^2 - delta1^2)/(2*0.25^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      delta0 <- delta1
      current <- prop
    }
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_apps, data = dat)
}

#### Try on Italy data ####
library(footBayes)
library(dplyr)
data(italy)
italy_2016_2022 <- italy %>%
  dplyr::select(Season, home, visitor, hgoal, vgoal, tier) %>%
  dplyr::filter(Season%in%c("2016", "2017", "2018", "2019", "2020", "2021", "2022") & tier == 1)
italy_2016_2022$goaldif <- italy_2016_2022$hgoal - italy_2016_2022$vgoal
# Need team IDs
italy_2016_2022$homeID <- as.numeric(as.factor(italy_2016_2022$home))
italy_2016_2022$awayID <- as.numeric(as.factor(italy_2016_2022$visitor))

#### Component-wise updates for team strengths ####
fit_mpcmp_skel_apps_2 <- function(data, niters = 1000){
  # Extract data
  dat <- data
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, -10:10)
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise
  beta01 <- log(mean(dat$hgoal))
  beta02 <- log(mean(dat$vgoal))
  beta_p <- c(beta01, beta02) # Prior mean for beta
  beta0 <- c(beta01, beta02)
  # Sum-to-zero constraint on attack and defence strengths
  team_att_str_h <- log(with(dat, tapply(hgoal, homeID, mean)) + 0.001) - beta01
  team_att_str_a <- log(with(dat, tapply(vgoal, awayID, mean)) + 0.001) - beta02
  team_att_str <- 0.5*(team_att_str_h + team_att_str_a)
  alpha0 <- team_att_str - mean(team_att_str)
  print(alpha0)
  team_def_str_h <- -log(with(dat, tapply(hgoal, awayID, mean)) + 0.001) + beta01
  team_def_str_a <- -log(with(dat, tapply(vgoal, homeID, mean)) + 0.001) + beta02
  team_def_str <- 0.5*(team_def_str_h + team_def_str_a)
  delta0 <- team_def_str - mean(team_def_str)
  print(delta0)
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams) # Initialise at equidispersion for each team
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta block update
    prop_covar <- 0.0001*matrix(c(0.15, 0.7*0.15, 0.7*0.15, 0.15), 2, 2)
    #beta1 <- rmvnorm(1, beta0, prop_covar)[1, ]
    beta1 <- rnorm(2, beta0, 0.005)
    lmu11 <- with(dat, beta1[1] + alpha0[homeID] - delta0[awayID])
    lmu21 <- with(dat, beta1[2] + alpha0[awayID] - delta0[homeID])
    if (max(lmu11) >= log(100)){lmu11 <- pmin(log(100), lmu11)}
    if (min(lmu11) < log(0.01)){lmu11 <- pmax(log(0.01), lmu11)}
    if (max(lmu21) >= log(100)){lmu21 <- pmin(log(100), lmu21)}
    if (min(lmu21) < log(0.01)){lmu21 <- pmax(log(0.01), lmu21)}
    lmu11_index <- floor((lmu11 - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((lmu21 - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((beta0 - beta_p)^2 - 
                              (beta1 - beta_p)^2)/(2*0.01^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      mu10 <- mu11
      mu20 <- mu21
      current <- prop
    }
    # nu update - on log scale (xi)
    xi1 <- rnorm(nteams, xi0, 0.35)
    if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
    if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
    nu11_lp <- exp(xi1[dat$homeID])
    nu21_lp <- exp(xi1[dat$awayID])
    nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
    nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
    prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (xi0^2 - xi1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    xi0[accept] <- xi1[accept]
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.003)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    lmu11 <- with(dat, beta0[1] + alpha1[homeID] - delta0[awayID])
    lmu21 <- with(dat, beta0[2] + alpha1[awayID] - delta0[homeID])
    if (max(lmu11) >= log(100)){lmu11 <- pmin(log(100), lmu11)}
    if (min(lmu11) < log(0.01)){lmu11 <- pmax(log(0.01), lmu11)}
    if (max(lmu21) >= log(100)){lmu21 <- pmin(log(100), lmu21)}
    if (min(lmu21) < log(0.01)){lmu21 <- pmax(log(0.01), lmu21)}
    lmu11_index <- floor((lmu11 - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((lmu21 - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum((alpha0^2 - alpha1^2)/(2*0.25^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      alpha0 <- alpha1
      current <- prop
    }
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.003)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    lmu11 <- with(dat, beta0[1] + alpha0[homeID] - delta1[awayID])
    lmu21 <- with(dat, beta0[2] + alpha0[awayID] - delta1[homeID])
    if (max(lmu11) >= log(100)){lmu11 <- pmin(log(100), lmu11)}
    if (min(lmu11) < log(0.01)){lmu11 <- pmax(log(0.01), lmu11)}
    if (max(lmu21) >= log(100)){lmu21 <- pmin(log(100), lmu21)}
    if (min(lmu21) < log(0.01)){lmu21 <- pmax(log(0.01), lmu21)}
    lmu11_index <- floor((lmu11 - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((lmu21 - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum((delta0^2 - delta1^2)/(2*0.25^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      delta0 <- delta1
      current <- prop
    }
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_apps, data = dat)
}

#### Code with larger grid for mu - component-wise updates ####
#### MLB fit ####
k_range <- min(mlb_dat2$goaldif):max(mlb_dat2$goaldif)
sum_index2_l <- pmax(0, -k_range)
sum_index2_u <- sum_index2_l + 10
sum_index1_l <- sum_index2_l + k_range
sum_index1_u <- sum_index2_u + k_range
m1 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index1_l), length(k_range), 11, byrow = T)
m2 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index2_l), length(k_range), 11, byrow = T)

fit_mpcmp_skel_MLB <- function(dat, niters = 1000, skel = FALSE){
  # Extract data
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  # An if statement here?
  mdiff <- match(diffs, min(dat$goaldif):max(dat$goaldif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals(dat)
  print(starts)
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  # Priors
  beta_m <- beta0
  beta_s <- 0.02
  alpha_m <- alpha0 #rep(0, nteams) #alpha0
  alpha_s <- 0.15
  delta_m <- delta0 #rep(0, nteams) #delta0 
  delta_s <- 0.15
  xi_s <- 0.50
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (one-at-a-time)
    beta1 <- rnorm(2, beta0, 0.05)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    prop <- log(fun4A(lam11, lam20, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[1] - beta_m[1])^2 - 
                          (beta1[1] - beta_m[1])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[1] <- beta1[1]
      lam10 <- lam11
      lmu10_index <- lmu11_index
      current <- prop
    }
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam10, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[2] - beta_m[2])^2 - 
                          (beta1[2] - beta_m[2])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[2] <- beta1[2]
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.40)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
      denrat <- with(dat, tapply(prop - current, homeID, sum) +
                       tapply(prop - current, awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.04)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.04)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_MLB, data = dat)
}


#### MLB fit with inflated one run home victories ####
k_range <- min(mlb2016_dat$goaldif):max(mlb2016_dat$goaldif)
sum_index2_l <- pmax(0, -k_range)
sum_index2_u <- sum_index2_l + 10
sum_index1_l <- sum_index2_l + k_range
sum_index1_u <- sum_index2_u + k_range
m1 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index1_l), length(k_range), 11, byrow = T)
m2 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index2_l), length(k_range), 11, byrow = T)

fit_mpcmp_skel_MLB_inf <- function(dat, niters = 1000, skel = FALSE){
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, min(dat$goaldif):max(dat$goaldif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Indicator for one run wins
  ii1R <- dat$hgoal == dat$vgoal + 1
  print(table(ii1R))
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  pis <- matrix(NA, nrow = K, ncol = 1) # Inflation parameter
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals(dat)
  print(starts)
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  pi0 <- 0.05
  # Priors
  beta_m <- beta0
  beta_s <- 0.01
  alpha_m <- alpha0 #rep(0, nteams) #alpha0
  alpha_s <- 0.05
  delta_m <- delta0 #rep(0, nteams) #delta0 
  delta_s <- 0.05
  xi_s <- 0.5
  # Normal prior for logit(pi)
  pil_m <- -3 
  pil_s <- 0.25
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
  current[ii1R] <- current[ii1R] + pi0
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (one-at-a-time)
    beta1 <- rnorm(2, beta0, 0.05)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    prop <- (1 - pi0)*fun4A(lam11, lam20, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + ((beta0[1] - beta_m[1])^2 - 
                          (beta1[1] - beta_m[1])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[1] <- beta1[1]
      lam10 <- lam11
      lmu10_index <- lmu11_index
      current <- prop
    }
    # Now home advantage term
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - pi0)*fun4A(lam10, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + ((beta0[2] - beta_m[2])^2 -
                          (beta1[2] - beta_m[2])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[2] <- beta1[2]
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.40)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- (1 - pi0)*fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs
      prop[ii1R] <- prop[ii1R] + pi0
      denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                       tapply(log(prop) - log(current), awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.02)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - pi0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                     tapply(log(prop) - log(current), awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.02)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - pi0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                     tapply(log(prop) - log(current), awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    # Inflation parameter update (on logit scale?)
    pil0 <- log(pi0/(1 - pi0)) 
    pil1 <- rnorm(1, pil0, 0.15)
    pi1 <- exp(pil1)/(1 + exp(pil1))
    prop <- (1 - pi1)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi1
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + sum(((pil0 - pil_m)^2 - 
                              (pil1 - pil_m)^2)/(2*pil_s^2))
    laccept <- min(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){pil0 <- pil1}
    pi0 <- exp(pil0)/(1 + exp(pil0))
    print(pi0)
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    lhood <- sum(log(current))
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    pis[k,] <- pi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis, pis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       pis = pis, ll = ll, cpu_time = cpu_time[3], 
       acc_rates = acc_rates, model = fit_mpcmp_skel_MLB_inf, 
       data = dat)
}

#### MLB fit with +1 inflation but no home advantage ####
fit_mpcmp_skel_MLB_inf2 <- function(dat, niters = 1000, skel = FALSE){
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, min(dat$goaldif):max(dat$goaldif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Indicator for one run wins
  ii1R <- dat$hgoal == dat$vgoal + 1
  print(table(ii1R))
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 1) # Overall mean
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  pis <- matrix(NA, nrow = K, ncol = 1) # Inflation parameter
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals(dat)
  print(starts)
  beta0 <- mean(starts$beta0)
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  pi0 <- 0.05
  # Priors
  beta_m <- beta0
  beta_s <- 0.01
  alpha_m <- alpha0 #rep(0, nteams) #alpha0
  alpha_s <- 0.05
  delta_m <- delta0 #rep(0, nteams) #delta0 
  delta_s <- 0.05
  xi_s <- 0.5
  # Normal prior for logit(pi)
  pil_m <- -3 
  pil_s <- 0.25
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
  current[ii1R] <- current[ii1R] + pi0
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (now represents overall mean for both H/A)
    beta1 <- rnorm(1, beta0, 0.01)
    mu11 <- with(dat, exp(beta1 + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1 + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - pi0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + ((beta0 - beta_m)^2 - 
                          (beta1 - beta_m)^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lam10 <- lam11
      lmu10_index <- lmu11_index
      lam20 <- lam21
      lmu20_index <- lmu21_index
      current <- prop
    }
    print(beta0)
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.40)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- (1 - pi0)*fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs
      prop[ii1R] <- prop[ii1R] + pi0
      denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                       tapply(log(prop) - log(current), awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.02)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0 + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0 + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - pi0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                     tapply(log(prop) - log(current), awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.02)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0 + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0 + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - pi0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0
    denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                     tapply(log(prop) - log(current), awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    # Inflation parameter update
    pil0 <- log(pi0/(1 - pi0)) 
    pil1 <- rnorm(1, pil0, 0.15)
    pi1 <- exp(pil1)/(1 + exp(pil1))
    prop <- (1 - pi1)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi1
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + sum(((pil0 - pil_m)^2 - 
                              (pil1 - pil_m)^2)/(2*pil_s^2))
    laccept <- min(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){pil0 <- pil1}
    pi0 <- exp(pil0)/(1 + exp(pil0))
    print(pi0)
    current <- (1 - pi0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0
    
    lhood <- sum(log(current))
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    pis[k,] <- pi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis, pis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       pis = pis, ll = ll, cpu_time = cpu_time[3], 
       acc_rates = acc_rates, model = fit_mpcmp_skel_MLB_inf, 
       data = dat)
}

#### MLB fit with +1, -1 inflation but no home advantage ####
fit_mpcmp_skel_MLB_inf3 <- function(dat, niters = 1000, skel = FALSE){
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, min(dat$goaldif):max(dat$goaldif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Indicator for one run home wins
  ii1R <- dat$hgoal == dat$vgoal + 1
  # Indicator for one run away wins
  ii1mR <- dat$hgoal == dat$vgoal - 1
  print(table(ii1R))
  print(table(ii1mR))
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 1) # Overall mean
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  pis <- matrix(NA, nrow = K, ncol = 2) # Inflation parameters
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals(dat)
  print(starts)
  beta0 <- mean(starts$beta0)
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  pi0 <- rep(0.05, 2)
  # Priors
  beta_m <- beta0
  beta_s <- 0.01
  alpha_m <- alpha0 #rep(0, nteams) #alpha0
  alpha_s <- 0.05
  delta_m <- delta0 #rep(0, nteams) #delta0 
  delta_s <- 0.05
  xi_s <- 0.5
  # Normal prior for logit(pi)
  pil_m <- -3 
  pil_s <- 0.25
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
  current[ii1R] <- current[ii1R] + pi0[1]
  current[ii1mR] <- current[ii1mR] + pi0[2]
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (now represents overall mean for both H/A)
    beta1 <- rnorm(1, beta0, 0.01)
    mu11 <- with(dat, exp(beta1 + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1 + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0[1]
    prop[ii1mR] <- prop[ii1mR] + pi0[2]
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + ((beta0 - beta_m)^2 - 
                          (beta1 - beta_m)^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lam10 <- lam11
      lmu10_index <- lmu11_index
      lam20 <- lam21
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.40)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs
      prop[ii1R] <- prop[ii1R] + pi0[1]
      prop[ii1mR] <- prop[ii1mR] + pi0[2]
      denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                       tapply(log(prop) - log(current), awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0[1]
    current[ii1mR] <- current[ii1mR] + pi0[2]
    
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.02)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0 + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0 + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0[1]
    prop[ii1mR] <- prop[ii1mR] + pi0[2]
    denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                     tapply(log(prop) - log(current), awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0[1]
    current[ii1mR] <- current[ii1mR] + pi0[2]
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.02)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0 + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0 + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0[1]
    prop[ii1mR] <- prop[ii1mR] + pi0[2]
    denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                     tapply(log(prop) - log(current), awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0[1]
    current[ii1mR] <- current[ii1mR] + pi0[2]
    
    # Inflation parameter update
    pil0 <- log(pi0/(1 - pi0)) 
    pil1 <- rnorm(2, pil0, 0.15)
    pi1 <- exp(pil1)/(1 + exp(pil1))
    prop <- (1 - sum(pi1))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi1[1]
    prop[ii1mR] <- prop[ii1mR] + pi1[2]
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + sum(((pil0 - pil_m)^2 - 
                              (pil1 - pil_m)^2)/(2*pil_s^2))
    laccept <- min(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){pil0 <- pil1}
    pi0 <- exp(pil0)/(1 + exp(pil0))
    print(pi0)
    current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0[1]
    current[ii1mR] <- current[ii1mR] + pi0[2]
    
    lhood <- sum(log(current))
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    pis[k,] <- pi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis, pis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       pis = pis, ll = ll, cpu_time = cpu_time[3], 
       acc_rates = acc_rates, model = fit_mpcmp_skel_MLB_inf, 
       data = dat)
}


#### Starting values ####
# For attack/defence model
start_vals <- function(data){
  pmp <- log(0.5*(with(data, tapply(hgoal, homeID, mean)) + 
                    with(data, tapply(vgoal, awayID, mean))))
  att_inits <- pmp - mean(pmp)
  omp <- -log(0.5*(with(data, tapply(hgoal, awayID, mean)) + 
                     with(data, tapply(vgoal, homeID, mean))))
  def_inits <- omp - mean(omp)
  hgoal_hat <- exp(att_inits[data$homeID] - def_inits[data$awayID])
  beta0hat <- log(mean(data$hgoal/hgoal_hat))
  vgoal_hat <- exp(att_inits[data$awayID] -
                     def_inits[data$homeID])
  beta1hat <- log(mean(data$vgoal/vgoal_hat))
  list(beta0 = c(beta0hat, beta1hat), alpha0 = att_inits,
       delta0 = def_inits)
}

# For corner constrained model
start_vals_cor <- function(data){
  pmp <- log(0.5*(with(data, tapply(hgoal, homeID, mean)) + 
                    with(data, tapply(vgoal, awayID, mean))))
  att_inits <- pmp - pmp[1]
  omp <- -log(0.5*(with(data, tapply(hgoal, awayID, mean)) + 
                     with(data, tapply(vgoal, homeID, mean))))
  def_inits <- omp - omp[1]
  hgoal_hat <- exp(att_inits[data$homeID] - def_inits[data$awayID])
  beta0hat <- log(mean(data$hgoal/hgoal_hat))
  vgoal_hat <- exp(att_inits[data$awayID] -
                     def_inits[data$homeID])
  beta1hat <- log(mean(data$vgoal/vgoal_hat))
  list(beta0 = c(beta0hat, beta1hat), alpha0 = att_inits,
       delta0 = def_inits)
}


# For overall strength model with yellow cards
# Want league-specific betas
start_vals2 <- function(data){
  omp <- log(0.5*(with(data, tapply(HY, homeID, mean)) + 
                    with(data, tapply(AY, awayID, mean))))
  PMP <- table(YC_skel$leagueID, YC_skel$homeID)
  nleague <- rowSums(PMP != 0)
  league_id <- rep(1, length(omp))
  league_id[(1:length(omp))[PMP[2,] != 0]] <- 2
  league_id[(1:length(omp))[PMP[3,] != 0]] <- 3
  league_id[(1:length(omp))[PMP[4,] != 0]] <- 4
  league_id[(1:length(omp))[PMP[5,] != 0]] <- 5
  league_m <- tapply(omp, league_id, mean)
  team_inits <- omp - league_m[league_id]
  hcards_hat <- exp(team_inits[data$homeID])
  ests_h <- log(with(data, tapply(HY/hcards_hat, leagueID, mean)))
  beta0hat <- ests_h[1]
  psihat_h <- ests_h[-1] - beta0hat
  acards_hat <- exp(team_inits[data$awayID])
  ests_a <- log(with(data, tapply(AY/acards_hat, leagueID, mean)))
  beta1hat <- ests_a[1]
  psihat_a <- ests_a[-1] - beta1hat
  league_inits <- 0.5*(psihat_h + psihat_a)
  list(beta0 = c(beta0hat, beta1hat), theta0 = team_inits,
       psi0 = league_inits)
}

#### Serie A fit ####
# Same as above but can tweak for acceptance rates, etc.
# Data first
library(footBayes)
library(dplyr)
italy_2016_2022 <- italy %>%
  dplyr::select(Season, home, visitor, hgoal, vgoal, tier) %>%
  dplyr::filter(Season%in%c("2016", "2017", "2018", "2019", "2020", "2021", "2022") & tier == 1)
italy_2016_2022$goaldif <- italy_2016_2022$hgoal - italy_2016_2022$vgoal
# Need team IDs
italy_2016_2022$homeID <- as.numeric(as.factor(italy_2016_2022$home))
italy_2016_2022$awayID <- as.numeric(as.factor(italy_2016_2022$visitor))

k_range <- (min(italy_2016_2022$goaldif) -3):(max(italy_2016_2022$goaldif) + 3)
sum_index2_l <- pmax(0, -k_range)
sum_index2_u <- sum_index2_l + 10
sum_index1_l <- sum_index2_l + k_range
sum_index1_u <- sum_index2_u + k_range
m1 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index1_l), length(k_range), 11, byrow = T)
m2 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index2_l), length(k_range), 11, byrow = T)

fit_mpcmp_skel_SA <- function(dat, niters = 1000, skel = FALSE){
  # Extract data
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  # An if statement here?
  mdiff <- match(diffs, (min(dat$goaldif) - 3):(max(dat$goaldif) + 3))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals(dat)
  print(starts)
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  # Priors
  beta_m <- beta0
  beta_s <- 0.02
  alpha_m <- alpha0
  alpha_s <- 0.15
  delta_m <- delta0
  delta_s <- 0.15
  xi_s <- 1.0
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (one-at-a-time)
    beta1 <- rnorm(2, beta0, 0.05)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    prop <- log(fun4A(lam11, lam20, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[1] - beta_m[1])^2 - 
                          (beta1[1] - beta_m[1])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[1] <- beta1[1]
      lam10 <- lam11
      lmu10_index <- lmu11_index
      current <- prop
    }
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam10, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[2] - beta_m[2])^2 - 
                          (beta1[2] - beta_m[2])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[2] <- beta1[2]
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.4)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
      denrat <- with(dat, tapply(prop - current, homeID, sum) +
                       tapply(prop - current, awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.03)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.03)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_SA, data = dat)
}

#### Draw-inflated Serie A fit ####
fit_zi_mpcmp_skel_SA <- function(dat, niters = 1000, skel = FALSE){
  # Extract data
  # Objects needed for fitting
  n <- nrow(dat)
  # Indicator for tied outcomes
  zzI <- dat$hgoal == dat$vgoal
  print(table(zzI))
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, min(dat$goaldif):max(dat$goaldif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  pis <- matrix(NA, nrow = K, ncol = 1) # (0, 0) inflation parameter
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals(dat)
  print(starts)
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  pi_z0 <- 0.02
  # Priors
  beta_m <- beta0
  beta_s <- 0.02
  alpha_m <- alpha0
  alpha_s <- 0.05
  delta_m <- delta0
  delta_s <- 0.05
  xi_s <- 0.5
  # Beta prior for pi_z0
  pi_a <- 0.1
  pi_b <- 1
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log((1 - pi_z0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  current[zzI] <- current[zzI] + log(pi_z0)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (one-at-a-time)
    beta1 <- rnorm(2, beta0, 0.05)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    prop <- log((1 - pi_z0)*fun4A(lam11, lam20, nu10_lp, nu20_lp)$probs)
    prop[zzI] <- prop[zzI] + log(pi_z0)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[1] - beta_m[1])^2 - 
                          (beta1[1] - beta_m[1])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[1] <- beta1[1]
      lam10 <- lam11
      lmu10_index <- lmu11_index
      current <- prop
    }
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log((1 - pi_z0)*fun4A(lam10, lam21, nu10_lp, nu20_lp)$probs)
    prop[zzI] <- prop[zzI] + log(pi_z0)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[2] - beta_m[2])^2 - 
                          (beta1[2] - beta_m[2])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[2] <- beta1[2]
      lmu20_index <- lmu21_index
      current <- prop
    }
    print(beta0)
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.4)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- log((1 - pi_z0)*fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
      prop[zzI] <- prop[zzI] + log(pi_z0)
      denrat <- with(dat, tapply(prop - current, homeID, sum) +
                       tapply(prop - current, awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log((1 - pi_z0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    current[zzI] <- current[zzI] + log(pi_z0)
    
    # Team attack updates (alphas)
    alpha1 <- rnorm(nteams, alpha0, 0.03)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log((1 - pi_z0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    prop[zzI] <- prop[zzI] + log(pi_z0)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log((1 - pi_z0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    current[zzI] <- current[zzI] + log(pi_z0)
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.03)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log((1 - pi_z0)*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    
    prop[zzI] <- prop[zzI] + log(pi_z0)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Update current log-likelihood based on alpha updates
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log((1 - pi_z0)*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    current[zzI] <- current[zzI] + log(pi_z0)
    
    # Draw inflation update
    
    
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    pis[k,] <- pi_z0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_SA, data = dat)
}


#### Yellow cards ####
load("Data/YC_skel.RData")
k_range <- min(YC_skel$card_dif):max(YC_skel$card_dif)
sum_index2_l <- pmax(0, -k_range)
sum_index2_u <- sum_index2_l + 10
sum_index1_l <- sum_index2_l + k_range
sum_index1_u <- sum_index2_u + k_range
m1 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index1_l), length(k_range), 11, byrow = T)
m2 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index2_l), length(k_range), 11, byrow = T)

#### Model to fit MPCMP-Skellam with league-wide nu ####
fit_mpcmp_skel_YC <- function(dat, niters = 1000, skel = FALSE){
  # Extract data
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$card_dif
  mdiff <- match(diffs, min(dat$card_dif):max(dat$card_dif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  # Get league-specific team IDs for league-specific mean-centering
  PMP <- table(dat$leagueID, dat$homeID)
  tl_id <- rep(1, nteams)
  tl_id[(1:nteams)[PMP[2,] != 0]] <- 2
  tl_id[(1:nteams)[PMP[3,] != 0]] <- 3
  tl_id[(1:nteams)[PMP[4,] != 0]] <- 4
  tl_id[(1:nteams)[PMP[5,] != 0]] <- 5
  nleagues <- length(table(dat$leagueID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  thetas <- matrix(NA, nrow = K, ncol = nteams) # Aggression parameter
  xis <- matrix(NA, nrow = K, ncol = nleagues) # League dispersions
  psis <- matrix(NA, nrow = K, ncol = nleagues) # League effects
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals2(dat)
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  theta0 <- starts$theta0
  psi0 <- c(0, starts$psi0) #rep(0, nleagues)
  # Priors
  beta_m <- beta0
  beta_s <- 0.03
  theta_m <- starts$theta0
  theta_s <- 0.15
  nu_s <- 0.50
  psi_m <- psi0
  psi_s <- 0.05
  # Initialise linear predictors
  mu10 <- with(dat, exp(beta01 + theta0[homeID] + psi0[leagueID]))
  mu20 <- with(dat, exp(beta02 + theta0[awayID] + psi0[leagueID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nleagues)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$leagueID]
  nu20_lp <- nu0[dat$leagueID]
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.01
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (one-at-a-time)
    beta1 <- rnorm(2, beta0, 0.03)
    mu11 <- with(dat, exp(beta1[1] + theta0[homeID] + psi0[leagueID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    prop <- log(fun4A(lam11, lam20, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[1] - beta_m[1])^2 - 
                          (beta1[1] - beta_m[1])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[1] <- beta1[1]
      lam10 <- lam11
      mu10 <- mu11
      current <- prop
    }
    mu21 <- with(dat, exp(beta1[2] + theta0[awayID] + psi0[leagueID]))
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam10, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + ((beta0[2] - beta_m[2])^2 - 
                          (beta1[2] - beta_m[2])^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0[2] <- beta1[2]
      lam20 <- lam21
      mu20 <- mu21
      current <- prop
    }
    # Update current log-likelihood based on beta updates
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    #current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nleagues, xi0, 0.25)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$leagueID])
      nu21_lp <- exp(xi1[dat$leagueID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
      denrat <- with(dat, tapply(prop - current, leagueID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*nu_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nleagues)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$leagueID]
    nu20_lp <- nu0[dat$leagueID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # theta updates - vector of teams, component-wise updates
    theta1 <- rnorm(nteams, theta0, 0.03)
    # Mean-centering by league
    #theta1 <- theta1 - mean(theta1)
    league_mean <- tapply(theta1, tl_id, mean)
    theta1 <- theta1 - league_mean[tl_id]
    mu11 <- with(dat, exp(beta0[1] + theta1[homeID] + psi0[leagueID]))
    mu21 <- with(dat, exp(beta0[2] + theta1[awayID] + psi0[leagueID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    # Why the line below in all the code? Hangover from block update
    #denrat <- sum(prop - current)
    denrat <- denrat + sum(((theta0 - theta_m)^2 -
                              (theta1 - theta_m)^2)/(2*theta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    theta0[accept] <- theta1[accept]
    #theta0 <- rep(0, nteams)
    # Update current log-likelihood based on theta updates
    mu10 <- with(dat, exp(beta0[1] + theta0[homeID] + psi0[leagueID]))
    mu20 <- with(dat, exp(beta0[2] + theta0[awayID] + psi0[leagueID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    # psi updates
    psi1 <- rnorm(nleagues, psi0, 0.05)
    mu11 <- with(dat, exp(beta0[1] + theta0[homeID] + psi1[leagueID]))
    mu21 <- with(dat, exp(beta0[2] + theta0[awayID] + psi1[leagueID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, leagueID, sum))
    denrat <- denrat + sum(((psi0 - psi_m)^2 - (psi1 - psi_m)^2)/(2*psi_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nleagues)) < laccept)
    psi0[accept] <- psi1[accept]
    # Corner constraint on league effects (for now)
    psi0[1] <- 0
    #print(psi0)
    # Update current log-likelihood based on psi updates
    mu10 <- with(dat, exp(beta0[1] + theta0[homeID] + psi0[leagueID]))
    mu20 <- with(dat, exp(beta0[2] + theta0[awayID] + psi0[leagueID]))
    if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
    if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
    if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
    if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
    lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
    lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    thetas[k, ] <- theta0
    xis[k,] <- xi0
    psis[k, ] <- psi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, thetas, xis, psis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, thetas = thetas, xis = xis, psis = psis,
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_YC, data = dat)
}


#### Model to fit Skellam ####
# Use existing code and fix nu = 1 as a patch
# Compare with footBayes for sanity check

#### Model to fit MPCMP-Skellam with fixed nu ####





#### Serie A fit (block updates) ####
fit_mpcmp_skel_SA_block <- function(dat, niters = 1000, skel = FALSE){
  # Extract data
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  # An if statement here?
  mdiff <- match(diffs, (min(dat$goaldif) - 3):(max(dat$goaldif) + 3))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 2) # Overall means for home/away
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals_cor(dat)
  print(starts)
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  # Priors
  beta_m <- beta0
  beta_s <- 0.02
  alpha_m <- alpha0
  alpha_s <- 0.15
  delta_m <- delta0
  delta_s <- 0.15
  xi_s <- 0.5
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.001
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update 
    beta1 <- rnorm(2, beta0, 0.035)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((beta0 - beta_m)^2 - 
                          (beta1 - beta_m)^2)/(2*beta_s^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lam10 <- lam11
      lam20 <- lam21
      lmu10_index <- lmu11_index
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu block update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.12)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
      denrat <- sum(prop - current)
        #with(dat, tapply(prop - current, homeID, sum) +
         #              tapply(prop - current, awayID, sum))
      denrat <- denrat + sum((xi0^2 - xi1^2)/(2*xi_s^2))
      laccept <- pmin(0, denrat)
      accept <- log(runif(1)) < laccept
      if (accept){
        xi0 <- xi1
        current <- prop
        lam10 <- lam11
        lam20 <- lam21
        nu10_lp <- nu11_lp
        nu20_lp <- nu21_lp
        nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
        nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
        }
    }
    nu0 <- exp(xi0)
    
    # alpha updates - vector of teams, block update and corner constraint
    alpha1 <- rnorm(nteams, alpha0, 0.035)
    alpha1[1] <- 0 # Corner constraint - check data for team 1
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){
      alpha0 <- alpha1
      # Update current log-likelihood based on alpha updates
      mu10 <- mu11
      mu20 <- mu21
      if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
      lmu10_index <- lmu11_index
      lmu20_index <- lmu21_index
      lam10 <- lam11
      lam20 <- lam21
      current <- prop
    }
    
    # Team defence block update (deltas)
    delta1 <- rnorm(nteams, delta0, 0.035)
    # Corner constraint
    delta1[1] <- 0
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){
      delta0 <- delta1
      # Update current log-likelihood based on alpha updates
      mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
      mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
      lmu10_index <- lmu11_index
      lmu20_index <- lmu21_index
      lam10 <- lam11
      lam20 <- lam21
      current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    }
    # Calculate log-likelihood
    lhood <- sum(current)
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates,
       model = fit_mpcmp_skel_SA, data = dat)
}

#### MLB fit (block updates) ####
# Assuming +1/-1 no home ad for for now
fit_mpcmp_skel_MLB_block <- function(dat, niters = 1000, skel = FALSE){
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, min(dat$goaldif):max(dat$goaldif))
  m11 <- m1[mdiff, ]
  m21 <- m2[mdiff, ]
  fm11 <- factorial(m11)
  fm21 <- factorial(m21)
  fM1 <- factorial(M1)
  fM2 <- factorial(M2)
  # Function to calculate probabilities
  fun4A <- function(lam1, lam2, nu1, nu2){
    PMP <- lam1^m11*lam2^m21/
      (factorial(m11)^nu1*factorial(m21)^nu2)
    GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
    list(probs = rowSums(PMP)/GG1)
  }
  # Indicator for one run home wins
  ii1R <- dat$hgoal == dat$vgoal + 1
  # Indicator for one run away wins
  ii1mR <- dat$hgoal == dat$vgoal - 1
  print(table(ii1R))
  print(table(ii1mR))
  # Some setting up
  K <- niters
  nteams <- length(table(dat$homeID))
  betas <- matrix(NA, nrow = K, ncol = 1) # Overall mean
  alphas <- matrix(NA, nrow = K, ncol = nteams) # Attacking strengths
  deltas <- matrix(NA, nrow = K, ncol = nteams) # Defensive strengths
  xis <- matrix(NA, nrow = K, ncol = nteams) # Team dispersions
  pis <- matrix(NA, nrow = K, ncol = 2) # Inflation parameters
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise via starting values
  starts <- start_vals_cor(dat)
  print(starts)
  beta0 <- mean(starts$beta0)
  alpha0 <- starts$alpha0
  delta0 <- starts$delta0
  pi0 <- rep(0.05, 2) # Could perhaps do better here
  # Priors
  beta_m <- beta0
  beta_s <- 0.02
  alpha_m <- alpha0
  alpha_s <- 0.1
  delta_m <- delta0 
  delta_s <- 0.1
  xi_s <- 0.5
  # Normal prior for logit(pi)
  pil_m <- -3 
  pil_s <- 0.25
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta0 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta0 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- rep(1, nteams)
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  # Lower bounds
  min_lmu <- log(0.01)
  min_nu <- 0.001
  step_lmu <- step_nu <- 0.001
  # Initalise lambdas
  lmu10_index <- floor((log(mu10) - min_lmu)/step_lmu + 1)
  nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
  lmu20_index <- floor((log(mu20) - min_lmu)/step_lmu + 1)
  nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
  lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
  lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
  current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
  current[ii1R] <- current[ii1R] + pi0[1]
  current[ii1mR] <- current[ii1mR] + pi0[2]
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update (now represents overall mean for both H/A)
    beta1 <- rnorm(1, beta0, 0.05)
    mu11 <- with(dat, exp(beta1 + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1 + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0[1]
    prop[ii1mR] <- prop[ii1mR] + pi0[2]
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + ((beta0 - beta_m)^2 - 
                          (beta1 - beta_m)^2)/(2*beta_s^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lam10 <- lam11
      lmu10_index <- lmu11_index
      lam20 <- lam21
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu update - on log scale (xi)
    if (!skel){
      xi1 <- rnorm(nteams, xi0, 0.50)
      if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
      if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
      nu11_lp <- exp(xi1[dat$homeID])
      nu21_lp <- exp(xi1[dat$awayID])
      nu11_index <- floor((nu11_lp - min_nu)/step_nu + 1)
      nu21_index <- floor((nu21_lp - min_nu)/step_nu + 1)
      lam11 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu11_index)])
      lam21 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu21_index)])
      prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs
      prop[ii1R] <- prop[ii1R] + pi0[1]
      prop[ii1mR] <- prop[ii1mR] + pi0[2]
      denrat <- with(dat, tapply(log(prop) - log(current), homeID, sum) +
                       tapply(log(prop) - log(current), awayID, sum))
      denrat <- denrat + (xi0^2 - xi1^2)/(2*xi_s^2)
      laccept <- pmin(0, denrat)
      accept <- (log(runif(nteams)) < laccept)
      xi0[accept] <- xi1[accept]
    }
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    nu10_index <- floor((nu10_lp - min_nu)/step_nu + 1)
    nu20_index <- floor((nu20_lp - min_nu)/step_nu + 1)
    lam10 <- exp(mpcmp_fast_grid[cbind(lmu10_index, nu10_index)])
    lam20 <- exp(mpcmp_fast_grid[cbind(lmu20_index, nu20_index)])
    current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0[1]
    current[ii1mR] <- current[ii1mR] + pi0[2]
    
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.025)
    alpha1[1] <- 0 # Corner constraint - check data for team 1
    mu11 <- with(dat, exp(beta0 + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0 + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0[1]
    prop[ii1mR] <- prop[ii1mR] + pi0[2]
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + sum(((alpha0 - alpha_m)^2 - 
                              (alpha1 - alpha_m)^2)/(2*alpha_s^2))
    laccept <- pmin(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){
      alpha0 <- alpha1
      # Update current log-likelihood based on alpha updates
      mu10 <- mu11
      mu20 <- mu21
      lmu10_index <- lmu11_index
      lmu20_index <- lmu21_index
      lam10 <- lam11
      lam20 <- lam21
      current <- prop
    }
    
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.025)
    # Corner constraint
    delta1[1] <- 0
    mu11 <- with(dat, exp(beta0 + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0 + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 100){mu11 <- pmin(100, mu11)}
    if (min(mu11) < 0.01){mu11 <- pmax(0.01, mu11)}
    if (max(mu21) >= 100){mu21 <- pmin(100, mu21)}
    if (min(mu21) < 0.01){mu21 <- pmax(0.01, mu21)}
    lmu11_index <- floor((log(mu11) - min_lmu)/step_lmu + 1)
    lmu21_index <- floor((log(mu21) - min_lmu)/step_lmu + 1)
    lam11 <- exp(mpcmp_fast_grid[cbind(lmu11_index, nu10_index)])
    lam21 <- exp(mpcmp_fast_grid[cbind(lmu21_index, nu20_index)])
    prop <- (1 - sum(pi0))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi0[1]
    prop[ii1mR] <- prop[ii1mR] + pi0[2]
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + sum(((delta0 - delta_m)^2 - 
                              (delta1 - delta_m)^2)/(2*delta_s^2))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      delta0 <- delta1
      # Update current log-likelihood based on alpha updates
      mu10 <- mu11
      mu20 <- mu21
      lmu10_index <- lmu11_index
      lmu20_index <- lmu21_index
      lam10 <- lam11
      lam20 <- lam21
      current <- prop
    }
    
    # Inflation parameter update
    pil0 <- log(pi0/(1 - pi0)) 
    pil1 <- rnorm(2, pil0, 0.125)
    pi1 <- exp(pil1)/(1 + exp(pil1))
    prop <- (1 - sum(pi1))*fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs
    prop[ii1R] <- prop[ii1R] + pi1[1]
    prop[ii1mR] <- prop[ii1mR] + pi1[2]
    denrat <- sum(log(prop) - log(current))
    denrat <- denrat + sum(((pil0 - pil_m)^2 - 
                              (pil1 - pil_m)^2)/(2*pil_s^2))
    laccept <- min(0, denrat)
    accept <- log(runif(1)) < laccept
    if (accept){pil0 <- pil1}
    pi0 <- exp(pil0)/(1 + exp(pil0))
    print(pi0)
    current <- (1 - sum(pi0))*fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs
    current[ii1R] <- current[ii1R] + pi0[1]
    current[ii1mR] <- current[ii1mR] + pi0[2]
    
    lhood <- sum(log(current))
    # Store parameters
    betas[k, ] <- beta0 
    alphas[k, ] <- alpha0
    deltas[k, ] <- delta0
    xis[k,] <- xi0
    pis[k,] <- pi0
    ll[k, ] <- lhood
    iter_times[k] <- proc.time()[3] - start
  }
  cpu_time <- proc.time() - ptm
  paras <- cbind(betas, alphas, deltas, xis, pis)
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/K
  }
  list(betas = betas, alphas = alphas, deltas = deltas, xis = xis, 
       pis = pis, ll = ll, cpu_time = cpu_time[3], 
       acc_rates = acc_rates, model = fit_mpcmp_skel_MLB_inf, 
       data = dat)
}
