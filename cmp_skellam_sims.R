# Simulations for CMP-Skellam model
devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)
# Three scenarios:
# 1. Both overdispersed
# 2. Both underdispersed
# 3. One over, one under
# Wrap them in one routine
sim_cmp_skellam <- function(nseasons = 5, nteam = 20, sdteam = 0.1,
                            beta = 0.4, gamma = 0.2, nu1 = 0.8, nu2 = 0.7){
  theta <- rnorm(nteam, 0, sdteam)
  # Corner constraint on first theta
  theta <- theta - theta[1]
  homeid <- rep(1:nteam, each = nteam)
  awayid <- rep(1:nteam, nteam)
  fixtures <- data.frame(hid = homeid, aid = awayid)
  fixtures <- fixtures[fixtures$hid != fixtures$aid, ]
  fixtures_big <- matrix(rep(t(fixtures), nseasons),
                         ncol = ncol(fixtures), byrow = TRUE)
  fixtures_big <- as.data.frame(fixtures_big)
  colnames(fixtures_big) <- c("hid", "aid")
  n <- nrow(fixtures_big)
  mu1 <- exp(beta + gamma + theta[fixtures_big$hid])
  mu2 <- exp(beta + theta[fixtures_big$aid])
  y1 <- rcomp(n, mu1, nu1)
  y2 <- rcomp(n, mu2, nu2)
  diff <- y1 - y2
  d <- data.frame(hid = fixtures_big$hid, aid = fixtures_big$aid, 
                  hg = y1, ag = y2, diff = diff)
  list(dat = d, team_str = theta)
}

#### Fit for overall strength parameter ####
# Dispersion for home/away
fit_sim_cmp_skellam <- function(data, niter = 1000){
  dat <- data$dat
  team_str <- data$team_str
  diffs <- dat$diff
  mdiff <- match(diffs, -10:10)
  n <- nrow(dat)
  lam1s <- NA
  lam2s <- NA
  nu1s <- NA
  nu2s <- NA
  ntheta <- length(table(dat$hid))
  K <- niter
  betas <- matrix(NA, nrow = K, ncol = 1) # Overall mean
  thetas <- matrix(NA, nrow = K, ncol = ntheta) # Team effects
  gammas <- matrix(NA, nrow = K, ncol = 1) # Home advantage
  nus <- matrix(NA, nrow = K, ncol = 2) # Dispersions for H/A
  iter_times <- rep(0, K)
  # Initialise (tried Skellam in footBayes with MLE - didn't work!)
  beta0 <- log(mean(dat$ag))
  theta0 <- team_str
  gamma0 <- log(mean(dat$hg)) - beta0
  lmu10 <- with(dat, beta0 + theta0[hid] + gamma0)
  lmu20 <- with(dat, beta0 + theta0[aid])
  nu10 <- nu20 <- 1.5
  # Find lambda (on log-scale!)
  lmu10_index <- round((lmu10 - log(0.01))/0.001 + 1)
  lmu20_index <- round((lmu20 - log(0.01))/0.001 + 1)
  lam10 <- mpcmp_fast_grid[cbind(lmu10_index, 
                                   round(nu10*1000))]
  lam20 <- mpcmp_fast_grid[cbind(lmu20_index, 
                                   round(nu20*1000))]
  current <- log(fun4B(lam10, lam20, nu10, nu20)$probs)
  date()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta update
    lmu11 <- with(dat, beta1 + theta0[hid] + gamma0)
    lmu21 <- with(dat, beta1 + theta0[aid])
    # Find lambda (on log-scale!)
    lmu11_index <- round((lmu11 - log(0.01))/0.001 + 1)
    lmu21_index <- round((lmu21 - log(0.01))/0.001 + 1)
    lam11 <- mpcmp_fast_grid[cbind(lmu11_index,
                                   round(nu10*1000))]
    lam21 <- mpcmp_fast_grid[cbind(lmu21_index,
                                   round(nu20*1000))]
    prop <- log(fun4B(lam11, lam21, nu10, nu20)$probs)
    denrat <- sum(prop - current)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      current <- prop
    }
    # gamma update
    gamma1 <- rnorm(1, gamma0, 0.10)
    lmu11 <- with(dat, beta0 + theta0[hid] + gamma1)
    lmu21 <- with(dat, beta0 + theta0[aid])
    # Find lambda (on log-scale!)
    lmu11_index <- round((lmu11 - log(0.01))/0.001 + 1)
    lmu21_index <- round((lmu21 - log(0.01))/0.001 + 1)
    lam11 <- mpcmp_fast_grid[cbind(lmu11_index,
                                   round(nu10*1000))]
    lam21 <- mpcmp_fast_grid[cbind(lmu21_index,
                                   round(nu20*1000))]
    prop <- log(fun4B(lam11, lam21, nu10, nu20)$probs)
    denrat <- sum(prop - current)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      gamma0 <- gamma1
      current <- prop
    }
    # theta updates - vector of teams, component-wise updates
    theta1 <- rnorm(ntheta, theta0, 0.15)
    theta1[1] <- 0
    lmu11 <- with(dat, beta0 + theta1[hid] + gamma0)
    lmu21 <- with(dat, beta0 + theta1[aid])
    # Find lambda (on log-scale!)
    lmu11_index <- round((lmu11 - log(0.01))/0.001 + 1)
    lmu21_index <- round((lmu21 - log(0.01))/0.001 + 1)
    lam11 <- mpcmp_fast_grid[cbind(lmu11_index,
                                   round(nu10*1000))]
    lam21 <- mpcmp_fast_grid[cbind(lmu21_index,
                                   round(nu20*1000))]
    prop <- log(fun4B(lam11, lam21, nu10, nu20)$probs)
    denrat <- with(dat, tapply(prop - current, hid, sum) +
                     tapply(prop - current, aid, sum))
    laccept <- pmin(0, denrat)
    accept <- (log(runif(ntheta)) < laccept)
    theta0[accept] <- theta1[accept]
    # Update linear predictor for nu updates
    lmu10 <- with(dat, beta0 + theta0[hid] + gamma0)
    lmu20 <- with(dat, beta0 + theta0[aid])
    # Find lambda (on log-scale!)
    lmu10_index <- round((lmu10 - log(0.01))/0.001 + 1)
    lmu20_index <- round((lmu20 - log(0.01))/0.001 + 1)
    lam10 <- mpcmp_fast_grid[cbind(lmu10_index,
                                   round(nu10*1000))]
    lam20 <- mpcmp_fast_grid[cbind(lmu20_index,
                                   round(nu20*1000))]
    current <- log(fun4B(lam10, lam20, nu10, nu20)$probs)
    # nu updates
    nuprop <- rlnorm(2, c(log(nu10), log(nu20)), 0.25)
    nu11 <- nuprop[1]; nu21 <- nuprop[2]
    # Find lambda (on log-scale!)
    # Use updated linear predictor from theta update above
    lam11 <- mpcmp_fast_grid[cbind(lmu10_index,
                                   round(nu11*1000))]
    lam21 <- mpcmp_fast_grid[cbind(lmu20_index,
                                   round(nu21*1000))]
    prop <- log(fun4B(lam11, lam21, nu11, nu21)$probs)
    denrat <- sum(prop - current) + sum(log(nuprop/c(nu10, nu20)))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    nuupdate <- c(nu10, nu20)
    if (accept){
      nuupdate <- nuprop
      nu10 <- nuupdate[1]; nu20 <- nuupdate[2]
      current <- prop
    }
    print(c(beta0, gamma0, theta0))
    print(c(nu10, nu20))
    # Store parameters
    betas[k] <- beta0 
    gammas[k] <- gamma0
    thetas[k, ] <- theta0
    nus[k,] <- c(nu10, nu20)
    iter_times[k] <- proc.time()[3] - start
  }
  list(betas = betas, gammas = gammas, thetas = thetas, nus = nus,
       iter_times = iter_times)
}

# Simulate data
simdat <- sim_cmp_skellam()
n <- nrow(simdat$dat)
# Some pre-fit formatting and preparation
M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
mdiff <- match(simdat$dat$diff, -10:10)
m11 <- m1[mdiff,]
m21 <- m2[mdiff,]
fm11 <- factorial(m11)
fm21 <- factorial(m21)
fM1 <- factorial(M1)
fM2 <- factorial(M2)
fun4A <- function(lam1, lam2, nu1, nu2){
  PMP <- lam1^m11*lam2^m21/
    (factorial(m11)^nu1*factorial(m21)^nu2)
  GG1 <- rowSums(lam1^M1*lam2^M2/(fM1^nu1*fM2^nu2))
  list(probs = rowSums(PMP)/GG1)
}
fun4B <- function(llam1, llam2, nu1, nu2){
  PMP <- exp(m11*llam1 + m21*llam2 -
               nu1*lfactorial(m11) - nu2*lfactorial(m21))
  GG1 <- rowSums(exp(M1*llam1 + M2*llam2 - 
                       (nu1*lfactorial(M1) + nu2*lfactorial(M2))))
  #GG1 <- rowSums(lam1^M1*lam2^M2/(factorial(M1)^nu1*factorial(M2)^nu2))
  list(probs = rowSums(PMP)/GG1)
}
# Now fit
date(); run <- fit_sim_cmp_skellam(simdat, niter = 1000); date()

# Try with just three thetas and lots of seasons
sim_simple <- sim_cmp_skellam(nseasons = 100, nteam = 3, 
                              nu1 = 1, nu2 = 1)
n <- nrow(sim_simple$dat)
# Some pre-fit formatting and preparation
M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
table(sim_simple$dat$diff)
mdiff <- match(sim_simple$dat$diff, -10:10)
m11 <- m1[mdiff,]
m21 <- m2[mdiff,]
fm11 <- factorial(m11)
fm21 <- factorial(m21)
fM1 <- factorial(M1)
fM2 <- factorial(M2)
date(); run <- fit_sim_cmp_skellam(sim_simple, niter = 1000); date()
# Compare with separate fits?
# Try underdisp too
# Frequentist Skellam anywhere else in R?
library(skellam)

# Fit to EPL data
n <- nrow(test$dat)
# Some pre-fit formatting and preparation
M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
table(test$dat$diff)
mdiff <- match(test$dat$diff, -10:10)
m11 <- m1[mdiff,]
m21 <- m2[mdiff,]
fm11 <- factorial(m11)
fm21 <- factorial(m21)
fM1 <- factorial(M1)
fM2 <- factorial(M2)
date(); run_epl <- fit_sim_cmp_skellam(test, niter = 10); date()

#### Simulation to mimic (small) league data with fixed nu ####
sim_cmp_skellam_epl <- function(nseasons = 5, nteam = 20, sdteam = 0,
                            beta1 = 0.6, beta2 = 0.4, nu = 1){
  alpha <- rnorm(nteam, 0, sdteam)
  delta <- rnorm(nteam, 0, sdteam)
  #xi <- rnorm(nteam, log(1), 0.2) # Tweak this... fix some?
  nu <- seq(nu, nu, length = nteam)
  # Sum-to-zero constraints on attack/defence
  alpha <- alpha - mean(alpha)
  delta <- delta - mean(delta)
  homeid <- rep(1:nteam, each = nteam)
  awayid <- rep(1:nteam, nteam)
  fixtures <- data.frame(hid = homeid, aid = awayid)
  fixtures <- fixtures[fixtures$hid != fixtures$aid, ]
  fixtures_big <- matrix(rep(t(fixtures), nseasons),
                         ncol = ncol(fixtures), byrow = TRUE)
  fixtures_big <- as.data.frame(fixtures_big)
  colnames(fixtures_big) <- c("hid", "aid")
  n <- nrow(fixtures_big)
  mu1 <- exp(beta1 + alpha[fixtures_big$hid] - delta[fixtures_big$aid])
  mu2 <- exp(beta2 + alpha[fixtures_big$aid] - delta[fixtures_big$hid])
  nu1_long <- nu[fixtures_big$hid]
  nu2_long <- nu[fixtures_big$aid]
  y1 <- rcomp(n, mu = mu1, nu = nu1_long)
  y2 <- rcomp(n, mu = mu2, nu = nu2_long)
  diff <- y1 - y2
  d <- data.frame(homeID = fixtures_big$hid, awayID = fixtures_big$aid, 
                  hgoal = y1, vgoal = y2, goaldif = diff)
  list(dat = d, beta1 = beta1, beta2 = beta2, 
       team_att = alpha, team_def = delta, true_nu = nu)
}
# Run it
test_dat <- sim_cmp_skellam_epl(nseasons = 500, nteam = 4, sdteam = 0.1)
# Pre-processing
n <- nrow(test_dat$dat)
M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
table(test_dat$dat$goaldif)
mdiff <- match(test_dat$dat$goaldif, -10:10)
m11 <- m1[mdiff,]
m21 <- m2[mdiff,]
fm11 <- factorial(m11)
fm21 <- factorial(m21)
fM1 <- factorial(M1)
fM2 <- factorial(M2)
fit_test_dat <- fit_epl_att_def3(test_dat$dat, niter = 1000)

#### Simulation studies ####
# 1. Fixed nu - compare with footBayes (nu = 1, nu < 1, nu > 1)
# This also gives effect of misspecified model
# Can also check if sim is correct?
# 2. Team-specific nu (nu  = 1, nu < 1, nu > 1, nu both. < & > 1)
# Compare with Skellam by fixing nu = 1 in code
# Actually, better to have separate fit for that to save fiddling

#### First sims ####
sim_skel_dat <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                    sdteam = 0.25, nu = 1)
# Test fit
fit_skel <- fit_att_def_skel(sim_skel_dat, 1000)

#### Sim study - equidispersion ####
res_const_nu <- NULL
true_vals <- NULL
for (i in 1:50){
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                   sdteam = 0.25, nu = 1)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel(info = sim_info, 1000)
  res_const_nu <- rbind(res_const_nu, colMeans(cbind(run$betas, run$alphas, run$deltas, 
                                               exp(run$xis))))
  true_vals <- rbind(true_vals, c(sim_info$beta1, sim_info$beta2,
                                  sim_info$team_att, sim_info$team_def, 
                                  sim_info$true_nu))
}

#### Sim study - underdispersion ####
res_const_nu_under <- NULL
true_vals_under <- NULL
for (i in 1:25){
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                  sdteam = 0.25, nu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel(info = sim_info, 1000)
  res_const_nu_under <- rbind(res_const_nu_under, 
                              colMeans(cbind(run$betas, run$alphas, 
                                             run$deltas, exp(run$xis))))
  true_vals_under <- rbind(true_vals_under, c(sim_info$beta1, sim_info$beta2,
                                  sim_info$team_att, sim_info$team_def, 
                                  sim_info$true_nu))
}

#### Sim study - overdispersion ####
res_const_nu_over <- NULL
true_vals_over <- NULL
for (i in 1:25){
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                  sdteam = 0.25, nu = 0.5)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel(info = sim_info, 1000)
  res_const_nu_over <- rbind(res_const_nu_over, 
                             colMeans(cbind(run$betas, run$alphas, 
                                            run$deltas,exp(run$xis))))
  true_vals_over <- rbind(true_vals_over, c(sim_info$beta1, sim_info$beta2,
                                  sim_info$team_att, sim_info$team_def, 
                                  sim_info$true_nu))
}


#### Fit to simulated data for attack and defence (Skellam) - fixed nu ####
# Can use Bessel functions to check when nu is fixed at 1
# Can compare with footBayes
fit_att_def_skel <- function(info, niters = 1000){
  dat <- info$dat
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
  xis <- matrix(NA, nrow = K, ncol = 1) # Overall (Poisson) dispersion
  ll <- matrix(NA, nrow = K, ncol = 1) # Posterior likelihood
  iter_times <- rep(0, K)
  # Initialise
  beta01 <- log(mean(dat$hgoal))
  beta02 <- log(mean(dat$vgoal))
  beta_p <- c(beta01, beta02) # To use as prior
  beta0 <- c(beta01, beta02)
  # Sum-to-zero constraint on attack and defence strengths
  alpha0 <- info$team_att
  #rnorm(nteams, 0, 0.1); alpha0 <- alpha0 - mean(alpha0)
  delta0 <- info$team_def
  #rnorm(nteams, 0, 0.1); delta0 <- delta0 - mean(delta0)
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
  if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
  if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
  if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
  nu0 <- 1
  xi0 <- log(nu0)
  # Get lambdas from grid for now
  # Find lambda
  nu0_lp <- rep(nu0, n)
  lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000), 
                                   round(nu0_lp*1000))] 
  lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000), 
                                   round(nu0_lp*1000))]
  # Need a bigger grid?
  current <- log(fun4A(lam10, lam20, nu0_lp, nu0_lp)$probs)
  ptm <- proc.time()
  for (k in 1:K){
    if(k %% 100 == 0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    start <- proc.time()[3]
    # beta block update
    prop_covar <- 0.0001*matrix(c(0.15, 0.7*0.15, 0.7*0.15, 0.15), 2, 2)
    #beta1 <- rmvnorm(1, beta0, prop_covar)[1, ]
    beta1 <- rnorm(2, beta0, 0.06)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu0_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu0_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu0_lp, nu0_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((beta0 - beta_p)^2 - 
                              (beta1 - beta_p)^2)/(2*0.25^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lam10 <- lam11
      lam20 <- lam21
      current <- prop
      # No mu in here?
    }
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.15)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID] ))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu0_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu0_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu0_lp, nu0_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (alpha0^2 - alpha1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    # Mean-centering
    alpha0 <- alpha0 - mean(alpha0)
    # Update linear predictors
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu0_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu0_lp*1000))]
    current <- log(fun4A(lam10, lam20, nu0_lp, nu0_lp)$probs)
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.15)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu0_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu0_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu0_lp, nu0_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (delta0^2 - delta1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Sum-to-zero
    delta0 <- delta0 - mean(delta0)
    # Update linear predictors
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu0_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu0_lp*1000))]
    current <- log(fun4A(lam10, lam20, nu0_lp, nu0_lp)$probs)
    # # nu update - on log scale (xi)
    #xi0 <- log(1)
    xi1 <- rnorm(1, xi0, 0.20)
    if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
    if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
    nu1_lp <- rep(exp(xi1), n)
    lam11 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu1_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu1_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu1_lp, nu1_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + (xi0^2 - xi1^2)/(2*0.25^2)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){xi0 <- xi1}
    nu0 <- exp(xi0)
    #print(nu0)
    # Update current log-likelihood based on xi(nu) updates
    nu0_lp <- rep(nu0, n)
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu0_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu0_lp*1000))]
    current <- log(fun4A(lam10, lam20, nu0_lp, nu0_lp)$probs)
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
       ll = ll, cpu_time = cpu_time[3], acc_rates = acc_rates)
}

#### Fit to simulated data for attack and defence (MPCMP-Skellam) - fixed nu ####
# Underdispersed
sim_skel_dat_mis1 <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                         sdteam = 0.1, nu = 2)
fit_skel_mis1 <- fit_att_def_skel(sim_skel_dat_mis1$dat, 1000)
# Overdispersed
sim_skel_dat_mis2 <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                         sdteam = 0.1, nu = 0.5)
fit_skel_mis2 <- fit_att_def_skel(sim_skel_dat_mis2$dat, 1000)
# Need to check the size of the differences in the below
sim_skel_dat_mis3 <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                         sdteam = 0.25, nu = 0.5)
fit_skel_mis3 <- fit_att_def_skel(sim_skel_dat_mis3$dat, 1000)

#### Fit to simulated data for attack and defence (Skellam) - nu by team ####
#### Fit to simulated data for attack and defence (MPCMP-Skellam) - nu by team ####
# Can use Bessel functions to check when nu fixed as 1
# Can only compare to footBayes when all nu are 1 (fixed, or estimated?)

#### Simulation to mimic (small) league data with team-varying nu ####
sim_cmp_skellam_epl2 <- function(nseasons = 5, nteam = 20, sdteam = 0,
                                beta1 = 0.6, beta2 = 0.4, 
                                nul = 1, nuu = 1){
  #alpha <- rnorm(nteam, 0, sdteam)
  #delta <- rnorm(nteam, 0, sdteam)
  alpha <- seq(-0.25, 0.25, length = nteam)
  delta <- seq(-0.25, 0.25, length = nteam)
  nu <- seq(nul, nuu, length = nteam)
  # Sum-to-zero constraints on attack/defence
  alpha <- alpha - mean(alpha)
  delta <- delta - mean(delta)
  homeid <- rep(1:nteam, each = nteam)
  awayid <- rep(1:nteam, nteam)
  fixtures <- data.frame(hid = homeid, aid = awayid)
  fixtures <- fixtures[fixtures$hid != fixtures$aid, ]
  fixtures_big <- matrix(rep(t(fixtures), nseasons),
                         ncol = ncol(fixtures), byrow = TRUE)
  fixtures_big <- as.data.frame(fixtures_big)
  colnames(fixtures_big) <- c("hid", "aid")
  n <- nrow(fixtures_big)
  mu1 <- exp(beta1 + alpha[fixtures_big$hid] - delta[fixtures_big$aid])
  mu2 <- exp(beta2 + alpha[fixtures_big$aid] - delta[fixtures_big$hid])
  nu1_long <- nu[fixtures_big$hid]
  nu2_long <- nu[fixtures_big$aid]
  y1 <- rcomp(n, mu = mu1, nu = nu1_long)
  y2 <- rcomp(n, mu = mu2, nu = nu2_long)
  diff <- y1 - y2
  d <- data.frame(homeID = fixtures_big$hid, awayID = fixtures_big$aid, 
                  hgoal = y1, vgoal = y2, goaldif = diff)
  list(dat = d, beta1 = beta1, beta2 = beta2, 
       team_att = alpha, team_def = delta, true_nu = nu)
}

sim_skel_mpcmp_vary_nu <- sim_cmp_skellam_epl2(nseasons = 100, nteam = 6, 
                                         sdteam = 0.25, 
                                         nul = 2, nuu = 2)
fit_att_def_skel2 <- function(info, niters = 1000){
  # Extract data
  dat <- info$dat
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
  beta_p <- c(beta01, beta02)
  beta0 <- c(beta01, beta02)
  # Sum-to-zero constraint on attack and defence strengths
  alpha0 <- info$team_att
  #rnorm(nteams, 0, 0.1); alpha0 <- alpha0 - mean(alpha0)
  delta0 <- info$team_def
  #rnorm(nteams, 0, 0.1); delta0 <- delta0 - mean(delta0)
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
  if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
  if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
  if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
  nu0 <- info$true_nu
  xi0 <- log(nu0)
  # Get lambdas from grid for now
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000), 
                                   round(nu10_lp*1000))] 
  lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000), 
                                   round(nu20_lp*1000))]
  # Need a bigger grid?
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
    beta1 <- rnorm(2, beta0, 0.02)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu10_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu20_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum(((beta0 - beta_p)^2 - 
                              (beta1 - beta_p)^2)/(2*0.01^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lam10 <- lam11
      lam20 <- lam21
      current <- prop
      # No mu in here as included in lambda updates
    }
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.20)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu10_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu20_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (alpha0^2 - alpha1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    alpha0[accept] <- alpha1[accept]
    alpha0 <- alpha0 - mean(alpha0)
    # Update linear predictors
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu10_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu20_lp*1000))]
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.20)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu10_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu20_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (delta0^2 - delta1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    delta0[accept] <- delta1[accept]
    # Sum-to-zero
    delta0 <- delta0 - mean(delta0)
    # Update linear predictors
    mu10 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta0[awayID]))
    mu20 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu10_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu20_lp*1000))]
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    # nu update - on log scale (xi)
    #xi0 <- log(info$true_nu)
    xi1 <- rnorm(nteams, xi0, 0.20)
    if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
    if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
    nu11_lp <- exp(xi1[dat$homeID])
    nu21_lp <- exp(xi1[dat$awayID])
    lam11 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                      round(nu11_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                      round(nu21_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (xi0^2 - xi1^2)/(2*0.25^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    xi0[accept] <- xi1[accept]
    nu0 <- exp(xi0)
    print(nu0)
    # # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                      round(nu10_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                    round(nu20_lp*1000))]
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
       model = fit_att_def_skel2, data = dat)
}

# Test fit
fit_cmp_skel <- fit_att_def_skel2(info = sim_skel_mpcmp_vary_nu, 1000)

# Small batch of sims
res <- NULL
true_vals <- NULL
for (i in 1:50){
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 100, nteam = 6, 
                                                 sdteam = 0.25, 
                                                 nul = 0.5, nuu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel2(info = sim_info, 1000)
  res <- rbind(res, colMeans(cbind(run$betas, run$alphas, run$deltas, 
                          exp(run$xis))))
  true_vals <- rbind(true_vals, c(sim_info$beta1, sim_info$beta2,
               sim_info$team_att, sim_info$team_def, sim_info$true_nu))
}

# Small batch of sims with constant (but estimated nu)
# Overdispersion
res_over <- NULL
true_vals <- NULL
for (i in 1:25){
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 120, nteam = 6, 
                                   sdteam = 0.25, 
                                   nul = 0.5, nuu = 0.5)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel2(info = sim_info, 1000)
  res_over <- rbind(res_over, colMeans(cbind(run$betas, run$alphas, run$deltas, 
                                   exp(run$xis))))
  true_vals <- rbind(true_vals, c(sim_info$beta1, sim_info$beta2,
                                  sim_info$team_att, sim_info$team_def, sim_info$true_nu))
}

# Small batch of sims with constant (but estimated nu)
# Underdispersion
res_under <- NULL
true_vals <- NULL
for (i in 1:25){
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 100, nteam = 6, 
                                   sdteam = 0.25, 
                                   nul = 2, nuu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel2(info = sim_info, 1000)
  res_under <- rbind(res_under, colMeans(cbind(run$betas, run$alphas, run$deltas, 
                                   exp(run$xis))))
  true_vals <- rbind(true_vals, c(sim_info$beta1, sim_info$beta2,
                                  sim_info$team_att, sim_info$team_def, sim_info$true_nu))
}

#### Fit with block updates for alpha & delta ####
k_range <- -25:25
sum_index2_l <- pmax(0, -k_range)
sum_index2_u <- sum_index2_l + 10
sum_index1_l <- sum_index2_l + k_range
sum_index1_u <- sum_index2_u + k_range
m1 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index1_l), length(k_range), 11, byrow = T)
m2 <- matrix(sequence(rep(11, length(k_range)), 
                      from = sum_index2_l), length(k_range), 11, byrow = T)

fit_att_def_skel3 <- function(info, niters = 1000){
  # Extract data
  dat <- info$dat
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, -30:30)
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
  alpha0 <- info$team_att
  #rnorm(nteams, 0, 0.1); alpha0 <- alpha0 - mean(alpha0)
  delta0 <- info$team_def
  #rnorm(nteams, 0, 0.1); delta0 <- delta0 - mean(delta0)
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
  if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
  if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
  if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
  nu0 <- info$true_nu
  xi0 <- log(nu0)
  nu10_lp <- nu0[dat$homeID]
  nu20_lp <- nu0[dat$awayID]
  # Initalise lambdas (from grid for now)
  lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000), 
                                   round(nu10_lp*1000))] 
  lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000), 
                                   round(nu20_lp*1000))]
  # Need a bigger grid?
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
    beta1 <- rnorm(2, beta0, 0.02)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu10_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu20_lp*1000))]
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
    xi1 <- rnorm(nteams, xi0, 0.25)
    if (max(xi1) >= log(9.9)){xi1 <- pmin(log(9.9), xi1)}
    if (min(xi1) < log(0.001)){xi1 <- pmax(log(0.001), xi1)}
    nu11_lp <- exp(xi1[dat$homeID])
    nu21_lp <- exp(xi1[dat$awayID])
    lam11 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu11_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu21_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu11_lp, nu21_lp)$probs)
    denrat <- with(dat, tapply(prop - current, homeID, sum) +
                     tapply(prop - current, awayID, sum))
    denrat <- denrat + (xi0^2 - xi1^2)/(2*0.5^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    xi0[accept] <- xi1[accept]
    nu0 <- exp(xi0)
    # Update current log-likelihood based on xi(nu) updates
    nu10_lp <- nu0[dat$homeID]
    nu20_lp <- nu0[dat$awayID]
    lam10 <- lam_grid_poly_10K[cbind(round(mu10*1000),
                                     round(nu10_lp*1000))]
    lam20 <- lam_grid_poly_10K[cbind(round(mu20*1000),
                                     round(nu20_lp*1000))]
    current <- log(fun4A(lam10, lam20, nu10_lp, nu20_lp)$probs)
    # alpha updates - vector of teams, component-wise updates
    alpha1 <- rnorm(nteams, alpha0, 0.03)
    # Mean-centering
    alpha1 <- alpha1 - mean(alpha1)
    mu11 <- with(dat, exp(beta0[1] + alpha1[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha1[awayID] - delta0[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu10_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu20_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum((alpha0^2 - alpha1^2)/(2*0.1^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      alpha0 <- alpha1
      current <- prop
      }
    # Team defence updates (deltas)
    delta1 <- rnorm(nteams, delta0, 0.03)
    # Sum-to-zero
    delta1 <- delta1 - mean(delta1)
    mu11 <- with(dat, exp(beta0[1] + alpha0[homeID] - delta1[awayID]))
    mu21 <- with(dat, exp(beta0[2] + alpha0[awayID] - delta1[homeID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    lam11 <- lam_grid_poly_10K[cbind(round(mu11*1000),
                                     round(nu10_lp*1000))]
    lam21 <- lam_grid_poly_10K[cbind(round(mu21*1000),
                                     round(nu20_lp*1000))]
    prop <- log(fun4A(lam11, lam21, nu10_lp, nu20_lp)$probs)
    denrat <- sum(prop - current)
    denrat <- denrat + sum((delta0^2 - delta1^2)/(2*0.1^2))
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
       model = fit_att_def_skel2, data = dat)
}

#### Block model sim study - underdispersion ####
res_const_nu_under <- NULL
true_vals_under <- NULL
for (i in 1:50){
  print("Underdispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                  sdteam = 0.25, nu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel4(info = sim_info, 1000)
  res_const_nu_under <- rbind(res_const_nu_under, 
                              colMeans(cbind(run$betas, run$alphas, 
                                             run$deltas, exp(run$xis))))
  true_vals_under <- rbind(true_vals_under, c(sim_info$beta1, sim_info$beta2,
                                              sim_info$team_att, sim_info$team_def, 
                                              sim_info$true_nu))
}

#### Block model sim study - overdispersion ####
res_const_nu_over <- NULL
true_vals_over <- NULL
for (i in 1:50){
  print("Overdispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 100, nteam = 6, 
                                  sdteam = 0.25, nu = 0.5, 
                                  beta1 = 2, beta2 = 1.6)
  ii_u <- sim_info$dat$goaldif > 30
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -30
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel4(info = sim_info, 1000)
  res_const_nu_over <- rbind(res_const_nu_over, 
                             colMeans(cbind(run$betas, run$alphas, 
                                            run$deltas,exp(run$xis))))
  true_vals_over <- rbind(true_vals_over, c(sim_info$beta1, sim_info$beta2,
                                            sim_info$team_att, sim_info$team_def, 
                                            sim_info$true_nu))
}

#### Block model sim study - bidispersion ####
res_biv_nu <- NULL
true_vals_biv <- NULL
for (i in 1:50){
  print("Bidispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 100, nteam = 6, 
                                  sdteam = 0.25, nul = 0.5, nuu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel4(info = sim_info, 1000)
  res_biv_nu <- rbind(res_biv_nu, 
                             colMeans(cbind(run$betas, run$alphas, 
                                            run$deltas,exp(run$xis))))
  true_vals_biv <- rbind(true_vals_biv, c(sim_info$beta1, sim_info$beta2,
                                            sim_info$team_att, sim_info$team_def, 
                                            sim_info$true_nu))
}


#### Code with larger grid for mu - block updates ####
# Seems more accurate for large mu and allows for larger values
# Can add interpolation later if required
# Otherwise straightforward, but just need to take log of mu value
fit_att_def_skel4 <- function(info, niters = 1000){
  # Extract data
  dat <- info$dat
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, -30:30)
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
  alpha0 <- info$team_att
  #rnorm(nteams, 0, 0.1); alpha0 <- alpha0 - mean(alpha0)
  delta0 <- info$team_def
  #rnorm(nteams, 0, 0.1); delta0 <- delta0 - mean(delta0)
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- info$true_nu
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
    beta1 <- rnorm(2, beta0, 0.02)
    mu11 <- with(dat, exp(beta1[1] + alpha0[homeID] - delta0[awayID]))
    mu21 <- with(dat, exp(beta1[2] + alpha0[awayID] - delta0[homeID]))
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
    denrat <- denrat + sum(((beta0 - beta_p)^2 - 
                              (beta1 - beta_p)^2)/(2*0.01^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if(accept){
      beta0 <- beta1
      lmu10_index <- lmu11_index
      lmu20_index <- lmu21_index
      current <- prop
    }
    # nu update - on log scale (xi)
    xi1 <- rnorm(nteams, xi0, 0.25)
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
    denrat <- denrat + (xi0^2 - xi1^2)/(2*0.5^2)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(nteams)) < laccept)
    xi0[accept] <- xi1[accept]
    nu0 <- exp(xi0)
    print(nu0)
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
    denrat <- sum(prop - current)
    denrat <- denrat + sum((alpha0^2 - alpha1^2)/(2*0.1^2))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      alpha0 <- alpha1
      current <- prop
    }
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
    denrat <- sum(prop - current)
    denrat <- denrat + sum((delta0^2 - delta1^2)/(2*0.1^2))
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
       model = fit_att_def_skel4, data = dat)
}

#### Code with larger grid for mu - component-wise updates ####
fit_att_def_skel5 <- function(info, niters = 1000, skel = FALSE){
  # Extract data
  dat <- info$dat
  # Objects needed for fitting
  n <- nrow(dat)
  M1 <- matrix(as.vector(m1), n, ncol(m1)*nrow(m1), byrow = T)
  M2 <- matrix(as.vector(m2), n, ncol(m2)*nrow(m2), byrow = T)
  diffs <- dat$goaldif
  mdiff <- match(diffs, min(k_range):max(k_range))
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
  beta0 <- starts$beta0
  beta01 <- beta0[1]
  beta02 <- beta0[2]
  # Priors
  beta_m <- beta0
  beta_s <- 0.01
  alpha_m <- starts$alpha0
  alpha_s <- 0.05
  delta_m <- starts$delta0
  delta_s <- 0.05
  xi_s <- 1.0
  # Sum-to-zero constraint on attack and defence strengths
  alpha0 <- info$team_att
  delta0 <- info$team_def
  # Initialise linear predictor
  mu10 <- with(dat, exp(beta01 + alpha0[homeID] - delta0[awayID]))
  mu20 <- with(dat, exp(beta02 + alpha0[awayID] - delta0[homeID]))
  if (max(mu10) >= 100){mu10 <- pmin(100, mu10)}
  if (min(mu10) < 0.01){mu10 <- pmax(0.01, mu10)}
  if (max(mu20) >= 100){mu20 <- pmin(100, mu20)}
  if (min(mu20) < 0.01){mu20 <- pmax(0.01, mu20)}
  nu0 <- info$true_nu
  if (skel){
    nu0 <- rep(1, nteams)
  }
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
    # beta update (one-at-a-time)
    beta1 <- rnorm(2, beta0, 0.03)
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
      #mu20 <- mu21 - NEEDED?
      current <- prop
    }
    # nu update - on log scale (xi)
    if (!skel){
    xi1 <- rnorm(nteams, xi0, 0.25)
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
       model = fit_att_def_skel4, data = dat)
}

#### Component-wise model sim study - underdispersion (same for each team) ####
res_const_nu_under <- NULL
true_vals_under <- NULL
print("Underdispersed case")
for (i in 1:50){
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 10, nteam = 20, 
                                  sdteam = 0.25, nu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel5(info = sim_info, 1000)
  res_const_nu_under <- rbind(res_const_nu_under, 
                              colMeans(cbind(run$betas, run$alphas, 
                                             run$deltas, exp(run$xis))))
  true_vals_under <- rbind(true_vals_under, c(sim_info$beta1, sim_info$beta2,
                                              sim_info$team_att, sim_info$team_def, 
                                              sim_info$true_nu))
}

#### Component-wise model sim study - overdispersion (same for each team) ####
# Need to update range
res_const_nu_over <- NULL
true_vals_over <- NULL
for (i in 1:25){
  print("Overdispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl(nseasons = 10, nteam = 20, 
                                  sdteam = 0.25, nu = 0.5, 
                                  beta1 = 2, beta2 = 1.6)
  ii_u <- sim_info$dat$goaldif > 30
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -30
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel5(info = sim_info, 1000)
  res_const_nu_over <- rbind(res_const_nu_over, 
                             colMeans(cbind(run$betas, run$alphas, 
                                            run$deltas,exp(run$xis))))
  true_vals_over <- rbind(true_vals_over, c(sim_info$beta1, sim_info$beta2,
                                            sim_info$team_att, sim_info$team_def, 
                                            sim_info$true_nu))
}


#### Component-wise model sim study - bidispersion ####
res_biv_nu <- NULL
true_vals_biv <- NULL
for (i in 1:50){
  print("Bidispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 100, nteam = 6, 
                                   sdteam = 0.25, nul = 0.5, nuu = 2)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel4(info = sim_info, 1000)
  res_biv_nu <- rbind(res_biv_nu, 
                      colMeans(cbind(run$betas, run$alphas, 
                                     run$deltas,exp(run$xis))))
  true_vals_biv <- rbind(true_vals_biv, c(sim_info$beta1, sim_info$beta2,
                                          sim_info$team_att, sim_info$team_def, 
                                          sim_info$true_nu))
}


res_biv_nu <- NULL
true_vals_biv <- NULL
for (i in 1:50){
  print("Bidispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 10, nteam = 20, 
                                   sdteam = 0.25, nul = 0.5, nuu = 2)
  ii_u <- sim_info$dat$goaldif > 20
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -20
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel5(info = sim_info, 1000)
  res_biv_nu <- rbind(res_biv_nu, 
                      colMeans(cbind(run$betas, run$alphas, 
                                     run$deltas,exp(run$xis))))
  true_vals_biv <- rbind(true_vals_biv, c(sim_info$beta1, sim_info$beta2,
                                          sim_info$team_att, sim_info$team_def, 
                                          sim_info$true_nu))
}



#### Component-wise model sim study - underdispersion (varies by team) ####
res_nu_under <- NULL
true_vals_under <- NULL
print("Underdispersed case")
for (i in 1:50){
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 10, nteam = 20, 
                                   nul = 1.5, nuu = 2.5)
  ii_u <- sim_info$dat$goaldif > 10
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -10
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel5(info = sim_info, 1000)
  res_nu_under <- rbind(res_nu_under, 
                        colMeans(cbind(run$betas, run$alphas, 
                                             run$deltas, exp(run$xis))))
  true_vals_under <- rbind(true_vals_under, c(sim_info$beta1, sim_info$beta2,
                                              sim_info$team_att, sim_info$team_def, 
                                              sim_info$true_nu))
}

#### Component-wise model sim study - overdispersion (varies by team) ####
# Need to update range
res_nu_over <- NULL
true_vals_over <- NULL
for (i in 1:50){
  print("Overdispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 10, nteam = 20, 
                                   nul = 0.25, nuu = 0.75, 
                                  beta1 = 2, beta2 = 1.6)
  ii_u <- sim_info$dat$goaldif > 30
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -30
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel5(info = sim_info, 1000)
  res_nu_over <- rbind(res_nu_over, 
                             colMeans(cbind(run$betas, run$alphas, 
                                            run$deltas,exp(run$xis))))
  true_vals_over <- rbind(true_vals_over, c(sim_info$beta1, sim_info$beta2,
                                            sim_info$team_att, sim_info$team_def, 
                                            sim_info$true_nu))
}

#### Component-wise model sim study - bidispersion ####
res_biv_nu <- NULL
true_vals_biv <- NULL
for (i in 1:50){
  print("Bidispersed case")
  print(i)
  sim_info <- sim_cmp_skellam_epl2(nseasons = 10, nteam = 20, 
                                   nul = 0.5, nuu = 2)
  ii_u <- sim_info$dat$goaldif > 25
  sim_info$dat <- sim_info$dat[!ii_u, ]
  ii_l <- sim_info$dat$goaldif < -25
  sim_info$dat <- sim_info$dat[!ii_l, ]
  run <- fit_att_def_skel5(info = sim_info, 1000)
  res_biv_nu <- rbind(res_biv_nu, 
                      colMeans(cbind(run$betas, run$alphas, 
                                     run$deltas,exp(run$xis))))
  true_vals_biv <- rbind(true_vals_biv, c(sim_info$beta1, sim_info$beta2,
                                          sim_info$team_att, sim_info$team_def, 
                                          sim_info$true_nu))
}




#### Plots ####
df_under <- data.frame(y = as.numeric(res_const_nu_under[, 43:62]), 
                       id = as.factor(rep(1:20, each = 50)))
plot_under <- ggplot(df_under, aes(x = id, y = y)) + 
  geom_boxplot(outlier.colour = "gray80", outlier.size = 0.5,
               colour = "black", fill = "gray50", show.legend = FALSE)
plot_under + geom_hline(yintercept = 2, linewidth = 1.2,
                        linetype = "dashed", colour = "darkgreen") + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", fill = "white"))
boxplot(res_const_nu_under[, 43:62], pch = 16, cex = 0.3, 
        xlab = "Team id", ylab = "Estimate")
# True values
points(1:20, rep(2, 20), pch = 19, col = "darkred")

#### Starting values ####
# log(mean) is biased upwards once team effects grow
# First attack and defensive strengths
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
(beta1hat <- log(mean(data$vgoal/vgoal_hat)))
list(beta0 = c(beta0hat, beta1hat), alpha0 = att_inits,
     delta0 = def_inits)
}

