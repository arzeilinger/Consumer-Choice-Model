#### ANALYSIS OF 2-CHOICE CONSUMER MOVEMENT DATA USING MLE
# Code includes maximum likelihood estimation of 4 model variants, model selection,
# estimation of variance using normal approximation method or jackknife method,
# model averaging, and calculation of equilibrial probabilities, Pj*

# Load required libraries
library(optimx)

## Input functions for P1 and P2 equations, 
## negative log likelihood (NLL) functions, 
## gradient functions, and AIC function
source("source_functions.R")


#### MLE Optimization Functions 
## Choice data must include the columns:
# t = time points of observations, must be greater than 0
# n1 = number of consumers on choice 1 for each time point
# n2 = number of consumers on choice 2 for each time point
# n3 = number of consumers in neutral space for each time point
# N = total sample size of experiment

## Example data set 
# Data for trials with Euschistus servus from experiment described in Zeilinger et al. (2014)
dat <- structure(list(t = c(0.167, 0.5, 1, 12, 18, 24, 36), n1 = c(1L, 
       2L, 2L, 2L, 4L, 4L, 3L), n2 = c(6L, 5L, 7L, 10L, 11L, 11L, 12L
       ), n3 = c(8L, 8L, 6L, 3L, 0L, 0L, 0L), N = c(15L, 15L, 15L, 15L, 
       15L, 15L, 15L)), .Names = c("t", "n1", "n2", "n3", "N"), class = "data.frame", row.names = c(NA, 
       7L))

## Fixed Model Optimization
op.fixed <- optimx(par = rep(0.01, 2), # define starting parameter values
                   fn = NLL.fixed, y = dat, # define NLL function and data set
                   gr = gr.fixed, # define gradient function
                   method = c("spg"), # define optimization method
                   lower = rep(0.0001, 2), # set lower inequality constraint
                   upper = rep(10, 2), # set upper inequality constraint
                   control = list(all.methods=FALSE, 
                                  ftol = 1e-20, maxit = 10000)) # set convergence tolerance and maximum number of iterations

## Free Attraction Model Optimization
op.p.choice <- optimx(par = rep(0.01, 3), 
                      fn = NLL.p.choice, y = dat, 
                      gr = gr.p.choice,
                      method = c("spg"), 
                      lower = rep(0.0001, 3),
                      upper = rep(10, 3),
                      control = list(all.methods=FALSE,
                                     ftol = 1e-20, maxit = 10000))

## Free Leaving Model Optimization
op.mu.choice <- optimx(par = rep(0.01, 3), 
                       fn = NLL.mu.choice, y = dat, 
                       gr = gr.mu.choice,
                       method = c("spg"), 
                       lower = rep(0.0001, 3),
                       upper = rep(10, 3),
                       control = list(all.methods=FALSE, 
                                      ftol = 1e-20, maxit = 10000))

## Free Model Optimization
op.choice <- optimx(par = rep(0.01, 4), 
                    fn = NLL.choice, y = dat, 
                    gr = gr.choice,
                    method = c("spg"),
                    lower = rep(0.0001, 4),
                    upper = rep(10, 4),
                    control = list(all.methods=FALSE, 
                                   ftol = 1e-20, maxit = 10000))


## AICc calculations for model variants
## Using aicc() function from 'source_functions.R' file
aic.fixed <- aicc(k = coef(op.fixed), L = op.fixed$value, n = dat$N[1])
aic.p.choice <- aicc(k = coef(op.p.choice), L = op.p.choice$value, n = dat$N[1])
aic.mu.choice <- aicc(k = coef(op.mu.choice), L = op.mu.choice$value, n = dat$N[1])
aic.choice <- aicc(k = coef(op.choice), L = op.choice$value, n = dat$N[1])
aic.comp <- data.frame("model" = c("fixed", "p.choice", "mu.choice", "choice"),
                       "AIC" = c(aic.fixed, aic.p.choice, aic.mu.choice, aic.choice))
aic.comp$dAIC <- abs(min(aic.comp$AIC) - aic.comp$AIC)
# Sort AIC results from lowest to highest dAIC
aic.ord <- aic.comp[order(aic.comp$dAIC),]
print(aic.ord)


#### Construct data frame with all parameter estimates and NLL values
# Extract parameter estimates
pars.fixed <- c(coef(op.fixed)[1], coef(op.fixed)[1], coef(op.fixed)[2], coef(op.fixed)[2])
pars.p.choice <- c(coef(op.p.choice)[1:2], coef(op.p.choice)[3], coef(op.p.choice)[3])
pars.mu.choice <- c(coef(op.mu.choice)[1], coef(op.mu.choice)[1], coef(op.mu.choice)[2:3])
pars.choice <- c(coef(op.choice))
# Extract NLL values
nll <- c(op.fixed$value, op.p.choice$value, op.mu.choice$value, op.choice$value)
# Add parameter estimates and NLL to AIC table
pars.hat <- data.frame(rbind(pars.fixed, pars.p.choice, pars.mu.choice, pars.choice), row.names = NULL)
names(pars.hat) <- c("p1", "p2", "mu1", "mu2")
mle.tab <- cbind(aic.comp, pars.hat, nll)


#########################################################################
#### Calculating Confidence Intervals using Quadratic approximation
#########################################################################
# Method based on description in Bolker (2008) pgs. 196 - 201.
# Method can only be used if MLE is close to global maximum.
# If parameter estimates are on or near an inequality constraint (e.g., ~0.0001),
# then the jackkife method must be used to estimate variance.

# Calculate SE and correlation matrices for each model
op.list <- list(op.fixed, op.p.choice, op.mu.choice, op.choice) # Place all optimx outputs in a list
se.list <- list(length(4)) # Create empty list for standard errors
cormat.list <- list(length(4)) # Create empty list for correlation matrices
for(i in 1:length(op.list)){
  mod.i <- op.list[[i]]
  atr.i <- attr(mod.i, "details") # Extract attributes from each model output
  hess.i <- atr.i[1,3][[1]] # Extract Hessian matrix from each output
  vcov.i <- solve(hess.i) # Calculate the Variance-Covariance Matrix
  cormat.list[[i]] <- cov2cor(vcov.i) # Calculate the the Corelation Matrix
  se.i <- sqrt(diag(vcov.i)) # Extract Variances and calculate standard errors
  se.list[[i]] <- se.i
}
# Extract variances for each model
se.fixed <- se.list[[1]]
se.p.choice <- se.list[[2]]
se.mu.choice <- se.list[[3]]
se.choice <- se.list[[4]]
# Create data frame of variances, with rows for each model and columns for each parameter
se.tab <- data.frame(rbind(c(se.fixed[1], se.fixed[1], se.fixed[2], se.fixed[2]),
                           c(se.p.choice[1], se.p.choice[2], se.p.choice[3], se.p.choice[3]),
                           c(se.mu.choice[1], se.mu.choice[1], se.mu.choice[2], se.mu.choice[3]),
                           c(se.choice)))
names(se.tab) <- c("se.p1", "se.p2", "se.mu1", "se.mu2") 
mle.tab <- cbind(mle.tab, se.tab) # Combine with parameter estimates and AIC results


###################################################################################
#### Jackknikfe standard error estimation
##########################################################################
# Choice data should be in "wide" format, with each row representing 1 trial,
# and the observation time points in separate columns

# Example data set from experiment described in Zeilinger et al. (2014)
dat.obs <- structure(list(hr_0.16667 = structure(c(2L, 2L, 3L, 1L, 3L, 3L, 
  2L, 2L, 2L, 3L, 2L, 3L, 3L, 2L, 2L), .Label = c("D", "NC", "U"
  ), class = "factor"), hr_0.5 = structure(c(2L, 3L, 3L, 2L, 3L, 
  3L, 2L, 2L, 2L, 1L, 2L, 3L, 2L, 2L, 1L), .Label = c("D", "NC", 
  "U"), class = "factor"), hr_1 = structure(c(3L, 4L, 4L, 3L, 4L, 
  4L, 3L, 3L, 3L, 2L, 3L, 4L, 4L, 4L, 2L), .Label = c("", "D", 
  "NC", "U"), class = "factor"), hr_12 = structure(c(4L, 4L, 4L, 
  3L, 4L, 4L, 4L, 4L, 3L, 2L, 2L, 4L, 4L, 4L, 3L), .Label = c("", 
  "D", "NC", "U"), class = "factor"), hr_18 = structure(c(4L, 4L, 
  4L, 2L, 4L, 4L, 4L, 4L, 4L, 2L, 2L, 4L, 4L, 4L, 2L), .Label = c("", 
  "D", "NC", "U"), class = "factor"), hr_24 = structure(c(4L, 4L, 
  4L, 2L, 4L, 4L, 4L, 4L, 4L, 2L, 2L, 4L, 4L, 4L, 2L), .Label = c("", 
  "D", "NC", "U"), class = "factor"), hr_36 = structure(c(4L, 4L, 
  4L, 2L, 4L, 4L, 4L, 4L, 4L, 2L, 2L, 4L, 4L, 4L, 4L), .Label = c("", 
  "D", "NC", "U"), class = "factor")), .Names = c("hr_0.16667", 
  "hr_0.5", "hr_1", "hr_12", "hr_18", "hr_24", "hr_36"), row.names = c(9L, 
  13L, 15L, 16L, 17L, 18L, 20L, 21L, 23L, 24L, 25L, 26L, 28L, 37L, 
  44L), class = "data.frame")

library(bootstrap)

#### Function to send to jackknife() function
# jkfunc() transforms "wide" observation data set (one trial per row)
# into data set for MLE then submits it to optimx function 
jkfunc <- function(x, dat.obs, par, model){
  # Reshape data set into "long" format
  datf <- reshape(dat.obs[x,], direction = "long", timevar = "time",
               times = names(dat.obs[,1:7]), varying = 1:7, v.names = "choice")
  # Create vector of observation time points, assuming columns of data set are in the form "hr_[time]"
  tvf <- as.numeric(unlist(strsplit(unique(datf$time), split = "_")))
  tvf <- tvf[!is.na(tvf)]
  datf$time <- factor(datf$time)
  # Function to transform data into table of frequencies of choices observed,
  # the format needed for maximum likelihood estimation
  ccfunc <- function(z){
    cc <- table(datf[datf$time == z,]$choice)[1:3]
    return(cc)
  }
  ccdat <- sapply(levels(datf$time), ccfunc)
  Cdat <- data.frame(tvf, t(ccdat))
  Cdat <- Cdat[order(Cdat$tvf),]
  names(Cdat) <- c("t", "n1", "n3", "n2")
  Cdat$N <- rowSums(Cdat[,-1])
  # Series of if/else statements to select which model variant to use and cycle through all of them
  if(model == "fixed")
    opf <- optimx(par = c(0.01, 0.01), 
                  fn = NLL.fixed, y = Cdat, 
                  gr = gr.fixed,
                  method = c("spg"), 
                  lower = rep(0.0001, 2),
                  upper = rep(10, 2),
                  control = list(all.methods=FALSE, 
                                 ftol = 1e-20, maxit = 10000)) 
  else
    if(model == "p.choice")
      opf <- optimx(par = rep(0.01, 3), 
                    fn = NLL.p.choice, y = Cdat, 
                    gr = gr.p.choice,
                    method = c("spg"), 
                    lower = rep(0.0001, 3),
                    upper = rep(10, 3),
                    control = list(all.methods=FALSE,
                                   ftol = 1e-20, maxit = 10000))
    else
      if(model == "mu.choice")
        opf <- optimx(par = rep(0.01, 3), 
                      fn = NLL.mu.choice, y = Cdat, 
                      gr = gr.mu.choice,
                      method = c("spg"), 
                      lower = rep(0.0001, 3),
                      upper = rep(10, 3),
                      control = list(all.methods=FALSE, 
                                     ftol = 1e-20, maxit = 10000))
  else
    opf <- optimx(par = rep(0.01, 4), 
                  fn = NLL.choice, y = Cdat, 
                  gr = gr.choice,
                  method = c("spg"),
                  lower = rep(0.0001, 4),
                  upper = rep(10, 4),
                  control = list(all.methods=FALSE, 
                                 ftol = 1e-20, maxit = 10000))
  parf <- coef(opf)[par] # Extract only the parameter estimates from optimx output
  return(parf)
}

#### jk.allmodels() calculates jackknife SE estimates 
#### for each parameter in each model variant using jkfunc()
jk.allmodels <- function(m){
  model <- m
  # Series of if/else statements to select correct model variant
  if(model == "fixed")
    p <- c(1,1,2,2)
  else
    if(model == "p.choice")
      p <- c(1,2,3,3)
  else
    if(model == "mu.choice")
      p <- c(1,1,2,3)
  else
    p <- 1:4
  # Function to apply jackknife method to each parameter in each model variant
  jk.allpars <- function(p){
    jk <- jackknife(x = 1:nrow(dat.obs), theta = jkfunc, dat = dat.obs, par = p, model = model)
    jk.par <- c(jk$jack.se)
    return(jk.par)
  }
  allpars <- sapply(p, jk.allpars)
  return(c(model, allpars))
}

models <- c("fixed", "p.choice", "mu.choice", "choice") # Vector of model names
jk.se <- sapply(models, jk.allmodels, simplify = TRUE) # Recursively run jk.allmodels()

#### Combine jackknife SD estimates and mle.tab data frame
## mle.tab includes AIC, NLL, and parameter estimates
jkse.tab <- data.frame(t(jk.se), row.names = NULL)
names(jkse.tab) <- c("model", "se.p1", "se.p2", "se.mu1", "se.mu2") 
mle.tab <- merge(x = mle.tab, y = jkse.tab, by = "model", sort = FALSE)
mle.tab$se.p1 <- as.numeric(levels(mle.tab$se.p1))[mle.tab$se.p1]
mle.tab$se.p2 <- as.numeric(levels(mle.tab$se.p2))[mle.tab$se.p2]
mle.tab$se.mu1 <- as.numeric(levels(mle.tab$se.mu1))[mle.tab$se.mu1]
mle.tab$se.mu2 <- as.numeric(levels(mle.tab$se.mu2))[mle.tab$se.mu2]


################################################################################
#### Model-averaged parameter estimates and variances ##########################
################################################################################
# Note: Averaging only models with some information, i.e., those with dAICc < 7

# All models likelihoods and weights
amd <- mle.tab$dAIC
alv <- exp(-amd/2)
tml <- sum(alv) # Total marginal likelihood
awv <- alv/tml

# Good model likelihoods and weights
gmd <- mle.tab[mle.tab$dAIC < 7,] # select only good models
gmd$gml <- exp(-gmd$dAIC/2) # Relative or marginal likelihoods for models with some information (within)
gmd$wts <- gmd$gml/tml # Calculate relative weights for each model

# Calculate averaged parameter estimates and unconditional variances for each parameter
# p1
p1av <- sum(gmd$wts*gmd$p1)
p1var <- gmd$se.p1^2
p1var.av <- sum(gmd$wts*(p1var + (gmd$p1 - p1av)^2))
# p2
p2av <- sum(gmd$wts*gmd$p2)
p2var <- gmd$se.p2^2
p2var.av <- sum(gmd$wts*(p2var + (gmd$p2 - p2av)^2))
# mu1
mu1av <- sum(gmd$wts*gmd$mu1)
mu1var <- gmd$se.mu1^2
mu1var.av <- sum(gmd$wts*(mu1var + (gmd$mu1 - mu1av)^2))
# mu2
mu2av <- sum(gmd$wts*gmd$mu2)
mu2var <- gmd$se.mu2^2
mu2var.av <- sum(gmd$wts*(mu2var + (gmd$mu2 - mu2av)^2))

# Table of model-averaged parameter estimates and variances
resav <- data.frame(params = c("p1", "p2", "mu1", "mu2"),
                    estimates = c(p1av, p2av, mu1av, mu2av), 
                     vars = c(p1var.av, p2var.av, mu1var.av, mu2var.av))
resav$se <- sqrt(resav$vars)




#####################################################################################
#### Equilibrial probabilities for stink bug data
#####################################################################################

## Functions to calculate equilibrium probabilities, P1* and P2* from parameter estimates
# Equilibrium equations from Equation A12 in Appendix A of Zeilinger et al. (2014)
P1eq.func <- function(p1, p2, mu1, mu2){
  P1eq <- (p1*mu2)/(mu1*p2 + mu2*p1 + mu1*mu2)
  return(P1eq)
}
P2eq.func <- function(p1, p2, mu1, mu2){
  P2eq <- (p2*mu1)/(mu1*p2 + mu2*p1 + mu1*mu2)
  return(P2eq)
}

p1hat <- resav$estimates[1]; p2hat <- resav$estimates[2]
mu1hat <- resav$estimates[3]; mu2hat <- resav$estimates[4]

P1.star <- P1eq.func(p1hat, p2hat, mu1hat, mu2hat) # P1 equilibrium
P2.star <- P2eq.func(p1hat, p2hat, mu1hat, mu2hat) # P2 equilibrium



