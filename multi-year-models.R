## R script for maximum likelihood estimation and model selection
## Load packages
require(bbmle); require(optimx)

## Example data set, for the Euschistus servus-Helicoverpa zea combination
dat <- structure(list(year = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("2009", "2010"), class = "factor"), 
    t = c(0.167, 0.5, 1, 12, 18, 24, 36, 0.167, 0.5, 1, 12, 18, 
    24, 36), n1 = c(1L, 2L, 2L, 2L, 4L, 4L, 3L, 3L, 4L, 6L, 7L, 
    7L, 4L, 4L), n2 = c(6L, 5L, 7L, 10L, 11L, 11L, 12L, 0L, 0L, 
    1L, 8L, 7L, 7L, 10L), n3 = c(8L, 8L, 6L, 3L, 0L, 0L, 0L, 
    12L, 11L, 8L, 0L, 1L, 4L, 1L), N = c(15L, 15L, 15L, 15L, 
    15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L)), .Names = c("year", 
"t", "n1", "n2", "n3", "N"), class = "data.frame", row.names = c(1L, 
2L, 3L, 4L, 5L, 6L, 7L, 22L, 23L, 24L, 25L, 26L, 27L, 28L))


###############################################################################
#### Probability models as functions

P1.func <- function(p1, p2, mu1, mu2, tau, N, n1m1, n2m1){
  P1f <- ((1/2)*p1*((2*mu2)/(mu2*p1 + mu1*(mu2 + p2)) + 
                      (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                               4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                         (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                       n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))/
                      exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                              4*(mu2*p1 + mu1*(mu2 + p2))))*tau)/sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                                                        4*(mu2*p1 + mu1*(mu2 + p2))) + 
                      (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                4*(mu2*p1 + mu1*(mu2 + p2))))*tau)*
                         (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                                 4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                            (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                          n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
                      sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))
  ifelse(P1f < 0, 0, P1f) 
}

P2.func <- function(p1, p2, mu1, mu2, tau, N, n1m1, n2m1){
  P2f <- ((mu1*p2)/(mu2*p1 + mu1*(mu2 + p2)) - 
            ((mu1 - mu2 + p1 - p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
               (-2 + (mu2*(mu1 + mu2 + p1 + p2 - sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                        4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                  (1/(N*p1))*(mu1*n1m1 - mu2*n1m1 + n1m1*p1 + 2*n2m1*p1 - n1m1*p2 + 
                                n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
            exp((1/2)*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                  tau)/(4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))) - 
            (exp((1/2)*(-mu1 - mu2 - p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
                   tau)*(mu1 - mu2 + p1 - p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))*
               (2 - (mu2*(mu1 + mu2 + p1 + p2 + sqrt((mu1 + mu2 + p1 + p2)^2 - 
                                                       4*(mu2*p1 + mu1*(mu2 + p2)))))/(mu2*p1 + mu1*(mu2 + p2)) + 
                  (1/(N*p1))*((-mu1)*n1m1 + mu2*n1m1 - n1m1*p1 - 2*n2m1*p1 + n1m1*p2 + 
                                n1m1*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2))))))/
            (4*sqrt((mu1 + mu2 + p1 + p2)^2 - 4*(mu2*p1 + mu1*(mu2 + p2)))))
  ifelse(P2f < 0, 0, P2f) 
}


###########################################################################
#### Negative log likelihood functions for mle2()
# Note: for the NLL functions with year effects, first number in parameter names
# correpsonds to the choice, while second number corresponds to year

# M1: Fixed Model
# Each rate is estimated over choice and year
NLL.fixed <- function(p1, p2 = p1, mu1, mu2 = mu1) {
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M2a: p by choice
# Attraction rates depend on choice
NLL.p.choice <- function(p1, p2, mu1, mu2 = mu1) {
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M2b: mu by choice
# Leaving rates depend on choice
NLL.mu.choice <- function(p1, p2 = p1, mu1, mu2) {
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M2c: p by year
# Attraction rates depend on year
NLL.p.year <- function(p1.1, p2.1 = p1.1, p1.2, p2.2 = p1.2, mu1, mu2 = mu1) {
  # Index p1 and p2 by year
  p1 <- c(p1.1, p1.2)[dat$year]
  p2 <- c(p2.1, p2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M2d: mu by year
# Leaving rates depend on year
NLL.mu.year <- function(p1, p2 = p1, mu1.1, mu2.1 = mu1.1, mu1.2, mu2.2 = mu1.2) {
  mu1 <- c(mu1.1, mu1.2)[dat$year]
  mu2 <- c(mu2.1, mu2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M3a: p and mu by choice
# Attraction and leaving rates depend on choice
NLL.choice <- function(p1, p2, mu1, mu2) {
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M3b: p and mu by year
# Attraction and leaving rates depend on year
NLL.year <- function(p1.1, p2.1 = p1.1, p1.2, p2.2 = p1.2, mu1.1, mu2.1 = mu1.1, mu1.2, mu2.2 = mu1.2) {
  p1 <- c(p1.1, p1.2)[dat$year]
  p2 <- c(p2.1, p2.2)[dat$year]
  mu1 <- c(mu1.1, mu1.2)[dat$year]
  mu2 <- c(mu2.1, mu2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M4a: p by choice and year; mu by choice
# Attraction rates depend on choice and year; leaving rates depend on choice
NLL.pcy.muc <- function(p1.1, p2.1, p1.2, p2.2, mu1, mu2) {
  p1 <- c(p1.1, p1.2)[dat$year]
  p2 <- c(p2.1, p2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M4b: p by choice; mu by choice and year
# Attraction rates depend on choice; leaving rates depend on choice and year
NLL.pc.mucy <- function(p1, p2, mu1.1, mu2.1, mu1.2, mu2.2) {
  mu1 <- c(mu1.1, mu1.2)[dat$year]
  mu2 <- c(mu2.1, mu2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M4c: p by choice and year; mu by year
# Attraction rates depend on choice and year; leaving rates depend on year
NLL.pcy.muy <- function(p1.1, p2.1, p1.2, p2.2, mu1.1, mu2.1 = mu1.1, mu1.2, mu2.2 = mu1.2) {
  p1 <- c(p1.1, p1.2)[dat$year]
  p2 <- c(p2.1, p2.2)[dat$year]
  mu1 <- c(mu1.1, mu1.2)[dat$year]
  mu2 <- c(mu2.1, mu2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


# M4d: p by year; mu by choice and year
# Attraction rates depend on year; and leaving rates depend on choice and year
NLL.py.mucy <- function(p1.1, p2.1 = p1.1, p1.2, p2.2 = p1.2, mu1.1, mu2.1, mu1.2, mu2.2) {
  p1 <- c(p1.1, p1.2)[dat$year]
  p2 <- c(p2.1, p2.2)[dat$year]
  mu1 <- c(mu1.1, mu1.2)[dat$year]
  mu2 <- c(mu2.1, mu2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}

# M5: Free Model
# Attraction and leaving rates depend on choice and year
NLL.free <- function(p1.1, p2.1, p1.2, p2.2, mu1.1, mu2.1, mu1.2, mu2.2) {
  p1 <- c(p1.1, p1.2)[dat$year]
  p2 <- c(p2.1, p2.2)[dat$year]
  mu1 <- c(mu1.1, mu1.2)[dat$year]
  mu2 <- c(mu2.1, mu2.2)[dat$year]
  # Specify parameters
  p1f <- p1
  p2f <- p2
  mu1f <- mu1
  mu2f <- mu2
  y <- dat
  tf <- y$t # Specify vector of time points from data set 'y'
  lt <- length(tf)/2
  tm1f <- c(0, tf[1:(lt-1)], 0, tf[(lt+1):(length(tf)-1)]) # specify tm1 ("t minus 1") vector from the tf vector
  tauf <- tf - tm1f # calculate tau = t - tm1
  # specify n1, n2, and n3 from data set 'y'
  n1f <- y$n1 
  n2f <- y$n2
  n3f <- y$n3
  n.all <- cbind(n1f, n2f, n3f)
  # calculate n1m1 ("n1 minus 1") and n2m1 from n1 vector 
  n1m1f <- c(0, n1f[1:(lt-1)], 0, n1f[(lt+1):(length(tf)-1)])
  n2m1f <- c(0, n2f[1:(lt-1)], 0, n2f[(lt+1):(length(tf)-1)])
  Nf <- y$N # specify total sample size, N
  # Calculate P1, P2, P3
  P1 <- P1.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P2 <- P2.func(p1 = p1f, p2 = p2f, mu1 = mu1f, mu2 = mu2f, tau = tauf, N = Nf, n1m1 = n1m1f, n2m1 = n2m1f)
  P3 <- 1 - P1 - P2
  p.all <- cbind(P1, P2, P3)
  # Calculate negative log likelihood and correct for possible undeflow or overflow of calculations
  ifelse(all(p.all > 0 & p.all < 1),
    -sum(sapply(1:nrow(dat), function(z) dmultinom(x = n.all[z,], prob = p.all[z,], log = TRUE), simplify = TRUE)),
    1e07)
}


#######################################################################
#### MLE Optimization Functions 
op.func <- function(dat){
  # Function to optimize each model and return list of results
  ## M1: Fixed Optimization
  op.fixed <- mle2(NLL.fixed, data = dat,
                    start = list(p1 = 0.01, mu1 = 0.01),
                    optimizer = "optimx", method = "spg",
                    lower = rep(0.0001, 2), upper = rep(10, 2),
                    control = list(ftol = 1e-20, maxit = 10000))
  ## M2a: p by choice Optimization
  op.p.choice <- mle2(NLL.p.choice, data = dat,
                       start = list(p1 = 0.01, p2 = 0.01, mu1 = 0.01),
                       optimizer = "optimx", method = "spg",
                       lower = rep(0.0001, 3), upper = rep(10, 3),
                       control = list(ftol = 1e-20, maxit = 10000))
  ## M2b: mu by choice Optimization
  op.mu.choice <- mle2(NLL.mu.choice, data = dat,
                        start = list(p1 = 0.01, mu1 = 0.01, mu2 = 0.01),
                        optimizer = "optimx", method = "spg",
                        lower = rep(0.0001, 3), upper = rep(10, 3),
                        control = list(ftol = 1e-20, maxit = 10000))
  ## M2c: p by year Optimization
  op.p.year <- mle2(NLL.p.year, data = dat,
                    start = list(p1.1 = 0.01, p1.2 = 0.01, mu1 = 0.01),
                    optimizer = "optimx", method = "spg",
                    lower = rep(0.0001, 3), upper = rep(10, 3),
                    control = list(ftol = 1e-20, maxit = 10000),
                    parameters = list(p1~year, p2~year))
  ## M2d: mu by year Optimization
  op.mu.year <- mle2(NLL.mu.year, data = dat,
                      start = list(p1 = 0.01, mu1.1 = 0.01, mu1.2 = 0.01),
                      optimizer = "optimx", method = "spg",
                      lower = rep(0.0001, 3), upper = rep(10, 3),
                      control = list(ftol = 1e-20, maxit = 10000),
                      parameters = list(mu1~year, mu2~year))
  ## M3a: p and mu by choice Optimization
  op.choice <- mle2(NLL.choice, data = dat,
                     start = list(p1 = 0.01, p2 = 0.01, mu1 = 0.01, mu2 = 0.01), 
                     optimizer = "optimx", method = "spg",
                     lower = rep(0.0001, 4), upper = rep(10, 4),
                     control = list(ftol = 1e-20, maxit = 10000))
  ## M3b: p and mu by year Optimization
  op.year <- mle2(NLL.year, data = dat,
                   start = list(p1.1 = 0.01, p1.2 = 0.01, mu1.1 = 0.01, mu1.2 = 0.01),
                   optimizer = "optimx", method = "spg",
                   lower = rep(0.0001, 4), upper = rep(10, 4),
                   control = list(ftol = 1e-20, maxit = 10000),
                   parameters = list(p1~year, p2~year, mu1~year, mu2~year))
  ## M4a: p by year and choice; mu by choice Optimization
  start.pcy.muc = list(p1.1 = 0.01, p2.1 = 0.01, p1.2 = 0.01, p2.2 = 0.01, mu1 = 0.01, mu2 = 0.01)
  op.pcy.muc <- mle2(NLL.pcy.muc, data = dat,
                      start = start.pcy.muc,
                      optimizer = "optimx", method = "spg",
                      lower = rep(0.0001, 6), upper = rep(10, 6),
                      control = list(ftol = 1e-20, maxit = 10000),
                      parameters = list(p1~year, p2~year))
  ## M4b: p by choice; mu by year and choice Optimization
  start.pc.mucy = list(p1 = 0.01, p2 = 0.01, mu1.1 = 0.01, mu2.1 = 0.01, mu1.2 = 0.01, mu2.2 = 0.01)
  op.pc.mucy <- mle2(NLL.pc.mucy, data = dat,
                      start = start.pc.mucy,
                      optimizer = "optimx", method = "spg",
                      lower = rep(0.0001, 6), upper = rep(10, 6),
                      control = list(ftol = 1e-20, maxit = 10000),
                      parameters = list(mu1~year, mu2~year))
  ## M4c: p by choice and year; mu by year Optimization
  start.pcy.muy = list(p1.1 = 0.01, p2.1 = 0.01, p1.2 = 0.01, p2.2 = 0.01, mu1.1 = 0.01, mu1.2 = 0.01)
  op.pcy.muy <- mle2(NLL.pcy.muy, data = dat,
                      start = start.pcy.muy,
                      optimizer = "optimx", method = "spg",
                      lower = rep(0.0001, 6), upper = rep(10, 6),
                      control = list(ftol = 1e-20, maxit = 10000),
                      parameters = list(p1~year, p2~year, mu1~year, mu2~year))
  ## M4d: p by year; mu by year and choice Optimization
  start.py.mucy = list(p1.1 = 0.01, p1.2 = 0.01, mu1.1 = 0.01, mu1.2 = 0.01, mu2.1 = 0.01, mu2.2 = 0.01)
  op.py.mucy <- mle2(NLL.py.mucy, data = dat,
                      start = start.py.mucy,
                      optimizer = "optimx", method = "spg",
                      lower = rep(0.0001, 6), upper = rep(10, 6),
                      control = list(ftol = 1e-20, maxit = 10000),
                      parameters = list(p1~year, p2~year, mu1~year, mu2~year))
  ## M5: Free model Optimization
  #   start.free = list(p1.1 = 0.3, p1.2 = 0.3, p2.1 = 0.05, p2.2 = 0.56,
  #                     mu1.1 = 0.008, mu1.2 = 0.001, mu2.1 = 0.017, mu2.2 = 0.07)
  start.free = list(p1.1 = 0.01, p1.2 = 0.01, p2.1 = 0.01, p2.2 = 0.01,
                    mu1.1 = 0.01, mu1.2 = 0.01, mu2.1 = 0.01, mu2.2 = 0.01)
  op.free <- mle2(NLL.free, data = dat,
                   start = start.free,
                   optimizer = "optimx", method = "spg",
                   lower = rep(0.0001, 8),
                   control = list(ftol = 1e-20, maxit = 10000),
                   parameters = list(p1~year, p2~year, mu1~year, mu2~year))
  ### Model Selection
  modsel <- ICtab(op.fixed, 
                  op.p.choice, op.mu.choice, op.p.year, op.mu.year,
                  op.choice, op.year,
                  op.pcy.muc, op.pc.mucy, op.pcy.muy, op.py.mucy,
                  op.free,
                  type = "AICc", sort = TRUE, delta = TRUE, base = TRUE, 
                  nobs = mean(dat[dat$t == dat$t[1],]$N))
  op.list <- list(op.fixed,
                op.p.choice, op.mu.choice, op.p.year, op.mu.year,
                op.choice, op.year,
                op.pcy.muc, op.pc.mucy, op.pcy.muy, op.py.mucy,
                op.free)
  names(op.list) <- c("op.fixed",
                      "op.p.choice", "op.mu.choice", "op.p.year", "op.mu.year",
                      "op.choice", "op.year",
                      "op.pcy.muc", "op.pc.mucy", "op.pcy.muy", "op.py.mucy",
                      "op.free")
  return(list(op.list, modsel))
}

# Run optimization and AIC selection for the E. servus-H. zea dataset
resultsList <- op.func(dat = dat)
# List of MLE results for each model
opList <- resultsList[[1]]
# AIC table
modelSelect <- resultsList[[2]]

