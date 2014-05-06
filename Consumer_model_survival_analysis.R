#### SURVIVAL ANALYSIS OF CHOICE DATA
## To determine constant attraction and leaving rates

## Load required libraries
library(survival); library(lattice)

## Example data set
dat <- structure(list(time0 = c(0.066666667, 0.016666667, 0.5, 0.5, 
  0.016666667, 0.15, 0.5, 12, 0.5, 12, 0.033333333, 0.5, 12, 12, 
  12, 18, 0.116666667, 1, 0.133333333), moved0 = c(1L, 1L, 1L, 
  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
  ), choice1 = structure(c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 
  2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L), .Label = c("D", "U"), class = "factor"), 
  moved1 = c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
  0L, 0L, 0L, 0L, 0L, 0L, 1L), time1 = c(0.266666667, 0.316666667, 
  6.166666667, 35.66666667, 35.98333333, 35.85, 35.66666667, 
  29.5, 35.66666667, 29.5, 35.96666667, 35.66666667, 29.5, 
  29.5, 29.5, 21, 35.88333333, 35.25, 0.2)), .Names = c("time0", 
  "moved0", "choice1", "moved1", "time1"), class = "data.frame", row.names = c(NA, 
  19L))

## Settling time
# Obtaining Kaplan-Meier estimates of (S(t)) and time by plant choice for attraction rates
# choice1 is a variable that records the first choice of each individual consumer
# time0 is the time that each consumer moved from the neutral space to its first choice
# moved0 is a binary variable recording whether each consumer moved (1) or did not move (0)
setl.fit = survfit(Surv(time0, moved0) ~ choice1, data = dat) # Calculate K-M survival
sumsetl <- summary(setl.fit)
# Extract time, S(t), and choice from the survfit object
setl.dat <- data.frame(cbind("time" = sumsetl$time,
                             "surv" = sumsetl$surv,
                             "choice" = as.character(sumsetl$strata)))
setl.dat$time <- as.numeric(levels(setl.dat$time))[setl.dat$time]
setl.dat$surv <- as.numeric(levels(setl.dat$surv))[setl.dat$surv]
setl.dat$rate <- rep("Attraction", nrow(setl.dat))

## Leaving time
# Obtaining Kaplan-Meier estimates of (S(t)) and time by plant choice for leaving rates
# time1 is the time that each consumer its first choice
# moved1 records whether each consumer left the first choice (1) or did not leave (0)
leav.fit = survfit(Surv(time1, moved1) ~ choice1, data = dat) # Calculate K-M survival
sumleav <- summary(leav.fit)
# Extract time, S(t), and choice from the survfit object
leav.dat <- data.frame(cbind("time" = sumleav$time,
                             "surv" = sumleav$surv,
                             "choice" = as.character(sumleav$strata)))
leav.dat$time <- as.numeric(levels(leav.dat$time))[leav.dat$time]
leav.dat$surv <- as.numeric(levels(leav.dat$surv))[leav.dat$surv]
leav.dat$rate <- rep("Leaving", nrow(leav.dat))


#### Combine settling data and leaving data
survdat <- rbind(setl.dat, leav.dat)
survdat <- survdat[survdat$surv > 0,] # Remove K-M estimates of 0
survdat$log.time <- log(survdat$time) 
survdat$log.surv <- log(-log(survdat$surv))


#### Plot K-M survival estimates vs. time
xyplot(log.surv ~ log.time|choice*rate, data = survdat, 
       scales = list(alternating = FALSE, tck = c(1, 0), cex = 1.1),
       xlab = list("ln(time)", cex = 1.3),
       ylab = list("ln{-ln[S(t)]}", cex = 1.3), aspect = 1,
       layout = c(2,2), as.table = TRUE, strip = TRUE,
       pch = 16, col = "black", type = c('p', 'r'))

