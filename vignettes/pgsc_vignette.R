## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------

library(pgsc)

### Parameters ###
NN <- 20                    # Number of units
TT <- 100                   # Number of periods
MM <- 2                     # Number of treatments
RR <- 3                     # Number of covariates
SS <- 2                     # Number of unit FEs
QQ <- 3                     # Number of time FEs
b <- c(1,2)                 # Treatment coefficients
sq <- matrix( c( 1, .2, .3, -.9, -.1, .2 ), nrow=SS, ncol=QQ )
                            # Weighting matrix on time-varying factors
p <- t( matrix(c( -1, 0, .2, .5, .2, 0), nrow=MM, ncol=RR ) )  
                            # The covariance of X and D
r <- .1 * c( .5, 1, 2)      # Coefficient on observed covariates
sig <- .2                   # Noise sd
sig.y <- 5                  # Unit FEs
sig.t <- 4                  # Time FE noise

### Data
set.seed(42)
fes <- matrix( rnorm(NN * SS, 0, sig.y), NN, SS )    # Unit fixed effects
tfes <- matrix( rnorm(TT * QQ, 0, sig.t ), TT, QQ )  # Time fixed effects
X <- array( rnorm(NN*RR*TT), c( NN, RR, TT ))        # Covariates
D <- array(NA, dim=c( NN, MM, TT ))
D[] <- apply(X, 3, function(x) x%*%p) 
D <- D + array( rnorm(NN*MM*TT), c( NN, MM, TT ))               # Treatments, correlated with X
Y <- sapply( 1:TT, function(i) D[,,i] %*% b + X[,,i] %*% r ) +  # Treatment & covariates
  fes[,1] + rep(1,NN) %*% t(tfes[,1]) +                         # FEs and TFE
  fes %*% sq %*% t(tfes) +                                      # Time-unit interaction
  rnorm( NN*TT, 0, sig )                                        # Noise
dta <- data.frame( n=state.abb[1:NN], t=rep(1:TT,each=NN), y=c(Y),
                   do.call(rbind,lapply(1:TT,function(i)D[,,i])),
                   do.call(rbind,lapply(1:TT,function(i)X[,,i])) )
names(dta) <- c('n','t','y', paste0('D',1:MM), paste0('X',1:RR) )
    # Bind into a data frame

## ------------------------------------------------------------------------
### Panel regression
library(plm)
pan <- plm( y ~ D1 + D2 + X1 + X2 + X3, dta, effect = 'twoways', index = c('n','t'))
summary(pan)

