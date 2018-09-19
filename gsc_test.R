rm(list=ls())
library(microbenchmark)
library(plm)
source('~/Dropbox/2017/research/AFG/conflict_revenue/code/gsc.R')
Rcpp::sourceCpp('src/gsc.cpp')

NN <- 34 # 35       # Number of units
TT <- 150           # Number of periods
MM <- 2             # Number of treatments
RR <- 3             # Number of covariates
SS <- 2             # Number of unit FEs
QQ <- 3             # Number of time FEs
set.seed(42)

### TO DO
# - Automate iteration/solution                     DONE
# - Compute GSC error                               DONE
# - Two-step estimaton weighted by GSC error        DONE
# - Generate an example where GSC beats panel       DONE
# - Reuse gsc.R code to estimate b_i for one-step   DONE
# - Apply to data: 
#     1) Figure out NaN handling
#     2) Compute estimates
# - Comparing weighed covariates in data
# - Do a monte carlo version of this to see if the estimation is ok.
# - Check what is going on with the second coefficient of b
# - Add asymptotic errors (need 100%, but only if estimates look good)
# - Check the GSC method when the interior assumption fails (?)
# - Charting for GSC
# - Proper computation of partial solution...
# - ...including derivatives

### Parameters ###
b <- c(1,2)         # Treatment coefficients
sq <- matrix( c( 1, .2, .3, -.9, -.1, .2 ), nrow=SS, ncol=QQ )
# Time-unit FE interactions for outcome
# v <- matrix( c( .2, .3, -.3, -.4, .7, 1, 
#                 -.4, -.6, -.5, -.7, .8, .2 ), nrow=SS, ncol=QQ*MM )
#                     # Time-unit FE interactions for treatment
r <- .1 * c( .5, 1, 2)
# Making this large will make estimation worse. Because the X are uncorrelated,
# then they are an extreme case of time-and-unit fixed effects.
p <- t( matrix(c( -1, 0, .2, .5, .2, 0), nrow=MM, ncol=RR ) )  
# Treatment affected by covariates
sig <- .2           # Noise sd
sig.y <- 5          # Unit FEs
sig.t <- 4          # Time FE noise

### Data ###
fes <- matrix( rnorm(NN * SS, 0, sig.y), NN, SS )    # Unit fixed effects
tfes <- matrix( rnorm(TT * QQ, 0, sig.t ), TT, QQ )  # Time fixed effects
X <- array( rnorm(NN*RR*TT), c( NN, RR, TT ))        # Covariates
D <- array(NA, dim=c( NN, MM, TT ))
D[] <- apply(X, 3, function(x) x%*%p) 
D <- D + array( rnorm(NN*MM*TT), c( NN, MM, TT ))               # Treatments
Y <- sapply( 1:TT, function(i) D[,,i] %*% b + X[,,i] %*% r ) +  # Treatment & covariates
  fes[,1] + rep(1,NN) %*% t(tfes[,1]) +                         # FEs and TFE
  fes %*% sq %*% t(tfes) +                                      # Time-unit interaction
  rnorm( NN*TT, 0, sig )                                        # Noise
# Outcomes.
plot( c(1,TT), range(Y), xlab='', ylab='', type='n' )
for( i in 1:NN ) lines( 1:TT, Y[i,], col=i )
# Plot
dta <- data.frame( n=1:NN, t=rep(1:TT,each=NN), y=c(Y),
                   do.call(rbind,lapply(1:TT,function(i)D[,,i])),
                   do.call(rbind,lapply(1:TT,function(i)X[,,i])) )
names(dta) <- c('n','t','y', paste0('D',1:MM), paste0('X',1:RR) )

### Panel regression ###
pan <- plm( y ~ D1 + D2 + X1 + X2 + X3, dta, effect = 'twoways', index = c('n','t'))
print(summary(pan))

### Test low-level evaluation functions ###
wt.init <- matrix( 1 / (NN-1), NN, NN-2 )
# Initial weights
gsc.init <- gsc_target( NN, TT, wt.init, Y, D, b, rep(1,NN) ) #, print_level = 2 )
n.obs <- gsc_target( NN, TT, wt.init, Y, D, b, rep(1,NN), return_n_obs = TRUE )
eps <- sapply( 1:TT, function(t) sapply( 1:NN, function(i) Y[i,t] - D[i,,t] %*% b ) )
gsc.init.check <- .5 * sum( ( eps - W_extract( wt.init, NN ) %*% eps ) ^ 2 ) / n.obs
gsc.grad.int.num <- gsc_target_grad_num( NN, TT, wt.init, Y, D, b, rep(1,NN) )
gsc.grad.int <- gsc_target_grad( NN, TT, wt.init, Y, D, b, rep(1,NN) )
cbind( gsc.grad.int.num, gsc.grad.int, gsc.grad.int.num - gsc.grad.int )
gsc.grad.int.num.2 <- gsc_target_grad_num( NN, TT, wt.init + 1, Y, D, b, rep(1,NN) )
gsc.grad.int.2 <- gsc_target_grad( NN, TT, wt.init + 1, Y, D, b, rep(1,NN) )
cbind( gsc.grad.int.num.2, gsc.grad.int.2, gsc.grad.int.num.2 - gsc.grad.int.2 )


### Naming of optimization functions ###
f <- function(x){
  gsc_target( NN, TT, matrix(x[1:(NN*(NN-2))], NN, NN-2), Y, D, tail(x,MM), rep(1,NN) )
}
f.grad <- function(x){
  gsc_target_grad( NN, TT, matrix(x[1:(NN*(NN-2))], NN, NN-2), Y, D, tail(x,MM), rep(1,NN) )
}
f.grad.b <- function(x){
  gsc_target_grad_b( NN, TT, matrix(x[1:(NN*(NN-2))], NN, NN-2), Y, D, tail(x,MM), rep(1,NN) )
}
f.grad.num <- function(x){
  gsc_target_grad_num( NN, TT, matrix(x[1:(NN*(NN-2))], NN, NN-2), Y, D, tail(x,MM), rep(1,NN) )
}
res <- microbenchmark(f(c(wt.init,b)), f.grad(c(wt.init,b)), f.grad.b(c(wt.init,b)), f.grad.num(c(wt.init,b)))
print(res)
# autoplot(res)
# print(tt <- system.time(f.grad.num(c(wt.init,b))))

### Compute iterative solution ###
sol.it <- gsc.iter(wt.init, Y, D, b)

### Compute direct optimization ###
ub <- c( rep(1,NN*(NN-2)), rep(Inf,MM) )
lb <- c( rep(0,NN*(NN-2)), -rep(Inf,MM) )
opts <- list( algorithm = "NLOPT_LD_SLSQP", maxeval=1000, print_level=0 ) #, 
# check_derivatives = TRUE ) #, check_derivatives_print = "all" )
x0 <- c(sol.it$wt,sol.it$b) 
# x0 <- c(wt.init, b ) 
# Initialize near the iterative solution
sol <- nloptr( x0, f, f.grad, ub=ub, lb=lb, opts=opts )
print( sol$iterations )
print( cbind( f.grad(x0) * (x0 > 1e-08),
              f.grad(sol$solution) * (sol$solution > 1e-08),
              f.grad.num(sol$solution) * (sol$solution > 1e-08) )  )
# Print the solution

### Extract the direct solution ###
wt.hat <- matrix( sol$solution[1:(NN*(NN-2))], NN, NN-2 )
W.hat <- W_extract( wt.hat, NN )
b.hat <- tail(sol$solution,MM)

### Compare the solutions ###
print( cbind( b.hat, sol.it$b, b.hat - sol.it$b ) )

### Compute the vector of unit-specific fits ###
sig.i.1 <- gsc_target_i( NN, TT, sol.it$wt, Y, D, 
                         matrix( sol.it$b, MM, NN ) )
print( c( mean(sig.i.1), gsc_target( NN, TT, sol.it$wt, Y, D, sol.it$b, rep(1,NN) ) ) )
# The average of the errors and the target function. Should be the same (close?)
sol.2.step.1 <- gsc.iter( wt.init, Y, D, b, sig.i.1)
# The two-step estimator without the best-fit weighting estimator
sol.i <- gsc.iter.i( wt.init, Y, D, b, print.level = 1 )
# Try again, this time fitting each unit alone
sig.i <- gsc_target_i( NN, TT, sol.i$wt, Y, D, t( sol.i$b ) )
# The resulting fit
barplot( t(cbind( 1/sig.i.1, 1/sig.i)), beside = TRUE )
sol.2.step <- gsc.iter( wt.init, Y, D, b, sig.i)
# The two-step estimator
sol.compare <- rbind( pan$coefficients[c('D1','D2')], sol.it$b, b.hat, sol.2.step.1$b, sol.2.step$b, b )
rownames(sol.compare) <- c( 'Panel FEs', 'GSC iter', 'GSC direct', 'Two-step GSC', 'Optimal two-step GSC', 'Truth')
print(sol.compare)

plot( c(1,NN), c(-.2,1), type='n' )
for( i in 1:NN ) lines( 1:NN, sol.2.step$W[i,], col=i)
abline(h=0)

### Testing ###
rest.val <- .85
g.1 <- function(b) b[1] - rest.val ; g.1.grad <- function(b) c(1,0)
sol.2.step.const <- gsc.iter( wt.init, Y, D, c(rest.val,b[2]), sig.i, g.i=g.1, g.i.grad=g.1.grad )
# The restricted solution
h.grad.1 <- gsc.grad.cl( sol.2.step.const$wt, Y, D, sol.2.step.const$b, 1, sig.i)
# The gradients of the objective function
mu.1 <- gsc.mu(h.grad.1)
s.i.1 <- gsc.s.i(h.grad.1)
S.1 <- gsc.wald(s.i.1)
v.S.1 <- gsc.wald.boot( s.i.1, 100000 )
print(mean( v.S.1 > S.1 ))

g.2 <- function(b) b[1] / ( 1 - b[2] ) + 1 ; g.2.grad <- function(b) c( 1 / ( 1 - b[2] ), b[1] / ( 1 - b[2] )^2 )
sol.2.step.const.2 <- gsc.iter( wt.init, Y, D, b, sig.i, g.i=g.2, g.i.grad=g.2.grad )
# The restricted solution
h.grad.2 <- gsc.grad.cl( sol.2.step.const.2$wt, Y, D, sol.2.step.const.2$b, 1, sig.i)
# The gradients of the objective function
mu.2 <- gsc.mu(h.grad.2)
s.i.2 <- gsc.s.i(h.grad.2)
S.2 <- gsc.wald(s.i.2)
v.S.2 <- gsc.wald.boot( s.i.2, 100000 )
print(mean( v.S.2 > S.2 ))






