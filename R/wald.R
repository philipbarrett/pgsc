#####################################################################################
# wald.R
# Provides functions for performing Wald tests on a restricted hypothesis
# Philip Barrett
# Washington, DC, 16sep2018
#####################################################################################



gsc.grad.cl <- function(W, Y, D, b, k.idx, sig.i=NULL ){
  # Returns the matrix of time-cluster-specific gradients for the kth treatment parameter
  NN <- dim(D)[1] ; MM <- length(b) ; TT <- dim(D)[3]   # Dimensions of the problem
  if( is.null(sig.i) ) sig.i <- rep( 1, NN )
  out <- matrix( NA, nrow=NN, ncol=TT)
  for( nn in 1:NN ){
    for( tt in 1:TT ){
      out[nn,tt] <- 1 / sig.i[nn] * ( ( D[nn, k.idx, tt] - sum( D[,k.idx,tt] * W[nn,]) ) * 
                                        ( Y[nn,tt] - sum( D[nn,,tt] * b ) - sum( W[nn,] * ( Y[,tt] - D[,,tt] %*% b ) ) ) )
    }
  }
  return(out)
}

gsc.mu <- function( grad.cl ){
  # Computes the mu residuals for each cluster after regressing on all prior clusters
  NN <- nrow(grad.cl) ; out <- 0 * grad.cl ; out[1,] <- grad.cl[1,]
  out[2,] <- lm( grad.cl[2,] ~ 0 + grad.cl[1,] )$residuals
  for(i in 3:NN){
    out[i,] <- lm( grad.cl[i,] ~ 0 + t(grad.cl[1:(i-1),]) )$residuals
  }
  return(out)
}

gsc.s.i <- function(grad.cl){
  # Computes the orthogonalized score from the gradients
  mu <- gsc.mu(grad.cl)
  return( s <- apply( mu, 1, mean ) )
}

gsc.s.i.k <- function( grad.cl ){
  # Computes the mean for each cluster after regressing on all prior clusters
  NN <- nrow(grad.cl) ; out <- rep( 0, NN ) ; out[1] <- mean(grad.cl[1,])
  out[2] <- lm( grad.cl[2,] ~ grad.cl[1,] )$coefficients[1]
  for(i in 3:NN){
    out[i] <- lm( grad.cl[i,] ~ t(grad.cl[1:(i-1),]) )$coefficients[1]
  }
  return(out)
}

gsc.wald <- function( s.i, W.d=NULL ){
  # Computes the Wald statistic for the weights W
  NN <- ncol(s.i)
  if(is.null(W.d)) W.d <- rep(1,NN)
  s.bar.0 <- apply(s.i, 1, mean)
  s.i <- s.i * ( rep( 1, nrow(s.i)) %*% t(W.d) )
  s.bar <- apply(s.i, 1, mean)
  sig.hat <- Reduce( '+', lapply( 1:NN, function(i) 
    ( s.i[,i] - s.bar ) %*% t( ( s.i[,i] - s.bar ) ) ) ) / (NN-1)
  # sig.hat <- sig.hat + 1e-09*diag(nrow(sig.hat))
  # Avoids numerical issues
  wald.s <- t(s.bar) %*% solve(sig.hat, s.bar)
  return(wald.s)
}

gsc.wald.boot <- function( s.i, n.it=1000 ){
  # Bootstraps the Wald statistic with appropriate weights
  NN <- ncol(s.i)
  W.d <- -1 + 2 * matrix( rbinom(NN*n.it, 1, .5), nrow=NN, ncol=n.it )
  v.s <- rep(NA, n.it)
  for( i in 1:n.it ) v.s[i] <- gsc.wald( s.i, W.d[,i])
  return(v.s)
}