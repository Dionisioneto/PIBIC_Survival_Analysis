
library(eha)

## Funcao que particiona os grids para censura intervalar

#--- Algumas fun??es importantes:

# Calculando as fun??es de taxa de falha e acumulada:

cal_ht_MEPP <- function( time.obs, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet){
  
  
  dens_MEPP <- alpha.par*(ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)^(alpha.par-1))*dpch(x =time.obs, cuts = grid.vet, levels = lambda.par)
  Cumu_MEPP <- ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)^alpha.par
  ht <- dens_MEPP/(1-Cumu_MEPP)
  
  return(ht)
}


cal_Ht_MEPP <- function(time.obs, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet){
  Cumu_MEPP <- ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)^alpha.par
  Ht  = -log(1-Cumu_MEPP)
  return(Ht)
}


loglikIC.mepp.int <- function(par, l, r, x.cov, grid){
  
  b <- length(grid)+1
  
  hazards = par[1:b] ## taxas de falha para os b intervalos
  alpha = par[b + 1] ## pparrametro de potencia
  
  n.cov = dim(x.cov)[2] 
  
  betas= par[(b + 2):(b + 1 + n.cov)]
  n.sample <- nrow(x.cov)
  
  delta <- ifelse(is.finite(r), 1, 0)
  
  lik <- rep(0, n.sample)
  
  s.l <- exp(-cal_Ht_MEPP(time.obs=l[delta==1], lambda.par=hazards,alpha.par=alpha, grid.vet=grid)*exp(x.cov[delta==1,] %*% as.matrix(betas)))

  
  s.r <- exp(-cal_Ht_MEPP(time.obs=r[delta==1], lambda.par=hazards,alpha.par=alpha, grid.vet=grid)*exp(x.cov[delta==1,] %*% as.matrix(betas)))
  
  lik[delta==1] <- s.l-s.r
  
  s.l2 <- exp(-cal_Ht_MEPP(time.obs=l[delta==0], lambda.par=hazards,alpha.par=alpha, grid.vet=grid)*exp(x.cov[delta==0,] %*% as.matrix(betas)))
  
  lik[delta==0] <- s.l2
  
  loglik = sum(log(lik))
  
  return(-1*loglik)
}
