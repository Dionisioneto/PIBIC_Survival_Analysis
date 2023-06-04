
## candidatos opostos
library(eha)

## Modelo Exponencial por partes

## taxa de falha
cal_ht_MEP <- function(time.obs, lambda.par=lambda.par, grid.vet=grid.vet){
  
  
  dens_MEP <- dpch(x =time.obs, cuts = grid.vet, levels = lambda.par)
  Cumu_MEP <- ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)
  ht <- dens_MEP/(1-Cumu_MEP)
  
  return(ht)
}


cal_Ht_MEP <- function(time.obs, lambda.par=lambda.par, grid.vet=grid.vet){
  Cumu_MEP <- ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)
  Ht  = -log(1-Cumu_MEP)
  return(Ht)
}


SpopMEP <- function(t=t, lambda.par=lambda.par, grid.vet=grid.vet, beta.par=beta.par, theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){
  
  S_MEP <- as.numeric(exp(-cal_Ht_MEP(time.obs=t, lambda.par=lambda.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEP
  return(spop)
}


loglikIC.MEP <- function(a, l=l, r=r, x.cure=x.cure, x.risk=x.risk, grid.vet=grid.vet){

  npar <- length(a)
  
  b <- length(grid.vet)+1
  
  hazards = a[1:b] ## taxas de falha para os b intervalos
  alpha = a[b + 1] ## parametro de potencia
  
  n.cov.cure = dim(x.cure)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  n.cov.risk = dim(x.risk)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  
  betas.cure = a[(b + 2):(b + 1 + n.cov.cure)]
  betas.risk = a[(b + 2 + n.cov.cure):(b + 1 + n.cov.cure + n.cov.risk)]
  
  
  n.sample <- nrow(x.cure)
  
  cens <- ifelse(is.finite(r), 1, 0)
  lik <- rep(0, n.sample)
  
  p2 <- SpopMEP(t=l[cens==1], lambda.par=hazards, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  p1 <- SpopMEP(t=r[cens==1], lambda.par=hazards, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  lik[cens==1] <- p2-p1
  
  p1 <- SpopMEP(t=l[cens==0], lambda.par=hazards, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==0,], x.risk=x.risk[cens==0,])
  lik[cens==0] <- p1
  return(sum(log(lik)))
  
}

particoes = 2

x.f <- cbind(hemo.icens$High, hemo.icens$Medium)
x.c <- cbind(1, hemo.icens$High)

left = hemo.icens$L
right = hemo.icens$R

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", bmax=particoes+1)
grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(rep(0.01, particoes+1), 1, rep(0.1, dim(x.f)[2] - 1 + dim(x.c)[2]))


test <- optim(par = chutes, fn=loglikIC.MEP, gr = NULL, method = "BFGS",
              control=list(fnscale=-1), hessian = TRUE, l=left,
              r=right, x.cure=x.c, x.risk=x.f, grid.vet=grid.obs)







