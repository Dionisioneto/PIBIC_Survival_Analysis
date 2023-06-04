
## estudo particular do MEPP

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

SpopMEPP <- function(t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){
  
  Ht_MEPP_cox = cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)
  S_MEPP <- as.numeric(exp(-Ht_MEPP_cox))
  
  elinpred <- as.numeric(exp(-1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEPP
  return(spop)
}


loglikIC <- function(a, l=l, r=r, x.cure=x.cure, x.risk=x.risk, grid.vet=grid.vet){
  
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
  p2 <- SpopMEPP(t=l[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  p1 <- SpopMEPP(t=r[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  
  lik[cens==1] <- p2-p1
  p1 <- SpopMEPP(t=l[cens==0], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==0,], x.risk=x.risk[cens==0,])
  lik[cens==0] <- p1
  
  return(-1*sum(log(lik)))
}

n = 500

alpha.f   <- 0.8 
lambda.f  <- c(1.1, 0.3, 0.9)
grid.time <- c(0.5, 2)
beta.f    <- c(-0.5, 0.8)
beta.c    <- c(1.2, 0.5, -0.5)

Theta = c(lambda.f, alpha.f,beta.c,beta.f)

dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f, 
                               grid.vet=grid.time, beta.par= beta.f, lambda.parc=1, 
                               theta.par = beta.c , A = 0.4, B =22)

x.f <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
x.c <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)

grid.obs=time.grid.interval(li=dadosIC$L, ri=dadosIC$R, type="OBS", bmax=length(lambda.f))
grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(rep(0.1, length(lambda.f)), 1, 1, 0.5, 0.5, 0.5, 0.5)

max.mepp = optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                 control=list(fnscale=-1), hessian = TRUE, l=dadosIC$L,
                 r=dadosIC$R, x.cure=x.c, x.risk=x.f, grid.vet=grid.obs)

cbind(max.mepp$par, Theta)

estimacao.int.cox = optim(par = chutes,
                          fn = loglik.int,
                          gr = NULL,
                          hessian = TRUE,
                          method = "BFGS",
                          time.l = dadosIC$L, 
                          time.r = dadosIC$R,
                          grid = grid.obs ,
                          delta = ifelse(dadosIC$R == Inf,0,1),
                          x.matrix = cbind(dadosIC$xi1, dadosIC$xi2))



















