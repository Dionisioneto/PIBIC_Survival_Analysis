################################################
## Code: Funcions for   suvival cure          ## 
## mixture model using MEP
################################################
## Author: Dionisio Neto                      ##
################################################
## Date: 03/05/2023                           ##
################################################

library(eha)
#--- Algumas funcoes importantes:

# Calculando as funcoes de taxa de falha e acumulada:

## funcao taxa de falha do mep

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


## modelo para populacoes mistas

SpopMEP <- function(t=t, lambda.par=lambda.par,
                     grid.vet=grid.vet, beta.par=beta.par,
                     theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){
  
  beta.par = as.matrix(beta.par)
  
  S_MEP <- as.numeric(exp(-cal_Ht_MEP(time.obs=t, lambda.par=lambda.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(-1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEP
  return(spop)
}


loglikIC.MEP <- function(a, l=l, r=r, x.cure=x.cure, x.risk=x.risk, grid.vet=grid.vet){
  
  npar <- length(a)
  b <- length(grid.vet)+1
  
  hazards = a[1:b] ## taxas de falha para os b intervalos
  
  n.cov.cure = dim(x.cure)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  n.cov.risk = dim(x.risk)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  
  betas.cure = a[(b + 1):(b + n.cov.cure)]
  betas.risk = a[(b + 1 + n.cov.cure):(b + n.cov.cure + n.cov.risk)]
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

fit.mepp.cf = function(L, R, n.int, cov.risco, cov.cura, start){
  
  ## extracao do grid observado
  
  grid.obs=time.grid.interval(li=L, ri=R, type="OBS", bmax= n.int)
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  
  est <- optim(par = start, fn=loglikIC, gr = NULL, method = "BFGS",
               control=list(fnscale=-1), hessian = TRUE, l=L, 
               r=R, x.cure=cov.cura, x.risk=cov.risco, grid.vet=grid.obs)
  
  estimated = est$par
  hessian = est$hessian
  loglik = est$value
  
  results = list(estimated = estimated, hessian = hessian, loglik = loglik)
  
  return(results)
}


fit.mep.cf = function(l, r, n.int, cov.risco, cov.cura, start){
  
  grid.obs=time.grid.interval(li=l, ri=r,
                              type="OBS", bmax=n.int)
  
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  
  max.mep = optim(par = start, fn=loglikIC.MEP, gr = NULL, method = "BFGS",
                  control=list(fnscale=-1), hessian = TRUE,
                  l=l, r=r,
                  x.cure=cov.cura, x.risk=cov.risco, 
                  grid.vet=grid.obs)
  
  estimated = max.mep$par
  hessian = max.mep$hessian
  loglik = max.mep$value
  
  results = list(estimated = estimated, hessian = hessian, loglik = loglik)
  
  return(results)
  
}

## lambda = 3
## beta.cura = 5
## beta.risco = 4

chute = c(0.1,0.2,0.1,
          1.2,0.1,0.1,0.1,0.1,
          0.1,0.1,0.1,0.1)

ajuste.mep = fit.mep.cf(l = smoke2009$Timept1, r = smoke2009$Timept2,
              n.int = 3, cov.risco = covariaveis, cov.cura = cbind(1,covariaveis),
            start = chute)

ajuste.mep$estimated
ajuste.mep$loglik

AIC.surv(ajuste.mep$loglik, n.param = length(ajuste.mep$estimated))

BIC.surv(ajuste.mep$loglik, n.param = length(ajuste.mep$estimated),
         n.sample = dim(smoke2009)[1])

HC.surv(ajuste.mep$loglik, n.param = length(ajuste.mep$estimated),
        n.sample = dim(smoke2009)[1])

## ---
## dados aneurysm
## ---

chute = c(0.1,0.1,
          1,
          1, 0.1,0.1,0.1,
          0.1,0.1,0.1)

gr_pad = (aneurysm$gr - mean(aneurysm$gr))/sd(aneurysm$gr) ## precisa normalizar para estimar

covariaveis.aneurysm = cbind(aneurysm$mo, gr_pad, aneurysm$lok)

ajuste.mep = fit.mep.cf(l = aneurysm$t.left, r = aneurysm$t.right,
                        n.int = 2, cov.risco = covariaveis.aneurysm, 
                        cov.cura = cbind(1,covariaveis.aneurysm),
                        start = chute)

ajuste.mep$estimated ## estimacao errada do parametro potencia
ajuste.mep$loglik

AIC.surv(ajuste.mep$loglik, n.param = length(ajuste.mep$estimated))

BIC.surv(ajuste.mep$loglik, n.param = length(ajuste.mep$estimated),
         n.sample = dim(smoke2009)[1])

HC.surv(ajuste.mep$loglik, n.param = length(ajuste.mep$estimated),
        n.sample = dim(smoke2009)[1])

## -----
## dados breast cancer
## -----

n.intervalos = 15

for (int in 2:15){
  
  n.intervalos = int
  
  chute = c(rep(0.1,n.intervalos),
            1.2,0.1,
            0.1)
  
  ajuste.mep.breast = fit.mep.cf(l = breast$left, r = breast$right,
                                 n.int = n.intervalos, cov.risco = cbind(breast$ther), 
                                 cov.cura = cbind(1, breast$ther),
                                 start = chute)
  
  ajuste.mep.breast$estimated ## estimacao errada do parametro potencia
  ajuste.mep.breast$loglik
  
  aicc = AIC.surv(ajuste.mep.breast$loglik,
           n.param = length(ajuste.mep.breast$estimated))
  
  bicc = BIC.surv(ajuste.mep.breast$loglik,
           n.param = length(ajuste.mep.breast$estimated),
           n.sample = dim(breast)[1])
  
  hcc = HC.surv(ajuste.mep.breast$loglik,
          n.param = length(ajuste.mep.breast$estimated),
          n.sample = dim(breast)[1])
  
  print(paste("n intervalo", int))
  print(paste("AIC: ", aicc))
  print(paste("BIC: ", bicc))
  print(paste("HC: ", hcc))
  
}









































