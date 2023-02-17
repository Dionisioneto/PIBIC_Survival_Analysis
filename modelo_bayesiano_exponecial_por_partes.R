## -----
## Implmentacao Bayesiana do Exponencial por partes
## -----

## Pacotes

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha)

## ajustar o script com as funcoes de sobrevivencia 
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')


########################################################
# Autor: Paulo Cerqueira Jr.                           #
# data: 18/06/2020                                     #
# Código: Calcula o MC para o modelo EPPP com          #
# censura intervalar Bayesiano                         #
########################################################


#--- Limpando o R:

rm(list=ls(all=TRUE))

#--- Load packages to parallelization:

library(snow)
library(snowfall)

#--- Definir o numero de CPU's (tenho no maximo 7):

sfInit(parallel=TRUE, cpus=6, type = "SOCK")


#--- Carregando alguns pacotes para cada CPU:
### Load packages
### -------------------

sfLibrary(runjags)
sfLibrary(eha)
sfLibrary(psych)
sfLibrary(coda)


parallel.bayes.ppe <- function(n.mc=n.mc){
  
  # set seed:
  set.seed(20000901)
  
  
  ### Data generation:
  ### -------------------
  
  #n.mc <- 1
  
  # important functions:
  
  source("functions_ajuste.R")
  
  # Configurations for MCMC:
  
  burnin  <- 4000 #
  lag <- 10
  Npost <- 1500
  Nsim <- burnin + (Npost)*lag
  #mcmc <- list(n.chain=n.chain, burnin=burnin, lag=lag, Nsim=Nsim)
  const <- 10000
  
  # Config. for generation:
  
  N      <- 50 # Tamanho amostral
  
  #--- Parâmetros falha:
  
  alpha.f = 0.8
  lambda.f <- c(1.1, 0.8, 0.5)
  grid.time <- c(0.5, 2)
  b <- length(lambda.f)
  
  beta.f <- c(-0.5, 0.5)
  xi1 <- rbinom(N, 1, 0.5)
  xi2 <- rnorm(N)
  X.real <- cbind(xi1, xi2)
  
  #--- Parâmetros cens:
  
  lambda.c <- 0.1
  
  #--- Definindo o diret?rio:
  
  direct.fun <- getwd()
  
  #--- Definindo o diret?rio para cada CPU:
  
  sub.direct <- paste(direct.fun,"/n", N, sep="")
  
  #--- Definindo a pasta a ser armazenada:
  
  dir.create(sub.direct)
  
  data.gen <- sim.ICdata(n=N, lambda.par=lambda.f, alpha.par=alpha.f, grid.vet=grid.time, beta.par=beta.f, x.mat=X.real,
                         lambda.parc=lambda.c)
  write.table(x = data.gen, file = paste(sub.direct, "/samples_2c_gen","_",N,"_",n.mc, ".txt", sep=""), sep=" ")
  
  
  #data.gen <- read.table(paste(sub.direct, "/samples_gen","_",N,"_",n.mc, ".txt", sep=""), h=T, sep="")
  
  
  ### Model in JAGS
  ### ----------------
  JAGS_PH_MEPP <- "
  model
  {  
 
    for(i in 1:N1)  # For left-censored  
    {
      eta.t1[i] <- inprod(X1[i,], beta.t[])
      Ht.up1[i] <- inprod(xi.falha.up1[i,],lambda.t[])
      St.up1[i] <- exp(-Ht.up1[i])
      Ft.up1[i] <- 1-exp(-Ht.up1[i])
      St.up1.ppe[i] <- pow( (1-pow((Ft.up1[i]),alpha)), exp(eta.t1[i]) )
      loglik1[i]  <- log(1-St.up1.ppe[i]) # funcao de vero.
     
      dummy1[i] ~ dpois(-loglik1[i] + const1)
     
    }
 
    for(i in 1:N2)  # For right-censored  
    {
      ht[i] <- tempo[i]
      St[i] <-
     
      loglik2[i] <- cens[i]*log(ht[i])+log(St[i])
     
     
      dummy2[i] ~ dpois(-loglik2[i] +const2)
    }
 
    for(i in 1:N3)  # For interval-censored  
    {
      eta.t3[i] <- inprod(X3[i,], beta.t[])
      Ht.up3[i] <- inprod(xi.falha.up3[i,],lambda.t[])
      St.up3[i] <- exp(-Ht.up3[i])
      Ft.up3[i] <- 1-exp(-Ht.up3[i])
      St.up3.ppe[i] <- pow( (1-pow((Ft.up3[i]),alpha)), exp(eta.t3[i]) )
     
     
      Ht.lo3[i] <- inprod(xi.falha.lo3[i,],lambda.t[])
      St.lo3[i] <- exp(-Ht.lo3[i])
      Ft.lo3[i] <- 1-exp(-Ht.lo3[i])
      St.lo3.ppe[i] <- pow( (1-pow((Ft.lo3[i]),alpha)), exp(eta.t3[i]) )

      loglik3[i] <- log(St.lo3.ppe[i]-St.up3.ppe[i])
     
      dummy3[i] ~ dpois(-loglik3[i]+const3)
    }
   
    alpha ~ dgamma(0.25, 0.25)
   
    for (j in 1:b)
    {  
          lambda.t[j] ~ dgamma(0.01, 0.01)
    }
   
    for (k in 1:p)
    {
      beta.t[k] ~ dnorm(0,0.001)      
    }
  }
  "
  
  ### Prepare the data such that
  ### they are suitable for JAGS
  ### --------------------------------------
  
  t.low <- data.gen$L
  t.upp <- data.gen$R
  X <- cbind(x1=data.gen$xi1, x2=data.gen$xi2)
  p <- ncol(X)
  
  # type = 1 (lcens), 2 (rcens), 3 (icens)
  type <- rep(3, N)
  type[t.low==0] <- 1
  type[is.infinite(t.upp)] <- 2
  lim <- cbind(t.low, t.upp)
  
  ## Auxiliary matrix for cumulative hazard function:
  
  #a <- time.grid.interval(li = ifelse(is.na(lim[,1]), 0, lim[,1]),
  #                        ri = ifelse(is.na(lim[,2]), Inf, lim[,2]),
  #                        type="OBS", bmax = b)
  
  a <- c(0, grid.time , Inf)
  
  xi.falha.lo <- matrix(0,N,b)
  xi.falha.up <- matrix(0,N,b)
  
  for(i in 1:N){                      
    if(type[i]==3){
      for(j in 1:b){
        xi.falha.lo[i,j] <- (min(lim[i,1], a[j+1])-a[j])*((lim[i,1]-a[j])>0)
        xi.falha.up[i,j] <- (min(lim[i,2], a[j+1])-a[j])*((lim[i,2]-a[j])>0)
      }  
    }else{
      if(type[i]==2){
        for(j in 1:b){
          xi.falha.lo[i,j] <- (min(lim[i,1], a[j+1])-a[j])*((lim[i,1]-a[j])>0)
          #xi.falha.up[i,j] <- (min(lim[i,2], a[j+1])-a[j])*((lim[i,2]-a[j])>0)
        }  
      }else{
        for(j in 1:b){
          #xi.falha.lo[i,j] <- (min(lim[i,1], a[j+1])-a[j])*((lim[i,1]-a[j])>0)
          xi.falha.up[i,j] <- (min(lim[i,2], a[j+1])-a[j])*((lim[i,2]-a[j])>0)
        }  
      }  
      
    }
  }    
  
  ### List with data to be passed to JAGS
  ### (data are divided by type of censoring: left, right, interval)
  ### ----------------------------------------------------------------
  data.gen.jags <- list(
    N1 = sum(type == 1), xi.falha.up1=xi.falha.up[type == 1,], X1 = as.matrix(X[type == 1,]),
    N2 = sum(type == 2), xi.falha.lo2=xi.falha.lo[type == 2,], X2 = as.matrix(X[type == 2,]),
    N3 = sum(type == 3), xi.falha.lo3=xi.falha.lo[type == 3,], xi.falha.up3=xi.falha.up[type == 3,], X3 = as.matrix(X[type == 3,]),
    p=p, b=b, const1=const, const2=const, const3=const ,
    dummy1=rep(0, sum(type==1)),
    dummy2=rep(0, sum(type==2)),
    dummy3=rep(0, sum(type==3)) )
  
  ### Initial parameters: 3 lists for 3 chains
  ### ------------------------------------------------
  
  #inits.gen <- list(lambda.t= rep(0.1,b), alpha=0.5,  beta.t=rep(0,p))  
  
  
  # inits.gen <- list(
  #   list(lambda.t= rep(0.1,b), alpha=0.5,  beta.t=rep(0,p)) ,
  #   list(lambda.t= rep(0.1,b), alpha=1,  beta.t=rep(0,p)))
  
  inits.gen <- list(
    list(lambda.t= rep(0.1,b), alpha=0.5,  beta.t=rnorm(p)) ,
    list(lambda.t= rep(0.5,b), alpha=1,  beta.t=rnorm(p,1)))
  
  
  #parameters.gen <- c("lambda.t", "alpha","beta.t", "loglik1", "loglik2", "loglik3")
  parameters.gen <- c("lambda.t", "alpha","beta.t")
  
  ### Test implementation of JAGS
  ### -----------------------------
  #testjags(silent = FALSE)
  
  
  ### Suppress admininstrative output (just try it)
  ### -----------------------------------------------
  #runjags.options(silent.runjags = TRUE, silent.jags = TRUE)
  
  
  ### JAGS option to run three chains in parallel
  ### ---------------------------------------------
  #runjags.options(method = "rjparallel")
  runjags.options(method = "rjags")       ## sets it back to run everything on just one core
  
  ### Run JAGS
  ### --------------------
  system.time(
    jagsfit <- run.jags(model = JAGS_PH_MEPP,
                        monitor = parameters.gen,
                        data = data.gen.jags, inits = inits.gen,
                        adapt = 1000, n.chains = 2, thin = lag,
                        burnin = burnin, sample = Npost))
  
  summary.jags <- jagsfit$summaries
  write.table(x = summary.jags, file = paste(sub.direct, "/table_2c_gen","_",N,"_",n.mc, ".txt", sep=""), sep=" ")
  
}

mc.val <- 1:500
#mc.val <- mc.val[-c(1:23, 50:56, 81:482)]
sfLapply(mc.val, fun=parallel.bayes.ppe)

sfStop()






