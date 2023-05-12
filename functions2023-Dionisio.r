################################################
## Code: Funcions for   suvival cure          ## 
## mixture model using 
################################################
## Author: Paulo Junior                       ##
################################################
## Date: 03/05/2023                           ##
################################################



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



time <- function( t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=u.unif ){
  exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.mat%*%beta.par)) - u.unif
}

gen.mepp <- function(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat){
  
  raiz <- uniroot(time, c(0, 10000), lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=runif(1))
  mepp.time <- raiz$root
  
  return(mepp.time)
}


sim.std.cure.ICdata <- function(n=n, lambda.par=lambda.par, alpha.par=alpha.par, 
                                grid.vet=grid.vet, beta.par=beta.par, lambda.parc=lambda.parc, 
                                theta.par = c(1, 0.5, 0), A = 5, B = 15){
  
  
  ########################################################
  ## Generating function for mixture cure model using   ##
  ## PPEM for noncured individuals.                     ##
  ##                                                    ##
  ## Arguments:
  ##
  ## n: Sample size;
  ## lambda.par: PPEM rates;
  ## alpha.par: Power PPEM parameter;
  ## grid.vet: Grid time for PPEM;
  ## beta.par: coef. for time to event;
  ## theta.par: coef. for cure;
  ## lambda.parc: Censoring rate;
  ## A: censoring parameter
  ## B: censoring parameter
  ########################################################
  
  
  ## Testing vals:
  
  # n = 100
  # lambda.par  = c(1.1, 0.8, 0.5)
  # alpha.par   = 0.8
  # grid.vet    = c(0.5, 2)
  # beta.par    = c(-0.5, 0.5)
  # lambda.parc = 1
  # theta.par       = c(1, 0.5, 0)
  # A = 5
  # B = 15
  # 
  ########################################
  
  intercept <- 1
  xi1 <- rbinom(n, 1, 0.5)
  xi2 <- rnorm(n)
  X_cure <- cbind(intercept, xi1, xi2)
  X <- cbind(xi1, xi2)
  
  
  #-- Cure probability mixture model:
  
  elinpred <- exp(-theta.par%*%t(X_cure))
  probY <- 1/(1+elinpred)
  Y <- rbinom(n, size=1, prob=probY)
  
  #-- Censoring times generation:
  
  C <- pmin(A, B*rexp(n, rate=lambda.parc))
  
  #-- Generating the survival times:
  
  t     <- rep(0,n) # tempos de falha
  t_obs <- rep(0,n)
  delta <- rep(0,n)
  
  for(i in 1:n){
    if(Y[i]==1){
      t[i] <- gen.mepp(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=X[i,])
      #t[i] <- gen.mepp(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=0)
      if(t[i]<C[i]){
        t_obs[i] <- t[i]
        delta[i] <- 1  
      }else{
        t_obs[i] <- C[i]
      }
    }else{
      t[i] <- C[i]
      t_obs[i] <- t[i] 
    }
  }
  
  tempo <- t_obs
  
  #-- Observed times:
  
  delta <- ifelse(t < C, 1, 0) 
  
  L <- R <- tempo * NA
  for (i in 1:n) {
    if (delta[i] == 0) {
      L[i] <- tempo[i]
      R[i] <- Inf
    }
    else {
      L[i] <- 0
      add <- stats::runif(1, 0.1, 0.5)
      R[i] <- add
      check <- (L[i] <= tempo[i] & tempo[i] < R[i])
      while (!check) {
        L[i] <- L[i] + add
        add <- stats::runif(1, 0.1, 0.5)
        R[i] <- R[i] + add
        check <- (L[i] <= tempo[i] & tempo[i] < R[i])
      }
    }
  }
  
  dados <- data.frame(L, R, tempo, Y, delta, xi1, xi2)
  
  return(dados)
}



SpopMEPP <- function(t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){
  
  # t          <- dadosIC$tempo
  # lambda.par <- lambda.f
  # alpha.par  <- alpha.f
  # grid.vet   <- grid.time
  # beta.par   <- beta.f
  # theta.par  <- beta.c
  # x.cure     <- x.c
  # x.risk     <- x.f
  
  S_MEPP <- as.numeric(exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEPP
  return(spop)
}


loglikIC <- function(a, l=l, r=r, x.cure=x.cure, x.risk=x.risk, grid.vet=grid.vet){
  
  
  
  # lambda.par  = c(1.1, 0.8, 0.5)
  # alpha.par   = 0.8
  # grid.vet    = c(0.5, 2)
  # beta.par    = c(-0.5, 0.5)
  # theta.par   = c(1, 0.5, 0)
  # 
  # 
  # 
  # l      <- dadosIC$L
  # r      <- dadosIC$R
  # x.cure <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)
  # x.risk <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
  # 
  # a<- c(lambda.par, alpha.par, theta.par,  beta.par)
  # 
  #a <- 1:length(param)
  
  npar <- length(a)
  
  b <- length(grid.vet)+1
  
  hazards = a[1:b] ## taxas de falha para os b intervalos
  alpha = a[b + 1] ## parametro de potencia
  
  n.cov.cure = dim(x.cure)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  n.cov.risk = dim(x.risk)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  
  betas.cure = a[(b + 2):(b + 1 + n.cov.cure)]
  betas.risk = a[(b + 5):(b + 4 + n.cov.risk)]
  

  n.sample <- nrow(x.cure)
  
  cens <- ifelse(is.finite(r), 1, 0)
  lik <- rep(0, n.sample)
  
  p2 <- SpopMEPP(t=l[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  p1 <- SpopMEPP(t=r[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  lik[cens==1] <- p2-p1
  
  p1 <- SpopMEPP(t=l[cens==0], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==0,], x.risk=x.risk[cens==0,])
  lik[cens==0] <- p1
  return(sum(log(lik)))
  
}


loglikIC2 <- function(a, l=l, r=r, x.cure=x.cure, grid.vet=grid.vet){
  
  
  
  # lambda.par  = c(1.1, 0.8, 0.5)
  # alpha.par   = 0.8
  # grid.vet    = c(0.5, 2)
  # beta.par    = c(-0.5, 0.5)
  # theta.par   = c(1, 0.5, 0)
  # 
  # 
  # 
  # l      <- dadosIC$L
  # r      <- dadosIC$R
  # x.cure <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)
  # x.risk <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
  # 
  # a<- c(lambda.par, alpha.par, theta.par,  beta.par)
  # 
  #a <- 1:length(param)
  
  npar <- length(a)
  
  b <- length(grid.vet)+1
  
  hazards = a[1:b] ## taxas de falha para os b intervalos
  alpha = a[b + 1] ## parametro de potencia
  
  n.cov.cure = dim(x.cure)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  #n.cov.risk = dim(x.risk)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  
  betas.cure = a[(b + 2):(b + 1 + n.cov.cure)]
  #betas.cure = a[(b + 5):(b + 4 + n.cov.cure)]
  
  
  
  betas.risk <- 0
  
  n.sample <- nrow(x.cure)
  
  cens <- ifelse(is.finite(r), 1, 0)
  lik <- rep(0, n.sample)
  
  p2 <- SpopMEPP(t=l[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=0)
  p1 <- SpopMEPP(t=r[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=0)
  lik[cens==1] <- p2-p1
  
  p1 <- SpopMEPP(t=l[cens==0], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==0,], x.risk=0)
  lik[cens==0] <- p1
  return(sum(log(lik)))
  
}

# testando o modelo:


library(eha)

n <- 10000 # Tamanho amostral

#--- Par?metros falha:
alpha.f   <- 0.8 
lambda.f  <- c(1.1, 0.8, 0.5)
n.intervals <- length(lambda.f)
grid.time <- c(0.5, 2)
beta.f    <- c(-0.5, 0.5)
#beta.f    <- 0

beta.c    <- c(1, 0.5, -0.5)
lambda.c <- 1



# Ajuste usando os intervalos de tempo:

ini.info <- c(lambda.f, alpha.f, beta.c,  beta.f)

npar <- length(c(lambda.f, alpha.f,beta.c, beta.f))

samp <- 100

est  <- matrix(NA, ncol=npar, nrow=samp)

for(i in 1:samp){
  
  dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f, 
                                 grid.vet=grid.time, beta.par= beta.f, lambda.parc=lambda.c, 
                                 theta.par = beta.c , A = 5, B = 15)
  
  
  x.f <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
  x.c <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)
 
  #ini.info <- rep(1, n.intervals)
  
  test <- optim(par = ini.info, fn=loglikIC, gr = NULL, method = "BFGS",
                control=list(fnscale=-1), hessian = TRUE, l=dadosIC$L, 
                r=dadosIC$R, x.cure=x.c, x.risk=x.f, grid.vet=grid.time)
  
  
  
  # ini.info <- c(lambda.f, alpha.f, beta.c)
  # test <- optim(par = ini.info, fn=loglikIC2, gr = NULL, method = "BFGS",
  #               control=list(fnscale=-1), hessian = TRUE, l=dadosIC$L, 
  #               r=dadosIC$R, x.cure=x.c, grid.vet=grid.time)
  # test$par
  
  est[i,] <- test$par
  cat("Iter=", i, "\n")
  
}

apply(est, 2, mean)



