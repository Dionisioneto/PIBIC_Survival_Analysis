################################################
## Code: Funcions for   suvival cure          ## 
## mixture model using 
################################################
## Author: Paulo Junior                       ##
################################################
## Date: 03/05/2023                           ##
################################################

time.grid.interval <- function(li=li, ri=ri, type=type, bmax=bmax)
{  
  ## Funcao que retorna os intervalos da partiÃƒÂ§ÃƒÂ£o mais fina
  ## baseada nos limites observados, distintos e finitos.
  ## Argumentos:
  ## li: limite inferior dos intervalos observados.
  ## ri: limite superior dos intervalos observados.
  ## bmax: numero mÃƒÂ¡ximo de intervalos.
  
  # li = dados$L
  # ri= dados$R
  
  #--- Inicio da funcao:
  
  #-- Construir uma grade tipo 1: baseando-se em tempos observaveis
  if(type=="OBS")
  {
    #grid.vet <- sort(unique(c(0, li, ri, Inf)))
    grid.vet <- sort(unique(c(0, li[is.finite(li)], ri[is.finite(ri)], Inf)))
    grid.size.vet <- length(grid.vet) # Grid time size
    
    if( isTRUE(bmax<grid.size.vet)==TRUE )
    {
      k        <- round((length(grid.vet)-1)/bmax,0)
      id.grid  <- round(seq(k,(length(grid.vet)-1), length.out=bmax),0)
      grid.vet <- c(0,grid.vet[-1][id.grid])
      return(grid.vet)
    }else{
      grid.vet <- sort(unique(c(0, li, ri, Inf)))
      return(grid.vet)
    }
  } #-- Construir uma grade tipo 2: espacos equiprovaveis
  if(type=="EQUI")
  {
    grade.vet <- seq(0, max(ri[ri!=Inf]), length.out=bmax)
    grid.vet <- c(grade.vet,Inf)
    return(grid.vet)
  }  
}



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


sim.std.cure.ICdata <- function(n, lambda.par, alpha.par, 
                                grid.vet, beta.par, lambda.parc, 
                                theta.par, A = 5, B = 15){
  
  
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
  
  
  elinpred <- exp(-1*(X_cure%*%theta.par))
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
  
  beta.par= as.matrix(beta.par)
  
  S_MEPP <- as.numeric(exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(-1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- (1-probY)+(probY*S_MEPP)
  return(spop)
}


loglikIC <- function(par, l, r, x.cure, x.risk, grid){
  
  b <- length(grid)+1
  
  hazards = par[1:b] ## taxas de falha para os b intervalos
  alpha = par[b + 1] ## pparrametro de potencia
  
  n.cov.cure = dim(x.cure)[2] 
  n.cov.risk = dim(x.risk)[2] 
  
  betas.cure = par[(b + 2):(b + 1 + n.cov.cure)]
  betas.risk = par[(b + 2 + n.cov.cure):(b + 1 + n.cov.cure + n.cov.risk)]
  n.sample <- nrow(x.cure)
  
  delta <- ifelse(is.finite(r), 1, 0)
  
  lik <- rep(0, n.sample)
  
  spop.l <- SpopMEPP(t=l[delta==1], lambda.par=hazards, alpha.par=alpha,
                     grid.vet=grid, beta.par=betas.risk, theta.par=betas.cure,
                     x.cure=x.cure[delta==1,], x.risk=x.risk[delta==1,])
  
  spop.r <- SpopMEPP(t=r[delta==1], lambda.par=hazards, alpha.par=alpha,
                     grid.vet=grid, beta.par=betas.risk, theta.par=betas.cure,
                     x.cure=x.cure[delta==1,], x.risk=x.risk[delta==1,])
  
  lik[delta==1] <- spop.l-spop.r
  
  spop.l2 <- SpopMEPP(t=l[delta==0], lambda.par=hazards, alpha.par=alpha,
                      grid.vet=grid, beta.par=betas.risk, theta.par=betas.cure,
                      x.cure=x.cure[delta==0,], x.risk=x.risk[delta==0,])
  lik[delta==0] <- spop.l2
  
  return(-1*sum(log(lik)))
}


## Funcao que particiona os grids para censura intervalar

time.grid.interval <- function(li=li, ri=ri, type=type, bmax=bmax){  
  ## Funcao que retorna os intervalos da partiÃƒÂ§ÃƒÂ£o mais fina
  ## baseada nos limites observados, distintos e finitos.
  ## Argumentos:
  ## li: limite inferior dos intervalos observados.
  ## ri: limite superior dos intervalos observados.
  ## bmax: numero mÃƒÂ¡ximo de intervalos.
  
  # li = dados$L
  # ri= dados$R
  
  #--- Inicio da funcao:
  
  #-- Construir uma grade tipo 1: baseando-se em tempos observaveis
  if(type=="OBS")
  {
    #grid.vet <- sort(unique(c(0, li, ri, Inf)))
    grid.vet <- sort(unique(c(0, li[is.finite(li)], ri[is.finite(ri)], Inf)))
    grid.size.vet <- length(grid.vet) # Grid time size
    
    if( isTRUE(bmax<grid.size.vet)==TRUE )
    {
      k        <- round((length(grid.vet)-1)/bmax,0)
      id.grid  <- round(seq(k,(length(grid.vet)-1), length.out=bmax),0)
      grid.vet <- c(0,grid.vet[-1][id.grid])
      return(grid.vet)
    }else{
      grid.vet <- sort(unique(c(0, li, ri, Inf)))
      return(grid.vet)
    }
  } #-- Construir uma grade tipo 2: espacos equiprovaveis
  if(type=="EQUI")
  {
    grade.vet <- seq(0, max(ri[ri!=Inf]), length.out=bmax)
    grid.vet <- c(grade.vet,Inf)
    return(grid.vet)
  }  
}


## ajuste do modelo exponencial por partes potencia 
## em censura intervalar e fracao de curados

fit.mepp.cf = function(L, R, n.int, cov.risco, cov.cura, start){
  
  ## extracao do grid observado
  
  grid.obs=time.grid.interval(li=L, ri=R, type="OBS", bmax= n.int)
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  
  est <- optim(par = start, fn=loglikIC, gr = NULL, method = "BFGS",
               control=list(fnscale=1), hessian = TRUE, l=L, 
               r=R, x.cure=cov.cura, x.risk=cov.risco, grid=grid.obs)
  
  estimated = est$par
  hessian = est$hessian
  loglik = est$value
  
  results = list(estimated = estimated, hessian = hessian, loglik = loglik)
  
  return(results)
}





