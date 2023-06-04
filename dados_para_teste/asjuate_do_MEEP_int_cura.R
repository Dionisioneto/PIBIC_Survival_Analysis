## ---------
## Proposta de ajuste para o
## Modelo Exponencial por Partes Potencia
## em censura intervalar e fracao de cura
## ---------



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


#--- Algumas funcoes importantes:

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


SpopMEPP <- function(t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, 
                     beta.par=beta.par, theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){

  
  S_MEPP <- as.numeric(exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEPP
  return(spop)
}


## funcao de verossimilhanca
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
  
  return(sum(log(lik)))
}

dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f, 
                               grid.vet=grid.time, beta.par= beta.f, lambda.parc=1, 
                               theta.par = beta.c , A = 5, B =10)

## dados da populcao de curados e nao-curados
particoes = 2

#left = dadosIC$L
#right = dadosIC$R

left = hemo.icens$L
right = hemo.icens$R

#x.f <- cbind(x1=dadosIC$xi1,x2=dadosIC$xi2)
#x.c <- cbind(1, x1=dadosIC$xi1,x2=dadosIC$xi2)

x.f = cbind(x1=hemo.icens$Low,x2=dadosIC$xi2[1:dim(hemo.icens)[1]])
x.c = cbind(1, x1=hemo.icens$Low,x2=dadosIC$xi2[1:dim(hemo.icens)[1]])

## algoritmo para escolha do grid

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", bmax=particoes+1)
grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(rep(0.01, particoes+1), 1, 0.9, rep(0.05, dim(x.f)[2] - 1 + dim(x.c)[2]))

test <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
              control=list(fnscale=-1), hessian = T, l=left, 
              r=right, x.cure=x.c, x.risk=x.f, grid.vet=grid.obs)


test2 <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "Nelder-Mead",
              control=list(fnscale=-1), hessian = F, l=left, 
              r=right, x.cure=x.c, x.risk=x.f, grid.vet=grid.obs)

test2$par


















