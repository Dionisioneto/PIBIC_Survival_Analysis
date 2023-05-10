## ------
## Funcoes PIBIC 2023
## Projeto: Modelos Semi-Parametricos para dados
## de sobrevivencia com censura intervalar
##
## Aplicacao a censura intervalar
## ------

## ---
## Pacotes 
## ---

source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(survival, ReIns)

## ------
## Funcao geradora de dados para censura intervalar e fracao de cura,
## com covariaveis
## ------

N =50

survival.time = rep(NA, N)
L = rep(NA, N)
R = rep(NA, N)
delta = rep(NA, N)

# # Colunas
# B: v.a Bernoulli para indicar curado e nao curado
# L: limite inferior do intervalo observado
# R: limite superior do intervalo observado
# delta: indicador de censura
# x1: valores para a primeira variavel continua (Normal padrao)
# x2: valores para a segunda variavel continua (Bernoulli)

#pi = exp(x.matriz %*% betas)/(1 + exp(x.matriz %*% betas))



cbind(survival.time, L , R, delta, Y)

dados.cens = survival::Surv(survival.time,delta)
km_fit <- survival::survfit(dados.cens ~ 1)

plot(km_fit)

sim.cure.icen = function(N, prob.ber, betas.cure, betas.risk,
                         c1, c2){
  
  ## Armazenameto de dados
  
  survival.time = rep(NA, N)
  L = rep(NA, N)
  R = rep(NA, N)
  delta = rep(NA, N)
  
  ## geracao dos dados
  x2 = rnorm(n = N, mean = 0, sd = 1)
  x1 = rbinom(n = N, size = 1, prob = prob.ber)
  matriz = cbind(rep(1, N), x1, x2)
  
  b = betas.cure
  
  z = b[1] + b[2]*x1 + b[3]*x2
  
  prob.z = 1/(1+exp(-z))
  
  Y = rbinom(N,1,prob.z)
  
  ## Proporcao de censuras
  
  A = rexp(n = N, rate = 1)
  
  C = cbind(c1,c2*A)
  C = pmin(C[,1], C[,2])
  
  ##------
  for(i in 1:N){
    if(Y[i] == 0){
      survival.time[i] = C[i]
      delta[i] = 0
    }
    
    else{
      survival.time[i] = exp(matriz[i,-1] %*% betas.risk)
      delta[i] = ifelse(survival.time[i] <= C[i], 1, 0)
    }
    
    if(delta[i] == 0){
      L[i] = C[i]
      R[i] = Inf
    }
    
    else{
      L[i] = 0
      Qj = runif(1, min = 0.1, max = 0.5)
      R[i] = Qj
      check = (L[i] <= survival.time[i] & survival.time[i] < R[i])
      
      while(!check){
        L[i] = L[i] + Qj
        Qj = runif(1, min = 0.1, max = 0.5)
        R[i] = R[i] + Qj
        check = (L[i] <= survival.time[i] & survival.time[i] < R[i])
      }
    }
  }
  data = data.frame(time = survival.time, L = L, R = R, 
                    delta = delta, risk = Y, X1 = matriz[,2], X2 = matriz[,3])
  return(data)
}

dados = sim.cure.icen(N = 100, prob.ber = 0.5, betas.cure = c(1.2,0.5,1.2),
              betas.risk = c(0.5,0.3), c1 = 12, c2 = 17)



dados.cens = Surv(dados$time, dados$delta)
km = survfit(dados.cens~1)

plot(km)

## -----
## Estudo da curva de sobrevivencia pelo estimador de
## Turnbull
## -----

tempo.aval = seq(0,12,0.1)
trnb.fit = Turnbull(x = tempo.aval, L = dados$L, R = dados$R,
                    censored = dados$delta)

plot(tempo.aval, trnb.fit$surv, type = "s", 
     ylab = "Estimador de Turnbull para S(t)")

## ------
## Funcao geradora de dados para censura intervalar e fracao de cura,
## com covariaveis: Versao do Modelo Exponencial por Partes potencia.
## ------

sim.cure.icen.mepp = function(N, prob.ber, betas.cure, betas.risk,
                         c1, c2, lambdas, grid, alpha){
  
  ## Armazenameto de dados
  
  survival.time = rep(NA, N)
  L = rep(NA, N)
  R = rep(NA, N)
  delta = rep(NA, N)
  
  ## geracao dos dados
  x2 = rnorm(n = N, mean = 0, sd = 1)
  x1 = rbinom(n = N, size = 1, prob = prob.ber)
  matriz = cbind(rep(1, N), x1, x2)
  
  b = betas.cure
  
  z = matriz %*% b
  
  prob.z = 1/(1+exp(-z))
  
  Y = rbinom(N,1,prob.z)
  
  ## Proporcao de censuras
  
  A = rexp(n = N, rate = 1)
  
  C = cbind(c1,c2*A)
  C = pmin(C[,1], C[,2])
  
  ##------
  for(i in 1:N){
    if(Y[i] == 0){
      survival.time[i] = C[i]
      delta[i] = 0
    }
    
    else{
      #survival.time[i] = exp(matriz[i,-1] %*% betas.risk)
      survival.time[i] = gen.mepp.cox(n = 1, lambda.par = lambdas, alpha.par = alpha,
                                      cuts = grid, beta.par = betas.risk, x.mat = matriz[,-1])
      
      delta[i] = ifelse(survival.time[i] <= C[i], 1, 0)
    }
    
    if(delta[i] == 0){
      L[i] = C[i]
      R[i] = Inf
    }
    
    else{
      L[i] = 0
      Qj = runif(1, min = 0.1, max = 0.5)
      R[i] = Qj
      check = (L[i] <= survival.time[i] & survival.time[i] < R[i])
      
      while(!check){
        L[i] = L[i] + Qj
        Qj = runif(1, min = 0.1, max = 0.5)
        R[i] = R[i] + Qj
        check = (L[i] <= survival.time[i] & survival.time[i] < R[i])
      }
    }
  }
  data = data.frame(time = survival.time, L = L, R = R, 
                    delta = delta, risk = Y, X1 = matriz[,2], X2 = matriz[,3])
  return(data)
}


dados = sim.cure.icen.mepp(N = 100, prob.ber = 0.5, betas.cure = c(1.2,0.5,1.2),
                      betas.risk = c(0.5,0.3), c1 = 12, c2 = 17,
                      lambdas = c(0.6, 0.8, 0.9), grid = c(0.3,0.6), alpha = 1.4)

head(dados)

dados.cens = Surv(dados$time, dados$delta)
km = survfit(dados.cens~1)

plot(km)


## -----
## Estudo da curva de sobrevivencia pelo estimador de
## Turnbull
## -----

tempo.aval = seq(0,12,0.1)
trnb.fit = Turnbull(x = tempo.aval, L = dados$L, R = dados$R,
                    censored = dados$delta)

plot(tempo.aval, trnb.fit$surv, type = "s", 
     ylab = "Estimador de Turnbull para S(t)")


## -----
## funcao de sobrevivencia para a populacao curada
## -----

spop.mepp = function(time, cuts, levels, alpha, cure.matrix, risk.matrix,
                     betas.cure, betas.risk){
  
  st = as.numeric(PPE(time = time, cuts = cuts, 
            levels = levels, alpha = alpha, type = "survival"))
  
  stcox = as.numeric(st^(exp(risk.matrix %*% betas.risk)))
  
  xb = as.numeric(exp(cure.matrix %*% betas.cure))
    
  pi = 1/(1+exp(-xb))
  
  spop = pi + ((1-pi)*stcox)
  
  return(spop)
}

# spop.mepp(time = dados$time, cuts = grid, levels = taxas.de.falha, alpha = potencia,
#           cure.matrix = as.matrix(cbind(1,dados$X1, dados$X2)),
#           risk.matrix = as.matrix(cbind(dados$X1, dados$X2)),
#           betas.cure = c(1, 0.3,0.6), betas.risk = c(-0.5,1.2))


sl = spop.mepp(time = dados$L, cuts = c(0.3,0.6), levels = c(0.6, 0.8, 0.9),
               alpha = 1.4,
              cure.matrix = as.matrix(cbind(1,dados$X1, dados$X2)),
              risk.matrix = as.matrix(cbind(dados$X1, dados$X2)),
              betas.cure = c(1, 0.3,0.6), betas.risk = c(-0.5,1.2))

sr = spop.mepp(time = dados$R, cuts = c(0.3,0.6), levels = c(0.6, 0.8, 0.9),
              alpha = 1.4,
              cure.matrix = as.matrix(cbind(1,dados$X1, dados$X2)),
              risk.matrix = as.matrix(cbind(dados$X1, dados$X2)),
              betas.cure = c(1, 0.3,0.6), betas.risk = c(-0.5,1.2))

## ------
## Funcao de verosimilhanca para censura intervalar e fracao de cura,
## com covariaveis
## ------


loglik.int.fc = function(par, time.l, time.r,
                         grid, delta, cure.matrix, risk.matrix){
  

  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  alpha = par[b + 1] ## parametro de potencia

  n.covars.cure = dim(cure.matrix)[2] ## numero de covariaveis com fracao de cura
  n.covars.risk = dim(risk.matrix)[2] ## numero de covariaveis com fracao de risco
  
  betas.cure = par[(b + 2):(b + 1 + n.covars.cure )]
  betas.risk = par[(b + 2 + n.covars.cure):(b + 1 +  n.covars.cure+ n.covars.risk)]
  
  lik <- rep(0, dim(cure.matrix)[1])
  
  ## informacoes exponencial por partes Potencia (PPE), sob fracao de cura, para esquerda.
  sl =  spop.mepp(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = alpha,
                  cure.matrix = cure.matrix[delta==1,],
                  risk.matrix = risk.matrix[delta==1,],
                  betas.cure = betas.cure, betas.risk = betas.risk)
  
  ## informacoes exponencial por partes Potencia (PPE), sob fracao de cura, para direita.
  sr =  spop.mepp(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = alpha,
                  cure.matrix = cure.matrix[delta==1,],
                  risk.matrix = risk.matrix[delta==1,],
                  betas.cure = betas.cure, betas.risk = betas.risk)
  
  ## contribuicao das covariaveis

  lik[delta==1] = sl-sr
  
  ## contribuicao da censura
  sl =  spop.mepp(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = alpha,
                    cure.matrix = cure.matrix[delta==0,],
                    risk.matrix = risk.matrix[delta==0,],
                    betas.cure = betas.cure, betas.risk = betas.risk)
   
  lik[delta==0] = sl
 
  log.vero = sum(log(lik))

 return(log.vero)

}


taxas.de.falha = c(0.2, 0.4, 0.9)
particoes = c(0.3,0.6)
potencia = 0.8

betas.cura = c(1.2,0.5,0.3)
betas.risco = c(1.3, 0.8)

parametros = c(taxas.de.falha, potencia, betas.cura, betas.risco)

dados = sim.cure.icen.mepp(N = 100, prob.ber = 0.5, betas.cure = betas.cura,
                           betas.risk = betas.risco, c1 = 4, c2 = 7,
                           lambdas = taxas.de.falha, grid = particoes, alpha = potencia)

## tentativa de estimacao por Maxima Verossimilhanca
grid.obs = time.grid.interval(li = dados$L, ri = dados$R, 
                          type = "OBS", bmax = length(taxas.de.falha))

grid.obs = grid.obs[-c(1, length(grid.obs))]

chutes = c(rep(0.1,length(grid)+1),1, 1, 0.1, 0.1, 0.1, 0.1)


## Metodo numerico BFGS
estimacao.intervalar = optim(par = parametros,
                             fn = loglik.int.fc,
                             gr = NULL,
                             hessian = TRUE,
                             method = "BFGS",
                             time.l =  dados$L, 
                             time.r = dados$R,
                             grid = particoes,
                             delta = dados$delta,
                             cure.matrix = cbind(1,dados$X1, dados$X2),
                             risk.matrix = cbind(dados$X1, dados$X2))

cbind(parametros, estimacao.intervalar$par)

estimacao.intervalar2 = optim(par = chutes,
                             fn = loglik.int.fc,
                             gr = NULL,
                             hessian = F,
                             method = "Nelder-Mead",
                             time.l =  dados$L, 
                             time.r = dados$R,
                             grid = grid.obs,
                             delta = dados$delta,
                             cure.matrix = cbind(1,dados$X1, dados$X2),
                             risk.matrix = cbind(dados$X1, dados$X2))
estimacao.intervalar2$par
parametros




## -----
## versao 2
## -----

spop.mepp = function(time, cuts, levels, alpha, cure.matrix, risk.matrix,
                     betas.cure, betas.risk){
  
  st = as.numeric(PPE(time = time, cuts = cuts, 
                      levels = levels, alpha = alpha, type = "survival"))
  
  stcox = as.numeric(st^(exp(risk.matrix %*% betas.risk)))
  
  xb = as.numeric(exp(cure.matrix %*% betas.cure))
  
  pi = 1/(1+exp(-xb))
  
  spop = pi + ((1-pi)*stcox)
  
  return(spop)
}

sl = spop.mepp(time = dados$L, cuts = c(0.3,0.6), levels = c(0.6, 0.8, 0.9),
               alpha = 1.4,
               cure.matrix = as.matrix(cbind(1,dados$X1, dados$X2)),
               risk.matrix = as.matrix(cbind(dados$X1, dados$X2)),
               betas.cure = c(1, 0.3,0.6), betas.risk = c(-0.5,1.2))

sr = spop.mepp(time = dados$R, cuts = c(0.3,0.6), levels = c(0.6, 0.8, 0.9),
               alpha = 1.4,
               cure.matrix = as.matrix(cbind(1,dados$X1, dados$X2)),
               risk.matrix = as.matrix(cbind(dados$X1, dados$X2)),
               betas.cure = c(1, 0.3,0.6), betas.risk = c(-0.5,1.2))


## funcao log-verossimilhanca

loglikIC <- function(a, l=l, r=r, x.cure=x.cure, x.risk=x.risk, grid.vet=grid.vet){
  
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
  
  p2 <- spop.mepp(time=l[cens==1], levels=hazards, alpha=alpha, 
                  cuts=grid.vet, betas.risk=betas.risk, betas.cure=betas.cure,
                  cure.matrix=x.cure[cens==1,], risk.matrix=x.risk[cens==1,])
  
  p1 <-  spop.mepp(time=r[cens==1], levels=hazards, alpha=alpha, 
                   cuts=grid.vet, betas.risk=betas.risk, betas.cure=betas.cure,
                   cure.matrix=x.cure[cens==1,], risk.matrix=x.risk[cens==1,])
  
  lik[cens==1] <- p2-p1
  
  p1 <- spop.mepp(time=l[cens==0], levels=hazards, alpha=alpha, 
                  cuts=grid.vet, betas.risk=betas.risk, betas.cure=betas.cure,
                  cure.matrix=x.cure[cens==1,], risk.matrix=x.risk[cens==1,])
  
  lik[cens==0] <- p1
  return(sum(log(lik)))
  
}

#--- Parametros falha:
alpha.f   <- 0.8 
lambda.f  <- c(1.1, 0.8, 0.5)
n.intervals <- length(lambda.f)
grid.time <- c(0.5, 2)
beta.f    <- c(-0.5, 0.5)
#beta.f    <- 0

beta.c    <- c(1, 0.5, -0.5)
lambda.c <- 1



taxas.de.falha = c(1.1, 0.8, 0.5)
particoes = c(0.5, 2)
potencia = 0.8

betas.cura = c(1, 0.5, -0.5)
betas.risco = c(-0.5, 0.5)

parametros = c(taxas.de.falha, potencia, betas.cura, betas.risco)

dados = sim.cure.icen.mepp(N = 100, prob.ber = 0.5, betas.cure = betas.cura,
                           betas.risk = betas.risco, c1 = 4, c2 = 7,
                           lambdas = taxas.de.falha, grid = particoes, alpha = potencia)



## Metodo numerico BFGS
estimacao.intervalar = optim(par = parametros, fn=loglikIC,
                             gr = NULL, method = "BFGS",
                             control=list(fnscale=-1),
                             hessian = TRUE,
                             l=dados$L, 
                             r=dados$R, 
                             x.cure=cbind(1, x1=dados$X1, x2=dados$X2),
                             x.risk=cbind(x1=dados$X1, x2=dados$X2),
                             grid.vet=particoes)














