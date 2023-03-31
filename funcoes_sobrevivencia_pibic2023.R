## ------
## Funcoes PIBIC 2023
## Projeto: Modelos Semi-Parametricos para dados
## de sobrevivencia com censura intervalar
## ------

## ---
## Pacotes 
## ---

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha)

# eha: pacote do exponencial por partes pot?ncia
# latex2exp: pacote para a inserir a linguagem latex em graficos


## ---
## Funcoes para a comparacao de modelos
## ---

## AIC

## Recomenda-se usar o AIC para selecionar modelos 
## quando o numero de observacoes
## (n) e maior que pelo menos 40 vezes o numero de 
## parametros (p)

AIC.surv = function(loglik, n.param) {
  aic_formula = 2*(n.param - loglik)
  return(aic_formula)
}


## BIC

## O BIC e obtido atraves de resultados assintoticos e da 
## suposicao de que os dados pertencem a familia exponencial

BIC.surv = function(loglik, n.param, n.sample){
  bayes_criterion = -2*loglik + (log(n.sample) * n.param)
  return(bayes_criterion)
}

## HC

HC.surv = function(loglik, n.param, n.sample){
  hc_criterion = -2*(loglik) + (2*n.param*log(log(n.sample)))
  return(hc_criterion)
}


## ---
## Funcao para trazer valores das
## funcoes do Exponencial por partes (EP, PEM)
## ---



PE = function(time, cuts, levels, type = "survival"){
  
  ## funcao de distribuicao e densidade do exponencial por partes (PE)
  F0.t =  ppch(q = time, cuts = cuts, levels = levels)
  
  ## sobrevivencia
  s0.t = 1 - F0.t
  
  ## densidade
  f0.t = (dpch(x = time, cuts = cuts, levels = levels))
  
  ## taxa de falha
  h0.t = f0.t/s0.t
  
  ## taxa de falha acumulada
  H0.t = -log(s0.t)
  
  if(type == "survival"){
    return(s0.t)
  } else {
    if(type == "density"){
      return(f0.t)
    } else {
      if(type == "hazard"){
        return(h0.t)
      } else{
        if(type == "cum_hazard"){
          return(H0.t)
        } else {
          if(type == "distribution"){
            return(F0.t) 
          } else{
            return("Please choose one type of function!!!")
          }
        }
      }
    }
  }
}

## ------------------------------------------------------------------------------

## ---
## Modelo Exponencial por partes por potencia (PPE)
## ---

PPE = function(time, cuts, levels, alpha, type = "survival"){
  
  ## ---
  ## Seccao do Exponencial por partes
  ## ---
  
  ## sobrevivencia
  F0.t = PE(time = time, cuts = cuts, levels = levels, type = 'distribution')
  
  ## densidade
  f0.t = PE(time = time, cuts = cuts, levels = levels, type = 'density')
  
  ## ---
  ## Atribuindo aos valores do Exponencial por partes potencia (PPE) 
  ## ---
  
  ## funcao distribuicao do PPE
  F1.t = (F0.t)^alpha
  
  ## funcao densidade do PPE
  f1.t = (alpha*((F0.t)^(1 - alpha))) * f0.t
  
  ## funcao de sobrevivencia do PPE
  s1.t = 1 - F1.t
  
  ## funcao taxa de falha do PPE
  h1.t = f1.t/s1.t
  
  ## taxa de falha acumulada do PPE
  H1.t = -log(s1.t)
  
  if(type == "survival"){
    return(s1.t)
  } else {
    if(type == "density"){
      return(f1.t)
    } else {
      if(type == "hazard"){
        return(h1.t)
      } else{
        if(type == "cum_hazard"){
          return(H1.t)
        } else {
          if(type == "distribution"){
            return(F1.t) 
          } else{
            return("Please choose one type of function!!!")
          }
        }
      }
    }
  }                                
}


## quantil na funcao de distribuicao, dado os grids, hazards e
## alpha no MEPP.

qppe = function(p, cuts, levels, alpha, lower.tail = T, log.p = F){
  if(log.p)
    p = exp(p)
  if (any(p >=1))
    stop("p must be < 1")	
  if(any(p <= 0))
    stop("p must be > 0")
  if(!lower.tail)
    p = 1 - p
  
  n = length(p)
  y = as.numeric(p)
  
  func = function(q, u){
    PPE(q, cuts = cuts, levels = levels, alpha = alpha, type = 'distribution') - u
  }
  
  for(i in 1:n){
    y[i] = uniroot(func, interval = c(0,1000), u = p[i])$root
  }	
  return(y)
} 

## Geracao de dados para MEPP, dado os grids, hazards e
## alpha.

rppe = function (n, cuts, levels, alpha) 
{
  x <- runif(n)
  qppe(p = x, cuts = cuts, levels = levels, alpha = alpha)
}



##------------------------------------------------------------------------------------------------------

## ---
## Funcao para realizar as particoes (grid) no suporte do tempo de falha
## --

time.grid.obs.t <- function(time, event, n.int=NULL)
{
  o <- order(time)  ## ordem dos tempos de falha 
  time <- time[o]   ## ordena o vetor do tempo de falha
  event <- event[o] ## ordena o vator do indicador de falha
  time.aux <- unique(time[event==1])
  if(is.null(n.int))
  {
    n.int <- length(time.aux)
  }
  
  m <- length(time.aux)
  if(n.int > m)
  {
    a <- c(0,unique(time[event==1]))
    a[length(a)] <- Inf
  }
  else
  {
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    a <- c(0,time.aux[idf])
    a[length(a)] <- Inf
  }
  return(a=a)
}





##------------------------------------------------------------------------------------------------------

## ------
## Funcao de verossiilhanca
## Exponencial por Partes (PE)
## (censura a direita)
## ------
## Funcao de Log-Verossimilhanca

loglik.PE = function(par,time,
                      delta,
                      cuts){
  
  b = length(cuts) + 1
  hazards = par[1:b]
  
  F.t =  ppch(q = time, cuts = cuts, levels = hazards)
  s.t = 1 - F.t
  f.t = dpch(x = time, cuts = cuts, levels = hazards)
    
  log.vero = sum((delta*log(f.t)) + ((1 - delta)*log(s.t)))
  
  #return(-1*log_verossimilhanca)
  return(-1*log.vero)
}

loglik.PE2 = function(par,time,
                     delta,
                     cuts){
  
  b = length(cuts) + 1
  hazards = par[1:b]
  
  f.t = PE(time = time, cuts = cuts, levels = hazards, type = 'distribution')
  s.t = PE(time = time, cuts = cuts, levels = hazards, type = 'survival')
  
  log.vero = sum((delta*log(f.t)) + ((1 - delta)*log(s.t)))
  
  #return(-1*log_verossimilhanca)
  return(-1*log.vero)
}





##------------------------------------------------------------------------------------------------------

## ------
## Funcao de verossiilhanca
## Exponencila por Partes Potencia (PPE)
## (censura a direita)
## ------

# loglik.PPE = function(par,time,
#                       delta,
#                       cuts
# ){
#   
#   b = length(cuts) + 1 ## numero de intervalos
#   hazards = par[1:b]             ## taxas de falha dos b intervalos
#   expoente = par[b +1] ## alpha: ultimo elemento dos parametros a ser estimado
#   
#   densidade = PPE(time = time, cuts = cuts,
#                   levels = hazards, alpha = expoente, type = 'density')
#   
#   sobrevivencia = PPE(time = time, cuts = cuts,
#                       levels = hazards, alpha = expoente, type = 'survival')
#   
#   log_vero = sum(delta*log(densidade) + 
#                    (1 - delta)*log(sobrevivencia))
#   
#   #return(-1*log_verossimilhanca)
#   return(-1*log_vero)
# }

loglik.PPE2 = function(par,time,
                      delta,
                      cuts
){
  
  b = length(cuts) + 1 ## numero de intervalos
  hazards = par[1:b]             ## taxas de falha dos b intervalos
  expoente = par[b +1] ## alpha: ultimo elemento dos parametrtos a ser estimado
  
  ## informacoes exponencial por partes (PE)
  F0.t =  ppch(q = time, cuts = cuts, levels = hazards) ## distribuicao
  
  s0.t = 1 - F0.t ## sobrevivencia
  
  f0.t = dpch(x = time, cuts = cuts, levels = hazards) ## densidade
  
  ## Distribuicao acumulada
  F.t = (F0.t)^expoente
  
  densidade = expoente*(F0.t^(expoente - 1)) * f0.t
  sobrevivencia = 1 - F.t
  
  log_vero = sum(delta*log(densidade) + (1 - delta)*log(sobrevivencia))
  
  #return(-1*log_verossimilhanca)
  return(-1*log_vero)
}


##------------------------------------------------------------------------------------------------------

## ------
## Geracao de dados do exponencial por partes potencia
## ------


time <- function(time.obs, lambda.par, alpha.par, cuts, u.unif){
  
  t <- exp(-PPE(time = time.obs, cuts = cuts, 
                levels = lambda.par, alpha = alpha.par, type = "cum_hazard")) - u.unif
  
  return(t)
  
}

## funcao para encontrar a raiz unitaria que e o tempo de sobrevivencia
gen.mepp <- function(n = n, lambda.par, alpha.par, cuts){
  
  values = c()
  
  for (ind in 1:n){
    raiz <- uniroot(time, c(0, 1000), lambda.par=lambda.par, alpha.par=alpha.par, 
                    cuts=cuts, u.unif=runif(1))
    mepp.time <- raiz$root
    values[ind] = mepp.time
  }
  
  
  return(values)
}

# ## teste
# 
# taxas_falha = c(0.8,1,0.2,0.6,0.3,0.6,0.4)
# 
# ## intervalos de escolha
# intervalos = seq(0.1,4, length = 6)



##------------------------------------------------------------------------------------------------------

## ------
## Geracao de dados do exponencial por partes potencial
## sob o efeito de covariaveis
## ------

## Geracao de dados para o Modelo Exponencial por Partes por Potencia com efeito de covariaveis


time.cox <- function(time.obs, lambda.par, alpha.par, cuts, beta.par, x.mat, u.unif){
    
    t <- exp(-PPE(time = time.obs, cuts = cuts, 
                  levels = lambda.par, alpha = alpha.par, type = "cum_hazard")*exp(x.mat%*%beta.par)) - u.unif
    
    return(t)
    
  }
  
## funcao para encontrar a raiz unitaria que e o tempo de sobrevivencia
gen.mepp.cox <- function(n = n, lambda.par, alpha.par, cuts, beta.par, x.mat){
    
    values = c()
    
    for (ind in 1:n){
      raiz <- uniroot(time.cox, c(0, 1000), lambda.par=lambda.par, alpha.par=alpha.par, 
                      cuts=cuts, beta.par=beta.par,
                      x.mat=x.mat[ind,], u.unif=runif(1))
      mepp.time <- raiz$root
      values[ind] = mepp.time
      }
    
    
    return(values)
  }


# taxas = c(1.1, 0.8, 0.5)
# intervalos = c(0.5, 2)
# t_PE = rpch(1000, intervalos, taxas)
# alpha = 2.5
# 
# beta_par = c(0.5, 2.3)
# x1 = rbinom(1000, 1, prob = 0.5)   ## continuo
# x2 = rnorm(1000, mean = 0, sd = 1) ## discreto
# covars = cbind(x1, x2)
# 
# valores = gen.mepp(n = 1000, 
#          lambda.par = taxas, 
#          alpha.par = alpha,
#          cuts = intervalos,
#          beta.par = beta_par,
#          x.mat = covars)

##------------------------------------------------------------------------------------------------------

## ------
## Funcao de maxima verossimilhanca do PPE
## para estimacao de parametros
## ------

## Goal: intervals, hazards, alpha, and betas.

loglik.cox = function(par,tempos,
                  censura,
                  intervalos,
                  covariaveis){
  
  b = length(intervalos) + 1 ## numero de intervalos
  hazards = par[1:b]             ## taxas de falha dos b intervalos
  exp = par[b +1] ## alpha: ultimo elemento dos parametrtos a ser estimado
  
  n.covars = dim(covariaveis)[2]
  betas = par[(b + 2): (b + 1 + n.covars)]
  
  ## informacoes exponencial por partes (PE)
  F0.t = ppch(q = tempos, cuts = intervalos, levels = hazards)
  f0.t = dpch(x = tempos, cuts = intervalos, levels = hazards)
  ## informacoes exponencial por partes Potencia (PPE)
  F1.t = (F0.t)^exp
  f1.t = exp*(F0.t)^(exp-1) * f0.t

  s1.t = 1 - F1.t
  
  h1.t = f1.t/s1.t 
  
  ## efeito das covariaveis
  pred.linear = covariaveis %*% betas
  
  harzard.cox = (h1.t*exp(pred.linear))
  sobrevivencia.cox = (s1.t^(exp(pred.linear)))
  
  ## Likelihood function
  log.verossimilhanca = sum((censura*log(harzard.cox)) + log(sobrevivencia.cox))

  return(-1*log.verossimilhanca)
}


##------------------------------------------------------------------------------------------------------

## ------
## Geracao de dados do exponencial por partes potencia
## para dados com censura intervalar
## ------

# sim.ICdata <- function(n=n, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat,
#                        lambda.parc=lambda.parc){

sim.ICdata <- function(n, lambda.param, alpha.param, grid.vector, beta.param, x.matrix,lambda.cens.param){
  
  t <- rexp(0,n) # tempos de falha
  c <- rexp(0,n) # tempos de censura
  cens <- rep(1,n) # variiavel indicadora
  tempo <- vector(length=n) # tempo observado
  
  
  t = gen.mepp.cox(n = n, lambda.par = lambda.param, alpha.par = alpha.param,
                   cuts = grid.vector, beta.par = beta.param, x.mat = x.matrix) # tempo de falha
  
  c = rexp(n, rate = lambda.cens.param) # censura
  
  tempo = pmin(t, c)
  
  #-- Observed times:
  
  delta <- ifelse(t<= c, 1, 0)
  
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
  
  dados <- data.frame(L, R, tempo, delta, x.matrix)
  
  return(dados)
}

##------------------------------------------------------------------------------------------------------

## ------
## Funcao de verosimilhanca para censura intervalar
## ------

loglik.int = function(par, time.r, time.l,
                          delta, gamma, cuts, x.matrix){
  
  b = length(cuts) +1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia
  
  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas = par[(b + 2): (b + 1 + n.covars)]

  ## informacoes exponencial por partes (PE) para esquerda
  F0.tl = ppch(q = time.l, cuts = cuts, levels = hazards)
  
  ## informacoes exponencial por partes Potencia (PPE)
  F1.tl = (F0.tl)^exp
  
  s1.tl = 1 - F1.tl
  
  ## informacoes exponencial por partes (PE) para direita
  F0.tr = ppch(q = time.r, cuts = cuts, levels = hazards)
  
  ## informacoes exponencial por partes Potencia (PPE)
  F1.tr = (F0.tr)^exp
  
  s1.tr = 1 - F1.tr

  ## contribuicao das covariaveis
  pred.linear =  x.matrix %*% betas 
  
  s1.tl.cox = s1.tl^(pred.linear)
  
  s1.tr.cox = s1.tr^(pred.linear)
  
  
  log.vero = delta(s1.tr^(pred.linear)) +
    gamma*(s1.tr^(pred.linear) - s1.tr^(pred.linear)) +
    (1-delta-gamma)*(1 - s1.tl^(pred.linear))
    
  return(-1*log.vero)
}

# 
# tamanho.amostral = 128
# 
# ## parametros do PPE
# taxas.falha = c(0.2, 0.4, 0.8)
# particoes = c(0.5, 0.9)
# potencia = 1.4
# lambda.cen = 0.9
# 
# ## pesos das covariaveis
# 
# beta = c(0.5, 2.3)
# x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
# x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
# x.matriz = as.matrix(cbind(x1, x2))
# 
# 
# 
# data.int.cen = sim.ICdata(n = tamanho.amostral, lambda.param = taxas.falha,
#            alpha.param = potencia, grid.vector = particoes,
#            beta.param = beta , x.matrix = x.matriz,
#            lambda.cens.param = lambda.cen)
# 
# ## criacao da coluna indicadora para a censura a direita
# delta = ifelse(data.int.cen$R == Inf, 1, 0)
# 
# 
# ## criacao da coluna indicadora para o intervalo
# gamma = ifelse(data.int.cen$R != Inf & data.int.cen$L != 0, 1, 0)
# 
# 
# estimacao.int.cox = optim(par = chutes,
#                             fn = loglik.cox,
#                             gr = NULL,
#                             hessian = TRUE,
#                             method = "BFGS",
#                             tempos = tempo,
#                             censura = delta,
#                             intervalos = grid,
#                             covariaveis = x.matriz)






