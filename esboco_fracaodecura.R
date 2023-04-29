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

## ------
## Funcao geradora de dados para censura intervalar e fracao de cura,
## com covariaveis
## ------

# # Colunas
# Y: v.a Bernoulli para indicar curado e nao curado
# L: limite inferior do intervalo observado
# R: limite superior do intervalo observado
# delta: indicador de censura
# x1: valores para a primeira variavel continua (Normal padrao)
# x2: valores para a segunda variavel continua (Bernoulli)

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
                    delta = delta, risk = Y,
                    X1 = matriz[,2], X2 = matriz[,3])
  return(data)
}

dados = sim.cure.icen(N = 100, prob.ber = 0.5, betas.cure = c(1,0.5,1.2),
              betas.risk = c(0.5,0.3), c1 = 4, c2 = 7)

#library(survival)

dados.cens = Surv(dados$time, dados$delta)
km = survfit(dados.cens~1)

plot(km)

## ---
## Geracao de dados com o tempo de falha sobre o MEPP
## ---

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
                    delta = delta, risk = Y,
                    X1 = matriz[,2], X2 = matriz[,3])
  return(data)
}

dados = sim.cure.icen(N = 100, prob.ber = 0.5, betas.cure = c(1,0.5,1.2),
                      betas.risk = c(0.5,0.3), c1 = 4, c2 = 7)



## ------
## Funcao de verosimilhanca para censura intervalar e fracao de cura,
## com covariaveis
## ------

## parametros a serem retornados:
## taxas de falha(lambdas), potencia, betas.cura, betas.risco

loglik.int.fc = function(par, time.l, time.r,
                         grid, delta, x.matriz){

  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia

  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas.cura = par[(b + 2):(b + 1 + n.covars)]
  betas.risco = par[(b + 5):(b + 4 +(n.covars-1))]
  
  ## calculo da populacao curada

  pi = exp(x.matriz %*% betas.cura)/(1 + exp(x.matriz %*% betas.cura))

  ## calculo da populacao nao-curada

  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")

  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")

  ## contribuicao das covariaveis

  sl = s0.tl^exp((x.matrix[delta==1,-1] %*% betas.risco))

  sr = s0.tr^exp((x.matrix[delta==1,-1] %*% betas.risco))

  p1 = delta*(1-pi)*(sl - sr)

  ## contribuicao da censura
  s0.tl = PPE(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^exp((x.matrix[delta==0,] %*% betas.risco))

  p2 = (1-delta)*(pi + ((1-pi)*sl))

  log.vero = log(p1 + p2)

  return(-1*log.vero)

}

## taxas de falha(lambdas), potencia, betas.cura, betas.risco

lambdas = c(0.5, 0.8, 0.9)
particoes = c(0.6, 1.2)




# 
# ## funcao de geracao de dados com censura intervalar e fracao de cura
# sim_bch <- function(N, theta = c(1, 0.5, 0),
#                     lambda = 1, A = 5, B = 15, prob = 0.5){
# 
#   ## geracao de dados fracao de cura
#   intercept <- 1
#   xi1 <- stats::rbinom(N, 1, prob)
#   xi2 <- stats::rnorm(N)
#   X <- data.frame(intercept, xi1, xi2)
# 
#   lambda_pois <- as.numeric(exp(theta %*% t(X)))
#   N_L <- stats::rpois(N, lambda_pois)
#   a <- stats::rexp(N)
#   C <- cbind(A,a * B)
#   C <- C[,1] * (C[,1] <= C[,2]) + C[,2] * (C[,1] > C[,2])
#   T <- c(1:N) * NA
# 
#   for (i in 1:N) {
#     T[i] <- ifelse(N_L[i] > 0, min(stats::rexp(N_L[i],
#                                                lambda)), Inf)
#   }
# 
#   delta <- ifelse(T <= C,1,0)
#   Z <- ifelse(T <= C, T, C)
#   L <- R <- Z * NA
#   for (i in 1:N){
#     if (delta[i] == 0){
#       L[i] <- Z[i]
#       R[i] <- Inf
#     }
#     else{
#       L[i] <- 0
#       add  <- stats::runif(1, 0.1, 0.5)
#       R[i] <- add
#       check <- (L[i] <= Z[i] & Z[i] < R[i])
#       while(!check){
#         L[i] <- L[i] + add
#         add <- stats::runif(1, 0.1, 0.5)
#         R[i] <- R[i] + add
#         check <- (L[i] <= Z[i] & Z[i] < R[i])
#       }
#     }
#   }
#   dados <- data.frame(Z, L, R, delta, xi1, xi2, N_latent = N_L)
#   return(dados)
# }
# 
# ## especificacoes amostrais
# N = 100
# prob.ber = 0.5
# 
# betas = c(1.2, 0.5, 1.5)
# A = 12
# B = 17
# lambda = 1
# 
# ## geracao de dados fracao de cura
# intercept <- 1
# xi1 <- stats::rbinom(N, 1, prob.ber)
# xi2 <- stats::rnorm(N)
# X <- data.frame(intercept, xi1, xi2)
# 
# lambda_pois <- as.numeric(exp(betas %*% t(X)))
# 
# N_L <- stats::rpois(N, lambda_pois)
# 
# a <- stats::rexp(N)
# C <- cbind(A,a * B)
# C <- pmin(C[,1], C[,2])
# 
# 
# tempo <- c(1:N) * NA
# 
# for(unidade in 1:N){
#   tempo[unidade] <- ifelse(N_L[unidade] >= 0, min(stats::rexp(N_L[unidade],lambda)), Inf)
# }
# 
# delta = ifelse(tempo <= C, 1,0)
# 
# Z.tempo = ifelse(tempo <= C, tempo, C)
# 
# L = R = Z.tempo*NA
# 
# for(unidade in 1:N){
#   if (delta[unidade] == 0){
#     L[unidade] = Z.tempo[unidade]
#     R[unidade] = Inf
#   }
#   else{
#     L[unidade] = 0
#     Qj = stats::runif(1, 0.1, 0.5)
#     R[unidade] = Qj
#     verifica = (L[unidade] <= Z.tempo[unidade] & Z.tempo[unidade] < R[unidade])
#   }
#   while(!check){
#     L[unidade] = L[unidade] + Qj
#     Qj = stats::runif(1, 0.1, 0.5)
#     R[unidade] = R[unidade] + Qj
#     verifica = (L[unidade] <= Z.tempo[unidade] & Z.tempo[unidade] < R[unidade])
#   }
# }
# 
# dados = data.frame()
# 
# for (i in 1:N){
#   if (delta[i] == 0){
#     L[i] <- Z[i]
#     R[i] <- Inf
#   }
#   else{
#     L[i] <- 0
#     add  <- stats::runif(1, 0.1, 0.5)
#     R[i] <- add
#     check <- (L[i] <= Z[i] & Z[i] < R[i])
#     while(!check){
#       L[i] <- L[i] + add
#       add <- stats::runif(1, 0.1, 0.5)
#       R[i] <- R[i] + add
#       check <- (L[i] <= Z[i] & Z[i] < R[i])
#     }
#   }
# }
# dados <- data.frame(Z, L, R, delta, xi1, xi2, N_latent = N_L)
# return(dados)
















