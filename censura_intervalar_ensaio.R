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

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha)

source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')

## ------
## Funcao geradora de dados para censura intervalar
## sem covariaveis
## ------

tamanho.amostral = 600
lambdas = c(0.3, 1.2, 0.6, 0.8)
alpha = 1.3
grid = c(0.3, 0.8, 1.2)

lambda.cens = 0.2

## funcao geradora de dados para censura intervalar, sem covariaveis
sim.IC <- function(n, lambda.param, alpha.param, grid.vector, lambda.cens.param){
  
  t <- rexp(0,n) # tempos de falha
  c <- rexp(0,n) # tempos de censura
  cens <- rep(1,n) # variiavel indicadora
  tempo <- vector(length=n) # tempo observado
  
  
  t = gen.mepp(n = n, lambda.par = lambda.param, alpha.par = alpha.param,
                   cuts = grid.vector) # tempo de falha
  
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
  
  dados <- data.frame(L, R, tempo, delta)
  
  return(dados)
}

dados.int = sim.IC(n = tamanho.amostral, lambda.param = lambdas, grid.vector = grid,
                       alpha.param = alpha,
                       lambda.cens.param = lambda.cens)

prop.table(table(dados.int$delta))
head(dados.int)

## ------
## Funcao de verosimilhanca para censura intervalar, sem covariaveis
## ------

loglik.int = function(par, time.r, time.l,
                      grid,
                      delta){
  
  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia
  
  like = rep(0, length(delta))  ## armazenamento de informacao da verossimilhanca
  
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  
  like[delta==1] = (s0.tl - s0.tr)
  
  ## contribuicao da censura
  s0.tl = PPE(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  like[delta==0] = s0.tl 
  
  log.vero = sum(log(like))
  
  return(-1*log.vero)
}

## escolha dos grids, com tempos observaveis
grid = time.grid.interval(li = dados.int$L, ri = dados.int$R, 
                          type = "OBS", bmax = length(lambdas))

grid = grid[-c(1, length(grid))]

grid2 = c(0.4, 0.7, 1.1)

## estudo dos grids
table(cut(dados.int$L,c(0.3, 0.8, 1.2))) ## com os grids parametros reais
table(cut(dados.int$L,grid)) ## com os grids estimados pela observacao
table(cut(dados.int$L,grid2)) ## com os grids estimados pela observacao

table(cut(dados.int$R,c(0.3, 0.8, 1.2))) ## com os grids parametros reais
table(cut(dados.int$R,grid)) ## com os grids estimados pela observacao
table(cut(dados.int$R,grid2)) ## com os grids estimados pela observacao


chutes = c(rep(0.5,length(grid)+1),1.3)

## Metodo numerico BFGS
estimacao.intervalar = optim(par = chutes,
                          fn = loglik.int,
                          gr = NULL,
                          hessian = TRUE,
                          method = "BFGS",
                          time.l = dados.int$L, 
                          time.r = dados.int$R,
                          grid = c(0.3, 0.8, 1.2),
                          delta = dados.int$delta)

estimacao.intervalar$par
c(lambdas, alpha)



## ------
## Funcao geradora de dados para censura intervalar
## com covariaveis
## ------

tamanho.amostral = 600
lambdas = c(0.3, 0.4, 0.9, 1.2)
alpha = 1.4
grid = c(0.4, 0.6, 1.2)

x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1) ## Normal, continua
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5) ## Bernoulli, discreta
x.matriz = as.matrix(cbind(x1, x2))
betas = c(2.2, 0.5)

lambda.cens = 1.5

dados.int = sim.ICdata(n = tamanho.amostral, lambda.param = lambdas, grid.vector = grid,
                       alpha.param = alpha, x.matrix = x.matriz, beta.param = betas,
                       lambda.cens.param = lambda.cens)

prop.table(table(dados.int$delta))
head(dados.int)


## ------
## Funcao de verosimilhanca para censura intervalar, com covariaveis
## ------

loglik.int = function(par, time.r, time.l,
                      grid,
                      delta, x.matrix){
  
  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia
  
  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas = par[(b + 2):(b + 1 + n.covars)]
  
  like = rep(0, dim(x.matrix)[1])  ## armazenamento de informacao da verossimilhanca
  
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## contribuicao do evento
  
  sl = s0.tl^(x.matrix[delta==1,] %*% betas)
  
  sr = s0.tr^(x.matrix[delta==1,] %*% betas)
  
  like[delta==1] = (sl - sr)
  
  ## contribuicao da censura
  s0.tl = PPE(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^(x.matrix[delta==0,] %*% betas)
  
  like[delta==0] = sl
  
  log.vero = sum(log(like))
  
  return(-1*log.vero)
}

head(dados.int)


grid = time.grid.interval(li = dados.int$L, ri = dados.int$R, 
                          type = "OBS", bmax = length(lambdas))

grid = grid[-c(1, length(grid))]

chutes = c(rep(1,length(grid)+1),1.5,1,2)

## Metodo numerico BFGS
estimacao.int.cox = optim(par = chutes,
                          fn = loglik.int,
                          gr = NULL,
                          hessian = TRUE,
                          method = "BFGS",
                          time.l = dados.int$L, 
                          time.r = dados.int$R,
                          grid = grid,
                          delta = dados.int$delta,
                          x.matrix = cbind(dados.int$x1, dados.int$x2))



## debug da funcao de log-verossimilhanca


  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas = par[(b + 2):(b + 1 + n.covars)]
  
  hazards = c(0.3, 0.4, 0.9, 1.2)
  exp = 1.2
  betas = c(2.2, 0.5)
  
  grid = time.grid.interval(li = dados.int$L, ri = dados.int$R, 
                            type = "OBS", bmax = length(lambdas))
  grid = grid[-c(1, length(grid))]
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = dados.int$L[dados.int$delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = dados.int$R[dados.int$delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  
  ## contribuicao do evento
  x.matrix = cbind(dados.int$x1, dados.int$x2)
  
  sl = s0.tl^(x.matrix[dados.int$delta==1,] %*% betas)
  
  sr = s0.tr^(x.matrix[dados.int$delta==1,] %*% betas)
  
  like = rep(0, dim(x.matrix)[1]) ## valores para a verossimilhances
  
  like[dados.int$delta==1] = (sl - sr)
  
  ## contribuicao da censura
  s0.tl = PPE(time = dados.int$L[dados.int$delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^(x.matrix[dados.int$delta==0,] %*% betas)
  
  like[dados.int$delta==0] = sl
  
  log.vero = sum(log(like))
  







for (c in seq(0.05,3, 0.05)){
  cat("lambda cens", c, "\n")
  
  tamanho.amostral = 500
  lambdas = c(0.3, 0.6, 0.9, 1.2)
  alpha = 1.2
  grid = c(0.4, 0.6, 1.2)
  
  x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1) ## Normal, continua
  x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5) ## Bernoulli, discreta
  x.matriz = as.matrix(cbind(x1, x2))
  betas = c(1.5, 0.5)
  
  lambda.cens = c
  
  dados.int = sim.ICdata(n = tamanho.amostral, lambda.param = lambdas, grid.vector = grid,
                         alpha.param = alpha, x.matrix = x.matriz, beta.param = betas,
                         lambda.cens.param = lambda.cens)
  
  print(prop.table(table(dados.int$delta)))
  
  grid = time.grid.interval(li = dados.int$L, ri = dados.int$R, 
                            type = "OBS", bmax = length(lambdas))
  grid = grid[-c(1, length(grid))]
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = dados.int$L[dados.int$delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = dados.int$R[dados.int$delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  
  ## contribuicao do evento
  x.matrix = cbind(dados.int$x1, dados.int$x2)
  
  sl = s0.tl^(x.matrix[dados.int$delta==1,] %*% betas)
  
  sr = s0.tr^(x.matrix[dados.int$delta==1,] %*% betas)
  
  like = rep(0, dim(x.matrix)[1]) ## valores para a verossimilhances
  
  like[dados.int$delta==1] = (sl - sr)
  
  ## contribuicao da censura
  s0.tl = PPE(time = dados.int$L[dados.int$delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^(x.matrix[dados.int$delta==0,] %*% betas)
  
  like[dados.int$delta==0] = sl
  
  log.vero = sum(log(like), na.rm = T)
  
  cat("quantidade de validos da loglik", length(log(like))-sum(is.na(log(like))), "\n")
  cat("quantidade de invalidos da loglik", sum(is.na(log(like))), "\n", "\n")
}

## valores razoaveis para o percentual de censura
## c(4, 4.25, 5, 5.25, 5.5)
## a partir de 6 ta muito bom
## em 8.5 temos zero quantidades de na


head(dados.int)



lpowpch_IC <- function(a, l=l, r=r, x.mat=x.mat, grid.vet=grid.vet){
  
  npar <- length(a)
  n.int <- length(grid.vet)+1
  elinpred <- as.numeric(exp(x.mat%*%a[(n.int+2):npar]))
  n.sample <- nrow(x.mat)
  cens <- ifelse(is.finite(r), 1, 0)
  lik <- rep(0, n.sample)
  
  lambda <- a[1:n.int]
  alpha.par <- a[n.int+1]
  
  p2 <- (1-(ppch(q=l[cens==1], cuts = grid.vet, levels = lambda )^alpha.par))^(elinpred[cens==1])
  p1 <- (1-(ppch(q=r[cens==1], cuts = grid.vet, levels = lambda )^alpha.par ))^(elinpred[cens==1])
  lik[cens==1] <- p2-p1
  p1 <- (1-(ppch(q=l[cens==0], cuts = grid.vet, levels = lambda )^alpha.par ))^(elinpred[cens==0])
  lik[cens==0] <- p1
  return(sum(log(lik), na.rm = T))
  
}


## ------
## Funcao geradora de dados para censura intervalar e fracao de cura,
## com covariaveis
## ------

## Colunas
# B: v.a Bernoulli para indicar curado e nao curado
# L: limite inferior do intervalo observado
# R: limite superior do intervalo observado
# delta: indicador de censura
# x1: valores para a primeira variavel continua (Normal padrao)
# x2: valores para a segunda variavel continua (Bernoulli)

pi = exp(x.matriz %*% betas)/(1 + exp(x.matriz %*% betas))

## ------
## Funcao de verosimilhanca para censura intervalar e fracao de cura, 
## com covariaveis
## ------

## informacoes exponencial por partes Potencia (PPE) para esquerda
s0.tl = PPE(time = time.l, cuts = grid.l, levels = hazards, alpha = exp, type = "survival")

## informacoes exponencial por partes Potencia (PPE) para direita
s0.tr = PPE(time = time.r, cuts = grid.r, levels = hazards, alpha = exp, type = "survival")

## sob a condicao de riscos proporcionais (Modelo de Cox)

sl = s0.tl^(exp(x.matriz %*% betas)) 
sr = s0.tr^(exp(x.matriz %*% betas)) 


loglik = sum(((1 - B) * log(1 - pi)) + (B * log(pi)) + (1 - B)*((delta)*(sl - sr) + (1 - delta)*(sl - sr)))



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




# loglik.int = function(par, time.r, time.l,
#                       delta, alpha, cuts, x.matrix){
#   
#   b = length(cuts) +1 ## numero de intervalos
#   hazards = par[1:b] ## taxas de falha para os b intervalos
#   exp = par[b + 1] ## parametro de potencia
#   
#   n.covars = dim(x.matrix)[2] ## numero de covariaveis
#   betas = par[(b + 2): (b + 1 + n.covars)]
#   
#   ## informacoes exponencial por partes (PE) para esquerda
#   F0.tl = ppch(q = time.l, cuts = cuts, levels = hazards)
#   
#   ## informacoes exponencial por partes Potencia (PPE) para esquerda
#   F1.tl = (F0.tl)^exp
#   
#   s1.tl = 1 - F1.tl
#   
#   ## informacoes exponencial por partes (PE) para direita
#   F0.tr = ppch(q = time.r, cuts = cuts, levels = hazards)
#   
#   ## informacoes exponencial por partes Potencia (PPE) para direita
#   F1.tr = (F0.tr)^exp
#   
#   s1.tr = 1 - F1.tr
#   
#   ## contribuicao das covariaveis
#   pred.linear =  x.matrix %*% betas 
#   
#   sl = s1.tl^(pred.linear)
#   
#   sr = s1.tr^(pred.linear)
#   
#   
#   log.vero = sum(delta*(sl - sr) + (1-delta)*(sl))
#   
#   return(-1*log.vero)
# }
