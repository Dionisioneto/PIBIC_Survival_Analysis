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
## com covariaveis
## ------

tamanho.amostral = 300
lambdas = c(0.3, 0.9, 1.2, 0.45)
alpha = 1.2
grid = c(0.8, 1.2, 4.5)

x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1) ## Normal, continua
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5) ## Bernoulli, discreta
x.matriz = as.matrix(cbind(x1, x2))
betas = c(0.5, 2.2)

lambda.cens = 0.9

dados.int = sim.ICdata(n = tamanho.amostral, lambda.param = lambdas, grid.vector = grid,
           alpha.param = alpha, x.matrix = x.matriz, beta.param = betas,
           lambda.cens.param = lambda.cens)

head(dados.int)

## ------
## Funcao de verosimilhanca para censura intervalar, com covariaveis
## ------

loglik.int = function(par, time.r, time.l,
                      grid.l, grid.r,
                      delta, cuts, x.matrix){
  
  b = length(grid.r) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia
  
  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas = par[(b + 2): (b + 1 + n.covars)]
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = time.l, cuts = grid.l, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = time.r, cuts = grid.r, levels = hazards, alpha = exp, type = "survival")
  
  ## contribuicao das covariaveis
  pred.linear =  x.matrix %*% betas 
  
  sl = s0.tl^(pred.linear)
  
  sr = s0.tr^(pred.linear)
  
  log.vero = sum((delta*log(sl - sr)) + ((1-delta)*log(sl)))
  
  return(-1*log.vero)
}

#b = length(part.r ) + 1 ## numero de intervalos
hazards = lambdas ## taxas de falha para os b intervalos
exp = alpha ## parametro de potencia
# 
betas = betas
# 
# ## informacoes exponencial por partes Potencia (PPE) para esquerda
s0.tl = PPE(time = dados.int$L[105:125], cuts = part.l, levels = hazards, alpha = exp, type = "survival")
# 
# ## informacoes exponencial por partes Potencia (PPE) para direita
s0.tr = PPE(time = dados.int$R[105:125], cuts = part.r, levels = hazards, alpha = exp, type = "survival")
# 
# 



#delta[105:125]

pred.linear =  x.matriz[105:125,] %*% betas 

sl = s0.tl^(pred.linear) 

sr = s0.tr^(pred.linear)


lv = sum((dados.int$delta[105:125]*log(sl - sr)) + ((1-dados.int$delta[105:125])*log(sl)))

## ------
## Otimizacao da log-verosimilhanca para censura intervalar, com covariaveis
## ------

# particao para a estimacao
# tentei juntar e tirar o unique, mas iria fica um vetor maior que o numero de taxas

part.l = time.grid.obs.t(time = dados.int$L, event = dados.int$delta, n.int = length(lambdas))
part.r = time.grid.obs.t(time = dados.int$R, event = dados.int$delta, n.int = length(lambdas))

part.l = part.l[-c(1, length(part.l))]
part.r = part.r[-c(1, length(part.r))]

chutes = c(rep(0.1,length(part.l)+1),1,0.5,1.5)

## Metodo numerico BFGS
estimacao.teste.cox = optim(par = chutes,
                            fn = loglik.int,
                            gr = NULL,
                            hessian = TRUE,
                            method = "BFGS",
                            time.l = dados.int$L, 
                            time.r = dados.int$R,
                            grid.l = part.l, 
                            grid.r = part.r,
                            delta = dados.int$delta,
                            cuts = grid,
                            x.matrix = x.matriz)



## Metodo numerico Nelder-Mead
estimacao.teste.cox = optim(par = chutes,
                            fn = loglik.int,
                            gr = NULL,
                            hessian = F,
                            method = "Nelder-Mead",
                            time.l = dados.int$L, 
                            time.r = dados.int$R,
                            grid.l = part.l, 
                            grid.r = part.r,
                            delta = dados.int$delta,
                            cuts = grid,
                            x.matrix = x.matriz)


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





























