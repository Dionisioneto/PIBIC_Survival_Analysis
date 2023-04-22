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

tamanho.amostral=1000

survival.time = rep(0, tamanho.amostral)
delta = rep(0, tamanho.amostral)

## Colunas
# B: v.a Bernoulli para indicar curado e nao curado
# L: limite inferior do intervalo observado
# R: limite superior do intervalo observado
# delta: indicador de censura
# x1: valores para a primeira variavel continua (Normal padrao)
# x2: valores para a segunda variavel continua (Bernoulli)

pi = exp(x.matriz %*% betas)/(1 + exp(x.matriz %*% betas))

## Passo 1: geracao de uma variavel aleatoria bernoulli, sob o modelo
## logistico com covariaveis

set.seed(10)

x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
matriz = cbind(rep(1, tamanho.amostral), x1, x2)

betas = c(1.2, 1.2, 0)

z = betas[1] + betas[2]*x1 + betas[3]*x2
  
prob.z = 1/(1+exp(-z))
Y = rbinom(tamanho.amostral,1,prob.z)


## Passo 2: Proporcao de censuras

c1 = 0.46
c2 = 0.32

A = rexp(n = tamanho.amostral, rate = 1)

C = pmin(rep(c1,tamanho.amostral), c2 * A)

## Passo 3: Se Y = 0, entao T = C e adota-se delta = 0.

survival.time[Y == 0] = C[Y==0]
delta[Y == 0] = rep(0, length(Y[Y==0]))


## Passo 4: Se Y = 1, gera-se T por meio de uma distribuicao de risco condicional
## Geracao de dados para o Modelo Exponencial por Partes por Potencia com efeito de covariaveis

betas = c(1.2, 0, 0.5)
t = exp(matriz%*%betas)

ifelse(t <= C,1,0)

delta[Y == 1] = rep(0, length(Y[Y==1]))





## ------
## Funcao de verosimilhanca para censura intervalar e fracao de cura, 
## com covariaveis
## ------


loglik.int.fc = function(par, time.l, time.r,
                         grid, delta, x.matriz){
  
  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia
  
  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas = par[(b + 2):(b + 1 + n.covars)]
  
  
  ## calculo da populacao curada
  
  pi = exp(x.matriz %*% betas)/(1 + exp(x.matriz %*% betas))
  
  ## calculo da populacao nao-curada
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## contribuicao das covariaveis
  
  sl = s0.tl^exp((x.matrix[delta==1,] %*% betas))
  
  sr = s0.tr^exp((x.matrix[delta==1,] %*% betas))
  
  
  ## contribuicao da censura
  s0.tl = PPE(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^exp((x.matrix[delta==0,] %*% betas))
  
  
  p1 = delta*(1-pi)*(sl - sr)
  p2 = (1-delta)*(pi + ((1-pi)*sl))
  
  log.vero = log(p1 + p2)
  
  return(-1*log.vero)
  
}





















