## ------
## Funcoes PIBIC 2023
## Projeto: Modelos Semi-Parametricos para dados
## de sobrevivencia com censura intervalar
##
## Geracao de dados para censura intervalar e fracao de cura
## ------


## ---
## Pacotes 
## ---

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha)

source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')

n.sample = 1000

##distribuicao dos erros (ruidos)
#u = rnorm(n = n.sample, mean = 0, sd = sqrt(5))
#v = rnorm(n = n.sample, mean = 0, sd = sqrt(1))

## geracao da variavel binaria sob o modelo de regressao logistica

logistic = function(betas, x.matrix){

  ones = rep(1, dim(x.matrix)[1])
  
  xbetas = cbind(ones,x.matrix) %*% betas 
  
  var.logistic = 1/(1 + exp(-xbetas))
  
  return(var.logistic)
}

## intercepto, bernoulli, normal 
betas.par = c(2.3,-0.5, 0)

X1 = rbinom(n = n.sample, size = 1, prob = 0.5)
X2 = rnorm(n = n.sample, mean = 0, sd = 1)
x.matriz = as.matrix(cbind(X1, X2))

prob.logistic = logistic(betas = betas.par, x.matrix = x.matriz)

prob.logistic[100:107]

## Geracao de uma variavel bernoulli, na qual cada probabilidade e 
## configurada pela proababilidade do modelo logistico

Y = rbinom(n = n.sample, size = 1, prob = prob.logistic)


## geracao dos vetores (L, R, delta)


geracao.dados.fc = function(n.size, var.ber, mean.exp, x.matrix, betas,
                            plato = 0.1){
  if(var.ber == 0){
    t = rexp(n =n.size, rate = 1/mean.exp) ## t  = C
    delta = rep(0, n.sample)               ## delta = 0
  }
  else{ ## var.ber = 1
    xbetas = x.matrix %*% betas
    t = exp(xbetas) 
    C = rexp(n =n.size, rate = 1/mean.exp)
    delta = ifelse(t < C, 1, 0)
  }
  
  if(delta == 0){
    L = L
    R = Inf 
  }
  else{
    m = n.size - length(delta == 0)
    Q = runif(n = m, min = plato, max = 1)
    L = sum(Q[-c(length(m))])
    R = sum(Q)
    }
  return(L, R, delta) 
  }

geracao.dados.fc(n.size = n.sample, var.ber = Y[1], mean.exp = 10,
                  x.matrix = x.matriz[1,], betas = betas.par[2:length(betas.par)])

<<<<<<< HEAD


=======
>>>>>>> 2b399fdf62e9d3c8ad74c6b842611c2f5d20ff87









