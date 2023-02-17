## -----
## Simulacao Monte Carlo para a validacao da estimacao do PPE
## -----

## Pacotes

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha)

## ajustar o script com as funcoes de sobrevivencia 
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')

set.seed(10)
## ajsutando a parametrizacao 
tamanho.amostral = 600

## parametros do PPE
taxas.falha = c(0.2, 0.4, 0.65)
particoes = c(0.5, 0.9)
potencia = 1.4

## pesos das covariaveis

beta = c(0.5, 2.3)
x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
x.matriz = as.matrix(cbind(x1, x2))

# numero de iteracoes a acontecer
n.iter = 500 


## armazenamento
taxas.falha.iter = matrix(data = 0, nrow = n.iter, ncol = length(taxas.falha))
potencia.iter = rep(0, n.iter)
betas.iter = matrix(data = 0, nrow = n.iter, ncol = length(beta))


for (i in 1:n.iter){
  
  cat("Realizando iteracao: ", i, "/", n.iter, "\n", sep = "")
  
  ## geracao de dados do tempo de falha e tempo com censura
  tempo.falha = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
                             cuts = particoes, x.mat = x.matriz, beta.par = beta)
  
  tempo.censura = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
                               cuts = particoes, x.mat = x.matriz, beta.par = beta)
  
  ## geracao do tempo e censura
  tempo = pmin(tempo.falha, tempo.censura)
  delta = ifelse(tempo.falha <= tempo.censura, 1, 0)
 
  
  ## particao para a estimacao
  grid = time.grid.obs.t(tempo, delta, n.int = 3)
  grid = grid[-c(1, length(grid))]
  
  chutes = c(rep(0.5,length(grid)+1),1,0.5,1)
  
  
  ## Metodo numerico BFGS
  estimacao.teste.cox = optim(par = chutes,
                              fn = loglik.cox,
                              gr = NULL,
                              hessian = TRUE,
                              method = "BFGS",
                              tempos = tempo,
                              censura = delta,
                              intervalos = grid,
                              covariaveis = x.matriz)
  
  ## salvar resultados

  taxas.falha.iter[i,] = estimacao.teste.cox$par[1:length(taxas.falha)]
  potencia.iter[i] = estimacao.teste.cox$par[length(taxas.falha) +1]
  betas.iter[i,] = estimacao.teste.cox$par[(length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))]
}



## 1. taxas de falha

taxas.falha
colMeans(taxas.falha.iter)
apply(taxas.falha.iter, MARGIN = 2, FUN = "sd")

## funcao para calcular o bias%
bias = function(est.matrix, param.matrix){
  bias.calculation = ((est.matrix - param.matrix)/param.matrix)*100
  return(bias.calculation)
}

bias(colMeans(taxas.falha.iter),taxas.falha)


## 2. potencia


potencia
mean(potencia.iter)
sd(potencia.iter)
#bias(mean(potencia.iter), potencia)

## 3. coeficientes

beta
colMeans(betas.iter)
apply(betas.iter, MARGIN = 2, FUN = "sd")
bias(colMeans(betas.iter), beta)







