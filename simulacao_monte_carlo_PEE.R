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
n.iter = 50 


## armazenamento
matrix.iter = matrix(data = 0, nrow = n.iter, 
                     ncol = length(taxas.falha) + length(potencia) + length(beta))

matrix.var = matrix(data = 0, nrow = n.iter, 
                    ncol = length(taxas.falha) + length(potencia) + length(beta))


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
  matrix.iter[i,1:length(taxas.falha)] = estimacao.teste.cox$par[1:length(taxas.falha)]
  matrix.iter[i,length(taxas.falha) + 1] = estimacao.teste.cox$par[length(taxas.falha) +1]
  matrix.iter[i, (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))] = estimacao.teste.cox$par[(length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))]
  
  matrix.var[i,] = diag(solve(-estimacao.teste.cox$hessian))
}

## funcao para calcular o bias%
bias = function(est.matrix, param.matrix){
  bias.calculation = ((est.matrix - param.matrix)/param.matrix)*100
  return(bias.calculation)
}


## 1. taxas de falha

taxas.falha
colMeans(matrix.iter[,1:length(taxas.falha)])
apply(matrix.iter[,1:length(taxas.falha)], MARGIN = 2, FUN = "sd")

bias(colMeans(matrix.iter[,1:length(taxas.falha)]),taxas.falha)

## 2. potencia

potencia
mean(matrix.iter[,length(taxas.falha) + 1])
sd(matrix.iter[,length(taxas.falha) + 1])
#bias(mean(potencia.iter), potencia)

## 3. coeficientes

beta
colMeans(matrix.iter[, (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))])
apply(matrix.iter[, (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))], MARGIN = 2, FUN = "sd")
bias(colMeans(matrix.iter[, (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))]), beta)

## codigo para informar o erro e tentar de novo a simulcao
## paranao quebrar o laco 

n.iter = 10

iter.error = c()

iteracao = 1

while (iteracao <= n.iter) {
  result = tryCatch({
    cat("Realizando iteracao: ", iteracao, "/", n.iter, "\n", sep = "")
    
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
    
    chutes = c(rep(4,length(grid)+1),1,0.5,1)
    
    
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
    matrix.iter[iteracao,1:length(taxas.falha)] = estimacao.teste.cox$par[1:length(taxas.falha)]
    matrix.iter[iteracao,length(taxas.falha) + 1] = estimacao.teste.cox$par[length(taxas.falha) +1]
    matrix.iter[iteracao, (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))] = estimacao.teste.cox$par[(length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))]
    
    matrix.var[iteracao,] = diag(solve(-estimacao.teste.cox$hessian))
    
    iteracao = iteracao + 1
    
    estimacao.teste.cox
    
  }, error = function(e){
    
    print(paste0("Erro na iteracao ", iteracao, ": ", conditionMessage(e)))
    # valor nulo nessa iteracao
    NULL
  }
  )
  
  # continua as iteracoes se tiver um erro
  if(is.null(result$convergence)){
    iter.error[iteracao] = as.character(iteracao)
    #iteracao = iteracao
  }
  
  # caso conttrario, com com os resultados
}

## retirando os NA
iter.error[!is.na(iter.error)]




