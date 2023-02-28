## -----
## Simulacao Monte Carlo para a validacao da estimacao do PPE
## -----

## Pacotes

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha, dplyr, maxLik)

## ajustar o script com as funcoes de sobrevivencia 
source('C:/Users/dionisio.neto/Desktop/Dionisio_Neto/Survival_Analysis/PIBIC_Survival_Analysis/funcoes_sobrevivencia_pibic2023.R')


## funcao para calcular o bias%
bias = function(est.matrix, param.matrix){
  bias.calculation = ((est.matrix - param.matrix)/param.matrix)*100
  return(bias.calculation)
}

## ------
## Realizacao do experimento
## ------

set.seed(10)
## ajustando a parametrizacao 
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
iteracao = 1 ## iniciador do laco while

## armazenamento
matrix.iter = matrix(data = 0, nrow = n.iter, 
                     ncol = length(taxas.falha) + length(potencia) + length(beta))

matrix.ep = matrix(data = 0, nrow = n.iter, 
                    ncol = length(taxas.falha) + length(potencia) + length(beta))

## codigo para informar o erro e tentar de novo a simulcao
## para nao quebrar o laco 


iter.error = c()

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
    grid = time.grid.obs.t(tempo, delta, n.int = length(taxas.falha))
    grid = grid[-c(1, length(grid))]
    
    chutes = c(rep(0.1,length(grid)+1),1,0.5,1)
    
    
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
    
    matrix.ep[iteracao,] = sqrt(diag(solve(estimacao.teste.cox$hessian)))
    
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
  
  # caso contrario, com os resultados
}

 ## retirando os NA
iter.error[!is.na(iter.error)]


## ------
## Verificacao do real com o predito
## ------

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

## ------
## Verificacao da matriz Hessiana
## ------

#*** pergunta: porque deu tudo negativo?
matrix.ep

## ------
## criacao do coverage probability
## ------


## O covergae probability e um conceito frequentista que busca
## avaliar dentro do total de repeticoes do seu experimento a porcentagem de
## intervalos de confiance que conseguem captar o seu real parametro

## IC(par) = par +- (quantil_t)*(diagonal_inverso_hessiano)

## numero de parametros estimados
n.param.est = length(taxas.falha) + length(potencia) + length(beta) 
  
## nivel de confianca
conf.level = 0.95

t.quant = qt(p = (1 - conf.level)/2, df = tamanho.amostral - n.param.est)*(-1)

sup.int = matrix.iter[,1:length(taxas.falha)] + (t.quant*matrix.ep[,1:length(taxas.falha)])
inf.int = matrix.iter[,1:length(taxas.falha)] - (t.quant*matrix.ep[,1:length(taxas.falha)])

matrix.taxas.par = matrix(rep(taxas.falha, n.iter), nrow = n.iter, ncol = length(taxas.falha),
       byrow = T)

## amplitude do intervalo de confianca
sup.int - inf.int 

## observacao dos 10 primeiros
sup.int[1:10, 1:3]
inf.int [1:10 , 1:3]

matrix.taxas.par[1,1:3]

## porcentagem de capturacao do intervalo de confianca para as taxas
colMeans(matrix.taxas.par >= inf.int & matrix.taxas.par <= sup.int)




##*** alguma coisa deve estar errada

