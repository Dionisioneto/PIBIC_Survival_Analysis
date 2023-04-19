## -----
## Simulacao Monte Carlo para a validacao da estimacao do PPE
## -----

## Pacotes

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha, dplyr, readxl)

## ajustar o script com as funcoes de sobrevivencia 
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')


## funcao para calcular o bias%
bias = function(est.matrix, param.matrix){
  bias.calculation = ((est.matrix - param.matrix)/param.matrix)
  return(colMeans(bias.calculation)*100)
}


## ------
## Realizacao do experimento
## ------

set.seed(10)
## ajustando a parametrizacao 
tamanho.amostral = 600

## parametros do PPE
taxas.falha = c(0.2, 0.4, 0.8)
particoes = c(0.5, 0.9)
potencia = 1.4

## pesos das covariaveis

beta = c(0.5, 2.3)
x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
x.matriz = as.matrix(cbind(x1, x2))

# numero de iteracoes a acontecer
n.iter = 500
iteracao = 1 ## iniciador do laco while

## armazenamento
matrix.iter = matrix(data = 0, nrow = n.iter, 
                     ncol = length(taxas.falha) + length(potencia) + length(beta))

matrix.ep = matrix(data = 0, nrow = n.iter, 
                    ncol = length(taxas.falha) + length(potencia) + length(beta))

## particao do banco de dados originais
col.taxas.de.falha = 1:length(taxas.falha)  
col.potencia = length(taxas.falha) + 1
col.betas = (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))

## codigo para informar o erro e tentar de novo a simulcao
## para nao quebrar o laco 

iter.error = c()

while (iteracao <= n.iter) {
  result = tryCatch({
    cat("Realizando iteracao: ", iteracao, "/", n.iter, "\n", sep = "")
    
    ## geracao de dados do tempo de falha e tempo com censura
    tempo.falha = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
                               cuts = particoes, x.mat = x.matriz, beta.par = beta)
    
    tempo.censura = rexp(n = tamanho.amostral, rate = 0.9)
    
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
    matrix.iter[iteracao,col.taxas.de.falha] = estimacao.teste.cox$par[col.taxas.de.falha]
    matrix.iter[iteracao,col.potencia] = estimacao.teste.cox$par[col.potencia]
    matrix.iter[iteracao,col.betas] = estimacao.teste.cox$par[col.betas]
    
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
ititer.error[!is.na(iter.error)]


## ------
## Verificacao da matriz Hessiana
## ------

#matrix.ep

## ------
## Verificacao do real com o predito
## ------

## 1. taxas de falha

taxas.falha
colMeans(matrix.iter[,col.taxas.de.falha])
apply(matrix.iter[,col.taxas.de.falha], MARGIN = 2, FUN = "sd")

matrix.taxas.par = matrix(rep(taxas.falha, n.iter), nrow = n.iter, ncol = length(taxas.falha),
                          byrow = T)

bias(matrix.iter[,col.taxas.de.falha],matrix.taxas.par)

## ------
## criacao do coverage probability
## ------

## O covergae probability e um conceito frequentista que busca
## avaliar dentro do total de repeticoes do seu experimento a porcentagem de
## intervalos de confiance que conseguem captar o seu real parametro

## IC(par) = par +- (quantil_normal-padrao)*(diagonal_inverso_hessiano)

## nivel de confianca
conf.level = 0.95

quant = qnorm(p = conf.level + ((1-conf.level)/2), mean = 0, sd = 1)

# ---
# taxas
# ---

sup.int = matrix.iter[,col.taxas.de.falha] + (quant*matrix.ep[,col.taxas.de.falha])
inf.int = matrix.iter[,col.taxas.de.falha] - (quant*matrix.ep[,col.taxas.de.falha])

matrix.taxas.par = matrix(rep(taxas.falha, n.iter), nrow = n.iter, ncol = length(taxas.falha),
                          byrow = T)

## amplitude do intervalo de confianca
# sup.int - inf.int 

## observacao dos 10 primeiros
#sup.int[1:10, 1:3]
#inf.int [1:10 , 1:3]

matrix.taxas.par[1,1:3]

## porcentagem de capturacao do intervalo de confianca para as taxas
colMeans(matrix.taxas.par >= inf.int & matrix.taxas.par <= sup.int)



## 2. potencia

potencia
mean(matrix.iter[,col.potencia])
sd(matrix.iter[,col.potencia])

bias(matrix.iter[,col.potencia], matrix(potencia, nrow = n.iter))

## 3. coeficientes

beta
colMeans(matrix.iter[, col.betas])
apply(matrix.iter[, col.betas], MARGIN = 2, FUN = "sd")

matrix.betas.par = matrix(rep(beta, n.iter), nrow = n.iter, ncol = length(beta),
                          byrow = T)

bias(matrix.iter[,col.betas], matrix.betas.par)

# ---
# potencia
# ---

sup.int.pot = matrix.iter[,col.potencia] + (quant*matrix.ep[,col.potencia])
inf.int.pot = matrix.iter[,col.potencia] - (quant*matrix.ep[,col.potencia])

## porcentagem de capturacao do intervalo de confianca para o potencia
mean(potencia <= sup.int.pot & potencia >= inf.int.pot)

# ---
# lambda
# ---

sup.int.beta = matrix.iter[,col.betas] + (quant*matrix.ep[,col.betas])
inf.int.beta = matrix.iter[,col.betas] - (quant*matrix.ep[,col.betas])

matrix.betas.par = matrix(rep(beta, n.iter), nrow = n.iter, ncol = length(beta),
                          byrow = T)

## porcentagem de capturacao do intervalo de confianca para as taxas
colMeans(matrix.betas.par >= inf.int.beta & matrix.betas.par <= sup.int.beta)


## organizacao em uma matrix

parametros = c(taxas.falha, potencia, beta)

valores.medios = c(colMeans(matrix.iter[,col.taxas.de.falha]),
                   mean(matrix.iter[,col.potencia]),
                   colMeans(matrix.iter[, col.betas]))

desvios.padroes = c(apply(matrix.iter[,col.taxas.de.falha], MARGIN = 2, FUN = "sd"),
                    sd(matrix.iter[,col.potencia]),
                    apply(matrix.iter[, col.betas], MARGIN = 2, FUN = "sd"))

bias.processo = c(bias(matrix.iter[,col.taxas.de.falha],matrix.taxas.par),
                  bias(matrix.iter[,col.potencia], matrix(potencia, nrow = n.iter)),
                  bias(matrix.iter[,col.betas], matrix.betas.par))

cp.processo = c(colMeans(matrix.taxas.par >= inf.int & matrix.taxas.par <= sup.int), 
                mean(potencia <= sup.int.pot & potencia >= inf.int.pot),
                colMeans(matrix.betas.par >= inf.int.beta & matrix.betas.par <= sup.int.beta))

#ensaio1 = cbind(parametros, valores.medios, desvios.padroes, bias.processo, cp.processo)
#ensaio100 = cbind(parametros, valores.medios, desvios.padroes, bias.processo, cp.processo)
#ensaio200 = cbind(parametros, valores.medios, desvios.padroes, bias.processo, cp.processo)
ensaio500 = cbind(parametros, valores.medios, desvios.padroes, bias.processo, cp.processo)
#ensaio1000 = cbind(parametros, valores.medios, desvios.padroes, bias.processo, cp.processo)

array.ensaios = array(0, dim = c(6, 6, 5))

array.ensaios[,,1] = ensaio1
array.ensaios[,,2] = ensaio100
array.ensaios[,,3] = ensaio200
array.ensaios[,,4] = ensaio500
array.ensaios[,,5] = ensaio1000


write.csv(ensaio1000, file = "ensaio1000.csv")











