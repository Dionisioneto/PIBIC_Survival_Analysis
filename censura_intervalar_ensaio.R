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
alpha = 1.4
grid = c(0.3, 0.8, 1.2)

lambda.cens = 0.45

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

chutes = c(rep(0.1,length(grid)+1),1)

## Metodo numerico BFGS
estimacao.intervalar = optim(par = chutes,
                          fn = loglik.int,
                          gr = NULL,
                          hessian = TRUE,
                          method = "BFGS",
                          time.l = dados.int$L, 
                          time.r = dados.int$R,
                          grid = grid,
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
betas = c(0, 0.5)

lambda.cens = 0.45

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
  
  sl = s0.tl^exp((x.matrix[delta==1,] %*% betas))
  
  sr = s0.tr^exp((x.matrix[delta==1,] %*% betas))
  
  like[delta==1] = (sl - sr)
  
  ## contribuicao da censura
  s0.tl = PPE(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^exp((x.matrix[delta==0,] %*% betas))
  
  like[delta==0] = sl
  
  log.vero = sum(log(like))
  
  return(-1*log.vero)
}

head(dados.int)

# realizacao de um experimento

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

estimacao.int.cox$par
c(lambdas, alpha, betas)

## ---
## Simulacao monte carlo
## ---

## tamanho amostral (n) = (50,100,200,500,1000,5000)

set.seed(10)

## ajustando a parametrizacao 
tamanho.amostral = 50

## parametros do PPE
taxas.falha = c(0.3, 0.4, 0.9, 1.2)
particoes = c(0.3, 0.4, 0.9, 1.2)
potencia = 1.4

## pesos das covariaveis

beta = c(0, 0.5)
x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
x.matriz = as.matrix(cbind(x1, x2))

# numero de iteracoes a acontecer
n.iter = 50
iteracao = 1 ## iniciador do laco while

## ----
## armazenamento
matrix.iter = matrix(data = 0, nrow = n.iter, 
                     ncol = length(taxas.falha) + length(potencia) + length(beta))

matrix.ep = matrix(data = 0, nrow = n.iter, 
                   ncol = length(taxas.falha) + length(potencia) + length(beta))

## particao do banco de dados originais
col.taxas.de.falha = 1:length(taxas.falha)  
col.potencia = length(taxas.falha) + 1
col.betas = (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))

## ----

## codigo para informar o erro e tentar de novo a simulcao
## para nao quebrar o laco 

iter.error = c()

while (iteracao <= n.iter) {
  result = tryCatch({
    cat("Realizando iteracao: ", iteracao, "/", n.iter, "\n", sep = "")
    
    ## geracao de dados do tempo de falha com censura intervalar
    dados.int = sim.ICdata(n = tamanho.amostral, lambda.param = lambdas, grid.vector = grid,
                           alpha.param = alpha, x.matrix = x.matriz, beta.param = betas,
                           lambda.cens.param = lambda.cens)
    
    ## particao para a estimacao
    grid = time.grid.interval(li = dados.int$L, ri = dados.int$R, 
                              type = "OBS", bmax = length(lambdas))
    
    grid = grid[-c(1, length(grid))]
    
    chutes = c(rep(1,length(grid)+1),1.5,1,2)
    

    ## Metodo numerico BFGS
    estimacao.teste.cox = optim(par = chutes,
                                fn = loglik.int,
                                gr = NULL,
                                hessian = TRUE,
                                method = "BFGS",
                                time.l = dados.int$L, 
                                time.r = dados.int$R,
                                grid = grid,
                                delta = dados.int$delta,
                                x.matrix = cbind(dados.int$x1, dados.int$x2))

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
  }
  
  # caso contrario, com os resultados
}

set.seed(123) 

# numero de iteracoes a acontecer
n.iter = 20
iteracao = 1 ## iniciador do laco while

## ----
## armazenamento
## ----

matrix.iter = matrix(data = 0, nrow = n.iter, 
                     ncol = length(taxas.falha) + length(potencia) + length(beta))

matrix.ep = matrix(data = 0, nrow = n.iter, 
                   ncol = length(taxas.falha) + length(potencia) + length(beta))

## particao do banco de dados originais
col.taxas.de.falha = 1:length(taxas.falha)  
col.potencia = length(taxas.falha) + 1
col.betas = (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))



for (iteracao in 1:n.iter) {
  error_occurred <- FALSE  # initialize the error flag
  repeat {
    tryCatch({
      cat("Realizando iteracao: ", iteracao, "/", n.iter, "\n", sep = "")
      
      ## geracao de dados do tempo de falha com censura intervalar
      dados.int = sim.ICdata(n = tamanho.amostral, lambda.param = lambdas, grid.vector = grid,
                             alpha.param = alpha, x.matrix = x.matriz, beta.param = betas,
                             lambda.cens.param = lambda.cens)
      
      ## particao para a estimacao
      grid = time.grid.interval(li = dados.int$L, ri = dados.int$R, 
                                type = "OBS", bmax = length(lambdas))
      
      grid = grid[-c(1, length(grid))]
      
      chutes = c(rep(1,length(grid)+1),1.5,1,2)
      
      
      ## Metodo numerico BFGS
      estimacao.teste.cox = optim(par = chutes,
                                  fn = loglik.int,
                                  gr = NULL,
                                  hessian = TRUE,
                                  method = "BFGS",
                                  time.l = dados.int$L, 
                                  time.r = dados.int$R,
                                  grid = grid,
                                  delta = dados.int$delta,
                                  x.matrix = cbind(dados.int$x1, dados.int$x2))
      
      if (is.null(estimacao.teste.cox)) {
        stop("A convergência frenquentista BFGS falhou.")
      }
    
      
      if (is.na(sqrt(diag(solve(estimacao.teste.cox$hessian))))) {
        warning("Alguma variância da matriz hessiana foi estimada negativamente")
        error_occurred <- TRUE  # set the error flag to true
      }
      
      ## salvar resultados das estimativas
      matrix.iter[iteracao,col.taxas.de.falha] = estimacao.teste.cox$par[col.taxas.de.falha]
      matrix.iter[iteracao,col.potencia] = estimacao.teste.cox$par[col.potencia]
      matrix.iter[iteracao,col.betas] = estimacao.teste.cox$par[col.betas]
      
      ## salvar resultados da matriz Hessiana (aproximacao da Hessiana de Fisher)
      matrix.ep[iteracao,] = sqrt(diag(solve(estimacao.teste.cox$hessian)))
      
      break ## quebra o loop se nao tiver nenhum erro
    
    }, warning = function(w){
      
      message("Warning: ", w$message)  # print a warning message
      error_occurred <- TRUE
    }, error = function(e) {
      message("Error: ", e$message)
      error_ocurred = TRUE
    }) 
    
    if (error_occurred) {
      break  # exit the inner loop if an error or warning occurred 
      }
    }
  }

      


    



## retirando os NA
iter.error[!is.na(iter.error)]


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

## O coverge probability e um conceito frequentista que busca
## avaliar dentro do total de repeticoes do seu experimento a porcentagem de
## intervalos de confiance que conseguem captar o seu real parametro

## IC(par) = par +- (quantil_normal-padrao)*(diagonal_inverso_hessiano)

## nivel de confianca
conf.level = 0.95

quant = qnorm(p = conf.level + ((1-conf.level)/2), mean = 0, sd = 1)

# ---
# taxas
# ---

sup.int.lambdas = matrix.iter[,col.taxas.de.falha] + (quant*matrix.ep[,col.taxas.de.falha])
inf.int.lambdas = matrix.iter[,col.taxas.de.falha] - (quant*matrix.ep[,col.taxas.de.falha])


colMeans(matrix.taxas.par >= inf.int & matrix.taxas.par <= sup.int.lambdas)

sup.int.pot = matrix.iter[,col.potencia] + (quant*matrix.ep[,col.potencia])
inf.int.pot = matrix.iter[,col.potencia] - (quant*matrix.ep[,col.potencia])

## porcentagem de capturacao do intervalo de confianca para o potencia
mean(potencia <= sup.int.pot & potencia >= inf.int.pot)



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

cbind(parametros, valores.medios, desvios.padroes, bias.processo, cp.processo)

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

