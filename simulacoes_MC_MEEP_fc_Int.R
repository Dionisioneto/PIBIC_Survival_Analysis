## -----------------
## Simulacoes Monte Carlo
## Modelo Exponencial por Partes Potencia (MEPP)
## para dados com censura intevalar 
## com a presenca de fracao de cura
## -----------------

## ---
## funcoes a serem utilizadas
## ---

source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')

## ---
## funcao para a geracao de dados do MEPP
## ---


time <- function( t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=u.unif ){
  exp(-PPE(time=t,cuts=grid.vet,levels = lambda.par, alpha = alpha.par,type = "cum_hazard")*exp(x.mat%*%beta.par)) - u.unif
}

gen.mepp <- function(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat){
  
  raiz <- uniroot(time, c(0, 10000), lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=runif(1))
  mepp.time <- raiz$root
  
  return(mepp.time)
}


## ---
## funcao para a geracao de dados sob fracao de cura
## ---

sim.std.cure.ICdata <- function(n=n, lambda.par=lambda.par, alpha.par=alpha.par, 
                                grid.vet=grid.vet, beta.par=beta.par, lambda.parc=lambda.parc, 
                                theta.par = c(1, 0.5, 0), A = 5, B = 15){
  
  intercept <- 1
  xi1 <- rbinom(n, 1, 0.5)
  xi2 <- rnorm(n)
  X_cure <- cbind(intercept, xi1, xi2)
  X <- cbind(xi1, xi2)
  
  
  #-- Cure probability mixture model:
  
  elinpred <- exp(-theta.par%*%t(X_cure))
  probY <- 1/(1+elinpred)
  Y <- rbinom(n, size=1, prob=probY)
  
  #-- Censoring times generation:
  
  C <- pmin(A, B*rexp(n, rate=lambda.parc))
  
  #-- Generating the survival times:
  
  t     <- rep(0,n) # tempos de falha
  t_obs <- rep(0,n)
  delta <- rep(0,n)
  
  for(i in 1:n){
    if(Y[i]==1){
      t[i] <- gen.mepp(lambda.par=lambda.par, alpha.par=alpha.par,
                       grid.vet=grid.vet, beta.par=beta.par, x.mat=X[i,])
      
      if(t[i]<C[i]){
        t_obs[i] <- t[i]
        delta[i] <- 1  
      }else{
        t_obs[i] <- C[i]
      }
    }else{
      t[i] <- C[i]
      t_obs[i] <- t[i] 
    }
  }
  
  tempo <- t_obs
  
  #-- Observed times:
  
  delta <- ifelse(t < C, 1, 0) 
  
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
  
  dados <- data.frame(L, R, tempo, Y, delta, xi1, xi2)
  
  return(dados)
}

## ---
## funcao o calculo da funcao de sobrevivencia populacional sob o MEPP
## presenca de dados censurados em intervalos
## ---

SpopMEPP <- function(t=t, lambda.par=lambda.par, alpha.par=alpha.par,
                     grid.vet=grid.vet, beta.par=beta.par, 
                     theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){
  
  
  S_MEPP <- as.numeric(exp(-PPE(time=t,cuts = grid.vet,levels = lambda.par,alpha = alpha.par,type = 'cum_hazard')*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEPP
  return(spop)
}


## ---
## funcao de verossimilanca para o MEPP, com censura intervalar 
## (utiliza-se de SpopMEPP)
## ---

loglikIC <- function(a, l=l, r=r, x.cure=x.cure, x.risk=x.risk, grid.vet=grid.vet){
 
  b <- length(grid.vet)+1
  
  hazards = a[1:b] ## taxas de falha para os b intervalos
  alpha = a[b + 1] ## parametro de potencia
  
  n.cov.cure = dim(x.cure)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  n.cov.risk = dim(x.risk)[2] ## numero de covariaveis com fracao de cura, para risco tiramos um (beta0)
  
  betas.cure = a[(b + 2):(b + 1 + n.cov.cure)]
  betas.risk = a[(b + 5):(b + 4 + n.cov.risk)]
  
  
  n.sample <- nrow(x.cure)
  
  cens <- ifelse(is.finite(r), 1, 0)
  lik <- rep(0, n.sample)
  
  p2 <- SpopMEPP(t=l[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  p1 <- SpopMEPP(t=r[cens==1], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==1,], x.risk=x.risk[cens==1,])
  lik[cens==1] <- p2-p1
  
  p1 <- SpopMEPP(t=l[cens==0], lambda.par=hazards, alpha.par=alpha, grid.vet=grid.vet, beta.par=betas.risk, theta.par=betas.cure, x.cure=x.cure[cens==0,], x.risk=x.risk[cens==0,])
  lik[cens==0] <- p1
  return(sum(log(lik)))
  
}


## ---
## Gerando a parametrizacao para o MEPP
## ---


taxas.de.falha =  c(1.1, 0.8, 0.5)
particoes = c(0.5, 2)
potencia = 0.8

## Pesos das Covariaveis 

betas.cura = c(1.2, 0.5, -0.5)
betas.risco = c(-0.5, 1.2)

## organizando em um dataframe
Theta = as.data.frame(c(taxas.de.falha, potencia, betas.cura, betas.risco))
colnames(Theta) = "Real"
rownames(Theta) = c("lambda1", "lambda2", "lambda3",
                    "potencia",
                    "b0", "b1", "b2",
                    "beta1", "beta2")

###################################################################################


###################################################################################

### ---
### Algotirmo Monte Carlo para as
### estimacoes de Maxima Verossimilhanca
### ---

# salvar resultados
## ensaio50, ensaio200, ensaio500, ensaio1000, ensaio5000


n.amostras = c(50,100,500,1000)
## replicacoes Monte Carlo

amostra = n.amostras[1]
 
iteracao = 1 # iniciador do laco while
n.iter = 10
  
## Organizador das iteracoes Monte Carlo

## armazenamento das iteracoes
matrix.iter = matrix(data = 0, nrow = n.iter, 
                     ncol = nrow(Theta))

matrix.ep = matrix(data = 0, nrow = n.iter, 
                   ncol = nrow(Theta))

## particao do banco de dados originais
n.intervals = length(particoes) + 1

col.taxas.de.falha = 1:n.intervals  
col.potencia = n.intervals + 1
col.betas.cura = (n.intervals + 2):(n.intervals + 1 + length(betas.cura))
col.betas.risco = (n.intervals + 2 + length(betas.cura)):(n.intervals + 1 +  length(betas.cura) + length(betas.risco))

## contagens nas quais nao houve a convergencia 
## pelo metodo numerico BFGS

iter.error = c()

while(iteracao <= n.iter) {
  result = tryCatch({
    
    data = sim.std.cure.ICdata(n=amostra, lambda.par=taxas.de.falha, alpha.par=potencia, 
                               grid.vet=particoes, beta.par=betas.risco, lambda.parc=1, 
                               theta.par=betas.cura, A = 5, B = 15)
    
    covariaveis.cura = cbind(1, data$xi1, data$xi2)
    covariaveis.risco = cbind(data$xi1, data$xi2)
    
    grid.observado = time.grid.obs.t(data$tempo, data$delta, n.int = length(taxas.de.falha))
    grid.observado = grid.observado[-c(1, length(grid.observado))]
    
    chutes = c(rep(0.1, length(taxas.de.falha)), 1, 1, 0.5, 0.5, 0.5, 0.5)
    
    ## Estimacao Metodo numerico BFGS
    
    estimacao = optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                      hessian = TRUE, l=data$L, control=list(fnscale=-1),
                      r=data$R, x.cure=covariaveis.cura, x.risk=covariaveis.risco,
                      grid.vet=grid.observado)
    
    ## salvar resultados
    matrix.iter[iteracao,] = estimacao$par

    vetor.ep = sqrt(diag(solve(-estimacao$hessian)))
    
    iteracao = iteracao + 1
    
    estimacao
    
  }, error = function(e){
    
    print(paste0("Erro na iteracao ", iteracao, ": ", conditionMessage(e)))
    # valor nulo nessa iteracao
    NULL
  }
  )
  
  if(anyNA(vetor.ep)){
    print("Raiz negativa gerada")
    iteracao = iteracao - 1
  } else{
    matrix.ep[iteracao,] = vetor.ep
  }
  
  # continua as iteracoes se tiver um erro
  if(is.null(result$convergence)){
    iter.error[iteracao] = as.character(iteracao)
  }
  
  # caso contrario, com os resultados
  cat("Realizando iteracao: ", iteracao, "/", n.iter, "\n", sep = "")
}

## retirando os NA
iter.error[!is.na(iter.error)]


## ------
## Verificando se a estimacao foi eficiente
## ------

## esperanca dos estimadores

Theta
esperanca.est = apply(matrix.iter, MARGIN = 2, FUN = mean)

cbind(Theta,esperanca.est)

## desvio-padrao dos estimadores
dp.est = apply(matrix.iter, MARGIN = 2, FUN = sd)

## calculo do vies

Theta.matrix = t(matrix(rep(as.vector(Theta)$Real, n.iter), nrow = 9))

bias.matrix = (matrix.iter-Theta.matrix)/Theta.matrix

## vies (Bias)
bias = colMeans(bias.matrix)

## ---
## coverage probability
## probabilidade de cobertura
## Nivel de 95% de confianca
## ---

## Thetaj +- (quantil_normal_padrao)*(erro-padrao)


nivel.conf = 0.95
quantil = qnorm(p = nivel.conf+((1-nivel.conf)/2), mean=0, sd=1)




limite.superior = matrix.iter[,] + (quantil*matrix.ep[,])
limite.inferior = matrix.iter[,] - (quantil*matrix.ep[,])

## porcentagem de capturacao do intervalo de confianca para as taxas
prob.cobertura = colMeans(Theta.matrix[,] >= limite.inferior & Theta.matrix[,] <= limite.superior)

matriz.resultados = cbind(Theta,esperanca.est, dp.est, bias, prob.cobertura)

objeto = paste("mc_fc_ci_mepp_", amostra, sep = "")
assign(objeto , matriz.resultados)




  



