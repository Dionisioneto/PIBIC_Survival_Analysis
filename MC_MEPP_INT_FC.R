################################################
## Code: Monte Carlo simualtion of the new    ##
## suvival cure model mixture model using     ##
################################################
## Author: Dionisio Neto                      ##
################################################
## Date: 27/05/2023                           ##
################################################

source("C:/Users/NetoDavi/Desktop/survival_pibic/supervisor_functions.r")

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg, eha)

n <- 500 # Tamanho amostral

#--- Parametros do modelo:

alpha.f   <- 0.8 
lambda.f  <- c(1.1, 0.3, 0.9)
n.intervals <- length(lambda.f)
grid.time <- c(0.5, 2)
beta.f    <- c(-0.5, 0.8)
#beta.f    <- 0

beta.c    <- c(1.2, 0.5, -0.5)
lambda.c <- 1

Theta = c(lambda.f, alpha.f,beta.c,beta.f)


# Ajuste usando os intervalos de tempo:

ini.info <- c(lambda.f, alpha.f, beta.c,  beta.f)

npar <- length(c(lambda.f, alpha.f,beta.c, beta.f))

iter.error = c()

samp <- 200
i = 1

est  <- matrix(NA, ncol=npar, nrow=samp)
matrix.ep = matrix(data = 0, nrow = samp, 
                   ncol = length(Theta))

prop.cens = matrix(data = 0, nrow = samp, 
                   ncol = 2)

prop.cura = matrix(data = 0, nrow = samp, 
                   ncol = 2)



while(i <= samp) {
  result = tryCatch({
    cat("Realizando iteracao: ", i, "/", samp, "\n", sep = "")
    dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f, 
                                   grid.vet=grid.time, beta.par= beta.f, lambda.parc=1, 
                                   theta.par = beta.c , A = 0.4, B =22)
    
    
    x.f <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
    x.c <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)
    
    grid.obs=time.grid.interval(li=dadosIC$L, ri=dadosIC$R, type="OBS", bmax=length(lambda.f ))
    grid.obs=grid.obs[-c(1, length(grid.obs))]
    chutes = c(rep(0.1, length(lambda.f)), 1, 1, 0.5, 0.5, 0.5, 0.5)
    
    test <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                  control=list(fnscale=-1), hessian = TRUE, l=dadosIC$L, 
                  r=dadosIC$R, x.cure=x.c, x.risk=x.f, grid.vet=grid.obs)
    
    est[i,] <- test$par
    prop.cens[i,] = prop.table(table(dadosIC$delta))
    prop.cura[i,] = prop.table(table(dadosIC$Y))
    
    
    vetor.ep = sqrt(diag(solve(-test$hessian)))
    
    
    i  = i+1
    
    test
    
  }, error = function(e){
    
    print(paste0("Erro na iteracao ", i, ": ", conditionMessage(e)))
    # valor nulo nessa iteracao
    NULL
  }
  )
  
  if(anyNA(vetor.ep)){
    print("Raiz negativa gerada")
    i = i - 1
  } else{
    matrix.ep[i,] = vetor.ep
  }
  
  # continua as iteracoes se tiver um erro
  if(is.null(test$convergence)){
    iter.error[i] = as.character(i)
  }
  
  # caso contrario, com os resultados
  
}


for(i in 1: samp){
  cat("Realizando iteracao: ", i, "/", samp, "\n", sep = "")
  dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f,
                                 grid.vet=grid.time, beta.par= beta.f, lambda.parc=1,
                                 theta.par = beta.c , A = 1.5, B = 22.5)
  
  
  x.f <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
  x.c <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)
  
  grid.obs=time.grid.interval(li=dadosIC$L, ri=dadosIC$R, type="OBS", bmax=length(lambda.f ))
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  chutes = c(rep(0.1, length(lambda.f)), 1, 1, 0.5, 0.5, 0.5, 0.5)
  
  test <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                control=list(fnscale=-1), hessian = TRUE, l=dadosIC$L,
                r=dadosIC$R, x.cure=x.c, x.risk=x.f, grid.vet=grid.obs)
  
  est[i,] <- test$par
  prop.cens[i,] = prop.table(table(dadosIC$delta))
  prop.cura[i,] = prop.table(table(dadosIC$Y))
  
  
  vetor.ep = sqrt(diag(solve(-test$hessian)))
  matrix.ep[i,] = vetor.ep
}


Theta.matrix = matrix(rep(Theta,samp), ncol = length(Theta), byrow = T)

## esperanca das iteracoes
esperanca.est = apply(est, MARGIN = 2, FUN = mean)

## desvio-padrao dos estimadores
dp.est = apply(est, MARGIN = 2, FUN = sd)

## calculo do vies

Theta.matrix = t(matrix(rep(as.vector(Theta),samp), nrow = 9))

bias.matrix = (est-Theta.matrix)/Theta.matrix*100

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


limite.superior = est[,] + (quantil*matrix.ep)
limite.inferior = est[,] - (quantil*matrix.ep)

## porcentagem de capturacao do intervalo de confianca para as taxas
prob.cobertura = colMeans(Theta.matrix >= limite.inferior & Theta.matrix <= limite.superior)

matriz.resultados = cbind(Theta,esperanca.est, dp.est, bias, prob.cobertura)

matriz.resultados

summary(prop.cura[,1])
boxplot(prop.cura[,1], ylim = c(0,1))

summary(prop.cens[,1])
boxplot(prop.cens[,1], ylim = c(0,1))

#setwd('C:\\Users\\NetoDavi\\Desktop\\survival_pibic')
#write.csv2(x = matriz.resultados, file = "resultado2_n1000.csv")

## histogramas
par(mfrow=c(3,3), mai = c(0.6, 0.6, 0.2, 0.1))
hist(est[,1], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(lambda)[1]))

hist(est[,2], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(lambda)[2]))
hist(est[,3], col = "steelblue", main = "",
     ylab = "Frequência", xlab = expression(hat(lambda)[3]))


hist(est[,4], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(alpha)))

hist(est[,5], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(b)[0]))
hist(est[,6], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(b)[1]))
hist(est[,7], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(b)[2]))

hist(est[,8], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(beta)[1]))
hist(est[,9], col = "steelblue", main = "", breaks = 15,
     ylab = "Frequência", xlab = expression(hat(beta)[2]))

# QQ-Plots

par(mfrow=c(3,3))

qqnorm(est[,1], pch = 1, frame = FALSE, main = expression(hat(lambda)[1]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,1], col = "steelblue", lwd = 2)

qqnorm(est[,2], pch = 1, frame = FALSE, main = expression(hat(lambda)[2]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,2], col = "steelblue", lwd = 2)

qqnorm(est[,3], pch = 1, frame = FALSE, main = expression(hat(lambda)[3]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,3], col = "steelblue", lwd = 2)


qqnorm(est[,4], pch = 1, frame = FALSE, main = expression(hat(alpha)),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,4], col = "steelblue", lwd = 2)


qqnorm(est[,5], pch = 1, frame = FALSE, main = expression(hat(b)[0]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,5], col = "steelblue", lwd = 2)

qqnorm(est[,6], pch = 1, frame = FALSE, main = expression(hat(b)[1]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,6], col = "steelblue", lwd = 2)

qqnorm(est[,7], pch = 1, frame = FALSE, main = expression(hat(b)[2]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,7], col = "steelblue", lwd = 2)

qqnorm(est[,8], pch = 1, frame = FALSE, main = expression(hat(beta)[1]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,8], col = "steelblue", lwd = 2)

qqnorm(est[,9], pch = 1, frame = FALSE, main = expression(hat(beta)[2]),
       xlab = "Quantis Teóricos", ylab = "Quantis Observados")
qqline(est[,9], col = "steelblue", lwd = 2)


## -----
## Estudo da curva de sobrevivencia pelo estimador de
## Turnbull
## -----
library(ReIns)

n = 500

dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f, 
                               grid.vet=grid.time, beta.par= beta.f, lambda.parc=lambda.c, 
                               theta.par = beta.c , A = 1.4, B =22)

prop.table(table(dadosIC$delta))

tempo.aval = seq(0,2.5,0.1)

trnb.fit = Turnbull(x = tempo.aval, L = dadosIC$L, R = dadosIC$R,
                    censored = dadosIC$delta)



par(mfrow=c(1,2))

plot(tempo.aval, trnb.fit$surv, type = "s", 
     ylab = "Estimador de Turnbull para S(t)")



turnbull_fit = ic_np(cbind(L, R)~0, data = dadosIC)

plot(turnbull_fit)



## ajuste do modelo exponencial por partes potencia 
## em censura intervalar e fracao de curados

fit.mepp.cf = function(L, R, n.int, cov.risco, cov.cura, start){
  
  ## extracao do grid observado
  
  grid.obs=time.grid.interval(li=L, ri=R, type="OBS", bmax= n.int)
  grid.obs=grid.obs[-c(1, length(grid.obs))]

  est <- optim(par = start, fn=loglikIC, gr = NULL, method = "BFGS",
                control=list(fnscale=-1), hessian = TRUE, l=L, 
                r=R, x.cure=cov.cura, x.risk=cov.risco, grid.vet=grid.obs)
  
  estimated = est$par
  hessian = est$hessian
  loglik = est$value
  
  results = list(estimated = estimated, hessian = hessian, loglik = loglik)
  
  return(results)
}

mepp.tent = fit.mepp.cf(L = dadosIC$L, R = dadosIC$R, n.int = 3,
            cov.risco = cbind(x1=dadosIC$xi1, x2=dadosIC$xi2),
            cov.cura = cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2),
            start = c(rep(0.1, length(lambda.f)), 1, 1, 0.5, 0.5, 0.5, 0.5))

# mepp.tent$estimated
# mepp.tent$hessian
# mepp.tent$loglik

## tentativa em aplicacao a dados reais


#setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')
setwd("C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste")
## hemofilia_icens
hemo.icens = read.table('hemofilia_icens.txt', header = T)

dim(hemo.icens)
head(hemo.icens)

colnames(hemo.icens)































