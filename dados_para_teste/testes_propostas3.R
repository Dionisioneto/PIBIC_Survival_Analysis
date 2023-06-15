## propostas para o modelo Weibull, MEP e MEPP
rm(list=ls())

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg, ReIns)

#setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')
source('C:/Users/NetoDavi/Desktop/survival_pibic/supervisor_functions.r')

#############
#############
## Analise em dados artificiais
#############
#############



n=500
alpha.f   <- 0.8 
lambda.f  <- c(1.1, 0.3, 0.9)
n.intervals <- length(lambda.f)
grid.time <- c(0.5, 2)
beta.f    <- c(-0.5, 0.8)


beta.c    <- c(1.2, 0.3, -0.5)
lambda.c <- 1

Theta = c(lambda.f, alpha.f,beta.c,beta.f)


dadosIC <- sim.std.cure.ICdata(n=n, lambda.par=lambda.f, alpha.par=alpha.f,
                               grid.vet=grid.time, beta.par= beta.f, lambda.parc=1,
                               theta.par = beta.c , A = 1.5, B = 22.5)


x.f <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
x.c <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)

grid.obs=time.grid.interval(li=dadosIC$L, ri=dadosIC$R, type="OBS", bmax=length(lambda.f ))
grid.obs=grid.obs[-c(1, length(grid.obs))]
chutes = c(rep(0.1, length(lambda.f)), 1, 1, 0.5, 0.5, 0.5, 0.5)

test <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
              control=list(fnscale=1), hessian = TRUE, l=dadosIC$L,
              r=dadosIC$R, x.cure=x.c, x.risk=x.f, grid=grid.obs)

for (int in 3:15){
  
  n.intervalos = int
  
  chute = c(rep(0.1,n.intervalos),
            0.8,
            1.2,0.1,
            0.1)
  
  x.f <- cbind(x1=dadosIC$xi1, x2=dadosIC$xi2)
  x.c <- cbind(1, x1=dadosIC$xi1, x2=dadosIC$xi2)
  
  grid.obs=time.grid.interval(li=dadosIC$L, ri=dadosIC$R, type="OBS", bmax=length(lambda.f ))
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  chutes = c(rep(0.1, length(lambda.f)), 1, 1, 0.5, 0.5, 0.5, 0.5)
  
  test <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                control=list(fnscale=1), hessian = TRUE, l=dadosIC$L,
                r=dadosIC$R, x.cure=x.c, x.risk=x.f, grid=grid.obs)
  
  
  aic = AIC.surv(loglik = test$value,
                 n.param = length(test$par))
  
  bic = BIC.surv(loglik = test$value,
                 n.param = length(test$par),
                 n.sample = dim(dadosIC)[1])
  
  hc = HC.surv(loglik = test$value,
               n.param = length(test$par),
               n.sample = dim(dadosIC)[1])
  
  print(paste("n intervalo", int))
  print(paste("AIC: ", aic))
  print(paste("BIC: ", bic))
  print(paste("HC: ", hc))
  
}

#############
#############
## Analise nos dados de hemofilia
#############
#############

setwd("C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste")

hemo.icens = read.table('hemofilia_icens.txt', header = T)

hemo.icens$L = hemo.icens$L/7
hemo.icens$R = hemo.icens$R/7

dim(hemo.icens)
head(hemo.icens)

colnames(hemo.icens)

#hemo.icens$L= hemo.icens$L/7
#hemo.icens$R = hemo.icens$R/7


n.intervalos = 4

grid.obs=time.grid.interval(li=hemo.icens$L, ri=hemo.icens$R,
                            type="OBS", bmax=n.intervalos)

#cut(hemo.icens$L,breaks = grid.obs, include.lowest = T)
#cut(hemo.icens$R,breaks = grid.obs, include.lowest = T)

#prop.table(table(cut(hemo.icens$L,breaks = grid.obs, include.lowest = T)))
#prop.table(table(cut(hemo.icens$R,breaks = grid.obs, include.lowest = T)))

grid.obs=grid.obs[-c(1, length(grid.obs))]

chute = c(c(0.5,0.05,0.05,0.5),
          3,
          1,1,
          -2)

x.c = cbind(1,hemo.icens$NoDose)
x.f = cbind(hemo.icens$NoDose)

test <- optim(par = chute, fn=loglikIC, gr = NULL, method = "BFGS",
              control=list(fnscale=1), hessian = T, l=hemo.icens$L,
              r=hemo.icens$R, x.cure=x.c, x.risk=x.f, grid=grid.obs)

test$par





AIC.surv(loglik = mepp.tent.hemo$loglik, n.param = length(mepp.tent.hemo$estimated))

BIC.surv(loglik = mepp.tent.hemo$loglik, n.param = length(mepp.tent.hemo$estimated), n.sample = dim(hemo.icens)[1])

HC.surv(loglik = mepp.tent.hemo$loglik, n.param = length(mepp.tent.hemo$estimated), n.sample = dim(hemo.icens)[1])



## debug dos dados



p1= SpopMEPP(t=hemo.icens$L, lambda.par=rep(0.01,n.intervalos), 
         alpha.par=1.2, grid.vet=grid.obs, 
         beta.par=c(0.1),
         theta.par=c(1,0.1), 
         x.cure=x.c, x.risk=x.f) - SpopMEPP(t=hemo.icens$R, lambda.par=rep(0.01,n.intervalos), 
                                            alpha.par=1.2, grid.vet=grid.obs, 
                                            beta.par=c(0.1),
                                            theta.par=c(1,0.1), 
                                            x.cure=x.c, x.risk=x.f)



#############
#############
## Analise nos dados smoke
#############
#############

smoke2009 = read.table('smoke_cessation_Bannerge2009.txt', header = T)

dim(smoke2009)
head(smoke2009)

## variaveis do estudo

tratamento = smoke2009$SIUC # Tipo de tratamento (SI/UC)
n_cigarros = smoke2009$F10Cigs_pad # Numero de cigarros fumados por dia, normalizado
duracao_dependente = smoke2009$Duration_pad # Duracao como dependente, normalizado
sexo = smoke2009$SexF

covariaveis = cbind(tratamento, sexo, n_cigarros, duracao_dependente)


n.intervalos = 7


grid.obs=time.grid.interval(li=smoke2009$Timept1, ri=smoke2009$Timept2,
                            type="OBS", bmax=n.intervalos)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chute = c(c(5,5,0.2,0.3,3,200,200),
          100,
          1.2,0.1,0.1,0.1,0.1,
          0.1,0.1,0.1,0.1)



est.mepp.smoke <- optim(par = chute, fn=loglikIC, gr = NULL, 
                        method = "BFGS",
                      control=list(fnscale=1), hessian = T,
                      l=smoke2009$Timept1,
                      r=smoke2009$Timept2, 
                      x.cure=covariaveis, x.risk=cbind(1, covariaveis),
                      grid=grid.obs)

est.mepp.smoke$par
est.mepp.smoke$value

AIC.surv(loglik = est.mepp.smoke$value*-1, 
         n.param = length(est.mepp.smoke$par))

BIC.surv(loglik = est.mepp.smoke$value*-1,
         n.param = length(est.mepp.smoke$par), 
         n.sample = dim(smoke2009)[1])

HC.surv(loglik = est.mepp.smoke$value*-1,
        n.param = length(est.mepp.smoke$par),
        n.sample = dim(smoke2009)[1])










