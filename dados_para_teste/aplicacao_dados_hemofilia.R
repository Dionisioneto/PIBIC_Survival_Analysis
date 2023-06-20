## propostas para o modelo Weibull, MEP e MEPP
rm(list=ls())

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg, ReIns)

#setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')
source('C:/Users/NetoDavi/Desktop/survival_pibic/source/funcoes_sobrevivencia_pibic2023.R')
source('C:/Users/NetoDavi/Desktop/survival_pibic/source/supervisor_functions.r')


setwd("C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste")
## hemofilia_icens
hemo.icens = read.table('hemofilia_icens.txt', header = T)

dim(hemo.icens)
head(hemo.icens)

colnames(hemo.icens)

hemo.icens$L= hemo.icens$L/7
hemo.icens$R = hemo.icens$R/7

## covariaveis: NoDose, Medium, High

## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)

turnbull_fit = ic_np(cbind(L, R)~0, data = hemo.icens, B = c(1,1))

plot(turnbull_fit, lwd = 2, bty = 'n', survRange = c(0,100), axes = F,
     ylab = "Probabilidade de Sobrevivência",
     xlab = "Tempo de Sobrevivência (Dias)")
#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(1,58,4), labels=seq(1,58,4), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), cex.axis=0.7)


## Ajuste do modelo Weibull aos dados
## ---

weibull = ic_par(formula = cbind(L, R) ~ NoDose, data = hemo.icens,
                 model = "ph", dist = "weibull")

weibull$coefficients

#plot(weibull)

AIC.surv(loglik = weibull$llk, n.param = 3, n.sample = dim(hemo.icens)[1])
BIC.surv(loglik = weibull$llk, n.param = 3, n.sample = dim(hemo.icens)[1])
HC.surv(loglik = weibull$llk, n.param = 3,n.sample = dim(hemo.icens)[1])


## ------
## Modelo Exponencial por partes com censura intervalar
## ------

source('C:/Users/NetoDavi/Desktop/survival_pibic/source/mep_interval_fc.R')

ll.hemo = hemo.icens$L/7
rr.hemo = hemo.icens$R/7

n.int = 7

x.f <- cbind(hemo.icens$NoDose)
x.c <- cbind(1, hemo.icens$NoDose)

grid.obs=time.grid.interval(li=ll.hemo, ri=rr.hemo, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,2,5,0.1,10,1,10,
           1,0.5,
           0.5)

max.mep.fc = optim(par = chutes, fn=loglikIC.MEP.fc,
                   gr = NULL, method = "BFGS",
                   control=list(fnscale=1),
                   hessian = TRUE,
                   l = ll.hemo, r = rr.hemo,
                   x.cure=x.c, x.risk=x.f, 
                   grid.vet=grid.obs)

max.mep.fc$par

AIC.surv(loglik = max.mep.fc$value*-1, n.param = length(max.mep.fc$par),
         n.sample = dim(hemo.icens)[1])

BIC.surv(loglik = max.mep.fc$value*-1, n.param = length(max.mep.fc$par),
         n.sample = nrow(hemo.icens))

HC.surv(loglik = max.mep.fc$value*-1, n.param = length(max.mep.fc$par),
        n.sample = nrow(hemo.icens))



## Modelo Exponencial por partes potencia sem racaode cura
source('C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste/proposta_MEPP_int_sem_cura.R')
ll.hemo = hemo.icens$L/7
rr.hemo = hemo.icens$R/7

n.int = 7

x.f <- cbind(hemo.icens$NoDose)
#x.c <- cbind(1, hemo.icens$NoDose)

grid.obs=time.grid.interval(li=ll.hemo, ri=rr.hemo, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(0.1,10,10,
           5,
           0.5)

max.mepp.sfc = optim(par = chutes, fn=loglikIC.mepp.int,
                     gr = NULL, method = "BFGS",
                     control=list(fnscale=1),
                     hessian = TRUE,
                     l = ll.hemo, r = rr.hemo,
                     x.cov=x.f, 
                     grid=grid.obs)

max.mepp.sfc$par

AIC.surv(loglik = max.mepp.sfc$value*-1,
         n.param = length(max.mepp.sfc$par),
         n.sample = dim(hemo.icens)[1])

BIC.surv(loglik = max.mepp.sfc$value*-1,
         n.param = length(max.mepp.sfc$par),
         n.sample = dim(hemo.icens)[1])

HC.surv(loglik = max.mepp.sfc$value*-1,
        n.param = length(max.mepp.sfc$par),
        n.sample = dim(hemo.icens)[1])




## Modelo Proposto (MEPP + fracao de cura)

ll.hemo = hemo.icens$L/7
rr.hemo = hemo.icens$R/7


n.int = 3

x.f <- cbind(hemo.icens$NoDose)
x.c <- cbind(1, hemo.icens$NoDose)

grid.obs=time.grid.interval(li=ll.hemo, ri=rr.hemo, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,2,5,
           2,
           1,0.5,
           0.5)

max.mepp = optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                 control=list(fnscale=1), hessian = TRUE,
                 l = ll.hemo, r = rr.hemo,
                 x.cure=x.c, x.risk=x.f, 
                 grid=grid.obs)

max.mepp$par

max.mepp$value

AIC.surv(loglik = max.mepp$value*-1,
         n.param = length(max.mepp$par),
         n.sample = dim(hemo.icens)[1])

BIC.surv(loglik = max.mepp$value*-1,
         n.param = length(max.mepp$par),
         n.sample = dim(hemo.icens)[1])

HC.surv(loglik = max.mepp$value*-1,
        n.param = length(max.mepp$par),
        n.sample = dim(hemo.icens)[1])






