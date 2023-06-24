## ---------------------------------------------------------------
## ajuste aos dados dados HIV_SUM_set2 
## ---------------------------------------------------------------

rm(list=ls())

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg, ReIns)

#setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')
source('C:/Users/NetoDavi/Desktop/survival_pibic/source/funcoes_sobrevivencia_pibic2023.R')
source('C:/Users/NetoDavi/Desktop/survival_pibic/source/supervisor_functions.r')


setwd("C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste")

hiv2 = read.table('HIV_SUM_set2.txt', header = T)

dim(hiv2)
head(hiv2) 


## ----
## Analise descritiva dos dados
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
## ----
par(mfrow = c(1,2))
npmle_fit_hiv2 <- ic_np(cbind(Li, Ri) ~ 0,
                          data = hiv2, B= c(1,1))

plot(npmle_fit_hiv2,
     xlab = "Tempo de Sobrevivência (em dias)",
     ylab = "Probabilidade de Sobrevivência",) ## temos o plato !!

## verificando outra curva

library(interval)

npmle_int_hiv <-icfit(Surv(Li,Ri,type="interval2")~DoseType, data=hiv2)

plot(npmle_int_hiv, shade = F,
     XLAB = "Tempo de Sobrevivência (em dias)",
     YLAB = "Probabilidade de Sobrevivência")

## can pick out just one group
# plot(npmle_int_hiv[1], shade = F, 
#      XLAB = "Tempo de Sobrevivência (em dias)",
#      YLAB = "Probabilidade de Sobrevivência")

par(mfrow = c(1,2))


## ------
## O modelo Weibull
## ------

hiv2 = read.table('HIV_SUM_set2.txt', header = T)

hiv2.weibull = ic_par(cbind(Li/7, Ri/7) ~ DoseType,
                      data = hiv2, model = "ph", dist = "weibull")

hiv2.weibull$coefficients

AIC.surv(loglik = hiv2.weibull$llk, 
         n.param = length(hiv2.weibull$coefficients),
         n.sample = dim(hiv2)[1])

BIC.surv(loglik = hiv2.weibull$llk, n.param = length(hiv2.weibull$coefficients),
         n.sample = dim(hiv2)[1])

HC.surv(loglik = hiv2.weibull$llk, n.param = length(hiv2.weibull$coefficients),
        n.sample = dim(hiv2)[1])

CAIC.surv(loglik = hiv2.weibull$llk, n.param = length(hiv2.weibull$coefficients),
          n.sample = dim(hiv2)[1])

## ---
## modelo exponencial por partes com fracao de cura
## ---

source('C:/Users/NetoDavi/Desktop/survival_pibic/source/mep_interval_fc.R')

hiv2 = read.table('HIV_SUM_set2.txt', header = T)

# Misturar as linhas

hiv2 = read.table('HIV_SUM_set2.txt', header = T)

hiv2 <- hiv2[sample(row.names(hiv2)), ]

left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = hiv2$DoseType

n.int = 5

x.f <- cbind(cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,0.1,1,2,2,
           1,0.5,
           0.5)

max.mep = optim(par = chutes, fn=loglikIC.MEP.fc, gr = NULL, method = "BFGS",
                control=list(fnscale=1), hessian = TRUE,
                l = left, r = right,
                x.cure=x.c, x.risk=x.f, 
                grid.vet=grid.obs)

max.mep$par



#max.mep$hessian
max.mep$value

AIC.surv(loglik = max.mep$value,
         n.param = length(max.mep$par),
         n.sample = dim(hiv2)[1])

BIC.surv(loglik = max.mep$value,
         n.param = length(max.mep$par),
         n.sample = dim(hiv2)[1])

HC.surv(loglik = max.mep$value,
        n.param = length(max.mep$par),
        n.sample = dim(hiv2)[1])

CAIC.surv(loglik = max.mep$value,
          n.param = length(max.mep$par),
          n.sample = dim(hiv2)[1])

parametros5_MEP_int = max.mep$par

## ---
## modelo exponencial por partes potencia com fracao de cura
## ---
source('C:/Users/NetoDavi/Desktop/survival_pibic/source/supervisor_functions.R')

hiv2 = read.table('HIV_SUM_set2.txt', header = T)


## Duas particoes

left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = hiv2$DoseType

n.int = 2

x.f <- cbind(x1=cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(0.1,1,
           5,
           1,0.5,
           0.5)


estimacao <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                   control=list(fnscale=1), hessian = TRUE, l=left,
                   r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)

estimacao$par

est.matrix.var.cov = solve(estimacao$hessian)
erro.padrao = sqrt(diag(est.matrix.var.cov))
erro.padrao

estimacao$value

AIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
         n.sample = nrow(hiv2))
BIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
         n.sample = nrow(hiv2))
HC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
        n.sample = nrow(hiv2))

CAIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
          n.sample = nrow(hiv2))


parametros2_int = estimacao$par


## cinco particoes

left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = hiv2$DoseType

n.int = 5

x.f <- cbind(x1=cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes =  c(0.1,1,1,5,0.1,
            5,
            1, 0.5,
            0.5)

estimacao <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                   control=list(fnscale=1), hessian = TRUE, l=left,
                   r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)

parametros5_int=estimacao$par

estimacao$value

AIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
         n.sample = nrow(hiv2))
BIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
         n.sample = nrow(hiv2))
HC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
        n.sample = nrow(hiv2))

CAIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
          n.sample = nrow(hiv2))

parametros5_int = estimacao$par

## ---
## Plot dos modelos escolhidos
## ---

setwd('C:\\Users\\NetoDavi\\Desktop\\survival_pibic\\imagens_tcc\\aplicacao_HIV')

dev.off()

npmle_int_hiv <-icfit(Surv(Li/7,Ri/7,type="interval2")~DoseType, data=hiv2)

plot(npmle_int_hiv, shade = F,
     XLAB = "Tempo de Sobrevivência (em semanas)",
     YLAB = "Probabilidade de Sobrevivência",
     LEGEND = F,
     XLEG = 0.2, YLEG = 0.2, axes = F,lwd = 2)

legend(x = 0.1, y = 0.2,           # Position
       legend = c("Sem o fator VIII",
                  "Dose baixa do fator VIII"),  # Legend texts
       lty = c(1, 2),           # Line types
       lwd = 2, bty = 'n')  

#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,8,1), labels=seq(0,8,1), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), cex.axis=0.7)


valores.tempo = seq(0,40,length.out = nrow(hiv2)) 

x.f <- cbind(x1=hiv2$DoseType)
x.c <- cbind(1, x1=hiv2$DoseType)

## ---
## MEPP com 2 particoes
## ---

## grupo Dosetype == 1
ind.dose1 = hiv2$DoseType == 1

valores.tempo1 = seq(0,40,length.out = nrow(hiv2[ind.dose1,])) 

grid.observavel2=time.grid.interval(li=hiv2$Li/7, ri=hiv2$Ri/7, type="OBS", 
                            bmax=2)

grid.observavel2=grid.observavel2[-c(1, length(grid.observavel2))]


surv2.estrat1 = SpopMEPP(t= valores.tempo1, lambda.par=parametros2_int[1:2],
                       alpha.par=parametros2_int[3],
                       grid.vet=grid.observavel2, 
                       theta.par=parametros2_int[4:5], 
                       beta.par=parametros2_int[6],
                       x.cure=x.c[ind.dose1,], x.risk=x.f[ind.dose1])

plot(npmle_int_hiv, shade = F,
     XLAB = "Tempo de Sobrevivência (em semanas)",
     YLAB = "Probabilidade de Sobrevivência",
     LEGEND = F,
     XLEG = 0.2, YLEG = 0.2, axes = F,lwd = 2)

#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,8,1), labels=seq(0,8,1), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), cex.axis=0.7)


lines(valores.tempo1,surv2.estrat1, col = 'red', lty=2, lwd = 2)

## ---
## grupo Dosetype == 0
## ---

valores.tempo0 = seq(0,40,length.out = nrow(hiv2[!ind.dose1,])) 

grid.observavel2=time.grid.interval(li=hiv2$Li/7, ri=hiv2$Ri/7, type="OBS", 
                                    bmax=2)

grid.observavel2=grid.observavel2[-c(1, length(grid.observavel2))]


surv2.estrat0 = SpopMEPP(t= valores.tempo0, lambda.par=parametros2_int[1:2],
                         alpha.par=parametros2_int[3],
                         grid.vet=grid.observavel2, 
                         theta.par=parametros2_int[4:5], 
                         beta.par=parametros2_int[6],
                         x.cure=x.c[!ind.dose1,], x.risk=x.f[!ind.dose1])

lines(valores.tempo0,surv2.estrat0, col = 'red', lwd = 2)


legend(x = 0.1, y = 0.25,           
       legend = c("Sem o fator VIII",
                  "Dose baixa do fator VIII",
                  "MEPP (2)"),  
       lty = c(1, 2),          
       lwd = 2, bty = 'n',
       col = c("black", "black", "red"))  


## ---
## MEPP com 5 particoes
## ---

grid.observavel5=time.grid.interval(li=hiv2$Li/7, ri=hiv2$Ri/7, type="OBS", 
                                    bmax=5)

grid.observavel5=grid.observavel5[-c(1, length(grid.observavel5))]

## grupo Dosetype == 1
valores.tempo1 = seq(0,40,length.out = nrow(hiv2[ind.dose1,])) 

surv5.estrat1 = SpopMEPP(t= valores.tempo1, 
                         lambda.par=parametros5_int[1:5],
                         alpha.par=parametros5_int[6],
                         grid.vet=grid.observavel5, 
                         theta.par=parametros5_int[7:8], 
                         beta.par=parametros5_int[9],
                         x.cure=x.c[ind.dose1,], x.risk=x.f[ind.dose1]) 

plot(npmle_int_hiv, shade = F,
     XLAB = "Tempo de Sobrevivência (em semanas)",
     YLAB = "Probabilidade de Sobrevivência",
     LEGEND = F,
     XLEG = 0.2, YLEG = 0.2, axes = F,lwd = 2)

#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,8,1), labels=seq(0,8,1), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), cex.axis=0.7)

lines(valores.tempo1,surv5.estrat1, col = 'red', lty=2, lwd = 2)

## grupo Dosetype == 0
valores.tempo0 = seq(0,40,length.out = nrow(hiv2[!ind.dose1,])) 

surv5.estrat0 = SpopMEPP(t= valores.tempo0, 
                        lambda.par=parametros5_int[1:5],
                        alpha.par=parametros5_int[6],
                        grid.vet=grid.observavel5, 
                        theta.par=parametros5_int[7:8], 
                        beta.par=parametros5_int[9],
                        x.cure=x.c[!ind.dose1,], x.risk=x.f[!ind.dose1]) 

lines( valores.tempo0,surv5.estrat0, col = 'red',  lwd = 2)

legend(x = 0.1, y = 0.25,           
       legend = c("Sem o fator VIII",
                  "Dose baixa do fator VIII",
                  "MEPP (5)"),  
       lty = c(1, 2),          
       lwd = 2, bty = 'n',
       col = c("black", "black", "red"))  



## ----
## Modelo exponencial por partes com fracao de cura
## versao com 5 particoes
## ----
grid.observavel5=time.grid.interval(li=hiv2$Li/7, ri=hiv2$Ri/7, type="OBS", 
                                    bmax=5)

grid.observavel5=grid.observavel5[-c(1, length(grid.observavel5))]

## grupo Dosetype == 1
valores.tempo1 = seq(0,40,length.out = nrow(hiv2[ind.dose1,])) 

parametros5_MEP_int

surv5.MEP.estrat1 = SpopMEP(t=valores.tempo1, 
                         lambda.par=parametros5_MEP_int[1:5],
                         grid.vet=grid.observavel5, 
                         theta.par=parametros5_MEP_int[6:7], 
                         beta.par=parametros5_MEP_int[8],
                         x.cure=x.c[ind.dose1,], x.risk=x.f[ind.dose1]) 

plot(npmle_int_hiv, shade = F,
     XLAB = "Tempo de Sobrevivência (em semanas)",
     YLAB = "Probabilidade de Sobrevivência",
     LEGEND = F,
     XLEG = 0.2, YLEG = 0.2, axes = F,lwd = 2)

#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,8,1), labels=seq(0,8,1), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), cex.axis=0.7)

lines(valores.tempo1,surv5.MEP.estrat1, col = 'red', lty=2, lwd = 2)

## grupo Dosetype == 0
valores.tempo0 = seq(0,40,length.out = nrow(hiv2[!ind.dose1,])) 

surv5.MEP.estrat0 = SpopMEP(t=valores.tempo0, 
                             lambda.par=parametros5_MEP_int[1:5],
                             grid.vet=grid.observavel5, 
                             theta.par=parametros5_MEP_int[6:7], 
                             beta.par=parametros5_MEP_int[8],
                             x.cure=x.c[!ind.dose1,], x.risk=x.f[!ind.dose1]) 

lines(valores.tempo0,surv5.MEP.estrat0, col = 'red',  lwd = 2)

legend(x = 0.1, y = 0.25,           
       legend = c("Sem o fator VIII",
                  "Dose baixa do fator VIII",
                  "MEP (5)"),  
       lty = c(1, 2),          
       lwd = 2, bty = 'n',
       col = c("black", "black", "red"))  


## ---
## Modelo Weibull
## ---

plot(npmle_int_hiv, shade = F,
     XLAB = "Tempo de Sobrevivência (em semanas)",
     YLAB = "Probabilidade de Sobrevivência",
     LEGEND = F,
     XLEG = 0.2, YLEG = 0.2, axes = F,lwd = 2)

#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,8,1), labels=seq(0,8,1), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,1,0.1), labels=seq(0,1,0.1), cex.axis=0.7)

hiv2.weibull$formula

surv.newdata = data.frame(DoseType = c(0, 1))

lines(hiv2.weibull,newdata = surv.newdata,
      col = c('red', 'red'), cis = F)

legend(x = 0.1, y = 0.25,           
       legend = c("Sem o fator VIII",
                  "Dose baixa do fator VIII",
                  "Weibull"),  
       lty = c(1, 2),          
       lwd = 2, bty = 'n',
       col = c("black", "black", "red")) 


### ---
### analise residual
### ---

### ---
## residuos de martingales
### ---
hiv2$Censind

spop_mepp2_cen1_l = SpopMEPP(t = hiv2$Li[hiv2$Censind==1]/7, lambda.par = parametros2_int[1:2], 
                       alpha.par = parametros2_int[3],
                       theta.par = parametros2_int[4:5],
                       beta.par = parametros2_int[6], 
                       x.cure = x.c[hiv2$Censind==1,],
                       x.risk = x.f[hiv2$Censind==1,],
                       grid.vet = grid.observavel2)

spop_mepp2_cen1_r = SpopMEPP(t = hiv2$Ri[hiv2$Censind==1]/7, lambda.par = parametros2_int[1:2], 
                        alpha.par = parametros2_int[3],
                        theta.par = parametros2_int[4:5],
                        beta.par = parametros2_int[6], 
                        x.cure = x.c[hiv2$Censind==1,],
                        x.risk = x.f[hiv2$Censind==1,],
                        grid.vet = grid.observavel2)

rm_hiv_cen1 = ((spop_mepp2_cen1_l*(1-log(spop_mepp2_cen1_l))) - (spop_mepp2_cen1_r*(1-log(spop_mepp2_cen1_r))))/(spop_mepp2_cen1_l-spop_mepp2_cen1_r)
rm_hiv_cen1prof = ((spop_mepp2_cen1_l*(log(spop_mepp2_cen1_l))) - (spop_mepp2_cen1_r*(log(spop_mepp2_cen1_r))))/(spop_mepp2_cen1_l-spop_mepp2_cen1_r)

###-----------------------------------
spop_mepp2_cen0_l = SpopMEPP(t = hiv2$Li[hiv2$Censind==0]/7, lambda.par = parametros2_int[1:2], 
                             alpha.par = parametros2_int[3],
                             theta.par = parametros2_int[4:5],
                             beta.par = parametros2_int[6], 
                             x.cure = x.c[hiv2$Censind==0,],
                             x.risk = x.f[hiv2$Censind==0,],
                             grid.vet = grid.observavel2)

spop_mepp2_cen0_r = SpopMEPP(t = hiv2$Ri[hiv2$Censind==0]/7, lambda.par = parametros2_int[1:2], 
                             alpha.par = parametros2_int[3],
                             theta.par = parametros2_int[4:5],
                             beta.par = parametros2_int[6], 
                             x.cure = x.c[hiv2$Censind==0,],
                             x.risk = x.f[hiv2$Censind==0,],
                             grid.vet = grid.observavel2)

rm_hiv_cen0 = ((spop_mepp2_cen0_l*(1-log(spop_mepp2_cen0_l))) - (spop_mepp2_cen0_r*(1-log(spop_mepp2_cen0_r))))/(spop_mepp2_cen0_l-spop_mepp2_cen0_r)
rm_hiv_cen0prof = ((spop_mepp2_cen0_l*(log(spop_mepp2_cen0_l))) - (spop_mepp2_cen0_r*(log(spop_mepp2_cen0_r))))/(spop_mepp2_cen0_l-spop_mepp2_cen0_r)


plot(1:length(rm_hiv_cen1), rm_hiv_cen1,
     xlim = c(0,270),
     ylim = c(0,1.1),
     col = 'steelblue',
     axes = F,
     xlab = "Observação",
     ylab = "Resíduos ajustados de Cox-Snell",
     pch = 19)

axis(side=1, at=seq(0,270,8),
     labels=seq(0,270,8), cex.axis=0.7)

axis(side=2, at=seq(0,1.1,0.1),
     labels=seq(0,1.1,0.1), cex.axis=0.7)

points(1:length(rm_hiv_cen0), rm_hiv_cen0,
       col = '#D6604D', pch = 16)

legend(x = 200, y = 0.3, bty = 'n',
       legend = c("Evento", "Censura"),
       col = c('steelblue', '#D6604D'),
       pch = c(19,19),
       cex = 1.3)

plot(1:length(rm_hiv_cen1prof), rm_hiv_cen1,
     xlim = c(0,270),
     ylim = c(0,1.1),
     col = 'blue')

points(1:length(rm_hiv_cen0prof), rm_hiv_cen0,
       col = 'red')




l.cr = -log(spop_mepp2_cen0_l)
r.cr = -log(spop_mepp2_cen0_r)

l.r.cr = as.data.frame(cbind(l.cr, r.cr))

turnbull.cr = ic_sp(cbind(l.cr , r.cr)~0, data = l.r.cr)

list.cr = seq(0,10, length.out = length(spop_mepp2_cen0_r))

Ht.rm_hiv_cen0 = -log(1-getFitEsts(turnbull.cr , q=list.cr))

plot(Ht.rm_hiv_cen0, rm_hiv_cen0)


l.cr1 = -log(spop_mepp2_cen1_l)
r.cr1 = -log(spop_mepp2_cen1_r)

l.r.cr = as.data.frame(cbind(l.cr, r.cr))

turnbull.cr = ic_sp(cbind(l.cr , r.cr)~0, data = l.r.cr)

list.cr = seq(0,10, length.out = length(spop_mepp2_cen0_r))

Ht.rm_hiv_cen0 = -log(1-getFitEsts(turnbull.cr , q=list.cr))

plot(Ht.rm_hiv_cen0, rm_hiv_cen0)

points(-log(spop_mepp2_cen1_l), rm_hiv_cen1, col = 'blue')
points(-log(spop_mepp2_cen1_r), rm_hiv_cen1, col = 'blue')


turnbull_fit_hiv = ic_sp(cbind(Li/7, Ri/7)~DoseType, data = hiv2)

list.rm_hiv_cen0 = seq(0,8, length.out = length(rm_hiv_cen0))

Ht.rm_hiv_cen0 = -log(1-getFitEsts(turnbull_fit_hiv, q=list.rm_hiv_cen0))

plot(Ht.rm_hiv_cen0, rm_hiv_cen0,
     col = 'red')


list.rm_hiv_cen1 = seq(0,3, length.out = length(rm_hiv_cen1))

Ht.rm_hiv_cen1 = -log(1-getFitEsts(turnbull_fit_hiv, q=list.rm_hiv_cen1))

plot(Ht.rm_hiv_cen1, rm_hiv_cen1, col = 'blue')





spop_mepp2l = SpopMEPP(t = hiv2$Li/7, lambda.par = parametros2_int[1:2], 
                             alpha.par = parametros2_int[3],
                             theta.par = parametros2_int[4:5],
                             beta.par = parametros2_int[6], 
                             x.cure = x.c,
                             x.risk = x.f,
                             grid.vet = grid.observavel2)

spop_mepp2r = SpopMEPP(t = hiv2$Ri/7, lambda.par = parametros2_int[1:2], 
                             alpha.par = parametros2_int[3],
                             theta.par = parametros2_int[4:5],
                             beta.par = parametros2_int[6], 
                             x.cure = x.c,
                             x.risk = x.f,
                             grid.vet = grid.observavel2)



rm_hiv = ((spop_mepp2_l*(1-log(spop_mepp2_l))) - (spop_mepp2_r*(1-log(spop_mepp2_r))))/(spop_mepp2_l-spop_mepp2_r)


l.r.cr = as.data.frame(cbind(l.cr, r.cr))

turnbull.cr = ic_sp(cbind(Li/7 , Ri/7)~DoseType, data = hiv2)

list.cr = seq(0,2, length.out = length(rm_hiv))

Ht.rm_hiv = -log(1-getFitEsts(turnbull.cr , q=list.cr))

plot(rm_hiv, Ht.rm_hiv)

plot(rm_hiv)



