### -------
## Validacao es estudo dos modelos de sobrevivencia
## em dados reais
## Data: 30/05/2023 
### -------

## propostas para o modelo Weibull, MEP e MEPP
rm(list=ls())

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg, ReIns)

#setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')
source('C:/Users/NetoDavi/Desktop/survival_pibic/supervisor_functions.r')


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

turnbull_fit = ic_np(cbind(L, R)~0, data = hemo.icens)

plot(turnbull_fit, lwd = 2, bty = 'n', survRange = c(0,100),
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

lines(weibull, col = 'red')

## ------
## Modelo Exponencial por partes com censura intervalar
## ------

source('C:/Users/NetoDavi/Desktop/survival_pibic/mep_interval_fc.R')

ll.hemo = hemo.icens$L/7
rr.hemo = hemo.icens$R/7

n.int = 7

x.f <- cbind(hemo.icens$NoDose)
x.c <- cbind(1, hemo.icens$NoDose)

grid.obs=time.grid.interval(li=ll.hemo, ri=rr.hemo, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,1,1,1,1,10,10,
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

l.turnbull = turnbull_fit$T_bull_Intervals[1,]
r.turnbull = turnbull_fit$T_bull_Intervals[2,]



spoprrh = SpopMEP(t = sort(rr.hemo,decreasing=T)/7, 
                  lambda.par = max.mep.fc$par[1:7],
                  grid.vet = grid.obs, theta.par = max.mep.fc$par[8:9],
                  beta.par = max.mep.fc$par[10], 
                  x.cure = x.c, x.risk = x.f)



## Modelo Exponencial por partes potencia sem racaode cura
source('C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste/proposta_MEPP_int_sem_cura.R')
ll.hemo = hemo.icens$L/7
rr.hemo = hemo.icens$R/7

n.int = 3

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










#mepp.tent.hemo$estimated

## breast cancer

breast = read.table('breast.txt', header = T)

dim(breast)
head(breast)


breast$left = breast$left/7
breast$right = breast$right/7

breast$right = ifelse(is.na(breast$right), Inf, breast$right)

## analise exploratoria
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)

turnbull_fit_breast = ic_np(cbind(left, right)~0, data = breast)

plot(turnbull_fit_breast)


## modelo weibull

weibull.breast = ic_par(formula = cbind(left, right) ~ ther, data = breast,
                 model = "ph", dist = "weibull")

weibull.breast$coefficients

plot(weibull.breast)


AIC.surv(loglik = weibull.breast$llk,
         n.param = length(weibull.breast$coefficients))

BIC.surv(loglik = weibull.breast$llk,
         n.param = length(weibull.breast$coefficients),
         n.sample = dim(breast)[1])

HC.surv(loglik = weibull.breast$llk,
        n.param = length(weibull.breast$coefficients),
        n.sample = dim(breast)[1])




n.intervalos = 15

for (int in 10:15){
  
  n.intervalos = 15
  
  chute = c(rep(0.1,n.intervalos),
            0.8,
            1.2,0.1,
            0.1)
  
  mepp.tent.breast = fit.mepp.cf(L = breast$left, R = breast$right, n.int = n.intervalos,
                                 cov.risco = cbind(breast$ther),
                                 cov.cura = cbind(1, breast$ther),
                                 start = chute)
  
  aic = AIC.surv(loglik = mepp.tent.breast$loglik,
           n.param = length(mepp.tent.breast$estimated))
  
  bic = BIC.surv(loglik = mepp.tent.breast$loglik,
           n.param = length(mepp.tent.breast$estimated),
           n.sample = dim(breast)[1])
  
  hc = HC.surv(loglik = mepp.tent.breast$loglik,
          n.param = length(mepp.tent.breast$estimated),
          n.sample = dim(breast)[1])
  
  print(paste("n intervalo", int))
  print(paste("AIC: ", aic))
  print(paste("BIC: ", bic))
  print(paste("HC: ", hc))
  
}

## o melhor modelo foi o com dois intervalos!!

n.intervalos = 2

chute = c(rep(0.1,n.intervalos),
          0.8,
          1.2,0.1,
          0.1)

mepp.tent.breast = fit.mepp.cf(L = breast$left, R = breast$right, n.int = n.intervalos,
                               cov.risco = cbind(breast$ther),
                               cov.cura = cbind(1, breast$ther),
                               start = chute)

## analise residual

parametros = mepp.tent.breast$estimated

lambdas.est.breast = parametros[1:n.intervalos]
alpha.est.breast = parametros[n.intervalos+1]
b.est.breast = parametros[(n.intervalos+2):(n.intervalos+3)]
betas.est.breast = parametros[n.intervalos+4]

grid.obs.breast=time.grid.interval(li=breast$left, ri=breast$right, type="OBS", bmax= n.intervalos)

grid.obs.breast=grid.obs.breast[-c(1, length(grid.obs.breast))]

## residuos de martingales

spop_mepp_l = SpopMEPP(t = breast$left, lambda.par = lambdas.est.breast, alpha.par = alpha.est.breast,
                       theta.par = b.est.breast, beta.par = betas.est.breast, 
                       x.cure = cbind(1, breast$ther), x.risk = cbind(breast$ther),
                       grid.vet = grid.obs.breast)

spop_mepp_r = SpopMEPP(t = breast$right, lambda.par = lambdas.est.breast, alpha.par = alpha.est.breast,
                       theta.par = b.est.breast, beta.par = betas.est.breast, 
                       x.cure = cbind(1, breast$ther), x.risk = cbind(breast$ther),
                       grid.vet = grid.obs.breast)

rm.breast = (spop_mepp_l*log(spop_mepp_l) - spop_mepp_r*log(spop_mepp_r))/(spop_mepp_l-spop_mepp_r)

r.Deviance.breast  = sign(rm.breast)*(-2*(rm.breast+log(rm.breast)))^(0.5)

## smoke_cessation_Bannerge2009

smoke2009 = read.table('smoke_cessation_Bannerge2009.txt', header = T)

dim(smoke2009)
head(smoke2009)

tempo.aval = seq(0,10,0.1)

delta = ifelse(smoke2009$Timept2 == Inf, 0, 1)

## ----
## Analise descritiva dos dados
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
## ----

npmle_fit <- ic_np(cbind(Timept1, Timept2) ~ 0, data = smoke2009)

plot(npmle_fit)

## ------
## O modelo Weibull
## ------

s.weibull = ic_par(formula = cbind(Timept1, Timept2) ~ (SIUC + F10Cigs + Duration + SexF),
                   data = smoke2009, model = "ph", dist = "weibull")

s.weibull$coefficients

#lines(s.weibull, col = 'red')


AIC.surv(loglik = s.weibull$llk, n.param = 6)
BIC.surv(loglik = s.weibull$llk, n.param = 6, n.sample = dim(smoke2009)[1])
HC.surv(loglik = s.weibull$llk, n.param = 6,n.sample = dim(smoke2009)[1])


## variaveis do estudo

tratamento = smoke2009$SIUC # Tipo de tratamento (SI/UC)
n_cigarros = smoke2009$F10Cigs_pad # Numero de cigarros fumados por dia, normalizado
duracao_dependente = smoke2009$Duration_pad # Duracao como dependente, normalizado
sexo = smoke2009$SexF

covariaveis = cbind(tratamento, sexo, n_cigarros, duracao_dependente)

## ---
## modelo exponencial por partes com fracao de cura
## ---

source('C:/Users/NetoDavi/Desktop/survival_pibic/mep_interval_fc.R')

covariaveis = cbind(smoke2009$SIUC, smoke2009$SexF,
                    smoke2009$Duration_pad, smoke2009$F10Cigs_pad)

n.int = 4

l.smoke= smoke2009$Timept1/7
r.smoke= smoke2009$Timept2/7

x.c = cbind(1,covariaveis)
x.f = covariaveis
  
grid.obs=time.grid.interval(l = l.smoke, r = r.smoke,
                            type="OBS", bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]


chutes = c(5,5,200,200,
           1,0.3,0.1,0.1,0.1,
           0.1,0.1,0.1,0.1)

max.mep = optim(par = chutes, fn=loglikIC.MEP.fc, gr = NULL, method = "BFGS",
                control=list(fnscale=1), hessian = TRUE,
                l = l.smoke, r = r.smoke,
                x.cure=x.c, x.risk= x.f, 
                grid.vet=grid.obs)
max.mep$par
#max.mep$hessian
max.mep$value

AIC.surv(loglik = max.mep$value*-1,
         n.param = length(max.mep$par))

BIC.surv(loglik = max.mep$value*-1,
         n.param = length(max.mep$par),
         n.sample = dim(smoke2009)[1])

HC.surv(loglik = max.mep$value*-1,
        n.param = length(max.mep$par),
        n.sample = dim(smoke2009)[1])

## ---
## modelo exponencial por partes potencia sem fracao de cura
## ---
n.int = 5

grid.obs=time.grid.interval(li=dadosIC$L, ri=dadosIC$R, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,1,1,1,1,
           1,
           0.5,0.5)

max.mepp.int = optim(par = chutes, fn=loglikIC.mepp.int,
                     gr = NULL, method = "BFGS",
                     control=list(fnscale=1), hessian = TRUE,
                     l = dadosIC$L, r = dadosIC$R,
                     x.cov=x.f, 
                     grid=grid.obs)

max.mepp.int$par

## ---
## modelo exponencial por partes potencia com fracao de cura
## ---

## lambda = 4
## alpha = 1
## beta.cura = 5
## beta.risco = 4

n.intervalos = 4

chute = c(rep(0.5, n.intervalos),
          10,
          1.2,0.1,0.1,0.1,0.1,
          0.1,0.1,0.1,0.1)
  

mepp.tent = fit.mepp.cf(L = smoke2009$Timept1/7,
                        R = smoke2009$Timept2/7,
                        n.int = n.intervalos,
                        cov.risco = covariaveis,
                        cov.cura = cbind(1, covariaveis),
                        start = chute)

mepp.tent$estimated
mepp.tent$loglik
 
AIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated))

BIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
         n.sample = dim(smoke2009)[1])

HC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
        n.sample = dim(smoke2009)[1])

## ---------------------------------------------------------------
## dados de dentes
## ---------------------------------------------------------------

## tooth

# tooth = read.table('tooth.txt', header = T)
# # 
# dim(tooth)
# head(tooth) ## para l = tooth$left, r = tooth$rightInf
# # 
# 
# delta.tooth = ifelse(tooth$rightInf == Inf, 0, 1)
# 
# ## ----
# ## Analise descritiva dos dados
# ## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
# ## ----
# 
# npmle_fit_tooth <- ic_np(cbind(left, rightInf) ~ 0, data = tooth)
# 
# plot(npmle_fit_tooth)

## ---------------------------------------------------------------
## dados tandmobAll_icensBKL
## ---------------------------------------------------------------

# tandmoball = read.table('tandmobAll_icensBKL.txt', header = T)
# # 
# dim(tandmoball)
# head(tandmoball) ## para l = tooth$left, r = tooth$rightInf
# 
# 
# delta.tooth = ifelse(tooth$rightInf == Inf, 0, 1)


## ---------------------------------------------------------------
## dados aneurysm
## ---------------------------------------------------------------

## dados do caso tipo 1

aneurysm = read.table('aneurysm.txt', header = T)

dim(aneurysm)
head(aneurysm) ## para l = tooth$left, r = tooth$rightInf

aneurysm$t.right = ifelse(is.na(aneurysm$t.right), Inf, aneurysm$t.right)

delta.aneurysm = ifelse(aneurysm$t.right == Inf, 0, 1)


## ----
## Analise descritiva dos dados
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
## ----

npmle_fit_aneurysm <- ic_np(cbind(t.left, t.right) ~ 0,
                            data = aneurysm)

plot(npmle_fit_aneurysm) ## temos o plato !!


## ------
## O modelo Weibull
## ------

# o z e indicador de censura 

aneurysm$gr_pad = (aneurysm$gr - mean(aneurysm$gr))/sd(aneurysm$gr)

a.weibull = ic_par(formula = cbind(t.left, t.right) ~ (mo+gr_pad+lok),
                   data = aneurysm, model = "ph", dist = "weibull")

a.weibull$coefficients

lines(a.weibull, col = 'red')


AIC.surv(loglik = a.weibull$llk, 
         n.param = length(a.weibull$coefficients))

BIC.surv(loglik = a.weibull$llk,
         n.param = length(a.weibull$coefficients),
         n.sample = dim(aneurysm)[1])

HC.surv(loglik = a.weibull$llk, 
        n.param = length(a.weibull$coefficients),
        n.sample = dim(aneurysm)[1])

## ---
## modelo exponencial por partes potencia com fracao de cura
## ---

## lambda = 4
## alpha = 1
## beta.cura = 4
## beta.risco = 3

gr_pad = (aneurysm$gr - mean(aneurysm$gr))/sd(aneurysm$gr) ## precisa normalizar para estimar

covariaveis.aneurysm = cbind(aneurysm$mo, gr_pad, aneurysm$lok)

a.l = aneurysm$t.left
a.r = aneurysm$t.right

n.int = 8

grid.obs=time.grid.interval(l = a.l, r = a.r,
                            type="OBS", bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]


chute = c(5,10,5,5,5,5,5,5,
          0.2,
          1.2,0.1,0.1,0.1,
          0.1,0.1,0.1)

max.mepp.fc = optim(par = chute, fn=loglikIC, gr = NULL, method = "BFGS",
                    control=list(fnscale=1), hessian = TRUE,
                    l = a.l, r = a.r,
                    x.cure= cbind(1,covariaveis.aneurysm),
                    x.risk= covariaveis.aneurysm, 
                    grid=grid.obs)

max.mepp.fc$par

AIC.surv(max.mepp.fc$value*-1,n.param = length(max.mepp.fc$par))

BIC.surv(max.mepp.fc$value*-1, n.param = length(max.mepp.fc$par),
         n.sample = dim(aneurysm)[1])

HC.surv(max.mepp.fc$value*-1, n.param = length(max.mepp.fc$par),
        n.sample = dim(aneurysm)[1])


## ---------------------------------------------------------------
## dados aidsCT_icensBKL
## ---------------------------------------------------------------


aids_CT = read.table('aidsCT_icensBKL.txt', header = T)

dim(aids_CT)
head(aids_CT) ## para l = tooth$left, r = tooth$rightInf


## ---------------------------------------------------------------
## dados aidsCohort_icensBKL
## ---------------------------------------------------------------

## parece um caso para censura intervalar bivariada
## com o caso Y e Z

aidscohort = read.table('aidsCohort_icensBKL.txt', header = T)

dim(aidscohort)
head(aidscohort) 

## troca dos limites inferior (L) e superior (R) de Y
aidscohort$L.Y = ifelse(is.na(aidscohort$L.Y), 0, aidscohort$L.Y)
aidscohort$R.Y = ifelse(is.na(aidscohort$R.Y), Inf, aidscohort$R.Y)

## troca dos limites inferior (L) e superior (R) de Z
aidscohort$L.Z = ifelse(is.na(aidscohort$L.Z), 0, aidscohort$L.Z)
aidscohort$R.Z = ifelse(is.na(aidscohort$R.Z), 999999, aidscohort$R.Z)

## codificando para zero e um a variavel idade
aidscohort$age = ifelse(aidscohort$age==2,1,0)

## ----
## Analise descritiva dos dados
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
## ----

npmle_fit_aids_Y <- ic_np(cbind(L.Y, R.Y) ~ 0,
                            data = aidscohort)

plot(npmle_fit_aids_Y) ## temos o plato !!

# ------------------------------------------------------------------------------

npmle_fit_aids_Z <- ic_np(cbind(L.Z, R.Z) ~ 0,
                          data = aidscohort, B = c(1,1))

plot(npmle_fit_aids_Z) ## temos o plato!! Bem alto por sinal!

## -----
## Analise e estimacao para os dados do caso Z 
## -----

## ------
## O modelo Weibull
## ------

aidsz.weibull = ic_par(formula = cbind(L.Z, R.Z) ~ (age+group),
                   data = aidscohort, model = "ph", dist = "weibull")

aidsz.weibull$coefficients

AIC.surv(loglik = aidsz.weibull$llk, n.param = length(aidsz.weibull$coefficients))
BIC.surv(loglik = aidsz.weibull$llk, n.param = length(aidsz.weibull$coefficients), n.sample = dim(aidscohort)[1])
HC.surv(loglik = aidsz.weibull$llk, n.param = length(aidsz.weibull$coefficients),n.sample = dim(aidscohort)[1])


## ---
## modelo exponencial por partes potencia com fracao de cura
## ---

ll.aidsz = aidscohort$L.Z/7
rr.aidsz = aidscohort$R.Z/7

n.int = 1

x.f <- cbind(aidscohort$age, aidscohort$group)
x.c <- cbind(1, aidscohort$age, aidscohort$group)

grid.obs=time.grid.interval(li=ll.aidsz, ri=rr.aidsz, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(10,
           200,
           1,0.1,0.1,
           0.1,0.1)

max.mepp = optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                 control=list(fnscale=1), hessian = TRUE,
                 l = ll.aidsz, r = rr.aidsz,
                 x.cure=x.c, x.risk=x.f, 
                 grid=grid.obs)


## -----
## Analise e estimacao para os dados do caso Y 
## -----

## ------
## O modelo Weibull
## ------

aidsy.weibull = ic_par(formula = cbind(L.Y, R.Y) ~ (age+group),
                   data = aidscohort, model = "ph", dist = "weibull")

aidsy.weibull$coefficients

AIC.surv(loglik = aidsy.weibull$llk, n.param = length(aidsy.weibull$coefficients))
BIC.surv(loglik = aidsy.weibull$llk, n.param = length(aidsy.weibull$coefficients), n.sample = dim(aidscohort)[1])
HC.surv(loglik = aidsy.weibull$llk, n.param = length(aidsy.weibull$coefficients),n.sample = dim(aidscohort)[1])


## ---
## modelo exponencial por partes potencia com fracao de cura
## ---

aidsy.l = aidscohort$L.Y
aidsy.r = aidscohort$R.Y
cov.aidsy = cbind(aidscohort$age,aidscohort$group)

n.int = 2

grid.obs=time.grid.interval(l = aidsy.l, r = aidsy.r,
                            type="OBS", bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(2,6,
          2,
          12,0.1,0.1,
          0.1,0.1)

max.mepp.fc = optim(par = chutes, fn=loglikIC, gr = NULL,
                    method = "Nelder-Mead",
                    control=list(fnscale=1), hessian = F,
                    l = aidsy.l, r = aidsy.r,
                    x.cure= cbind(1,cov.aidsy),
                    x.risk= cov.aidsy, 
                    grid=grid.obs)


## ---------------------------------------------------------------
## dados HIV_SUM_set2 
## ---------------------------------------------------------------


hiv2 = read.table('HIV_SUM_set2.txt', header = T)

dim(hiv2)
head(hiv2) 


## ----
## Analise descritiva dos dados
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
## ----

#npmle_fit_hiv2 <- ic_np(cbind(Li, Ri) ~ 0,
#                          data = hiv2)

#plot(npmle_fit_hiv2) ## temos o plato !!


## ------
## O modelo Weibull
## ------

hiv2 = read.table('HIV_SUM_set2.txt', header = T)

hiv2.weibull = ic_par(cbind(Li/7, Ri/7) ~ DoseType,
                      data = hiv2, model = "ph", dist = "weibull")

hiv2.weibull$coefficients

AIC.surv(loglik = hiv2.weibull$llk, n.param = length(hiv2.weibull$coefficients))

BIC.surv(loglik = hiv2.weibull$llk, n.param = length(hiv2.weibull$coefficients),
         n.sample = dim(hiv2)[1])

HC.surv(loglik = hiv2.weibull$llk, n.param = length(hiv2.weibull$coefficients),
        n.sample = dim(hiv2)[1])

## ---
## modelo exponencial por partes com fracao de cura
## ---

source('C:/Users/NetoDavi/Desktop/survival_pibic/mep_interval_fc.R')

hiv2 = read.table('HIV_SUM_set2.txt', header = T)

# Misturar as linhas

hiv2 <- hiv2[sample(row.names(hiv2)), ]


left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = hiv2$DoseType

n.int = 8

x.f <- cbind(x1=cov)
#x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,0.1,1,2,2,2,2,5,
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

AIC.surv(loglik = max.mep$value*-1,
         n.param = length(max.mep$par))

BIC.surv(loglik = max.mep$value*-1,
         n.param = length(max.mep$par),
         n.sample = dim(hiv2)[1])

HC.surv(loglik = max.mep$value*-1,
        n.param = length(max.mep$par),
        n.sample = dim(hiv2)[1])

## ---
## modelo exponencial por partes potencia sem a fracao de cura
## ---

source('C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste/proposta_MEPP_int_sem_cura.R')

hiv2 = read.table('HIV_SUM_set2.txt', header = T)
hiv2 <- hiv2[sample(row.names(hiv2)), ]

left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = cbind(x1=hiv2$DoseType)

n.int = 5

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,0.1,1,5,1,
           5,
           0.5)

max.mepp.int = optim(par = chutes, fn=loglikIC.mepp.int,
                     gr = NULL, method = "BFGS",
                     control=list(fnscale=1), hessian = TRUE,
                     l = left, r = right,
                     x.cov=cov, 
                     grid=grid.obs)

max.mepp.int$par
max.mepp.int$par



AIC.surv(loglik = max.mepp.int$value*-1,
         n.param = length(max.mepp.int$par))

BIC.surv(loglik = max.mepp.int$value*-1,
         n.param = length(max.mepp.int$par),
         n.sample = dim(hiv2)[1])

HC.surv(loglik = max.mepp.int$value*-1,
        n.param = length(max.mepp.int$par),
        n.sample = dim(hiv2)[1])





## ---
## modelo exponencial por partes potencia com fracao de cura
## ---
source('C:/Users/NetoDavi/Desktop/survival_pibic/mep_interval_fc.R')

left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = hiv2$DoseType

n.int = 5

x.f <- cbind(x1=cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,0.2,1,5,4.5,
           5, 
           1, 0.5,
           0.5)

estimacao <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
              control=list(fnscale=1), hessian = TRUE, l=left,
              r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)

estimacao$par

AIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par))
BIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
        n.sample = nrow(hiv2))
HC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
        n.sample = nrow(hiv2))

## ---------------------------------------------------------------
## dados yoghurt_icensBKL 
## ---------------------------------------------------------------


yoghurt = read.table('yoghurt_icensBKL.txt', header = T)

dim(yoghurt)
head(yoghurt) 

yoghurt$left = ifelse(is.na(yoghurt$left), 0, yoghurt$left)
yoghurt$right = ifelse(is.na(yoghurt$right), Inf, yoghurt$right)

## ----
## Analise descritiva dos dados
## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)
## ----

npmle_fit_yoghurt <- ic_np(cbind(left, right) ~ 0,
                        data = yoghurt, B = c(1,1))

plot(npmle_fit_yoghurt) ## temos o plato !!

## ------
## O modelo Weibull
## ------

yoghurt.weibull = ic_par(formula = cbind(left, right) ~ (adult),
                       data = yoghurt, model = "ph", dist = "weibull")

yoghurt.weibull$coefficients

AIC.surv(loglik = yoghurt.weibull$llk,
         n.param = length(yoghurt.weibull$coefficients))

BIC.surv(loglik = yoghurt.weibull$llk,
         n.param = length(yoghurt.weibull$coefficients),
         n.sample = dim(yoghurt)[1])

HC.surv(loglik = yoghurt.weibull$llk,
        n.param = length(yoghurt.weibull$coefficients),
        n.sample = dim(yoghurt)[1])


intervalos = 7

chute = c(rep(0.1,intervalos),
          1,
          1.2,0.1,
          0.1)


mepp.tentativa.yoghurt = fit.mepp.cf(L = yoghurt$left, R = yoghurt$right ,
                                  n.int = intervalos,
                                  cov.risco = cbind(yoghurt$adult),
                                  cov.cura = cbind(1, yoghurt$adult),
                                  start = chute)

mepp.tentativa.yoghurt$estimated

AIC.surv(loglik = mepp.tentativa.yoghurt$loglik,
         n.param = length(mepp.tentativa.yoghurt$estimated))

BIC.surv(loglik = mepp.tentativa.yoghurt$loglik,
         n.param = length(mepp.tentativa.yoghurt$estimated),
         n.sample = dim(yoghurt)[1])

HC.surv(loglik = mepp.tentativa.yoghurt$loglik,
        n.param = length(mepp.tentativa.yoghurt$estimated),
        n.sample = dim(yoghurt)[1])







