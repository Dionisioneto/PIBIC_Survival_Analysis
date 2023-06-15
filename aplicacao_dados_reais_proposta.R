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

hemo.icens$L= hemo.icens$L/30
hemo.icens$R = hemo.icens$R/30

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

AIC.surv(loglik = weibull$llk, n.param = 3)
BIC.surv(loglik = weibull$llk, n.param = 3, n.sample = dim(hemo.icens)[1])
HC.surv(loglik = weibull$llk, n.param = 3,n.sample = dim(hemo.icens)[1])



## Modelo Exponencial por partes com censura intervalar


## Modelo Exponencial por partes com censura intervalar e fracao de cura

## Modelo Proposto (MEPP + fracao de cura)


n.intervalos = 3

chute = c(rep(0.1,n.intervalos),
          1,
          1,0.1,
          0.1)

mepp.tent.hemo = fit.mepp.cf(L = hemo.icens$L , R = hemo.icens$R, 
                             n.int = n.intervalos,
                            cov.risco = cbind(hemo.icens$NoDose),
                            cov.cura = cbind(1,hemo.icens$NoDose),
                            start = chute)

AIC.surv(loglik = mepp.tent.hemo$loglik, n.param = length(mepp.tent.hemo$estimated))

BIC.surv(loglik = mepp.tent.hemo$loglik, n.param = length(mepp.tent.hemo$estimated), n.sample = dim(hemo.icens)[1])

HC.surv(loglik = mepp.tent.hemo$loglik, n.param = length(mepp.tent.hemo$estimated), n.sample = dim(hemo.icens)[1])

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

lines(s.weibull, col = 'red')


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
  

mepp.tent = fit.mepp.cf(L = smoke2009$Timept1,
                        R = smoke2009$Timept2,
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

a.weibull = ic_par(formula = cbind(t.left, t.right) ~ (mo+gr+lok),
                   data = aneurysm, model = "ph", dist = "weibull")

a.weibull$coefficients

lines(a.weibull, col = 'red')


AIC.surv(loglik = a.weibull$llk, n.param = length(a.weibull$coefficients))
BIC.surv(loglik = a.weibull$llk, n.param = length(a.weibull$coefficients), n.sample = dim(aneurysm)[1])
HC.surv(loglik = a.weibull$llk, n.param = length(a.weibull$coefficients),n.sample = dim(aneurysm)[1])

## ---
## modelo exponencial por partes potencia com fracao de cura
## ---

## lambda = 4
## alpha = 1
## beta.cura = 4
## beta.risco = 3

numero.intervalos = 2

chute = c(rep(0.1,numero.intervalos),
          1,
          1.2,0.1,0.1,
          0.1,0.1)

gr_pad = (aneurysm$gr - mean(aneurysm$gr))/sd(aneurysm$gr) ## precisa normalizar para estimar

covariaveis.aneurysm = cbind(aneurysm$mo, gr_pad, aneurysm$lok)

a.l = aneurysm$t.left
a.r = aneurysm$t.right

mepp.tent = fit.mepp.cf(L = a.l, R = a.r, n.int = numero.intervalos,
                        cov.risco = covariaveis.aneurysm,
                        cov.cura = cbind(1, covariaveis.aneurysm),
                        start = chute)

mepp.tent$estimated
mepp.tent$loglik

AIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated))

BIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
         n.sample = dim(aneurysm)[1])

HC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
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
aidscohort$R.Z = ifelse(is.na(aidscohort$R.Z), Inf, aidscohort$R.Z)

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
                          data = aidscohort)

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

## lambda = 3
## alpha = 1
## beta.cura = 3
## beta.risco = 2

numero.intervalos = 2

chute = c(rep(0.1,numero.intervalos ),
          1,
          1.2,0.1,0.1,
          0.1,0.1)

covariaveis.Z = cbind(aidscohort$age, aidscohort$group)

llz = as.numeric(aidscohort$L.Z)/7 
rrz = as.numeric(aidscohort$R.Z)/7 

mepp.tent.z = fit.mepp.cf(L = llz, R = rrz, n.int = numero.intervalos ,
                        cov.risco = covariaveis.Z,
                        cov.cura = cbind(1, covariaveis.Z ),
                        start = chute)

mepp.tent.z$estimated
mepp.tent.z$loglik

AIC.surv(mepp.tent.z$loglik, n.param = length(mepp.tent.z$estimated))

BIC.surv(mepp.tent.z$loglik, n.param = length(mepp.tent.z$estimated),
         n.sample = dim(aidscohort)[1])

HC.surv(mepp.tent.z$loglik, n.param = length(mepp.tent.z$estimated),
        n.sample = dim(aidscohort)[1])

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

## lambda = 2
## alpha = 1
## beta.cura = 4
## beta.risco = 3

numero.intervalos = 12

chute = c(rep(0.1,numero.intervalos ),
          1,
          1.2,0.1,0.1,
          0.1,0.1)

covariaveis.y = cbind(aidscohort$age, aidscohort$group)

lly = as.numeric(aidscohort$L.Y) + 0.0000010
rry = as.numeric(aidscohort$R.Y) + 0.0000011

mepp.tent.y = fit.mepp.cf(L = lly , R = rry, n.int = numero.intervalos,
                        cov.risco = covariaveis.y ,
                        cov.cura = cbind(1, covariaveis.y),
                        start = chute)
mepp.tent.y$estimated
mepp.tent.y$loglik

AIC.surv(mepp.tent.y$loglik, n.param = length(mepp.tent.y$estimated))

BIC.surv(mepp.tent.y$loglik, n.param = length(mepp.tent.y$estimated),
         n.sample = dim(aidscohort)[1])

HC.surv(mepp.tent.y$loglik, n.param = length(mepp.tent.y$estimated),
        n.sample = dim(aidscohort)[1])


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

npmle_fit_hiv2 <- ic_np(cbind(Li, Ri) ~ 0,
                          data = hiv2)

plot(npmle_fit_hemo_teste) ## temos o plato !!

left = hiv2$Li + 0.000001
rigth = hiv2$Ri + 0.000002
cov = hiv2$DoseType
  
intervalos = 7

chute = c(rep(0.1,intervalos),
          1,
          1.2,0.1,
          0.1)

mepp.tentativa.hemo = fit.mepp.cf(L = left, R = rigth , n.int = intervalos,
                          cov.risco = cbind(cov), cov.cura = cbind(1, cov),
                          start = chute)

mepp.tentativa.hemo$estimated


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







