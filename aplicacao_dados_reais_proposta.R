### -------
## Validacao es estudo dos modelos de sobrevivencia
## em dados reais
## Data: 30/05/2023 
### -------

## propostas para o modelo Weibull, MEP e MEPP


if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg)

#setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')
setwd("C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste")
## hemofilia_icens
hemo.icens = read.table('hemofilia_icens.txt', header = T)

dim(hemo.icens)
head(hemo.icens)

colnames(hemo.icens)

## covariaveis: NoDose, Medium, High

## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)

turnbull_fit = ic_np(cbind(L, R)~0, data = hemo.icens)

plot(turnbull_fit)

## ---
## Ajuste do modelo Weibull aos dados
## ---

weibull = ic_par(formula = cbind(L, R) ~ High, data = hemo.icens,
       model = "ph", dist = "weibull")

weibull$coefficients

plot(weibull)


AIC.surv(loglik = weibull$llk, n.param = 3)
BIC.surv(loglik = weibull$llk, n.param = 3, n.sample = dim(hemo.icens)[1])
HC.surv(loglik = weibull$llk, n.param = 3,n.sample = dim(hemo.icens)[1])



## Modelo Exponencial por partes com censura intervalar


## Modelo Exponencial por partes com censura intervalar e fracao de cura

## Modelo Proposto (MEPP + fracao de cura)

hemo.icens$L
hemo.icens$R
hemo.icens$Low

chute = c(0.1,0.1,0.1, 0.8, 1.2,0.3,0.3)

x2_falso = rnorm(n = dim(hemo.icens)[1])

# nao convergiu
mepp.tent = fit.mepp.cf(L = hemo.icens$L, R = hemo.icens$R, n.int = 2,
                        cov.risco = cbind(hemo.icens$Low, x2_falso),
                        cov.cura = cbind(1, hemo.icens$Low, x2_falso),
                        start = chute)

# mepp.tent$estimated
# mepp.tent$hessian
# mepp.tent$loglik

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

covariaveis = cbind(tratamento, n_cigarros, duracao_dependente, sexo)

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

chute = c(0.1,0.1,0.1,
          8,
          1.2,-0.5,0.2,-0.3,0.4,
          0.4,-0.5,0.3,0.3)
  
mepp.tent = fit.mepp.cf(L = smoke2009$Timept1, R = smoke2009$Timept2, n.int = 3,
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
# 
# dim(tooth)
# head(tooth) ## para l = tooth$left, r = tooth$rightInf
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
# 
# dim(tandmoball)
# head(tandmoball) ## para l = tooth$left, r = tooth$rightInf


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

chute = c(0.1,0.1,0.1,0.1,
          1,
          1.2,0.1,0.1,0.1,
          0.1,0.1,0.1)

gr_pad = (aneurysm$gr - mean(aneurysm$gr))/sd(aneurysm$gr) ## precisa normalizar para estimar

covariaveis.aneurysm = cbind(aneurysm$mo, gr_pad, aneurysm$lok)

mepp.tent = fit.mepp.cf(L = aneurysm$t.left, R = aneurysm$t.right, n.int = 4,
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

## lambda = 2
## alpha = 1
## beta.cura = 4
## beta.risco = 3

chute = c(0.1,0.1,0.1,
          1,
          1.2,0.1,0.1,
          0.1,0.1)

covariaveis.Z = cbind(aidscohort$age, aidscohort$group)
  
mepp.tent = fit.mepp.cf(L = aidscohort$L.Z, R = aidscohort$R.Z, n.int = 3,
                        cov.risco = covariaveis.Z ,
                        cov.cura = cbind(1, covariaveis.Z ),
                        start = chute)

mepp.tent$estimated
mepp.tent$loglik

AIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated))

BIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
         n.sample = dim(aneurysm)[1])

HC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
        n.sample = dim(aneurysm)[1])

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

## lambda = 2
## alpha = 1
## beta.cura = 4
## beta.risco = 3

chute = c(0.1,0.1,0.1,
          1,
          1.2,0.1,0.1,
          0.1,0.1)

covariaveis.Z = cbind(aidscohort$age, aidscohort$group)
  
mepp.tent = fit.mepp.cf(L = aidscohort$L.Z, R = aidscohort$R.Z, n.int = 3,
                        cov.risco = covariaveis.Z ,
                        cov.cura = cbind(1, covariaveis.Z ),
                        start = chute)

mepp.tent$estimated
mepp.tent$loglik

AIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated))

BIC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
         n.sample = dim(aneurysm)[1])

HC.surv(mepp.tent$loglik, n.param = length(mepp.tent$estimated),
        n.sample = dim(aneurysm)[1])











