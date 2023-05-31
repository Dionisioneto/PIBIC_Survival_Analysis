### -------
## Validacao es estudo dos modelos de sobrevivencia
## em dados reais
## Data: 30/05/2023 
### -------

## propostas para o modelo Weibull, MEP e MEPP


if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(icenReg)

setwd('C:\\Users\\Dionisio\\Desktop\\Dionisio_Neto\\PIBIC_Survival_Analysis\\dados_para_teste')

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




## Modelo Exponencial por partes



## Modelo Exponencial por partes

## Modelo Proposto (MEPP + fracao de cura)


## smoke_cessation_Bannerge2009

smoke2009 = read.table('smoke_cessation_Bannerge2009.txt', header = T)

dim(smoke2009)
head(smoke2009)

tempo.aval = seq(0,10,0.1)

delta = ifelse(smoke2009$Timept2 == Inf, 0, 1)

## estimador de turnbull ou Non-Parametric Maximum Likelihood Estimator (NPMLE)

npmle_fit <- ic_np(cbind(Timept1, Timept2) ~ 0, data = smoke2009)

plot(npmle_fit)





















