## ------
## Aplicacao de modelos probabilisticos em analise de sobrevivencia
## Apresentacao 1
## ------


if(!require(pacman)) install.packages('pacman'); library(pacman)

p_load(flexsurv,survival, eha, ICglm, latex2exp)

## carregar o script de funcoes
source('/cloud/project/funcoes_sobrevivencia_pibic2023.R')

head(lung)
dim(lung)

unique(lung$status)


tempo = lung$time
censura = ifelse(lung$status == 1, 0, 1)

## ------
## 2. Ajuste do modelo nao-parametrico de Kaplan-Meier
## ------

km = survfit(Surv(tempo, censura)~1)
km

# intervalos de tempo do estimado km
int_tempo = km$time

# Funcao de sobrevivencia para cada intervalo do Kaplan-Meier
sobrevivencia_km = km$surv



# Visualizacao
jpeg("estimador_km.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo", lwd = 2)
dev.off()

## ------
## 1. Modelos probabilisticos para o ajuste dados
## ------

## modelo exponencial

fit_exp = survreg(Surv(tempo, censura)~1, dist = 'exponential')
fit_exp

# valor do parametro theta do modelo exponencial
lambda = exp(fit_exp$coefficients[1]) ## lambda
lambda  # valor do parametro para a funcao 


st_exp = exp(-int_tempo/lambda)

# Visualizacao

jpeg("ajuste_exponencial.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_exp, col = 'red', lwd = 2)

legend('topright', "Ajuste Exponencial",col = "red", lty = 1, bty = 'n')
dev.off()

## modelo weibull
fit_weibull = survreg(Surv(tempo, censura)~1, dist = "weibull")
fit_weibull

alpha = exp(coef(fit_weibull))
gamma = 1/fit_weibull$scale

st_wb = exp(-(int_tempo/alpha)^gamma)
st_wb 

# Visualizacao

jpeg("ajuste_weibull.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_wb, col = 'purple', lwd = 2)

legend('topright', "Ajuste Weibull",col = "purple", lty = 1, bty = 'n')
dev.off()

## modelo log-normal
fit_lognorm = survreg(Surv(tempo, censura)~1, dist = "lognormal")
fit_lognorm

media = fit_lognorm$coefficients
dp = fit_lognorm$scale

st_lognorm = pnorm(-(log(int_tempo) - media)/dp)


jpeg("ajuste_lognorm.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_lognorm, col = 'forestgreen', lwd = 2)
legend('topright', "Ajuste Log-Normal",col = "forestgreen", lty = 1, bty = 'n')
dev.off()



## modelo log-logistico

fit_loglog = survreg(Surv(tempo, censura)~1, dist = "loglogistic")
fit_loglog

a = -(fit_loglog$coefficients/fit_loglog$scale)
b = 1/fit_loglog$scale

cbind(a,b)

st_loglog = 1/(1 + (exp(a) * (int_tempo)^b))


jpeg("ajuste_loglogistico.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_loglog, col = 'magenta', lwd = 2)

legend('topright', "Ajuste Log-Log?stico",col = "magenta", lty = 1, bty = 'n')
dev.off()

## modelo gamma
fit_gamma = flexsurvreg(Surv(tempo, censura)~1, dist = "gamma")
fit_gamma

alpha_gamma = exp(fit_gamma$coefficients[1])
beta_gamma = exp(fit_gamma$coefficients[2])

st_gamma = 1 - pgamma(int_tempo, alpha_gamma, beta_gamma)

jpeg("ajuste_gamma.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_gamma, col = 'dodgerblue4', lwd = 2)
legend('topright', "Ajuste Gama",
       col = "dodgerblue4", lty = 1, bty = 'n')
dev.off()

## modelo gama-generalizado
fit_gammagen = flexsurvreg(Surv(tempo, censura)~1, dist = "gengamma")
fit_gammagen

mu = fit_gammagen$coefficients[1]
sigma = exp(fit_gammagen$coefficients[2])
k = fit_gammagen$coefficients[3]

mu;sigma;k


st_gammagen = 1 - pgengamma(int_tempo, mu, sigma, k)

jpeg("ajuste_gamma_generalizado.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_gammagen, col = 'goldenrod3', lwd = 2)
legend('topright', "Ajuste Gama generalizado",col = "goldenrod3", lty = 1, bty = 'n')
dev.off()


jpeg("ajuste_modelos.jpeg", width = 1000, height = 700)
plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo")

lines(int_tempo, st_exp, col = 'red', lwd = 2)
lines(int_tempo, st_wb, col = 'purple', lwd = 2)
lines(int_tempo, st_lognorm, col = 'forestgreen', lwd = 2)
lines(int_tempo, st_loglog, col = 'magenta', lwd = 2)
lines(int_tempo, st_gamma, col = 'dodgerblue4', lwd = 2)
lines(int_tempo, st_gammagen, col = 'goldenrod3', lwd = 2)

legend('topright', c("Exponencial", "Weibull", "Log-Normal", 
                     "Log-Logistico", "Gamma", "Gamma generalizado")
       ,col = c('red','purple', 'forestgreen',  'magenta',
                'dodgerblue4', "goldenrod3"), lty = 1)
dev.off()


## Calculos do AIC, BIC e HC

# ---
# AIC
# ---

# Exponencial
AIC.surv(loglik = fit_exp$loglik[1], n.param = 1) 

# Weibull
AIC.surv(loglik = fit_weibull$loglik[1], n.param = 2)

# Log-Normal
AIC.surv(loglik = fit_lognorm$loglik[1], n.param = 2)

# Log-Logistico
AIC.surv(loglik = fit_loglog$loglik[1], n.param = 2)

# Gama
AIC.surv(loglik = fit_gamma$loglik[1], n.param = 2)

# Gama generalizado
AIC.surv(loglik = fit_gammagen$loglik[1], n.param = 3)


# ---
# BIC
# ---

# Exponencial
BIC(loglik = fit_exp$loglik[1],
    n.param = 1,
    n.sample= summary(fit_exp)$n)

# Weibull
BIC(loglik = fit_weibull$loglik[1],
    n.param = 2,
    n.sample = summary(fit_weibull)$n)

# Log-Normal
BIC(loglik = fit_lognorm$loglik[1],
    n.param = 2,
    n.sample = summary(fit_lognorm)$n)

# Log-Logistico
BIC(loglik = fit_loglog$loglik[1],
    n.param = 2,
    n.sample = summary(fit_loglog)$n)

# Gama
BIC(loglik = fit_gamma$loglik,
    n.param = 2,
    n.sample = fit_gamma$N)

# Gama generalizado
BIC(loglik =fit_gammagen$loglik,
    n.param = 3,
    n.sample = fit_gammagen$N)

# HC

# Exponencial
HC(loglik = fit_exp$loglik[1],
   n.param= length(fit_exp$coefficients),
   n.sample = summary(fit_exp)$n)

# Weibull
HC(loglik = fit_weibull$loglik[1],
   n.param= 2,
   n.sample = summary(fit_weibull)$n)

# Log-Normal
HC(loglik = fit_lognorm$loglik[1],
   n.param = 2,
   n.sample = summary(fit_lognorm)$n)

# Log-Logistico
HC(loglik = fit_loglog$loglik[1],
   n.param = 2,
   n.sample = summary(fit_loglog)$n)

# Gama
HC(loglik = fit_gamma$loglik,
   n.param = fit_gamma$npars,
   n.sample = fit_gamma$N)

# Gama generalizado
HC(loglik = fit_gammagen$loglik,
   n.param = fit_gammagen$npars,
   n.sample = fit_gammagen$N)


## ---
## modelo exponencial por partes
## ---

## visualizar o o estimador KM para escolher os intervalos

plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo", lwd = 2)

#tempo; censura

## ---
## caso 1
## ---

## tres grids, quatro taxas de falha

grids1 = time.grid.obs.t(time = tempo, event = censura, n.int = 5)
table(cut(tempo, grids1))

grids1 = grids1[-c(1,length(grids1))]

# quantas taxas eu preciso tentar
length(grids1) + 1

est = optim(par = rep(0.01,length(grids1) + 1),
            fn = loglik.PE,
            gr = NULL,
            hessian = F,
            method = "Nelder-Mead",
            time = tempo,
            delta = censura,
            cuts = grids1)

est$par

par(mfrow = c(1,2))
surv.pe = PE(time = int_tempo, cuts = grids1, levels = est2$par,
             type = 'survival')

plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = 'S(t) estimada',
     xlab = "Tempo", lwd = 2)

lines(int_tempo, surv.pe, col = '#D885A3', lwd = 2)

legend('topright', legend = 'PE(4)',
       lwd=2, bty = 'n', col = '#D885A3', cex = 0.8)

AIC.surv(loglik = -est$value, n.param = length(grids1) + 1)

BIC.surv(loglik = -est$value, n.param = length(grids1) + 1,
         n.sample = length(tempo))

HC.surv(loglik = -est$value, n.param = length(grids1) + 1,
         n.sample = length(tempo))


## testar para ate 6 particoes no tempo de falha

particoes = 1:6

for (particao in particoes){
  value.particao = time.grid.obs.t(time = tempo, event = censura, n.int = particao)
  
  grids1 = value.particao[-c(1,length(value.particao))]
  
  # quantas taxas eu preciso tentar
  length(grids1) + 1
  
  estimacao = optim(par = rep(0.01,length(grids1) + 1),
              fn = loglik.PE,
              gr = NULL,
              hessian = F,
              method = "BFGS",
              time = tempo,
              delta = censura,
              cuts = grids1)
  
  aic.res = AIC.surv(loglik = -estimacao$value, n.param = length(grids1) + 1)
  
  bic.res = BIC.surv(loglik = -estimacao$value, n.param = length(grids1) + 1,
           n.sample = length(tempo))
  
  hc.res = HC.surv(loglik = -estimacao$value, n.param = length(grids1) + 1,
          n.sample = length(tempo))
  
  cat("\n", "Quantidade de particoes: ", particao, "\n", "\n")
  cat("Valores para as particoes: ", value.particao, "\n" , "\n")
  
  cat("valor AIC: ", aic.res, "\n",
      "valor BIC: ", bic.res, "\n",
      "valor HC: ", hc.res, "\n")
}

## Laco pelo Nelder-Mead

for (particao in particoes){
  value.particao = time.grid.obs.t(time = tempo, event = censura, n.int = particao)
  
  grids1 = value.particao[-c(1,length(value.particao))]
  
  # quantas taxas eu preciso tentar
  length(grids1) + 1
  
  estimacao = optim(par = rep(0.01,length(grids1) + 1),
                    fn = loglik.PE,
                    gr = NULL,
                    hessian = F,
                    method = "Nelder-Mead",
                    time = tempo,
                    delta = censura,
                    cuts = grids1)
  
  aic.res = AIC.surv(loglik = -estimacao$value, n.param = length(grids1) + 1)
  
  bic.res = BIC.surv(loglik = -estimacao$value, n.param = length(grids1) + 1,
                     n.sample = length(tempo))
  
  hc.res = HC.surv(loglik = -estimacao$value, n.param = length(grids1) + 1,
                   n.sample = length(tempo))
  
  cat("\n", "Quantidade de particoes: ", particao, "\n", "\n")
  cat("Valores para as particoes: ", value.particao, "\n" , "\n")
  
  cat("valor AIC: ", aic.res, "\n",
      "valor BIC: ", bic.res, "\n",
      "valor HC: ", hc.res, "\n")
}

## ----
## Modelo Exponencial por Partes Potencia
## ----

## ---
## caso 1
## ---

## tres grids, quatro taxas de falha

grids1 = time.grid.obs.t(time = tempo, event = censura, n.int = 3)
table(cut(tempo, grids1))

grids1 = grids1[-c(1,length(grids1))]

# quantas taxas eu preciso tentar
length(grids1) + 1

est.ppe = optim(par = c(rep(0.01,length(grids1) + 1),1.5),
                fn = loglik.PPE2,
                gr = NULL,
                hessian = T,
                method = "BFGS",
                time = tempo,
                delta = censura,
                cuts = grids1)

#est.ppe

taxas.ppe.est1 = est.ppe$par[1:(length(grids1)+1)]
alpha.ppe.est1 = est.ppe$par[(length(grids1)+2)]

taxas.ppe.est1; alpha.ppe.est1

par(mfrow = c(1,2))

surv.ppe1 = PPE(time = int_tempo, cuts = grids1, levels = taxas.ppe.est1,
    alpha = alpha.ppe.est1, type = 'survival')

plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = 'S(t) estimada',
     xlab = "Tempo", lwd = 2)


lines(int_tempo, surv.ppe1 , col = '#FF87CA', lwd = 2)

legend('topright', legend= 'PPE(4)',
       lwd=2, bty = 'n', col = '#FF87CA', cex = 0.8)

## Medidas de qualidade de ajuste

AIC.surv(loglik = -est.ppe$value, n.param = length(grids1) + 1 + 1)
BIC.surv(loglik = -est.ppe$value, n.param = length(grids1) + 1 + 1,
          n.sample = length(tempo))
HC.surv(loglik = -est.ppe$value, n.param = length(grids1) + 1 + 1,
         n.sample = length(tempo))


## testar para ate 6 particoes no tempo de falha

for (particao in particoes){
  value.particao = time.grid.obs.t(time = tempo, event = censura, n.int = particao)
  
  grids1 = value.particao[-c(1,length(value.particao))]
  
  
  estimacao.ppe = optim(par = c(rep(0.01,length(grids1) + 1),1.5),
                  fn = loglik.PPE2,
                  gr = NULL,
                  hessian = T,
                  method = "BFGS",
                  time = tempo,
                  delta = censura,
                  cuts = grids1)

  
  aic.res = AIC.surv(loglik = -estimacao.ppe$value, n.param = length(grids1) + 1)
  
  bic.res = BIC.surv(loglik = -estimacao.ppe$value, n.param = length(grids1) + 1,
                     n.sample = length(tempo))
  
  hc.res = HC.surv(loglik = -estimacao.ppe$value, n.param = length(grids1) + 1,
                   n.sample = length(tempo))
  
  cat("\n", "Quantidade de particoes: ", particao, "\n", "\n")
  cat("Valores para as particoes: ", value.particao, "\n" , "\n")
  
  cat("valor AIC: ", aic.res, "\n",
      "valor BIC: ", bic.res, "\n",
      "valor HC: ", hc.res, "\n")
}

## verificando com o Nelder-Mead

for (particao in particoes){
  value.particao = time.grid.obs.t(time = tempo, event = censura, n.int = particao)
  
  grids1 = value.particao[-c(1,length(value.particao))]
  
  
  estimacao.ppe = optim(par = c(rep(0.01,length(grids1) + 1),1.5),
                        fn = loglik.PPE2,
                        gr = NULL,
                        hessian = F,
                        method = "Nelder-Mead",
                        time = tempo,
                        delta = censura,
                        cuts = grids1)
  
  
  aic.res = AIC.surv(loglik = -estimacao.ppe$value, n.param = length(grids1) + 1)
  
  bic.res = BIC.surv(loglik = -estimacao.ppe$value, n.param = length(grids1) + 1,
                     n.sample = length(tempo))
  
  hc.res = HC.surv(loglik = -estimacao.ppe$value, n.param = length(grids1) + 1,
                   n.sample = length(tempo))
  
  cat("\n", "Quantidade de particoes: ", particao, "\n", "\n")
  cat("Valores para as particoes: ", value.particao, "\n" , "\n")
  
  cat("valor AIC: ", aic.res, "\n",
      "valor BIC: ", bic.res, "\n",
      "valor HC: ", hc.res, "\n")
}


## Teste de Razao de Verrosimilhancas

# Modelo = c("Gama Generalizado", "Exponencial", "Log-Logistico",
#            "Log-Normal", "Weibull", "gamma")
# 
# Verossimilhanca = c(fit_gammagen$loglik, fit_exp$loglik,
#                     fit_loglog$loglik, fit_lognorm$loglik,
#                     fit_weibull$loglik, fit_gamma$loglik)
# 
# TRV = 2*(fit_gammagen$loglik-Verossimilhanca)
# 
# valor_p = round(pchisq(TRV,df=2,lower.tail=FALSE), 2)
# 
# resultado = data.frame(Modelo=Modelo, 
#                        Verossimilhanca = Verossimilhanca, 
#                        TRV=TRV, 
#                        valor_p=valor_p)
# 
# 
# int_tempo

## fitar o modelo exponencial por parte
## colocar o valor do KM dentro da fun??o de sobrevivencia com a funcao












