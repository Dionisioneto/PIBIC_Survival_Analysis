## ------
## Aplicacao de modelos probabilisticos em analise de sobrevivencia
## Apresentacao 1
## ------


if(!require(pacman)) install.packages('pacman'); library(pacman)

p_load(flexsurv,survival, eha, ICglm, latex2exp)



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

legend('topright', "Ajuste Log-Logístico",col = "magenta", lty = 1, bty = 'n')
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

# AIC
AIC(fit_exp) # Exponencial
AIC(fit_weibull) # Weibull
AIC(fit_lognorm) # Log-Normal
AIC(fit_loglog) # Log-Logistico
AIC(fit_gamma) # Gama
AIC(fit_gammagen) # Gama generalizado




# BIC

BIC = function(log_verossimilhanca, n, k){
  criterio_bayes = -2*log_verossimilhanca + (log(n) * k)
  return(criterio_bayes)
}

# Exponencial
BIC(log_verossimilhanca = fit_exp$loglik[1],
    k = 1,
    n = summary(fit_exp)$n)

# Weibull
BIC(log_verossimilhanca = fit_weibull$loglik[1],
    k = 2,
    n = summary(fit_weibull)$n)

# Log-Normal
BIC(log_verossimilhanca = fit_lognorm$loglik[1],
    k = 2,
    n = summary(fit_lognorm)$n)

# Log-Logistico
BIC(log_verossimilhanca = fit_loglog$loglik[1],
    k = 2,
    n = summary(fit_loglog)$n)

# Gama
BIC(log_verossimilhanca = fit_gamma$loglik,
    k = 2,
    n = fit_gamma$N)

# Gama generalizado
BIC(log_verossimilhanca =fit_gammagen$loglik,
    k = 3,
    n = fit_gammagen$N)

# HC

HC = function(log_vero, n, k){
  criterio = -2*(log_vero) + 2*k*log(log(n))
  return(criterio)
  }

# Exponencial
HC(log_vero = fit_exp$loglik[1],
   k = length(fit_exp$coefficients),
   n = summary(fit_exp)$n)


# Weibull
HC(log_vero = fit_weibull$loglik[1],
   k = 2,
   n = summary(fit_weibull)$n)

# Log-Normal
HC(log_vero = fit_lognorm$loglik[1],
   k = 2,
   n = summary(fit_lognorm)$n)

# Log-Logistico
HC(log_vero = fit_loglog$loglik[1],
   k = 2,
   n = summary(fit_loglog)$n)

# Gama
HC(log_vero = fit_gamma$loglik,
   k = fit_gamma$npars,
   n = fit_gamma$N)

# Gama generalizado
HC(log_vero = fit_gammagen$loglik,
   k = fit_gammagen$npars,
   n = fit_gammagen$N)

source('C:/Users/NetoDavi/Downloads/funcoes_sobrevivencia_pibic2023 (1).R')


## modelo exponencial por partes

## visualizar o o estimador KM para escolher os intervalos

plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = expression(hat(S(t))),
     xlab = "Tempo", lwd = 2)

#tempo; censura

## ---
## caso 1
## ---

## tres grids, quatro taxas de falha

grids1 = c(100, 200, 305.2)

table(cut(tempo, c(0,grids1, Inf)))

# quantas taxas eu preciso tentar
length(grids1) + 1

# est = optim(par = rep(0.005,4),
#                   fn = loglik.PE2,
#                   gr = NULL,
#                   hessian = T,
#                   method = "BFGS",
#                   time = tempo,
#                   delta = censura,
#                   cuts = grids)

est2 = optim(par = rep(0.001,4),
            fn = loglik.PE,
            gr = NULL,
            hessian = F,
            method = "Nelder-Mead",
            time = tempo,
            delta = censura,
            cuts = grids1)

est2$par

par(mfrow = c(1,2))
surv.pe = PE(time = int_tempo, cuts = grids1, levels = est2$par,
             type = 'survival')

plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = 'S(t) estimada',
     xlab = "Tempo", lwd = 2)

lines(int_tempo, surv.pe, col = '#D885A3', lwd = 2)

legend('topright', legend = 'PE(4)',
       lwd=2, bty = 'n', col = '#D885A3', cex = 0.8)


## AIC

# (2*2) - (2*est2$value)
# 
# function (model) 
# {
#   LL <- logLik(object = model)
#   df <- attr(LL, "df")
#   c(-2 * LL + 2 * df)
# }
# 
# (-2 * -est2$value) + (2 * 3)
# 
# BIC(-est2$value, n = 228, k = 3)
# HC(-est2$value, n = 228, k = 3)

## AIC
(-2 * -est2$value) + (2 * 4)

BIC(-est2$value, n = 228, k = 4)
HC(-est2$value, n = 228, k = 4)

## ---
## caso 2
## ---



## dois grids, tres taxas de falha
grids2 = c(100, 400)

table(cut(tempo, c(0,grids2, Inf)))

# quantas taxas eu preciso tentar
length(grids2) + 1

est3 = optim(par = rep(0.001,length(grids2) + 1),
             fn = loglik.PE,
             gr = NULL,
             hessian = F,
             method = "Nelder-Mead",
             time = tempo,
             delta = censura,
             cuts = grids2)

est3$par


surv.pe3 = PE(time = int_tempo, cuts = grids2, levels = est3$par,
             type = 'survival')


plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = 'S(t) estimada',
     xlab = "Tempo", lwd = 2)

lines(int_tempo, surv.pe3, col = '#8CC0DE', lwd = 2)

legend('topright', legend= 'PE(3)',
       lwd=2, bty = 'n', col = '#8CC0DE', cex = 0.8)

(-2 * -est3$value) + (2 * 3)

BIC(-est3$value, n = 228, k = 3)
HC(-est3$value, n = 228, k = 3)



## ----
## modelo exponencial por partes potencia
## ----

## ---
## caso 1
## ---

## tres grids, quatro taxas de falha

# grids1 = c(100, 200, 305.2)

# table(cut(tempo, c(0,grids, Inf)))

est.ppe = optim(par = c(rep(0.01,length(grids1) + 1),5),
            fn = loglik.PPE2,
            gr = NULL,
            hessian = F,
            method = "BFGS",
            time = tempo,
            delta = censura,
            cuts = grids1)

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

(-2 * -est.ppe$value) + (2 * 4)

BIC(-est.ppe$value, n = 228, k = 4)
HC(-est.ppe$value, n = 228, k = 4)

## ---
## caso 2
## ---

## dois grids, tres taxas de falha
#grids2 = c(100, 400)

#table(cut(tempo, c(0,grids2, Inf)))

# quantas taxas eu preciso tentar
#length(grids2) + 1

est.ppe2 = optim(par = c(rep(0.01,length(grids2) + 1),5),
                fn = loglik.PPE2,
                gr = NULL,
                hessian = F,
                method = "BFGS",
                time = tempo,
                delta = censura,
                cuts = grids2)

taxas.ppe.est2 = est.ppe2$par[1:(length(grids2)+1)]
alpha.ppe.est2 = est.ppe2$par[(length(grids2)+2)]

taxas.ppe.est2; alpha.ppe.est2

surv.ppe2 = PPE(time = int_tempo, cuts = grids2, levels = taxas.ppe.est2,
                alpha = alpha.ppe.est2, type = 'survival')

plot(int_tempo, sobrevivencia_km, type = 's',
     ylab = 'S(t) estimada',
     xlab = "Tempo", lwd = 2)


lines(int_tempo, surv.ppe2 , col = 'steelblue', lwd = 2)

legend('topright', legend= 'PPE(3)',
       lwd=2, bty = 'n', col = 'steelblue', cex = 0.8)

(-2 * -est.ppe2$value) + (2 * 3) ## AIC

BIC(-est.ppe2$value, n = 228, k = 3)
HC(-est.ppe2$value, n = 228, k = 3)




## Teste de Razao de Verrosimilhancas

Modelo = c("Gama Generalizado", "Exponencial", "Log-Logistico",
           "Log-Normal", "Weibull", "gamma")

Verossimilhanca = c(fit_gammagen$loglik, fit_exp$loglik,
                    fit_loglog$loglik, fit_lognorm$loglik,
                    fit_weibull$loglik, fit_gamma$loglik)

TRV = 2*(fit_gammagen$loglik-Verossimilhanca)

valor_p = round(pchisq(TRV,df=2,lower.tail=FALSE), 2)

resultado = data.frame(Modelo=Modelo, 
                       Verossimilhanca = Verossimilhanca, 
                       TRV=TRV, 
                       valor_p=valor_p)


int_tempo

## fitar o modelo exponencial por parte
## colocar o valor do KM dentro da função de sobrevivencia com a funcao












