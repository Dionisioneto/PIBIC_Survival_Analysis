## ---------------------------------------------------------------
## ajuste aos dados dados HIV_SUM_set2 
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

AIC.surv(loglik = hiv2.weibull$llk, 
         n.param = length(hiv2.weibull$coefficients),
         n.sample = dim(hiv2)[1])

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

x.f <- cbind(cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,0.1,1,2,2,2,5,2,
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
         n.param = length(max.mep$par),
         n.sample = dim(hiv2)[1])

BIC.surv(loglik = max.mep$value*-1,
         n.param = length(max.mep$par),
         n.sample = dim(hiv2)[1])

HC.surv(loglik = max.mep$value*-1,
        n.param = length(max.mep$par),
        n.sample = dim(hiv2)[1])


SpopMEP()





# ## ---
# ## modelo exponencial por partes potencia sem a fracao de cura
# ## ---
# 
# source('C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste/proposta_MEPP_int_sem_cura.R')
# 
# hiv2 = read.table('HIV_SUM_set2.txt', header = T)
# hiv2 <- hiv2[sample(row.names(hiv2)), ]
# 
# left = (hiv2$Li)/7
# right = (hiv2$Ri)/7
# cov = cbind(x1=hiv2$DoseType)
# 
# n.int = 2
# 
# grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
#                             bmax=n.int)
# 
# grid.obs=grid.obs[-c(1, length(grid.obs))]
# 
# chutes = c(1,0.1,
#            5,
#            0.5)
# 
# max.mepp.int = optim(par = chutes, fn=loglikIC.mepp.int,
#                      gr = NULL, method = "BFGS",
#                      control=list(fnscale=1), hessian = TRUE,
#                      l = left, r = right,
#                      x.cov=cov, 
#                      grid=grid.obs)
# 
# max.mepp.int$par
# max.mepp.int$par



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

n.int = 8

x.f <- cbind(x1=cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(1,0.2,1,10,1,1,1,5,
           9,
           1, 0.5,
           0.5)

estimacao <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                   control=list(fnscale=1), hessian = TRUE, l=left,
                   r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)

estimacao$par

AIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
         n.sample = nrow(hiv2))
BIC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
         n.sample = nrow(hiv2))
HC.surv(loglik = estimacao$value*-1, n.param = length(estimacao$par),
        n.sample = nrow(hiv2))














