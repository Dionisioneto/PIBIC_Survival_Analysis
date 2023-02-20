## ------
## Projeto: Modelos Semi-Parametricos para dados
## de sobrevivencia com censura intervalar
##
## Aplicacao do exponencial por partes sob censura intervalar com covariaveis
## ------

## ---
## Pacotes 
## ---

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(eha, survival)

## ajustar o script com as funcoes de sobrevivencia 
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')


## ---
## Presenca de covariaveis
## ---

set.seed(10)

tamanho.amostral = 500
taxas.falha = c(0.2, 0.4, 0.95)
particoes = c(0.5, 0.9)
potencia = 1.3


beta = c(0.5, 2.3)
x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
x.matriz = as.matrix(cbind(x1, x2))

tempo.falha = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
                           cuts = particoes, x.mat = x.matriz, beta.par = beta)

tempo.censura = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
                             cuts = particoes, x.mat = x.matriz, beta.par = beta)

tempo = pmin(tempo.falha, tempo.censura)
delta = ifelse(tempo.falha <= tempo.censura, 1, 0)

#table(cut(tempo, c(0,particoes, Inf)))


#chutes = rep(0.5, 6)

grid = time.grid.obs.t(tempo, delta, n.int = 3)
grid = grid[-c(1, length(grid))]
chutes = c(rep(0.5,length(grid)+1),1,0.5,1)

## Metodo numerico BFGS
estimacao.teste.cox = optim(par = chutes,
                            fn = loglik.cox,
                            gr = NULL,
                            hessian = TRUE,
                            method = "BFGS",
                            tempos = tempo,
                            censura = delta,
                            intervalos = grid,
                            covariaveis = x.matriz)
taxas.falha
estimacao.teste.cox$par[1:3]

potencia
estimacao.teste.cox$par[4]

beta
estimacao.teste.cox$par[5:6]

## ---
## Aplicacoes aos dados de cancer de pulmao
## ---

data(lung)

head(lung)

tempo = lung$time
censura = ifelse(lung$status == 1, 0, 1)

## variaveis a serem utilizadas

# 1. sex
# 2. age
# 3. meal.cal: quantidade de calorias consumidas nas refeicoes;
# 4. wt.loss: quantidade de peso perdido no  ́ultimos sei meses.

lung$sex
lung$age
#lung$meal.cal
#lung$wt.loss

## padronizacao da idade
age.z = (lung$age - mean(lung$age))/sd(lung$age)


## categorizacao bernoulli do sexo
sex.b = ifelse(lung$sex == 1, 0, 1)

## categorizacao da qualidade em cuidado proprio (acima/abaixo da mediana)
#lung$ph.karno[is.na(lung$ph.karno)] = median(lung$ph.karno, na.rm = T)

#phkarno.z =  ifelse(lung$ph.karno >= median(lung$ph.karno, na.rm = T), 1, 0)

x.matriz = as.matrix(cbind(age.z, sex.b))

#table(cut(tempo, c(0,particoes, Inf)))

#chutes = rep(0.5, 6)

grid = time.grid.obs.t(tempo, censura, n.int = 3)
grid = grid[-c(1, length(grid))]
chutes = c(rep(0.01,length(grid)+1),1,0.5,1)


## Metodo numerico BFGS
est.cov.cox = optim(par = chutes,
                          fn = loglik.cox,
                          gr = NULL,
                          hessian = T,
                          method = "BFGS",
                          tempos = tempo,
                          censura = censura,
                          intervalos = grid,
                          covariaveis = x.matriz)

est.cov.cox$par


## ---
## Censura intervalar com covariaveis
## ---

tamanho.amostral = 200
taxas.falha = c(0.2, 0.4, 0.95)
particoes = c(1.3, 2.3)
potencia = 1.3

set.seed(10)
beta = c(0.5, 2.3)
x1 = rnorm(n = tamanho.amostral, mean = 0, sd = 1)
x2 = rbinom(n = tamanho.amostral, size = 1, prob = 0.5)
x.matriz = as.matrix(cbind(x1, x2))

taxa.censura = 0.8

tempo.int.cox = sim.ICdata(n = tamanho.amostral,
           lambda.param = taxas.falha,
           alpha.param = potencia, 
           grid.vector = particoes,
           beta.param = beta,
           x.matrix = x.matriz,
           lambda.cens.param = taxa.censura)

head(tempo.int.cox)

## ----
## Estimacao de maxima verossmilhanca dos parametros
## ----

chutes = c(rep(0.1,3),
           1,
           0.25,2)

est.int.cox = optim(par = chutes,
                    fn = loglik.int,
                    gr = NULL,
                    hessian = F,
                    method = "Nelder-Mead",
                    time.r = tempo.int.cox$R,
                    time.l = tempo.int.cox$L,
                    delta = tempo.int.cox$delta,
                    cuts = particoes,
                    x.matrix = x.matriz)



est.int.cox.2 = optim(par = chutes,
                    fn = loglik.int,
                    gr = NULL,
                    hessian = T,
                    method = "BFGS",
                    time.r = tempo.int.cox$R,
                    time.l = tempo.int.cox$L,
                    delta = tempo.int.cox$delta,
                    cuts = particoes,
                    x.matrix = x.matriz)







