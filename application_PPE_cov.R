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
p_load(eha)

## ajustar o script com as funcoes de sobrevivencia 
source('C:/Users/NetoDavi/Desktop/survival_pibic/funcoes_sobrevivencia_pibic2023.R')


## ---
## Presenca de covariaveis
## ---

tamanho.amostral = 250
taxas.falha = c(0.2, 0.4, 0.95)
particoes = c(0.3, 1.9)
potencia = 1.3

set.seed(10)
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

table(cut(tempo, c(0,particoes, Inf)))

grids.estimacao = time.grid.obs.t(time = tempo, event = delta, n.int = 2)
grids.estimacao = grids.estimacao[-c(1, length(grids.estimacao))]

chutes = c(rep(1,3),1,0.5,1)
# chutes = rep(0.5, 6)

## Metodo numerico BFGS
estimacao.teste.cox = optim(par = chutes,
                          fn = loglik.cox,
                          gr = NULL,
                          hessian = TRUE,
                          method = "BFGS",
                          tempos = tempo,
                          censura = delta,
                          intervalos = grids.estimacao,
                         covariaveis = x.matriz)

taxas.falha
estimacao.teste.cox$par[1:3]

potencia
estimacao.teste.cox$par[4]

beta
estimacao.teste.cox$par[5:6]



## Metodo numerico de Nelder-Mead
estimacao.teste.cox2 = optim(par = chutes,
                            fn = loglik.cox,
                            gr = NULL,
                            hessian = FALSE,
                            method = "Nelder-Mead",
                            tempos = tempo,
                            censura = delta,
                            intervalos = grids.estimacao,
                            covariaveis = x.matriz)

taxas.falha
estimacao.teste.cox2$par[1:3] 

potencia


beta
estimacao.teste.cox2$par[5:6] 


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







