## ------
## Aplicacao da funcao de verossimilhanca de modelos 
## ------

## ---
## Modelo Exponencial por partes por potencia (PPE)
## ---

## pacote para a insercao dos comandos latex
install.packages('latex2exp')
 
library(latex2exp)
require(eha)


source('C:/Users/NetoDavi/Downloads/funcoes_sobrevivencia_pibic2023 (1).R')

## Geracao de dados

set.seed(100)

n = 100
intervalos = c(0.1, 0.15, 0.8)
taxas = c(0.3, 0.5, 0.6, 0.9)
alpha = 1.8

intervalos; taxas; alpha

time.fail = rppe(n = n, cuts = intervalos, levels = taxas, alpha = alpha)
time.cens = rppe(n = n, cuts = intervalos, levels = taxas, alpha = alpha)

#time.fail = rpch(n, cuts = intervalos, levels = taxas)
#time.cens = rpch(n, cuts = intervalos, levels = taxas)

tempos = pmin(time.fail,time.cens)
cens = ifelse(time.fail < time.cens, 1, 0)
mean(cens) ## proporcao de censura


hist(time.fail)
hist(time.cens)
hist(tempos)

table(cut(tempos, c(0,intervalos, Inf)))

estimacao = optim(par = c(rep(0.1,length(taxas)),2),
                  fn = loglik.PPE,
                  gr = NULL,
                  hessian = TRUE,
                  method = "BFGS",
                  time = tempos,
                  delta = cens,
                  cuts = intervalos)

estimacao$par[1:length(taxas)]
taxas

estimacao$par[length(taxas) + 1]
alpha

estimacao.teste2 = optim(par = c(rep(0.1,length(taxas)),2),
                  fn = loglik.PPE2,
                  gr = NULL,
                  hessian = TRUE,
                  method = "BFGS",
                  time = tempos,
                  delta = cens,
                  cuts = intervalos)

estimacao.teste2$par[1:length(taxas)]
taxas

estimacao.teste2$par[length(taxas) + 1]
alpha


## --- 
## Algumas visualizacoes
## ---
plot(t,estimacao$survival, type = 's')


## Replicacao dos graficos da Yollanda

# Para funcao densidade

## taxas de falha
taxas = list(c(0.2, 0.7),
              c(0.4,0.1),
              c(0.4,0.6),
              c(0.1,0.4))

alphas = c(0.8, 0.8, 1.0, 1.5)

# geracao de tempos de sobrevivencia por uma uniforme
tempo.surv = round(runif(1000, min = 0, max = 5),2)
tempo.surv = sort(tempo.surv)


estimacao1 = PPE(tempo = tempo.surv, 
                 intervalos = median(tempo.surv), 
                 taxas = taxas[[1]],
                 alpha = alphas[1])

for(i in 2011:2015){
  list[[paste0("A_",i)]] <- df[, year := as.numeric(i)]
} 

## For loop para todos os casos

for (exemplo in 1:4){
  variavel =  paste0("estimacao",exemplo)
  assign(variavel, PPE(tempo = tempo.surv, 
                      intervalos = median(tempo.surv), 
                      taxas = taxas[[exemplo]],
                      alpha = alphas[exemplo]))
  
}



estimacao1$density

## Funcao de densidade
plot(tempo.surv, estimacao1$density, type = 'l', lwd = 3,
     ylab = "Função Densidade", xlab ="Tempo")
lines(tempo.surv, estimacao2$density, col = 'red', lwd = 3)
lines(tempo.surv, estimacao3$density, col = 'green', lwd = 3)
lines(tempo.surv, estimacao4$density, col = 'blue', lwd = 3)

lambda1 = c(0.2, 0.4, 0.4, 0.1)
lambda2 = c(0.7, 0.1, 0.6, 0.4)
alpha =   c(0.8, 0.8, 1.0, 1.5)


legend('topright', 
       legend=TeX(sprintf(r'($\lambda_1 = %f, lambda_2 = %f, alpha = %f$)', 
                          lambda1, lambda2, alpha)), 
       lwd=c(3,3,3,3), 
       col=c('black', 'red', 'green', 'blue'), bty = 'n')

legend("topright", legend = expression(\lambda_1 = 0,2, \lambda_1 = 0,7, \alpha = 0,8 ),
       bty = 'n')

## Funcao de hazard
estimacao1$hazard

plot(tempo.surv, estimacao1$hazard, type = 'l', lwd = 3,
     ylab = "Função Taxa de Falha", xlab ="Tempo", ylim = c(0,1))
lines(tempo.surv, estimacao2$hazard, col = 'red', lwd = 3)
lines(tempo.surv, estimacao3$hazard, col = 'green', lwd = 3)
lines(tempo.surv, estimacao4$hazard, col = 'blue', lwd = 3)


## Funcao de Log-Verossimilhanca para censura a direita
## utilizando os conceitos ja establecdos


loglik = function(par,tempos,
                  censura,
                  intervalos
                  ){
  
  b = length(intervalos) + 1 ## numero de intervalos
  hazards = par[1:b]             ## taxas de falha dos b intervalos
  expoente = par[b +1] ## alpha: ultimo elemento dos parametrtos a ser estimado
  
  ## informacoes exponencial por partes (PE)
  F0.t =  ppch(q = tempos, cuts = intervalos, levels = hazards) ## distribuicao
   
  s0.t = 1 - F0.t ## sobrevivencia
  
  f0.t = (dpch(x = tempos, cuts = intervalos, levels = hazards)) ## densidade
   
  ## Distribuicao acumulada
  F.t = (1 - s0.t)^expoente
  
  densidade = expoente*(F0.t^(expoente - 1)) * f0.t
  sobrevivencia = 1 - F.t
  
  log_vero = sum(censura*log(densidade) + 
                   (1 - censura)*log(sobrevivencia))
  
  #return(-1*log_verossimilhanca)
  return(-1*log_vero)
}










