
## Estudo da Deviance Estatistica


# D =  Dinc + Dlat = 2*sum(dinc) + 2*sum(dlat)

## ssat: estimação do modelo saturado
## s: estimação com o modelo escolhido

## ---
## modelo exponencial por partes potencia com fracao de cura
## ---
source('C:/Users/NetoDavi/Desktop/survival_pibic/source/supervisor_functions.R')

hiv2 = read.table('HIV_SUM_set2.txt', header = T)


## Duas particoes

left = (hiv2$Li)/7
right = (hiv2$Ri)/7
cov = hiv2$DoseType

n.int = 2

x.f <- cbind(x1=cov)
x.c <- cbind(1, x1=cov)

grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                            bmax=n.int)

grid.obs=grid.obs[-c(1, length(grid.obs))]

chutes = c(0.1,1,
           1,
           1,0.5,
           0.5)


estimacao <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                   control=list(fnscale=1), hessian = TRUE, l=left,
                   r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)

parametros2_int = estimacao$par

## estimação de 5 particoes

grid.observavel5=time.grid.interval(li=hiv2$Li/7, ri=hiv2$Ri/7, type="OBS", 
                                    bmax=5)

grid.observavel5=grid.observavel5[-c(1, length(grid.observavel5))]

## grupo Dosetype == 1
valores.tempo1 = seq(0,40,length.out = nrow(hiv2[ind.dose1,])) 

parametros5_MEP_int

surv5.MEP.estrat1 = SpopMEP(t=valores.tempo1, 
                            lambda.par=parametros5_MEP_int[1:5],
                            grid.vet=grid.observavel5, 
                            theta.par=parametros5_MEP_int[6:7], 
                            beta.par=parametros5_MEP_int[8],
                            x.cure=x.c[ind.dose1,], x.risk=x.f[ind.dose1]) 

## verossimilhanca

Dev.loglik.IC.MEPP = function(par, 
                              l, r, 
                              x.cure, x.risk, 
                              grid){
  b <- length(grid)+1
  hazards = par[1:b] ## taxas de falha para os b intervalos
  alpha = par[b + 1] ## pparrametro de potencia
  
  n.cov.cure = dim(x.cure)[2] 
  n.cov.risk = dim(x.risk)[2] 
  
  betas.cure = par[(b + 2):(b + 1 + n.cov.cure)]
  betas.risk = par[(b + 2 + n.cov.cure):(b + 1 + n.cov.cure + n.cov.risk)]
  n.sample <- nrow(x.cure)
  
  delta <- ifelse(is.finite(r), 1, 0)
  
  lik <- rep(0, n.sample)
  
  spop.l <- SpopMEPP(t=l[delta==1], lambda.par=hazards, alpha.par=alpha,
                     grid.vet=grid, beta.par=betas.risk, theta.par=betas.cure,
                     x.cure=x.cure[delta==1,], x.risk=x.risk[delta==1,])
  
  spop.r <- SpopMEPP(t=r[delta==1], lambda.par=hazards, alpha.par=alpha,
                     grid.vet=grid, beta.par=betas.risk, theta.par=betas.cure,
                     x.cure=x.cure[delta==1,], x.risk=x.risk[delta==1,])
  
  lik[delta==1] <- spop.l-spop.r
  
  spop.l2 <- SpopMEPP(t=l[delta==0], lambda.par=hazards, alpha.par=alpha,
                      grid.vet=grid, beta.par=betas.risk, theta.par=betas.cure,
                      x.cure=x.cure[delta==0,], x.risk=x.risk[delta==0,])
  lik[delta==0] <- spop.l2
  
  return(log(lik))
}

loglik2 = Dev.loglik.IC.MEPP(par = parametros2_int,
                   l = hiv2$Li/7,
                   r = hiv2$Ri/7,
                   x.cure = x.c,
                   x.risk = x.f,
                   grid = grid.observavel2)
  

loglik5 = Dev.loglik.IC.MEPP(par = parametros5_int,
                   l = hiv2$Li/7,
                   r = hiv2$Ri/7,
                   x.cure = x.c,
                   x.risk = x.f,
                   grid = grid.observavel5)

dev = loglik5 - loglik2

D.O = sum(dev)*2


dev0 = dev[hiv2$Censind==0]
dev1 = dev[hiv2$Censind==1]



summary(dev)  

plot(1:length(dev1),dev1, xlim = c(1,270), ylim = c(-2.5,3.5), col = 'blue')  
points(1:length(dev0),dev0, col = 'red')
abline(h=c(-2,2), lty=c(2,2), lwd=c(2, 2))

plot(x = rep(parametros2_int[6],length(dev)), y = dev )  
plot(x = rep(parametros2_int[5],length(dev)), y = dev )  
plot(x = rep(parametros2_int[4],length(dev)), y = dev )  

 



