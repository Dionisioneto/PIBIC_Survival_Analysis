
## estudos de medidas de influencia global
## e distancia de verossimilhancas

dim(hiv2)

## estimacao completa do MEPP com 2 particoes

setwd("C:/Users/NetoDavi/Desktop/survival_pibic/dados_para_teste")

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
           2,
           1,0.5,
           0.5)


estimacao <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                   control=list(fnscale=1), hessian = TRUE, l=left,
                   r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)

estimacao$par

est.matrix.var.cov = solve(estimacao$hessian)
erro.padrao = sqrt(diag(est.matrix.var.cov))
erro.padrao


parametros2_int = estimacao$par
loglik_estimacao2_int = estimacao$value

## matriz para o armazenamento dos parametros
matriz_est_global = matrix(0, nrow = nrow(hiv2),
                           ncol = length(parametros2_int))

loglik_est_diff = rep(0, length.out = nrow(hiv2))


for(individuo in 1:nrow(hiv2)){
  
  ## Duas particoes
  data = hiv2[-individuo,]
  
  left = (data$Li)/7
  right = (data$Ri)/7
  cov = data$DoseType
  
  n.int = 2
  
  x.f <- cbind(x1=cov)
  x.c <- cbind(1, x1=cov)
  
  grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                              bmax=n.int)
  
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  
  chutes = c(1,0.1,
             2,
             1,0.5,
             0.5)
  
  
  est <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
               control=list(fnscale=1), hessian = TRUE, l=left,
               r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)
  
  matriz_est_global[individuo,] = est$par
  loglik_est_diff[individuo] = est$value
  
  print(paste("retirando a observacao: ", individuo))
}


## pontos que necessitaram de uma mudanca no chute

for(individuo in c(154,185,272,280,292,296,356)){
  
  ## Duas particoes
  data = hiv2[-individuo,]
  
  left = (data$Li)/7
  right = (data$Ri)/7
  cov = data$DoseType
  
  n.int = 2
  
  x.f <- cbind(x1=cov)
  x.c <- cbind(1, x1=cov)
  
  grid.obs=time.grid.interval(li=left, ri=right, type="OBS", 
                              bmax=n.int)
  
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  
  chutes = c(1,0.1,
             2,
             1,0.5,
             0.5)
  
  
  est <- optim(par = chutes, fn=loglikIC, gr = NULL, method = "BFGS",
                     control=list(fnscale=1), hessian = TRUE, l=left,
                     r=right, x.cure=x.c, x.risk=x.f, grid=grid.obs)
  
  matriz_est_global[individuo,] = est$par
  loglik_est_diff[individuo] = est$value
  
  print(paste("retirando a observacao: ", individuo))
}

## 154,185,272,280,292,296,356

matriz_est_global[c(154,185,272,280,292,296,356),]
loglik_est_diff


loglik_est_diff.matrix = as.matrix(loglik_est_diff)

dim(matriz_est_global)
dim(loglik_est_diff.matrix)



loglik_estimacao2_int

matrix_loglik_2_int = matrix(data = loglik_estimacao2_int,
                          nrow = nrow(hiv2))

matrix_est_2_int = matrix(data=parametros2_int,
                         nrow = nrow(hiv2),
                         ncol = length(parametros2_int),
                         byrow = T)


GDi_matrix = matrix(NA, nrow = nrow(hiv2),
             ncol = length(parametros2_int))

GDi_matrix

for(obs in 1:nrow(hiv2)){
  diff = matriz_est_global[obs,]-matrix_est_2_int[obs,]
  GDi_matrix[obs,] = abs((diff)^2*(-loglik_est_diff.matrix[1,1]))
}

plot(1:length(GDi_matrix[,1]),GDi_matrix[,1], 
     type = 'h', ylim = c(0,90),
     xlab = "Observação",
     ylab = "|Distância Generalizada de Cook|",
     main = expression(hat(lambda)[1]),
     axes = F) ## lambda1

#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,nrow(hiv2),8), labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,max(GDi_matrix[,1]+5),5),
     labels=seq(0,max(GDi_matrix[,1]+5),5), cex.axis=0.7)

pontos_influentes1 <- which(GDi_matrix[,1] > 20)
points(pontos_influentes1, GDi_matrix[,1][pontos_influentes1],
       col = "red", pch = 16)
text(pontos_influentes1, GDi_matrix[,1][pontos_influentes1], 
     labels = pontos_influentes1, pos = 3)


## lambda2

plot(1:length(GDi_matrix[,2]),GDi_matrix[,2], type = 'h',
     ylim = c(0,max(GDi_matrix[,2]+1)),
     xlab = "Observação",
     ylab = "|Distância Generalizada de Cook|",
     main = expression(hat(lambda)[2]),
     axes = F)


#Criando o novo eixo-x contendo valores incrementados de 1 em 1
axis(side=1, at=seq(0,nrow(hiv2),8),
     labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

#Criando o novo eixo-y contendo valores incrementados de 10 em 10
axis(side=2, at=seq(0,max(GDi_matrix[,2]+1),1),
     labels=seq(0,max(GDi_matrix[,2]+1),1), cex.axis=0.7)

pontos_influentes2 <- which(GDi_matrix[,2] > 10)
points(pontos_influentes2, GDi_matrix[,2][pontos_influentes2],
       col = "red", pch = 16)

text(pontos_influentes2, GDi_matrix[,2][pontos_influentes2], 
     labels = pontos_influentes2, pos = 3)

## alpha
plot(1:length(GDi_matrix[,3]),GDi_matrix[,3], type = 'h',
     ylim = c(0,max(GDi_matrix[,3]+1000)),
     xlab = "Observação",
     ylab = "|Distância Generalizada de Cook|",
     main = expression(hat(alpha)),
     axes = F) 

axis(side=1, at=seq(0,nrow(hiv2),8),
     labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

axis(side=2, at=seq(0,max(GDi_matrix[,3]),1000),
     labels=seq(0,max(GDi_matrix[,3]),1000), cex.axis=0.7)

pontos_influentes3 <- which(GDi_matrix[,3] > 2000)
points(pontos_influentes3, GDi_matrix[,3][pontos_influentes3],
       col = "red", pch = 16)

text(pontos_influentes3, GDi_matrix[,3][pontos_influentes3], 
     labels = pontos_influentes3, pos = 3)


## beta0
plot(1:length(GDi_matrix[,4]),GDi_matrix[,4], 
     type = 'h',
     ylim = c(0,max(GDi_matrix[,4]+1000)),
     xlab = "Observação",
     ylab = "|Distância Generalizada de Cook|",
     main = expression(hat(b)[0]),
     axes = F) 

axis(side=1, at=seq(0,nrow(hiv2),8),
     labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

axis(side=2, at=seq(0,max(GDi_matrix[,4]+500),500),
     labels=seq(0,max(GDi_matrix[,4]+500),500), cex.axis=0.7)

pontos_influentes4 <- which(GDi_matrix[,4] > 5000)
points(pontos_influentes4, GDi_matrix[,4][pontos_influentes4],
       col = "red", pch = 16)

text(pontos_influentes4, GDi_matrix[,4][pontos_influentes4], 
     labels = pontos_influentes4, pos = 3)

## beta1
plot(1:length(GDi_matrix[,5]),GDi_matrix[,5], type = 'h',
     ylim = c(0,max(GDi_matrix[,5]+1000)),
     xlab = "Observação",
     ylab = "|Distância Generalizada de Cook|",
     main = expression(hat(b)[0]),
     axes = F) 


axis(side=1, at=seq(0,nrow(hiv2),8),
     labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

axis(side=2, at=seq(0,max(GDi_matrix[,5]),1000),
     labels=seq(0,max(GDi_matrix[,5]),1000), cex.axis=0.7)


## nao houve distancia muitos grandes, ou um padrao a ser permitido
 
# pontos_influentes5 <- which(GDi_matrix[,5] > 30000)
# points(pontos_influentes5, GDi_matrix[,5][pontos_influentes5],
#        col = "red", pch = 16)
# 
# text(pontos_influentes5, GDi_matrix[,5][pontos_influentes5], 
#      labels = pontos_influentes5, pos = 3)

## Beta1
plot(1:length(GDi_matrix[,6]),GDi_matrix[,6], type = 'h',
     ylim = c(0,max(GDi_matrix[,6])+150),
     xlab = "Observação",
     ylab = "|Distância Generalizada de Cook|",
     main = expression(hat(beta)[1]),
     axes = F) 

axis(side=1, at=seq(0,nrow(hiv2),8),
     labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

axis(side=2, at=seq(0,max(GDi_matrix[,6]),100),
     labels=seq(0,max(GDi_matrix[,6]),100), cex.axis=0.7)


pontos_influentes6 <- which(GDi_matrix[,6] > 1000)
adj_positions <- ifelse(duplicated(GDi_matrix[,6][pontos_influentes6]),
                        3, 0.5)

points(pontos_influentes6, GDi_matrix[,6][pontos_influentes6],
       col = "red", pch = 16)

text(pontos_influentes6[-2], GDi_matrix[,6][pontos_influentes6[-2]], 
     labels = pontos_influentes6[-2], 
     pos = 3)

text(pontos_influentes6[2], GDi_matrix[,6][pontos_influentes6[2]], 
     labels = pontos_influentes6[2], 
     pos = 4)

## distancia de verossimilhancas

LD_theta = abs(2*(loglik_est_diff-(matrix_loglik_2_int)))

plot(1:length(LD_theta),LD_theta, type = 'h',
     ylim = c(0,max(LD_theta)+5),
     xlab = "Observação",
     ylab = "|Distância de log-Verossimilhanças|",
     main = "",
     axes = F)

axis(side=1, at=seq(0,nrow(hiv2),8),
     labels=seq(0,nrow(hiv2),8), cex.axis=0.7)

axis(side=2, at=seq(0,max(LD_theta)+5,1),
     labels=seq(0,max(LD_theta)+5,1), cex.axis=0.7)


pontos_influentes_log <- which(LD_theta > 50)
points(pontos_influentes_log, LD_theta[pontos_influentes_log],
       col = "red", pch = 16)

text(pontos_influentes_log, LD_theta[pontos_influentes_log], 
     labels = pontos_influentes_log, pos = 3)



