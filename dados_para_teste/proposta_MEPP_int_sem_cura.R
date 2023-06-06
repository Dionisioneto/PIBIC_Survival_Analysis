
library(eha)

## Funcao que particiona os grids para censura intervalar

time.grid.interval <- function(li=li, ri=ri, type=type, bmax=bmax)
{  
  ## Funcao que retorna os intervalos da partiÃƒÂ§ÃƒÂ£o mais fina
  ## baseada nos limites observados, distintos e finitos.
  ## Argumentos:
  ## li: limite inferior dos intervalos observados.
  ## ri: limite superior dos intervalos observados.
  ## bmax: numero mÃƒÂ¡ximo de intervalos.
  
  # li = dados$L
  # ri= dados$R
  
  #--- Inicio da funcao:
  
  #-- Construir uma grade tipo 1: baseando-se em tempos observaveis
  if(type=="OBS")
  {
    #grid.vet <- sort(unique(c(0, li, ri, Inf)))
    grid.vet <- sort(unique(c(0, li[is.finite(li)], ri[is.finite(ri)], Inf)))
    grid.size.vet <- length(grid.vet) # Grid time size
    
    if( isTRUE(bmax<grid.size.vet)==TRUE )
    {
      k        <- round((length(grid.vet)-1)/bmax,0)
      id.grid  <- round(seq(k,(length(grid.vet)-1), length.out=bmax),0)
      grid.vet <- c(0,grid.vet[-1][id.grid])
      return(grid.vet)
    }else{
      grid.vet <- sort(unique(c(0, li, ri, Inf)))
      return(grid.vet)
    }
  } #-- Construir uma grade tipo 2: espacos equiprovaveis
  if(type=="EQUI")
  {
    grade.vet <- seq(0, max(ri[ri!=Inf]), length.out=bmax)
    grid.vet <- c(grade.vet,Inf)
    return(grid.vet)
  }  
}

## Funcoes do exponencial por partes potencia


## ---
## Modelo Exponencial por partes por potencia (PPE)
## ---

PPE = function(time, cuts, levels, alpha, type = "survival"){
  
  ## ---
  ## Seccao do Exponencial por partes
  ## ---
  
  ## sobrevivencia
  F0.t = PE(time = time, cuts = cuts, levels = levels, type = 'distribution')
  
  ## densidade
  f0.t = PE(time = time, cuts = cuts, levels = levels, type = 'density')
  
  ## ---
  ## Atribuindo aos valores do Exponencial por partes potencia (PPE) 
  ## ---
  
  ## funcao distribuicao do PPE
  F1.t = (F0.t)^alpha
  
  ## funcao densidade do PPE
  f1.t = (alpha*((F0.t)^(1 - alpha))) * f0.t
  
  ## funcao de sobrevivencia do PPE
  s1.t = 1 - F1.t
  
  ## funcao taxa de falha do PPE
  h1.t = f1.t/s1.t
  
  ## taxa de falha acumulada do PPE
  H1.t = -log(s1.t)
  
  if(type == "survival"){
    return(s1.t)
  } else {
    if(type == "density"){
      return(f1.t)
    } else {
      if(type == "hazard"){
        return(h1.t)
      } else{
        if(type == "cum_hazard"){
          return(H1.t)
        } else {
          if(type == "distribution"){
            return(F1.t) 
          } else{
            return("Please choose one type of function!!!")
          }
        }
      }
    }
  }                                
}


## ------
## Funcao de verosimilhanca para censura intervalar, com covariaveis
## ------

loglik.int.mepp = function(par, time.r, time.l,
                      grid, x.matrix){
  
  x.matrix = as.matrix(x.matrix)
  
  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  exp = par[b + 1] ## parametro de potencia
  
  n.covars = dim(x.matrix)[2] ## numero de covariaveis
  betas = par[(b + 2):(b + 1 + n.covars)]
  betas = as.matrix(betas)
  
  delta = ifelse(time.r == Inf, 0,1)
  
  like = rep(0, dim(x.matrix)[1])  ## armazenamento de informacao da verossimilhanca
  
  
  ## informacoes exponencial por partes Potencia (PPE) para esquerda
  s0.tl = PPE(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## informacoes exponencial por partes Potencia (PPE) para direita
  s0.tr = PPE(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  
  ## contribuicao do evento
  
  sl = s0.tl^(exp(x.matrix[delta==1,] %*% betas))
  
  sr = s0.tr^(exp(x.matrix[delta==1,] %*% betas))
  
  like[delta==1] = (sl - sr)
  
  ## contribuicao da censura
  s0.tl = PPE(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = exp, type = "survival")
  sl = s0.tl^(exp(x.matrix[delta==0,] %*% betas))
  
  like[delta==0] = sl
  
  log.vero = sum(log(like))
  
  return(log.vero)
}


fit.mepp = function(L, R, n.int, cov, start){
  
  ## extracao do grid observado
  
  grid.obs=time.grid.interval(li=L, ri=R, type="OBS", bmax= n.int)
  grid.obs=grid.obs[-c(1, length(grid.obs))]
  
  mepp.est <- optim(par = start, fn=loglik.int.mepp, gr = NULL, method = "BFGS",
               control=list(fnscale=-1), hessian = TRUE, time.l=L, 
               time.r=R, x.matrix = cov, grid=grid.obs)
  
  estimated = mepp.est$par
  hessian = mepp.est$hessian
  loglik = mepp.est$value
  
  results = list(estimated = estimated, hessian = hessian, loglik = loglik)
  
  return(results)
}


## lambda = 2
## alpha = 1
## betas = 4

chute = c(0.1,0.2,0.1,
          1.2,
          0.1,0.1,0.1,0.1)

ajuste.mepp.nr = fit.mepp(L = smoke2009$Timept1, R = smoke2009$Timept2,
                            n.int = 3, cov = covariaveis, start = chute)

ajuste.mepp.nr$estimated
ajuste.mepp.nr$loglik

AIC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated))

BIC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated),
         n.sample = dim(smoke2009)[1])

HC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated),
        n.sample = dim(smoke2009)[1])

## ---
## dados aneurysm
## ---

chute = c(0.1,0.1,0.1,
          1,
          0.1,0.1,0.1)

gr_pad = (aneurysm$gr - mean(aneurysm$gr))/sd(aneurysm$gr) ## precisa normalizar para estimar

covariaveis.aneurysm = cbind(aneurysm$mo, gr_pad, aneurysm$lok)

ajuste.mepp.nr = fit.mepp(L = aneurysm$t.left, R = aneurysm$t.right,
                          n.int = 3, cov = covariaveis.aneurysm, start = chute)

ajuste.mepp.nr$estimated
ajuste.mepp.nr$loglik

AIC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated))

BIC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated),
         n.sample = dim(aneurysm)[1])

HC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated),
        n.sample = dim(aneurysm)[1])

## ---
## dados breast
## ---

n.intervalos = 15

for (int in 2:15){

  n.intervalos = int
  
  chute = c(rep(0.1,int),
            1,
            0.1)
  
  
  ajuste.mepp.nr = fit.mepp(L = breast$left, R = breast$right,
                            n.int = int, cov = cbind(breast$ther),
                            start = chute)
  
  aic.sf = AIC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated))
  
  bic.sf = BIC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated),
           n.sample = dim(breast)[1])
  
  hc.sf = HC.surv(ajuste.mepp.nr$loglik, n.param = length(ajuste.mepp.nr$estimated),
          n.sample = dim(breast)[1])
  
  print(paste("n intervalo", int))
  print(paste("AIC: ", aic.sf))
  print(paste("BIC: ", bic.sf))
  print(paste("HC: ", hc.sf))

}












