
library(eha)
## ---
## Calculo para a funcao de taxa de falha acumulada do MEPP
## ---

PPE = function(time, cuts, levels, alpha, type = "survival"){
  ## Seccao do Exponencial por partes
  
  ## sobrevivencia
  F0.t = PE(time = time, cuts = cuts, levels = levels, type = 'distribution')
  
  ## densidade
  f0.t = PE(time = time, cuts = cuts, levels = levels, type = 'density')
  
  ## Atribuindo aos valores do Exponencial por partes potencia (PPE) 
  
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

cal_ht_MEPP <- function( time.obs, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet){
  
  
  dens_MEPP <- alpha.par*(ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)^(alpha.par-1))*dpch(x =time.obs, cuts = grid.vet, levels = lambda.par)
  Cumu_MEPP <- ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)^alpha.par
  ht <- dens_MEPP/(1-Cumu_MEPP)
  
  return(ht)
}

cal_Ht_MEPP <- function(time.obs, lambda.par=lambda.par,
                        alpha.par=alpha.par, grid.vet=grid.vet){
  Cumu_MEPP <- ppch(q=time.obs, cuts = grid.vet, levels = lambda.par)^alpha.par
  Ht  = -log(1-Cumu_MEPP)
  return(Ht)
}

#time, cuts, levels, alpha,

time = runif(n = 10, min = 1.2, max = 2.4)
cuts = c(0.9, 1.7)
taxas = c(2, 0.3, 0.6)
alpha = 0.8

exp(-cal_Ht_MEPP(time.obs = sort(time), lambda.par = taxas, 
           alpha.par = alpha, grid.vet = cuts))

PPE(time=sort(time), cuts=cuts, levels=taxas,
    alpha=alpha, type = "survival")

####### FUNCAO PARA CHMAR A TAXA DE FALHA E OUTRAS FUNCOES OK!!!!



time <- function( t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=u.unif ){
  exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.mat%*%beta.par)) - u.unif
}

gen.mepp <- function(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat){
  
  raiz <- uniroot(time, c(0, 10000), lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=runif(1))
  mepp.time <- raiz$root
  
  return(mepp.time)
}



N = 10
cuts = c(0.9, 1.7)
taxas = c(2, 0.3, 0.6)
alpha = 0.8
betas = c(-0.5,0.2)
matriz.dados = cbind(x1 = rnorm(n = N, mean = 0, sd = 1),
                    x2 = rbinom(n =N, size = 1, prob = 0.5))

set.seed(10)
for (c in 1:10){
  a = gen.mepp(lambda.par=taxas, alpha.par=alpha, grid.vet=cuts, beta.par=betas, x.mat=matriz.dados[c,])
  print(a)
  }

time.cox <- function(time.obs, lambda.par, alpha.par, cuts, beta.par, x.mat, u.unif){
  
  t <- exp(-PPE(time = time.obs, cuts = cuts, 
                levels = lambda.par, alpha = alpha.par, type = "cum_hazard")*exp(x.mat%*%beta.par)) - u.unif
  
  return(t)
  
}


## funcao para encontrar a raiz unitaria que e o tempo de sobrevivencia
gen.mepp.cox <- function(n = n, lambda.par, alpha.par, cuts, beta.par, x.mat){
  
  values = c()
  
  for (ind in 1:n){
    raiz <- uniroot(time.cox, c(0, 1000), lambda.par=lambda.par, alpha.par=alpha.par, 
                    cuts=cuts, beta.par=beta.par,
                    x.mat=x.mat[ind,], u.unif=runif(1))
    mepp.time <- raiz$root
    values[ind] = mepp.time
  }
  
  
  return(values)
}

set.seed(10)
gen.mepp.cox(n = 10, lambda.par = taxas, alpha.par = alpha,
             cuts = cuts, beta.par= betas, x.mat = matriz.dados)



####### FUNCAO PARA CHMAR A TAXA DE FALHA E OUTRAS FUNCOES APARETEMENTE OK!!!!

#--- Algumas fun??es importantes:

# Calculando as fun??es de taxa de falha e acumulada:

sim.std.cure.ICdata <- function(n=n, lambda.par=lambda.par, alpha.par=alpha.par, 
                                grid.vet=grid.vet, beta.par=beta.par, lambda.parc=lambda.parc, 
                                theta.par = c(1, 0.5, 0), A = 5, B = 15){
  
  
  ########################################################
  ## Generating function for mixture cure model using   ##
  ## PPEM for noncured individuals.                     ##
  ##                                                    ##
  ## Arguments:
  ##
  ## n: Sample size;
  ## lambda.par: PPEM rates;
  ## alpha.par: Power PPEM parameter;
  ## grid.vet: Grid time for PPEM;
  ## beta.par: coef. for time to event;
  ## theta.par: coef. for cure;
  ## lambda.parc: Censoring rate;
  ## A: censoring parameter
  ## B: censoring parameter
  ########################################################
  
  
  ## Testing vals:
  
  # n = 100
  # lambda.par  = c(1.1, 0.8, 0.5)
  # alpha.par   = 0.8
  # grid.vet    = c(0.5, 2)
  # beta.par    = c(-0.5, 0.5)
  # lambda.parc = 1
  # theta.par       = c(1, 0.5, 0)
  # A = 5
  # B = 15
  # 
  ########################################
  
  intercept <- 1
  xi1 <- rbinom(n, 1, 0.5)
  xi2 <- rnorm(n)
  X_cure <- cbind(intercept, xi1, xi2)
  X <- cbind(xi1, xi2)
  
  
  #-- Cure probability mixture model:
  
  elinpred <- exp(-theta.par%*%t(X_cure))
  probY <- 1/(1+elinpred)
  Y <- rbinom(n, size=1, prob=probY)
  
  #-- Censoring times generation:
  
  C <- pmin(A, B*rexp(n, rate=lambda.parc))
  
  #-- Generating the survival times:
  
  t     <- rep(0,n) # tempos de falha
  t_obs <- rep(0,n)
  delta <- rep(0,n)
  
  for(i in 1:n){
    if(Y[i]==1){
      t[i] <- gen.mepp(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=X[i,])
      #t[i] <- gen.mepp(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=0)
      if(t[i]<C[i]){
        t_obs[i] <- t[i]
        delta[i] <- 1  
      }else{
        t_obs[i] <- C[i]
      }
    }else{
      t[i] <- C[i]
      t_obs[i] <- t[i] 
    }
  }
  
  tempo <- t_obs
  
  #-- Observed times:
  
  delta <- ifelse(t < C, 1, 0) 
  
  L <- R <- tempo * NA
  for (i in 1:n) {
    if (delta[i] == 0) {
      L[i] <- tempo[i]
      R[i] <- Inf
    }
    else {
      L[i] <- 0
      add <- stats::runif(1, 0.1, 0.5)
      R[i] <- add
      check <- (L[i] <= tempo[i] & tempo[i] < R[i])
      while (!check) {
        L[i] <- L[i] + add
        add <- stats::runif(1, 0.1, 0.5)
        R[i] <- R[i] + add
        check <- (L[i] <= tempo[i] & tempo[i] < R[i])
      }
    }
  }
  
  dados <- data.frame(L, R, tempo, Y, delta, xi1, xi2)
  
  return(dados)
}



cuts = c(0.9, 1.7)
taxas = c(2, 0.3, 0.6)
alpha = 0.8
betas = c(-0.5,0.2)
betas.cura = c(1.2,0.3,0.5)
matriz.dados = cbind(x1 = rnorm(n = N, mean = 0, sd = 1),
                     x2 = rbinom(n =N, size = 1, prob = 0.5))


set.seed(10)
dadosIC <- sim.std.cure.ICdata(n=N, lambda.par=taxas, alpha.par=alpha, 
                               grid.vet=cuts, beta.par= betas, lambda.parc=0.9, 
                               theta.par =betas.cura  , A = 5, B = 15)




sim.cure.icen.mepp = function(N, prob.ber, betas.cure, betas.risk,
                              c1, c2, lambdas, grid, alpha){
  
  ## Armazenameto de dados
  
  survival.time = rep(NA, N)
  L = rep(NA, N)
  R = rep(NA, N)
  delta = rep(NA, N)
  
  ## geracao dos dados
  x2 = rnorm(n = N, mean = 0, sd = 1)
  x1 = rbinom(n = N, size = 1, prob = prob.ber)
  matriz = cbind(rep(1, N), x1, x2)
  
  b = betas.cure
  
  z = matriz %*% b
  
  prob.z = 1/(1+exp(-z))
  
  Y = rbinom(N,1,prob.z)
  
  ## Proporcao de censuras
  
  A = rexp(n = N, rate = 1)
  
  C = cbind(c1,c2*A)
  C = pmin(C[,1], C[,2])
  
  ##------
  for(i in 1:N){
    if(Y[i] == 0){
      survival.time[i] = C[i]
      delta[i] = 0
    }
    
    else{
      #survival.time[i] = exp(matriz[i,-1] %*% betas.risk)
      survival.time[i] = gen.mepp.cox(n = 1, lambda.par = lambdas, alpha.par = alpha,
                                      cuts = grid, beta.par = betas.risk, x.mat = matriz[,-1])
      
      delta[i] = ifelse(survival.time[i] <= C[i], 1, 0)
    }
    
    if(delta[i] == 0){
      L[i] = C[i]
      R[i] = Inf
    }
    
    else{
      L[i] = 0
      Qj = runif(1, min = 0.1, max = 0.5)
      R[i] = Qj
      check = (L[i] <= survival.time[i] & survival.time[i] < R[i])
      
      while(!check){
        L[i] = L[i] + Qj
        Qj = runif(1, min = 0.1, max = 0.5)
        R[i] = R[i] + Qj
        check = (L[i] <= survival.time[i] & survival.time[i] < R[i])
      }
    }
  }
  data = data.frame(time = survival.time, L = L, R = R, 
                    delta = delta, risk = Y, X1 = matriz[,2], X2 = matriz[,3])
  return(data)
}

set.seed(10)
dadosic2 = sim.cure.icen.mepp(N = N, prob.ber = 0.5, betas.cure = betas.cura,
                           betas.risk =betas, c1 = 5, c2 = 15,
                           lambdas = taxas, grid = cuts, alpha = alpha)

####### Parece estarem condizendes as duas geracoes de dados para censura intervalar com fracao de cura


time <- function( t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=u.unif ){
  exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.mat%*%beta.par)) - u.unif
}

gen.mepp <- function(lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat){
  
  raiz <- uniroot(time, c(0, 10000), lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, x.mat=x.mat, u.unif=runif(1))
  mepp.time <- raiz$root
  
  return(mepp.time)
}


SpopMEPP <- function(t=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet, beta.par=beta.par, theta.par=theta.par, x.cure=x.cure, x.risk=x.risk){
  
  
  S_MEPP <- as.numeric(exp(-cal_Ht_MEPP(time.obs=t, lambda.par=lambda.par, alpha.par=alpha.par, grid.vet=grid.vet)*exp(x.risk%*%beta.par)))
  
  elinpred <- as.numeric(exp(1*(x.cure%*%theta.par)))
  probY    <- 1/(1+elinpred)
  spop     <- probY+(1-probY)*S_MEPP
  return(spop)
}


s1pop = SpopMEPP(t=dadosIC$tempo, lambda.par=c(0.6, 0.8, 0.9), alpha.par=1.4,
         grid.vet=c(0.3,0.6), beta.par=c(-0.5,1.2),
         theta.par=c(1, 0.3,0.6),
         x.cure=as.matrix(cbind(1,dadosIC$xi1, dadosIC$xi2)),
         x.risk=as.matrix(cbind(dadosIC$xi1, dadosIC$xi2)))

## -----
## funcao de sobrevivencia para a populacao curada
## -----

spop.mepp = function(time, cuts, levels, alpha, cure.matrix, risk.matrix,
                     betas.cure, betas.risk){
  
  st = as.numeric(PPE(time = time, cuts = cuts, 
                      levels = levels, alpha = alpha, type = "survival"))
  
  stcox = as.numeric(st^(exp(risk.matrix %*% betas.risk)))
  
  xb = as.numeric(exp(cure.matrix %*% betas.cure))
  
  pi = 1/(1+exp(xb))
  
  spop = pi + ((1-pi)*stcox)
  
  return(spop)
}

s2pop = spop.mepp(time = dadosIC$tempo, cuts = c(0.3,0.6), levels = c(0.6, 0.8, 0.9),
               alpha = 1.4,
               cure.matrix = as.matrix(cbind(1,dadosIC$xi1, dadosIC$xi2)),
               risk.matrix = as.matrix(cbind(dadosIC$xi1, dadosIC$xi2)),
               betas.cure = c(1, 0.3,0.6), betas.risk = c(-0.5,1.2))

loglik.int.fc = function(par, time.l, time.r,
                         grid, delta, cure.matrix, risk.matrix){
  
  
  b = length(grid) + 1 ## numero de intervalos
  hazards = par[1:b] ## taxas de falha para os b intervalos
  alpha = par[b + 1] ## parametro de potencia
  
  n.covars.cure = dim(cure.matrix)[2] ## numero de covariaveis com fracao de cura
  n.covars.risk = dim(risk.matrix)[2] ## numero de covariaveis com fracao de risco
  
  betas.cure = par[(b + 2):(b + 1 + n.covars.cure )]
  betas.risk = par[(b + 2 + n.covars.cure):(b + 1 +  n.covars.cure + n.covars.risk)]
  
  lik <- rep(0, dim(cure.matrix)[1])
  
  ## informacoes exponencial por partes Potencia (PPE), sob fracao de cura, para esquerda.
  sl =  spop.mepp(time = time.l[delta==1], cuts = grid, levels = hazards, alpha = alpha,
                  cure.matrix = cure.matrix[delta==1,],
                  risk.matrix = risk.matrix[delta==1,],
                  betas.cure = betas.cure, betas.risk = betas.risk)
  
  ## informacoes exponencial por partes Potencia (PPE), sob fracao de cura, para direita.
  sr =  spop.mepp(time = time.r[delta==1], cuts = grid, levels = hazards, alpha = alpha,
                  cure.matrix = cure.matrix[delta==1,],
                  risk.matrix = risk.matrix[delta==1,],
                  betas.cure = betas.cure, betas.risk = betas.risk)

  ## contribuicao das covariaveis
  
  lik[delta==1] = sl-sr
  
  ## contribuicao da censura
  sl =  spop.mepp(time = time.l[delta==0], cuts = grid, levels = hazards, alpha = alpha,
                  cure.matrix = cure.matrix[delta==0,],
                  risk.matrix = risk.matrix[delta==0,],
                  betas.cure = betas.cure, betas.risk = betas.risk)
  
  lik[delta==0] = sl
  
  log.vero = sum(log(lik))
  
  return(-1*log.vero)
  
}




