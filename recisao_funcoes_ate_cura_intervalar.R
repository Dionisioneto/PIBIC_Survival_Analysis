
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



















