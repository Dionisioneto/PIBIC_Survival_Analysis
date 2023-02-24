
## versao 1 do monte carlo

# 
# for (i in 1:n.iter){
#   
#   cat("Realizando iteracao: ", i, "/", n.iter, "\n", sep = "")
#   
#   ## geracao de dados do tempo de falha e tempo com censura
#   tempo.falha = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
#                              cuts = particoes, x.mat = x.matriz, beta.par = beta)
#   
#   tempo.censura = gen.mepp.cox(n = tamanho.amostral, lambda.par = taxas.falha, alpha.par = potencia,
#                                cuts = particoes, x.mat = x.matriz, beta.par = beta)
#   
#   ## geracao do tempo e censura
#   tempo = pmin(tempo.falha, tempo.censura)
#   delta = ifelse(tempo.falha <= tempo.censura, 1, 0)
#  
#   
#   ## particao para a estimacao
#   grid = time.grid.obs.t(tempo, delta, n.int = 3)
#   grid = grid[-c(1, length(grid))]
#   
#   chutes = c(rep(0.5,length(grid)+1),1,0.5,1)
#   
#   
#   ## Metodo numerico BFGS
#   estimacao.teste.cox = optim(par = chutes,
#                               fn = loglik.cox,
#                               gr = NULL,
#                               hessian = TRUE,
#                               method = "BFGS",
#                               tempos = tempo,
#                               censura = delta,
#                               intervalos = grid,
#                               covariaveis = x.matriz)
#   
#   ## salvar resultados
#   matrix.iter[i,1:length(taxas.falha)] = estimacao.teste.cox$par[1:length(taxas.falha)]
#   matrix.iter[i,length(taxas.falha) + 1] = estimacao.teste.cox$par[length(taxas.falha) +1]
#   matrix.iter[i, (length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))] = estimacao.teste.cox$par[(length(taxas.falha)+2):(length(taxas.falha)+1+length(beta))]
#   
#   matrix.var[i,] = diag(solve(-estimacao.teste.cox$hessian))
# }
