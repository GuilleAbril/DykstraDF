#### function for including users prjection

dykstra_proyec = function(f_proyec = list(), arguments = list(),
                          W = diag(length(arguments[[1]][[1]]))){
  x0 = arguments[[1]][[1]]
  M = length(f_proyec) ## nº de subconjuntos convexos (restricciones)
  N = length(x0) ## dimension del vector a proyectar = nº de variables
  E = matrix(rep(0,N*M),nrow = M, ncol = N)
  n = 1
  # variables auxiliares para calcular el error/criterio de parada
  error_vec = c()
  error = 1
  while ((error>0.0001 || n<=M) & n<100) {
    i = 1+(n-1)%%M
    x0e = as.numeric(x0) + E[i,]
    arguments[[i]][[1]] = x0e
    x = do.call(f_proyec[[i]],arguments[[i]])
    e = as.numeric(x0) + E[i,] - x
    error_vec[i] = calc_error(e,E[i,],x,x0,W = W)
    error = sum(error_vec)
    E[i,] = e
    x0 = x
    n = n + 1
  }
  return(x)
}

dykstra_proyecDF = function(DF, f_proyec = list(), arguments = list(),
                            W = diag(length(arguments[[1]][[1]]))){
  # en principio al no conocer como son las restricciones a priori
  # no puedo checkear que filas las cumplen o no
  df_final = data.frame()
  df_proyec = data.frame()
  for (i in 1:length(DF[,1])) {
    arguments[[1]][[1]] = DF[i,]# el primer argumento de la función es el vector a proyectar
    x = dykstra_proyec(f_proyec = f_proyec, arguments = arguments, W = W)
    df_final = rbind(df_final,x)
    if (any(x != DF[i,], na.rm = TRUE)){df_proyec = rbind(df_proyec,x)}
  }
  colnames(df_final) = colnames(DF)
  colnames(df_proyec) = colnames(DF)
  return(list(df_proyec = df_proyec,
              df_original = DF,
              df_final = df_final))
}
