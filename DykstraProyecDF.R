
#### función para incluir la proyección del usuario
f_p = function(x0,a,b){
  x = x0 - (a + b)
  return(as.vector(x))
}

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
dykstra_proyec(f_proyec = list("f_p", "P_Half_space"),
               arguments = list(list(6,0,0),list(x= 1, a = 1, b = 3)))

dykstra_proyec(f_proyec = list("f_p"),
               arguments = list(list(6,7,8)))


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

dykstra_proyecDF(DF = DF, f_proyec = list("f_p", "P_Half_space"),
                 arguments = list(list(6,0,0),list(x= c(1,0,0), a = c(1,1,1), b = 1)))

dykstra_proyecDF(DF = DF, f_proyec = list("P_Half_space"),
                 arguments = list(list(x= c(1,0,0), a = c(1,1,1), b = 1)))

dykstra_proyecDF(DF = DF, f_proyec = list("f_p"),
               arguments = list(list(6,7,8)))

dykstra_DF = function(DF,
                      typeConstraints = 'list() or c() of characters',# "ball", "lineal" "other"
                      # lo puedo pedir como un vector de caracteres(o una lista)
                      f_proyec = "list of functions returning vec",
                      arguments = "a list of arguments for the function f_proyec",
                      A = 'matrix',
                      b = 'vector',
                      r = rep(1,1),
                      centers = matrix(rep(0,length(DF[1,])*length(r)),
                                      nrow = length(r),
                                      ncol = length(DF[1,])),
                      eq = c(), #"vector of characteres: '<=', '==' ",
                      I = matrix(rep(1,dim(DF)[1]*dim(DF)[2]),
                                 ncol = dim(DF)[2],
                                 nrow = dim(DF)[1]),
                      W = diag(length(DF[1,]))){
  
  if (all(list("lineal", "ball","other") %in% typeConstraints)) {
    
    ###### ball and lineal constraints    
  } else if (all(list("lineal", "ball") %in% typeConstraints) & !("other" %in% typeConstraints)){
    if (length(eq) == 0) eq = rep("<=",length(b)+length(r)) 
    results = dykstra_linealBallDF(DF, A = A, b = b,
                                   r = r, centers = centers,
                                   eq = eq, I = I, W = W)
    
    ######## only type constraints define by user ############
  } else if ("other" %in% typeConstraints) {
    results = dykstra_proyecDF(DF = DF, f_proyec = f_proyec, arguments = arguments, W = W)
    
    
    ####### only type ball constraints 
  } else if ("ball" %in% typeConstraints){
    if (length(eq) == 0) eq = rep("<=",length(r))
    results = dykstra_ballDF(DF = DF, r = r, centers = centers, eq = eq, I = I, W = W)
    
    ######## only lineal constraints ############
  } else if ("lineal" %in% typeConstraints){
    if (length(eq) == 0) eq = rep("<=",length(b))
    results = dykstra_linealDF(DF = DF, A = A, b = b, eq = eq, I = I, W = W)
    
  }
  return(results)
}

dykstra_DF(DF = DF, typeConstraints = c("other"), 
           f_proyec = list("f_p", "P_Half_space"),
           arguments = list(list(6,0,0),list(x= c(1,0,0), a = c(1,1,1), b = 1)))

dykstra_DF(DF, typeConstraints = "lineal", A = A, b = b)
dykstra_DF(DF, typeConstraints = "lineal", A = A, b = b, I = I)
dykstra_DF(DF, typeConstraints = "ball")
#dykstra_DF(DF, typeConstraints = "ball", I = I)
dykstra_DF(DF, typeConstraints = "ball",
           r = c(1,2), 
           centers = matrix(c(0,0,0,2,0,0),
                            nrow = 2,
                            ncol = 3,byrow = TRUE))
#dykstra_DF(DF, typeConstraints = c("lineal","ball"), A = A, b = b)
#dykstra_DF(DF, typeConstraints = c("lineal","ball"), A = A, b = b, I = I)

