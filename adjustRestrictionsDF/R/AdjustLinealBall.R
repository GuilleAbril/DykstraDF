

dykstra_linealBall = function(x0,
                              A = 'matrix',
                              b = 'vector',
                              r = rep(1,1),
                              centers = matrix(rep(0,length(x0)*length(r)),
                                               nrow = length(r), 
                                               ncol = length(x0)),
                              eq = rep("<=",M),
                              W = diag(length(x0))){
    M = length(r) + length(b) ## nº de bolas + conj. lin. (subconjuntos convexos = restricciones)
    N = length(x0) ## dimension del vector a proyectar = nº de variables
    E = matrix(rep(0,N*M),nrow = M, ncol = N)
    n = 1
    # variables auxiliares para calcular el error/criterio de parada
    error_vec = c()
    error = 1
    row.names(centers) = c((length(b)+1):M)
    names(r) = c((length(b)+1):M)
    while ( (error>0.0001 || n<=M) & n<1000) { # criterio de parada, obligamos a recorrer todos los 
                                            # subconjuntos antes de que tenga en cuenta el error
      i = 1+(n-1)%%M # esta i recorre todos los conjuntos ciclicamente
      if (i <= length(b)) { # primero aplicamos Dykstra sobre los subespacios
        x0e = as.numeric(x0) + E[i,]
        if (eq[i]== "==") {
          x = P_subespace(x0e,as.matrix(A)[i,],as.vector(b)[i],W)
        } else {
          x = P_Half_space(x0e,as.matrix(A)[i,],as.vector(b)[i],W)
        }
        e = as.numeric(x0) + E[i,] - x
        error_vec[i] = calc_error(e,E[i,],x,x0,W = W)
        error = sum(error_vec)
        E[i,] = e
        x0 = x
      } else { # ahora sobre las restricciones tipo bola
        i_char = as.character(i)
        x0_prima = (x0 - centers[i_char,])/r[i_char] # reescalamos (cambiamos el sist. de coord) para proyectar sobre bola unidad centrada en origen
        x0_primae0 = as.numeric(x0_prima) + E[i,]
        if (eq[i] == "<="){
          x = P_ball(x0_primae0,W) # proyectamos sobre la "nueva" bola unidad (en el nuevo sist.)
        } else {x = P_circunference(x0_primae0,W)}
        x = x * r[i_char] + centers[i_char,]# recuperamos el sistema de coordenadas anterior
        e = as.numeric(x0) + E[i,] - x
        error_vec[i] = calc_error(e,E[i,],x,x0,W = W)
        error = sum(error_vec)
        E[i,] = e
        x0 = x
      }
      n = n+1
      #print(n)
      #print(error)
      #print(x)
    }
  return(x)
}
# tengo calcular el error cuando ha pasado al menos una vez por todos los subespacios,
# porque puede ser como en este caso que el error para el primer subespacio sea cero, porque ya cumple
# esa restricción y el algoritmo ya no proyecta más
                       
dykstra_linealBallDF = function(DF,
                                A = 'matrix',
                                b = 'vector',
                                r = c(1),
                                centers = matrix(rep(0,length(DF[1,])*length(r)),
                                                 nrow = length(r),
                                                 ncol = length(DF[1,])),
                                eq = rep("<=",M),
                                I = matrix(rep(1,dim(DF)[1]*dim(DF)[2]), 
                                           ncol = dim(DF)[2],
                                           nrow = dim(DF)[1]),
                                W = diag(length(DF[1,]))){
  M = length(r) + length(b)
  df_checkLineal = check_rows_lineal(DF,A,b,eq[1:length(b)])
  df_checkBall = check_rows_ball(DF,r,centers,eq[(length(b)+1):M])
  df_check = merge(df_checkLineal,df_checkBall, all = TRUE)
  df_final = data.frame()
  df_proyec = data.frame()
  for (i in 1:length(DF[,1])) {
    xi = DF[i,]
    ii = which(I[i,][1:length(I[i,])]==1)
    if (any(I[i,]==1)) {
      if (all(I[i,]==1)) {
        x = dykstra_linealBall(xi,A,b,r,centers,eq,W)
      } else {
        xi_prima = xi[ii]
        Wi = t(as.matrix(W[,ii]))[,ii]
        A_prima = as.matrix(A[,ii])
        bi_prima = b - as.matrix(A[,-ii])%*%as.numeric(xi[-ii])
        centers_prima = t(centers[,ii])
        r_prima = c()
        for (k in 1:length(r)) { # bucle para ir recorriendo las bolas (cada k un subespacio)
          r_prima[k] = sqrt(abs(r^2 - sum((xi[-ii]-centers[k,][-ii])^2)))#
        }
        x = dykstra_linealBall(xi_prima,
                               A = A_prima, b = bi_prima, 
                               r = r_prima, centers = centers_prima, eq, Wi)
      }
      xi[ii] = x
      df_final = rbind(df_final,xi)
      #df_proyec = rbind(df_proyec,xi)
      if (any(xi!= DF[i,], na.rm = TRUE)){df_proyec = rbind(df_proyec,xi)}
      
    } else {df_final = rbind(df_final,xi)}
  }
  return(
    list(df_proyec = df_proyec, 
         df_check = df_check,
         df_original = DF,
         df_final= df_final
    ))
}

  