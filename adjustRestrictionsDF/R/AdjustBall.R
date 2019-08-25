
# algoritmo de Dykstra
# esto aplicaría el algoritmo de dykstra sobre un vector x0 de N variables 
# de radio ri (elemento i-ésimo del vector r) 
# y con centro ci que será un vector de la dim de x0 (fila i-ésima de la matriz centers),
# si en el DF para las distintas variables hay que proyectar sobre bolas distintas este filtro
# se haría en la función dykstraDF_bola
dykstra_ball = function (x0,
                         r = rep(1,1),
                         centers = matrix(rep(0,length(x0)*length(r)),nrow = length(r), ncol = length(x0)),
                         eq = rep("<=",length(r)),
                         W= diag(length(x0))){
  M = length(r) ## nº de bolas(subconjuntos convexos = restricciones)
  N = length(x0) ## dimension del vecotr a proyectar = nº de variables
  E = matrix(rep(0,N*M),nrow = M, ncol = N)
  n = 1
  # variables auxiliares para calcular el error/criterio de parada
  error_vec = c()
  error = 1
  #print(x0)
  while ((error>0.0001 || n<=M) & n<100) {
    i = 1+(n-1)%%M
    x0e = as.numeric(x0) + E[i,]
    x0e_prima = (x0e - centers[i,])/r[i]
    if (eq[i] == "<="){
      x = P_ball(x0e_prima,W) # proyectamos sobre la "nueva" bola unidad (en el nuevo sist.)
    } else {x = P_circunference(x0_primae,W)}
    x = x * r[i] + centers[i,]# recuperamos el sistema de coordenadas anterior
    e = as.numeric(x0) + E[i,] - x
    error_vec[i] = calc_error(e,E[i,],x,x0,W = W)
    error = sum(error_vec)
    E[i,] = e
    x0 = x
    n = n+1
    #cat(x ,n, "\n")
    #print(e0)
    #print(error)
  }
  return(x)
}

check_rows_ball = function(DF,
                          r = rep(1,1),
                          centers = matrix(rep(0,length(DF[1,])*length(r)),nrow = length(r), ncol = length(DF[1,])),
                          eq = rep("<=",length(r)),
                          W= diag(length(DF[1,]))){
  
  wrong_rows = data.frame()
  for (i in 1:length(DF[,1])) {
    xi = DF[i,]
    checks_vec = c()
      for(k in 1:length(r)){
      r_prima = (as.numeric(xi)-centers[k,])%*%W%*%(as.numeric(xi)-centers[k,])
        if(eq[k] == "<="){
          checks_vec[k] = as.vector(r_prima)<=r[k]} else {
            checks_vec[k] = as.vector(r_prima)==r[k]
          }
        }
      if(any(checks_vec == FALSE, na.rm = TRUE)){
        wrong_rows = rbind(wrong_rows,xi)#data frame con los vectores que no cumplen las restricciones(y su índice original)
      }
    }
  return(wrong_rows)
}

dykstra_ballDF = function(DF,
                          r = rep(1,1),
                          centers = matrix(rep(0,length(DF[1,])*length(r)),
                                           nrow = length(r),
                                           ncol = length(DF[1,])),
                          eq = rep("<=",length(r)),
                          I = matrix(rep(1,length(DF[1,])*length(DF[,1])),
                                     ncol = length(DF[1,]),
                                     nrow = length(DF[,1])),
                          W= diag(length(DF[1,]))){
  df_check = check_rows_ball(DF = DF,r = r,centers = centers,eq = eq,W=W)
  df_final = data.frame()
  df_proyec = data.frame()
  for (i in 1:length(DF[,1])) {
    xi = DF[i,]
    ii = which(I[i,][1:length(I[i,])]==1)
    if (any(I[i,]==1)) {
      if (all(I[i,]==1)) {
        x = dykstra_ball(xi,r,centers,eq,W)
      } else {
        xi_prima = xi[ii]
        centers_prima = t(centers[,ii])# traspuesta para que lo ponga como una matriz por filas
        r_prima = c()
        for (k in 1:length(r)) { # bucle para ir recorriendo las bolas (cada k un subespacio)
          r_prima[k] = sqrt(abs(r^2 - as.numeric(xi[-ii])%*%as.matrix(W[,-ii])[-ii,]%*%as.numeric(xi[-ii])))
          # no debería hacer falta el abs(si es dato xi[-ii] no imputado es correcto la resta respecto de su centro no debería dar negativo)
        }
        Wi = t(as.matrix(W[,ii]))[,ii]
        x = dykstra_ball(xi_prima,r_prima,centers_prima,eq,Wi)# tengo que añadir a este vector x las coordenadas que no se cambiaron, porque no se imputaron
        }
      xi[ii] = x
      df_final = rbind(df_final,xi)
      if (any(xi!= DF[i,], na.rm = TRUE)){ df_proyec = rbind(df_proyec,xi)}
      
    } else {
      df_final = rbind(df_final,xi)}
  }
  return(
    list(df_proyec = df_proyec, 
         df_check = df_check,
         df_original = DF,
        df_final= df_final
    ))
  }

