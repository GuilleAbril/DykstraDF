
# primero porbaré con un dykstra para una bola de radio 1 y centrada en el origen
# proyección bola r=1 centrada en el origen
P_ball = function(x, W = diag(length(x))){
  xx = x%*%W%*%x
  Pb = 1/max(sqrt(xx),1)*x
  return(Pb)
}

P_circunference = function(x, W = diag(length(x))){
  xx = as.numeric(x%*%W%*%x)
  Pb = 1/(sqrt(xx))*x
  return(Pb)
}
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

x0 = c(2,3)
dykstra_ball(x0)
dykstra_ball(x0, r = c(2), centers = matrix(c(2,0),ncol = 2, nrow = 1))
dykstra_ball(c(0,2), r = c(1,2), centers = matrix(c(0,0,2,0),nrow = 2, ncol = 2,byrow = TRUE)) 
# esta última, proyección del algoritmo de dykstra para:
# una bola de radio 1 centrada en el origen y una de radio 2 centrada en (0,2)
# converge en la iteración 44 sol.(0.25,0.9682) 

col1b = c(0,0.3,0.3,0.6)
col2b = c(0,1,0.4,0.7)
col3b = c(-0.5,0.5,1,3)
df_bola = data.frame(col1b,col2b,col3b)

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

check_rows_ball(df_bola, r = rep(1,1))

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

dykstra_ballDF(df_bola)
dykstra_ballDF(df_bola, I = I)

######## si escalamos o no los vectores E da el mismo resultado
dykstra_ball2 = function (x0,
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
    x0_prima = (x0 - centers[i,])/r[i] # reescalamos (cambiamos el sist. de coord) para proyectar sobre bola unidad centrada en origen
    x0_primae = as.numeric(x0_prima) + E[i,]
    if (eq[i] == "<="){
      x = P_ball(x0_primae,W) # proyectamos sobre la "nueva" bola unidad (en el nuevo sist.)
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
x0 = c(2,3)
dykstra_ball2(x0)
dykstra_ball2(x0, r = c(2), centers = matrix(c(2,0),ncol = 2, nrow = 1))
dykstra_ball2(c(0,2), r = c(1,2), centers = matrix(c(0,0,2,0),nrow = 2, ncol = 2,byrow = TRUE)) 
#> x0 = c(2,3)
#> dykstra_ball2(x0)
#[1] 0.5547002 0.8320503
#> dykstra_ball2(x0, r = c(2), centers = matrix(c(2,0),ncol = 2, nrow = 1))
#[1] 2 2
#> dykstra_ball2(c(0,2), r = c(1,2), centers = matrix(c(0,0,2,0),nrow = 2, ncol = 2,byrow = TRUE))
#[1] 0.2499944 0.9682357
#> dykstra_ball(x0)
#[1] 0.5547002 0.8320503
#> dykstra_ball(x0, r = c(2), centers = matrix(c(2,0),ncol = 2, nrow = 1))
#[1] 2 2
#> dykstra_ball(c(0,2), r = c(1,2), centers = matrix(c(0,0,2,0),nrow = 2, ncol = 2,byrow = TRUE))
#[1] 0.2499699 0.9682536


