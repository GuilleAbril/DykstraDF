
## Funciones Dykstra

# DF ejemplo para ir probando
col1 = c(0,1,0,2)
col2 = c(0,7,1,1)
col3 = c(-3,0,1,3)
DF = data.frame(col1,col2,col3) #  4 individuos 3 variables
#W = matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, ncol = 3)

# Input de restricciones:
# 1. se los podemos pedir como una matriz A, con los coeficientes que acompañan a cada variable
  # en el orden que aparezcan en su DF
  # y que la columna de los coeficientes indepen. nos los de con un vector a parte
  a1 = c(1,1,1)
  a2 = c(1,2,3)
  a3 = c(1,0,0)
  A = matrix(c(a1, a2, a3), nrow = 3, ncol = length(a1), byrow = TRUE)
  b = c(1,-6,3)
  # o como una matriz donde la última columna sean los coef. indepen.
  Ab = matrix(c(a1, a2, a3, b), nrow = 3, ncol = length(a1)+1, byrow = FALSE)

  
##### Función que comprueba que filas del DF no cumplen las restricciones
  # devuelve estas solo
## cosas: estamos suponiendo todo desigualdades, hay que meter la opción para igualdades
## opción para hacer esto que el usuario diga a que filas de la matriz A le corresponde una igualdad
## y cuales una desigualdad (si es por orden mejor):
## Una manera, que nos diga un vector de la misma longitud que b con caracteres de la forma 
## "==" o "<=" en la posición correspondiente al b_i.
## por ejemplo para la 1ª y 3ª desigualdades y para la 2ª igualdad sería
check_rows_lineal = function(df,A,b,eq = rep('<=',length(b))){
  wrong_rows = data.frame()
  for (i in 1:length(df[,1])) {
    xi = df[i,]
    b_prima = A%*%as.numeric(xi)
    checks_vec = c()
    for(k in 1:length(b)){
      if(eq[k] == "<="){
      checks_vec[k] = as.vector(b_prima)[k]<=b[k]} else {
        checks_vec[k] = as.vector(b_prima)[k]==b[k]
        }
      }
      if(any(checks_vec == FALSE)){
        wrong_rows = rbind(wrong_rows,xi)#data frame con los vectores que no cumplen las restricciones(y su índice original)
      }
    }
  return(wrong_rows)
}

igualdades = c("<=", '<=', '<=')
df_wrong_rows = check_rows_lineal(DF,A,b,igualdades)
igualdades1 = c("<=", '==', '<=')
df_wrong_rows_eq = check_rows_lineal(DF,A,b,igualdades1)
# en este caso ninguna fila satisface las restricciones, la 1 primera fila ahora no cumple
# la segunda igualdad (probar al retailers)

#proyección sobre semiespacio
P_Half_space = function (x,a,b,W = diag(length(x))){
  xa = as.numeric(as.numeric(x)%*%W%*%a)
  aa = as.numeric(a%*%W%*%a)
  if (aa == 0) {  # evitar la división entre cero...
    pki = x
  } else {
    alpha = max(xa - as.numeric(b),0)
    pki = x-(1/aa)*alpha*a
  }
  return(pki)
}
# proyección sobre subespacios
P_subspace = function(x,a,b,W = diag(length(x))){
  xa = as.numeric(as.numeric(x)%*%W%*%a)
  aa = as.numeric(a%*%W%*%a)
  if (aa == 0) {
    pk = x
  } else {
  pk = x-(1/aa)*(xa - as.numeric(b))*a
  }
  return(pk)
}

## Función que implementa el algoritmo de Dykstra a un punto, 
## de momento con la proyección lineal

dykstra_lineal = function (x0,A,b,eq,W = diag(length(x0))){
  M = length(as.vector(b)) ## nº de subconjuntos convexos (restricciones)
  N = length(x0) ## dimension del vector a proyectar = nº de variables
  E = matrix(rep(0,N*M),nrow = M, ncol = N)
  n = 1
  # variables auxiliares para calcular el error/criterio de parada
  error_vec = c()
  error = 1
  #print(x0)
  while ((error>0.0001 || n<=M) & n<100) {# criterio de parada, obligamos a recorrer todos los 
    # subconjuntos antes de que tenga en cuenta el error
    i = 1+(n-1)%%M
    x0e = as.numeric(x0) + E[i,]
    if (eq[i]== "==") {
      x = P_subspace(x0e,as.matrix(A)[i,],as.vector(b)[i],W)
    } else {
      x = P_Half_space(x0e,as.matrix(A)[i,],as.vector(b)[i],W)
    }
    e = as.numeric(x0) + E[i,] - x
    ## calcular aquí el criterio de parada (el error):
    error_vec[i] = calc_error(e,E[i,],x,x0,W = W)
    error = sum(error_vec)
    # redefinnimos la varibles iniciales
    E[i,] = e
    x0 = x
    n = n+1
    #print(error)
  }
  return(x)
}

dykstra_lineal(DF[2,],A,b,eq = rep("<=",length(b)))

# índces que marcan el lugar donde se ha imputado una variable
i1 = c(0,0,1) # 1ª fila: todos datos originales
i2 = c(0,1,0) # 2ª fila: 2ª variable imputada, 1ª y 3ª originales
i3 = c(1,0,1) # 3ª fila: 1ª y 3ª variables imputadas
i4 = c(0,0,1) # 4ª fila: la última variable imputada
I = rbind(i1,i2,i3,i4) #i2 = I[2,]

# función que dado el data frame que el usuario quiere ajustar a las restricciones
# ajusta solo las variables que diga el usuario para cada vector, medinate la matriz I de índices
dykstra_linealDF = function(DF, A, b, 
                                           eq = rep("<=",length(b)), 
                                           I = matrix(rep(1,dim(DF)[1]*dim(DF)[2]), 
                                                      ncol = dim(DF)[2],
                                                      nrow = dim(DF)[1]),
                                           W = diag(length(DF[1,]))){
  df_check = check_rows_lineal(DF,A,b,eq) #comprobamos que filas no cumplen las restricciones
  # asi podemos comparar si había alguna que no tenía valores imputados que supiese el ususario
  # pero que no las cumpliese. hacer esto
  df_final = data.frame()
  df_proyec = data.frame()
    for (i in 1:length(DF[,1])) {
    xi = DF[i,]
    ii = which(I[i,][1:length(I[i,])]==1)# índice que marca la posición del dato que queremos cambiar en el vector i
    xi_prima = xi[ii] #'vector' con los elementos imputados para el individuo i, nuevo vector a proyectar
    A_prima = as.matrix(A[,ii]) # nueva matriz A de las restricciones, nuevos subconjuntos donde proyectar
    Wi = t(as.matrix(W[,ii]))[,ii]
    bi_prima = b - as.matrix(A[,-ii])%*%as.numeric(xi[-ii]) # términos indep. que definen los nuevos subconjuntos donde proyect.
    x = dykstra_lineal(xi_prima,A_prima,bi_prima,eq,Wi)# tengo que añadir a este vector x las coordenadas que no se cambiaron, porque no se imputaron
    xi[ii] = x
    df_final = rbind(df_final,xi)
    if (any(xi!= DF[i,])){ df_proyec = rbind(df_proyec,xi)}
  }
  return(
    list(df_proyec = df_proyec, 
         df_check = df_check,
         df_original = DF,
         df_final= df_final
         ))
}

dykstra_linealDF(DF, A, b, igualdades, I = I)
dykstra_linealDF(DF, A = A,  b = b, eq = igualdades)
    
# Try to change bucles in DFs to DataTable options or apply
# in functions Check_rows and Dykstra_lineal


# first lets proof apply

##### Función que comprueba que filas del DF no cumplen las restricciones
# devuelve estas solo
## cosas: estamos suponiendo todo desigualdades, hay que meter la opción para igualdades
## opción para hacer esto que el usuario diga a que filas de la matriz A le corresponde una igualdad
## y cuales una desigualdad (si es por orden mejor):
## Una manera, que nos diga un vector de la misma longitud que b con caracteres de la forma 
## "==" o "<=" en la posición correspondiente al b_i.
## por ejemplo para la 1ª y 3ª desigualdades y para la 2ª igualdad sería

## function which check if a vector satisfy the linear restrictions
## and return this vector if is not the case
check_vector_lineal = function(x, A, b, eq = '<='){
  b_prima = A%*%x
  logical_vect = b <= b_prima
  if(any(logical_vect == FALSE)){
    return(x)
  }
}

check_rows_lineal = function(df,A,b,eq = rep('<=',length(b))){
  df_checks = apply(df, 1, check_vector_lineal, A = A, b = b, eq)
  wrong_rows = data.frame()
  for (i in 1:length(df[,1])) {
    xi = df[i,]
    b_prima = A%*%as.numeric(xi)
    checks_vec = c()
    for(k in 1:length(b)){
      if(eq[k] == "<="){
        checks_vec[k] = as.vector(b_prima)[k]<=b[k]} else {
          checks_vec[k] = as.vector(b_prima)[k]==b[k]
        }
    }
    if(any(checks_vec == FALSE)){
      wrong_rows = rbind(wrong_rows,xi)#data frame con los vectores que no cumplen las restricciones(y su índice original)
    }
  }
  return(wrong_rows)
}

igualdades = c("<=", '<=', '<=')
df_wrong_rows = check_rows_lineal(DF,A,b,igualdades)
igualdades1 = c("<=", '==', '<=')
df_wrong_rows_eq = check_rows_lineal(DF,A,b,igualdades1)
# en este caso ninguna fila satisface las restricciones, la 1 primera fila ahora no cumple
# la segunda igualdad (probar al retailers)
