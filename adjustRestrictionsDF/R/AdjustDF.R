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
