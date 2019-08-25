## Projection functions

#projection on Half-Spaces
P_Half_space = function (x,a,b,W = diag(length(x))){
  xa = as.numeric(as.numeric(x)%*%W%*%a)
  aa = as.numeric(a%*%W%*%a)
  if (aa == 0) {  # avoid divide by 0
    pki = x
  } else {
    alpha = max(xa - as.numeric(b),0)
    pki = x-(1/aa)*alpha*a
  }
  return(pki)
}
# projection on Subspaces
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

# Prjections on ball
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

