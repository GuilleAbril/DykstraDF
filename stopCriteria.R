# stop criteria
## 2.4 c_E = sum((E-E0)W(E-E0)), suma sobre los subconjuntos
## 2.5 C_E + 2sum(E_0*W*(x-x0))
calc_error = function(e,e0,x,x0,W){
  c_i = (e-e0)%*%W%*%(e-e0) # fórmula 2.4 para un subconjunto
  error_i = abs(c_i + 2*as.numeric(e0)%*%W%*%as.numeric(x-x0)) # fórmula 2.5 para un subconjunto
  # La salida de esta función se va almacenando en un vector en la función principal
  # y el sumatorio a todos los subconjuntos se realiza en esta
  return(error_i)
}



