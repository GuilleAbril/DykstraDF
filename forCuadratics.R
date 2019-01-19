######### for cuadratics #########
#### que es lo que queremos
### 1/2xWAx + xWu <= beta

dykstra_cuadratics = function()
{
  #lo que sea
}
### para ejemplo
# Términos cuadráticos:
A = c(1,2,3)
# términos lineales
u = c(2,4,5)
# termino indep. beta = (inf(f(R^n)),f(x0))
# con f(x) = 1/2*(x*A*x+x*u)
# P(x0)(=x) es la sol. del sist.:
###   f(x) = beta
###   x = ((Id+t*A)^-1)(z-t*u) => (Id+t*A)x = x0-t*u => (Id+t*A)x - (x0-t*u) = 0
###   t > 0

## Usar nleqslv para resolver el sistema
## (futura posibilidad: diagonalizar antes A(base ortogonal para esta), resolver el sistema
## y deshacer el cambio de base...)
library(nleqslv)
Ac = matrix(c(1,0,0,1), nrow = 2, byrow = TRUE)
u = c(3,0)
beta = -11
Ac2 = matrix(c(1,2,3,1),nrow = 2) 
project_cuadratics = function(x0,A,u,beta){

  fn = function(xt){
    x = xt[1:(length(xt)-1)]
    t = xt[length(xt)]
    y = numeric(length(xt))
    
    f0 = 1/2*(x%*%A%*%x)+x%*%u-beta
    fi = (diag(length(x0))+t*A)%*%x - (x0-t*u)
    f = c(f0, fi)
    
    return(f)
  }
  
  Jac <- function(x) {
    J <- matrix(0,nrow=length(x),ncol=length(x))
    for (i in 1:(length(x)-1)) {
      J[,1][i] <- c(a_1*x1+u_1,a_2*x_2+u_2,...,a_N*x_N+u_N,
                    +a_1)  # parciales de todas las f respecto de la primera variable
    }
    J[,2] <- c(1-x[1],0,1)
    J[,3] <- c(0,1-x[1],1)
    J
  }
  
  return(x)
}

#### para el dibujo, caso 2D
r = 2
a = 2
b = 3
c = 3
x = seq(-5,5,0.1)
y = a*x^2 + b*x + c
plot(x, y,type = "l")
z1 = x ^2
z2 = y
plot(z1,z2,type = "l")
# a mano
#z2 = a*z1 + b*sqrt(z1) + c igual pero no considera la rama negativa de la raíz
plot(x,r-x^2)
