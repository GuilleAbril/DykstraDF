### probar a representar gráficamente en 2D el procedimiento (las p es de plot)
#####################################
# proyección sobre semiplanos
#a1p = c(-2,1)
a1p = c(0.01,1)
a2p = c(2,1)
Ap = rbind(a1p, a2p)
bp = c(0.5,-1)
# x0 = c(0,0)
x0 = c(2,1)
W = rbind(c(1,0),c(0,1))
plot(x0[1],x0[2], xlim = c(-1,2), ylim = c(-1.5,1.5), asp = 1,
     xlab = expression('x'[1]),
     ylab = expression('x'[2]))#,xlim = c(-0.6,0.5), ylim = c(-1.3,0.3)

abline(bp[1]/a1p[2],-a1p[1]/a1p[2])
abline(bp[2]/a2p[2],-a2p[1]/a2p[2])
#polygon(x= seq(-2,-0.5,0.01),y = -0.01*seq(-2,0.2,0.01)+1,col= "grey",border = "red")
M = length(as.vector(bp)) ## nº de subconjuntos convexos (restricciones)
N = length(x0) ## dimension del vecotr a proyectar = nº de variables
e0 = matrix(rep(0,N*M),nrow = M, ncol = N)
#print(x0)
for (n in 1:10) {
  i = 1+(n-1)%%M
  x0e0 = as.numeric(x0) + e0[i,] # vector a proyectar, el cambio respecto a los subespacios
  lines(c(x0[1],x0e0[1]),c(x0[2],x0e0[2]), col='blue')
  x = P_Half_space(x0e0,as.matrix(Ap)[i,],as.vector(bp)[i],W)
  lines(c(x0e0[1],x[1]),c(x0e0[2],x[2]), col='red')
  #lines(c(x0[1],x[1]),c(x0[2],x[2]), col='red')
  #e = as.numeric(x0) + e0[i,] - x 
  #e0[i,] = e
  x0 = x
  #print(x)
}
#########################
# dibujo proyección sobre un semiespacio y una bola (sb)
#x0 = c(0.8,2)
x0 = c(3,0)
asb = c(0,1)
bsb = -3/4
plot(3,0,xlim = c(-1.5,3), ylim = c(-1.5,2.5), asp = 1,
     xlab = expression('x'[1]),
     ylab = expression('x'[2]))
# Bola unidad
#plot(0,0,cex=0.1,xlim=c(0.5,3),ylim=c(-1,0),xlab="",ylab="",asp=1)
r <- 1
xbol <- seq(-r,r,0.001)
ybol <- sqrt(r^2-xbol^2)
polygon(c(xbol,-xbol),c(ybol,-ybol),cex=0.1)
# Bola y semiespacio
h <- sqrt(7)/4
xbolH <- seq(-h,h,0.001)
ybolH <- -sqrt(r^2-xbolH^2)
polygon(xbolH,ybolH,col="gray")
library(plotrix)
draw.circle(0,0,1)
abline(bsb/asb[2],asb[1]/asb[2])
e0 = matrix(rep(0,2*2),nrow = 2, ncol = 2)
for (n in 1:6) {
  i = 1 #circulo
  x0e0 = as.numeric(x0) + e0[i,]
  lines(c(x0[1],x0e0[1]),c(x0[2],x0e0[2]), col='blue')
  x = P_ball(x0e0)
  lines(c(x0e0[1],x[1]),c(x0e0[2],x[2]), col='red')
  e = as.numeric(x0) + e0[i,] - x 
  e0[i,] = e
  x0 = x
  i = 2 #semiespacio
  x0e0 = as.numeric(x0) + e0[i,]
  lines(c(x0[1],x0e0[1]),c(x0[2],x0e0[2]), col='blue')
  x = P_Half_space(x0e0,asb,bsb,W=diag(length(x0)))
  lines(c(x0e0[1],x[1]),c(x0e0[2],x[2]), col='red')
  #lines(c(x0[1],x[1]),c(x0[2],x[2]), col='red')
  e = as.numeric(x0) + e0[i,] - x 
  e0[i,] = e
  x0 = x
}


# dibujo proyección sobre dos bolas
x0 =c(0,2.5)
r = c(2,1)
centers = matrix(c(2,0,0,0),nrow = 2, ncol = 2,byrow = TRUE)
library(plotrix)
plot(0,2.5, xlim = c(-0.5,1.5), ylim = c(0.5,2.5), asp = 1,
     xlab = expression('x'[1]),
     ylab = expression('x'[2]))
draw.circle(0,0,1)
draw.circle(2,0,2)
# Bola unidad
#plot(0,0,cex=0.1,xlim=c(0.5,3),ylim=c(-1,0),xlab="",ylab="",asp=1)
r1 <- 1
xbol <- seq(-r1,r1,0.001)
ybol <- sqrt(r1^2-xbol^2)
polygon(c(xbol,-xbol),c(ybol,-ybol),cex=0.1)
# Bola y semiespacio
r2 <- 2
xbolH <- seq(0,2,0.001)
ybolH <- -sqrt(r2^2-xbolH^2)
polygon(xbolH,ybolH,col="gray")
M = length(r) ## nº de bolas(subconjuntos convexos = restricciones)
N = length(x0) ## dimension del vecotr a proyectar = nº de variables
E = matrix(rep(0,N*M),nrow = M, ncol = N)
for (n in 1:6) {
  i = 1+(n-1)%%M
  x0e = as.numeric(x0) + E[i,]
  lines(c(x0[1],x0e[1]),c(x0[2],x0e[2]), col='blue')
  x0e_prima = (x0e - centers[i,])/r[i]# reescalamos 
  x = P_ball(x0e_prima,W) # proyectamos sobre la "nueva" bola unidad (en el nuevo sist.)
  x = x * r[i] + centers[i,]# recuperamos el sistema de coordenadas anterior
  lines(c(x0e[1],x[1]),c(x0e[2],x[2]), col='red')
  #lines(c(x0[1],x[1]),c(x0[2],x[2]), col='red')
  e = as.numeric(x0) + E[i,] - x
  E[i,] = e
  x0 = x
  #cat(x ,n, "\n")
  print(e)
  #print(error)
}

######################### pa la función
dykstra_plot = function(x0,A,b,W = diag(length(x0))){
  x0 = as.numeric(x0)
  plot(x0[1],x0[2], asp = 1, xlim = c(-20,20), ylim = c(-20,20),
       xlab = expression('x'[1]),
       ylab = expression('x'[2]))#,xlim = c(-0.6,0.5), ylim = c(-1.3,0.3)
  #axis(4)
  for (i in 1:length(b)) {
    abline(b[i]/A[i,][2],A[i,][1]/A[i,][2])
  }
  #abline(b[1]/A[1,][2],A[1,][1]/A[1,][2])
  #abline(b[2]/A[2,][2],A[2,][1]/A[2,][2])
  
  #abline(bp[1]/a1p[2],a1p[1]/a1p[2])
  #abline(bp[2]/a2p[2],a2p[1]/a2p[2])
  M = length(as.vector(b)) ## nº de subconjuntos convexos (restricciones)
  N = length(x0) ## dimension del vecotr a proyectar = nº de variables
  e0 = matrix(rep(0,N*M),nrow = M, ncol = N)
  #print(x0)
  for (n in 1:10) {
    i = 1+(n-1)%%M
    x0e0 = as.numeric(x0) + e0[i,] # vector a proyectar, el cambio respecto a los subespacios
    x = P_Half_space(x0e0,as.matrix(A)[i,],as.vector(b)[i],W)
    lines(c(x0[1],x[1]),c(x0[2],x[2]), col='red')
    e = as.numeric(x0) + e0[i,] - x 
    e0[i,] = e
    x0 = x
    #print(x)
  }
}

dykstra_plot(DF[2,],A,b)
