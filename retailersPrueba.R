## Siguiendo el ejemplo que se muestra en el tema 10 de la referencia [jongeloo], probaremos la función
## para corregir los datos sobre el conjunto de datos Retailers.\\
## Este conjunto se encuentra en el paquete validate de R y consiste en un data frame con las siguientes
## columnas:
# 1ª: size
# 2ª: incl.prob
# 3ª: staff
# 4ª: turnover
# 5ª: other.rev
# 6ª: total.rev
# 7ª: staff.costs
# 8ª: total.costs
# 9ª: profit
# 10ª: bat
# Que están sometidos a las siguientes restricciones
#staff >= 0
# turnover >= 0
# other.rev >= 0
# turnover + other.rev == total.rev
# total.rev - total.costs == profit
# total.costs >= 0
### Prueba con Retailers
library(validate)##########
library(simputation)
data(retailers)

## Restricciones de la página 271(Jogen_loo), matricialmente serían:
a1r = -c(0,0,1,0,0,0,0,0,0,0)
a2r = -c(0,0,0,1,0,0,0,0,0,0)
a3r = -c(0,0,0,0,1,0,0,0,0,0)
a4r = -c(0,0,0,1,1,-1,0,0,0,0)
a5r = -c(0,0,0,0,0,1,0,-1,-1,0)
a6r = -c(0,0,0,0,0,0,0,1,0,0)
Ar = matrix(c(a1r, a2r, a3r, a4r, a5r, a6r), nrow = 6, ncol = length(a1r), byrow = TRUE)
br = rep(0,6)

# imputación sobre el DF original con la función impute_mf que usa vanderloo
#library(simputation)
retailers_imputed = impute_mf(retailers, . ~ .)
eq = c(rep('<=',3), '==', '==', '<=')
## Aplicamos ahora el algoritmo de Dykstra
#la función tiene que recibir un dataframe con valores numéricos
# en este caso la primera columna no nos interesa 
# 1º cambiando todos los datos
#Wr = diag(length(df_imputRetail[1,][-1])) 
# 2ºcambiando solo los imputados
Ir = is.na(retailers)# DF de TRUE o FALSE indicando si el dato ha sido imputado(TRUE) o no 
retailers_lineal = dykstra_DF(DF = retailers_imputed[,-1],
                           typeConstraints = "lineal",
                           A = Ar[,-1], 
                           b = br, 
                           eq = eq)
#> head(retailers_lineal$df_final,5)
#incl.prob staff turnover  other.rev total.rev staff.costs total.costs    profit     vat
#1      0.02 75.00 1985.132 15079.1687 17064.301     12021.3    7967.150 9097.1504 3486.63
#2      0.14  9.00 1589.353    29.4125  1618.765       131.0    1549.882   68.8825 1349.95
#3      0.14 22.64 6905.800     0.0000  6905.800       324.0    6486.400  419.4000 5802.69
#4      0.14 18.87 3861.000    13.0000  3874.000       290.0    3600.000  274.0000 3540.53
#5      0.14 20.73 5412.450   128.5300  5540.980       314.0    5499.490   41.4900 4371.47
retailers_linealI = dykstra_DF(DF = retailers_imputed[,-1],
                            typeConstraints = "lineal",
                            A = Ar[,-1], 
                            b = br, 
                            eq = eq,
                            I = Ir[,-1])
#> head(retailers_linealI$df_final,5)
#incl.prob staff turnover other.rev total.rev staff.costs total.costs profit     vat
#1      0.02 75.00        0  1130.183      1130     12021.3       18915  20045 3486.63
#2      0.14  9.00     1607     0.000      1607       131.0        1544     63 1349.95
#3      0.14 22.64     6886   -33.000      6919       324.0        6493    426 5802.69
#4      0.14 18.87     3861    13.000      3874       290.0        3600    274 3540.53
#5      0.14 20.73     5565    37.000      5602       314.0        5530     72 4371.47
# En este caso no todas las filas que no cumplen las restricciones han sido modificadas,
# ya que los índices que marcan que datos cambiar están dados según los datos faltantes
# del conjunto de datos original. En concreto las filas que no cumplen las restricciones
# pero no han sido cambiadas son:
#which(!(row.names(retail_linealI$df_check) %in% row.names(retail_linealI$df_proyec)))
#[1]  3 22 25 26 36

v <- validator(
  staff >= 0
  , turnover >= 0
  , other.rev >= 0
  , turnover + other.rev == total.rev
  , total.rev - total.costs == profit
  , total.costs >= 0
)
retailers_linealRspa <- match_restrictions(retailers_imputed, v)
#> head(retailers_lineal2,5)
#size incl.prob staff turnover    other.rev total.rev staff.costs total.costs     profit     vat
#1  sc0      0.02 75.00 1985.132 1.507917e+04 17064.301     12021.3    7967.150 9097.15048 3486.63
#2  sc3      0.14  9.00 1589.353 2.941280e+01  1618.765       131.0    1549.882   68.88240 1349.95
#3  sc3      0.14 22.64 6905.801 5.553284e-16  6905.799       324.0    6486.400  419.39950 5802.69
#4  sc3      0.14 18.87 3861.000 1.300000e+01  3874.000       290.0    3600.000  274.00000 3540.53
#5  sc3      0.14 20.73 5412.450 1.285298e+02  5540.980       314.0    5499.490   41.49006 4371.47
retailers_linealIRspa <- match_restrictions(retailers_imputed, v, adjust=Ir)
#> head(retailers_linealIRspa)
#size incl.prob staff turnover other.rev total.rev staff.costs total.costs profit     vat
#1  sc0      0.02 75.00      NaN       NaN      1130  12021.3000       18915  20045 3486.63
#2  sc3      0.14  9.00     1607         0      1607    131.0000        1544     63 1349.95
#3  sc3      0.14   NaN     6886       -33      6919    324.0000        6493    426 5802.69
#4  sc3      0.14 18.87     3861        13      3874    290.0000        3600    274 3540.53
#5  sc3      0.14 20.73     5565        37      5602    314.0000        5530     72 4371.47
#6  sc0      0.02  1.00       25         0        25    125.2533          22      3  187.61
#Warning messages:
 # 1: no convergence in record number 1: could not allocate enough memory
 # 2: no convergence in record number 3: could not allocate enough memory
####################################################

# Haciendo uso de la función replace_errors del paquete errorlocate
# para sustituir los datos que no cumplen as restricciones por NA antes de realizar
# la imputación

library(validate) 
library(simputation) 
library(errorlocate) 
library(rspa) 

# siguiendo la notación de Vanderloo
d1 <- replace_errors(retailers,v)
#miss <- is.na(d1)
# imputa los NA con un algoritmo de arbol de decisión
d2 <- impute_mf(d1, . ~ .)
d3 <- match_restrictions(d2, v)
d4 <- match_restrictions(d2, v, adjust=is.na(d1))
d_lineal = dykstra_DF(DF = d2[,-1],
                              typeConstraints = "lineal",
                              A = Ar[,-1], 
                              b = br, 
                              eq = eq)

d_linealI = dykstra_DF(DF = d2[,-1],
                               typeConstraints = "lineal",
                               A = Ar[,-1], 
                               b = br, 
                               eq = eq,
                               I = is.na(d1)[,-1])

# En este caso obtenemos los siguientes resultados:
#> head(d3)
#size incl.prob staff     turnover other.rev    total.rev staff.costs total.costs       profit     vat
#1  sc0      0.02 75.00 105168.14776 47.027757 105215.17483    23550.52 52042.58741 53172.587414 4532.71
#2  sc3      0.14  9.00   1603.22006  6.300064   1609.51996      131.00  1545.25998    64.259979 1346.68
#3  sc3      0.14 25.77   6879.42636 43.956361   6923.38243      324.00  6495.19121   428.191213 6291.12
#4  sc3      0.14 21.24   3861.00000 13.000000   3874.00000      290.00  3600.00000   274.000000 3817.87
#5  sc3      0.14 22.65   5508.61869 70.828686   5579.44754      314.00  5518.72377    60.723771 4763.30
#6  sc0      0.02  1.00     11.83772 21.937723     33.77485      435.32    26.38743     7.387426  150.63

#> head(d_lineal$df_final)
#incl.prob staff    turnover other.rev  total.rev staff.costs total.costs      profit     vat
#1      0.02 75.00 105168.1475  47.02750 105215.175    23550.52  52042.5875 53172.58750 4532.71
#2      0.14  9.00   1603.2200   6.30000   1609.520      131.00   1545.2600    64.26000 1346.68
#3      0.14 25.77   6879.4263  43.95625   6923.383      324.00   6495.1912   428.19125 6291.12
#4      0.14 21.24   3861.0000  13.00000   3874.000      290.00   3600.0000   274.00000 3817.87
#5      0.14 22.65   5508.6187  70.82875   5579.447      314.00   5518.7238    60.72375 4763.30
#6      0.02  1.00     11.8375  21.93750     33.775      435.32     26.3875     7.38750  150.63

#> head(d4)
#size incl.prob staff turnover     other.rev total.rev staff.costs total.costs profit     vat
#1  sc0      0.02 75.00 38960.01 -1.629876e-12     38960    23550.52       18915  20045 4532.71
#2  sc3      0.14  9.00  1607.00  0.000000e+00      1607      131.00        1544     63 1346.68
#3  sc3      0.14 25.77  6886.00  3.300000e+01      6919      324.00        6493    426 6291.12
#4  sc3      0.14 21.24  3861.00  1.300000e+01      3874      290.00        3600    274 3817.87
#5  sc3      0.14 22.65  5565.00  3.700000e+01      5602      314.00        5530     72 4763.30
#6  sc0      0.02  1.00    25.00  0.000000e+00        25      435.32          22      3  150.63

#> head(d_linealI$df_final)
#incl.prob staff turnover other.rev total.rev staff.costs total.costs profit     vat
#1      0.02 75.00 39123.75         0     38960    23550.52       18915  20045 4532.71
#2      0.14  9.00  1607.00         0      1607      131.00        1544     63 1346.68
#3      0.14 25.77  6886.00        33      6919      324.00        6493    426 6291.12
#4      0.14 21.24  3861.00        13      3874      290.00        3600    274 3817.87
#5      0.14 22.65  5565.00        37      5602      314.00        5530     72 4763.30
#6      0.02  1.00    25.00         0        25      435.32          22      3  150.63

#####################

# Prueba con una restricción de tipo bola con las variables profit y vat por ejemplo:
# (profit - mean(profit))^2 + (vat-mean(vat)) <= mean(colMeans(retailers_imputed[,c("vat","profit")])
#colMeans(retailers_imputed[,c("vat","profit")])
#vat   profit 
#2079.995 6939.221 

retailers_ball  = dykstra_DF(retailers_imputed[,c("vat","profit")],
                             typeConstraints = "ball",
                             r = mean(colMeans(retailers_imputed[,c("vat","profit")])), 
                             centers =  t(colMeans(retailers_imputed[,c("vat","profit")])))

retailers_ballI  = dykstra_DF(retailers_imputed[,c("vat","profit")],
                             typeConstraints = "ball",
                             r = mean(colMeans(retailers_imputed[,c("vat","profit")])), 
                             centers =  t(colMeans(retailers_imputed[,c("vat","profit")])),
                             I =Ir[,c("vat","profit")])

# o imponemos que 



### si pruebo con la métrica diagonal con la media de cada columna
Wr = diag(colMeans(retailers_imputed[,-1]))
retail_linealIW = dykstra_DF(DF = retailers_imputed[,-1],
                            typeConstraints = "lineal",
                            A = Ar[,-1], 
                            b = br, 
                            eq = eq,  
                            #I = Ir[,-1],
                            W = Wr)



