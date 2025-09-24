########################################################################################
# Este script permite seleccionar aleatoriamente 100 de filas de la base de datos ######
# Las variables que son montos de dinero se van a construir en miles de pesos ##########
# Variable respuesta: Ingresos Totales del Hogar (ING_TOTAL) ###########################
# Variables regresoras: Total ingresos sueldos y salarios (ING_T_D), ###################
####################### Total ingresos trabajo independiente (ING_T_I), ################
####################### Total ingresos por jubilaciones (ING_JUB), #####################
####################### Total ingresos por pensiones (ING_PEN), ########################
# Funcion link: identity ###############################################################
########################################################################################

## Para cargar los datos y definir parámetros de la función bsreg.fit()
source("/home/helton/Dropbox/quantile_reg_bs/Codigos/BSreg.fit.R")

library(haven)
HCB_SinCerosING_TOTAL <- read_sav("/home/helton/Dropbox/quantile_reg_bs/Codigos/HCB_SinCerosING_TOTAL.sav")

set.seed(20)

i <- sample(1:nrow(HCB_SinCerosING_TOTAL),size=100)

Banco_Novo <- HCB_SinCerosING_TOTAL[i,]

head(Banco_Novo)

write.csv(Banco_Novo, "/home/helton/Dropbox/quantile_reg_bs/Codigos/datos.csv")


# Se define la matriz de diseño con las covariables ya acordadas
# Se divide por 1000 para tener números más tratables
COL1 <- matrix(rep(1,33501),33501,1)
matrixX.aux <- data.matrix(cbind(COL1,
                                 HCB_SinCerosING_TOTAL[,57]/1000,
                                 HCB_SinCerosING_TOTAL[,63]/1000,
                                 HCB_SinCerosING_TOTAL[,81]/1000,
                                 HCB_SinCerosING_TOTAL[,87]/1000))


# Se define el vector respuesta: variable ING_TOTAL
vectorT.aux <- data.matrix(HCB_SinCerosING_TOTAL[,101]/1000)

# Lo que sigue es para extraer filas aleatoriamente


S       = sample(1:33501, 100) # Genera 100 números aleatorios de 1 a 33501
matrixX = matrix(rep(0,500),100,5) # Genera una matriz nula de 100 filas y 5 columnas
vectorT = matrix(rep(0,100),100,1) # Genera un vetor nulo con 100 elementos

for(i in 1:100){
  matrixX[i,] = matrixX.aux[S[i],]  # Extrae las filas de matrixX (con 33501 casos) de acuerdo a S
  vectorT[i,] = vectorT.aux[S[i],]  # Extrae los elementos de vectorT de acuerdo a S
}

colnames(matrixX) = colnames(cbind(COL1, 
                                   HCB_SinCerosING_TOTAL[,57],
                                   HCB_SinCerosING_TOTAL[,63],
                                   HCB_SinCerosING_TOTAL[,81],
                                   HCB_SinCerosING_TOTAL[,87]))

colnames(vectorT) = colnames(HCB_SinCerosING_TOTAL[,101])



#######################

dados <- read.csv("C:/Users/User/Documents/GitHub/MSc-est/MODELAGEM COM APOIO COMPUTACIONAL/article/referencias/datos.csv")

COL1 <- matrix(rep(1,100),100,1)

matrixX <- data.matrix(cbind(COL1,
                              dados[,57]/1000,
                              dados[,63]/1000,
                              dados[,81]/1000,
                              dados[,87]/1000))

vectorT <- data.matrix(dados[,101]/1000)


model = bsreg.fit(x=matrixX, y=vectorT, link="identity")

summary(model)
