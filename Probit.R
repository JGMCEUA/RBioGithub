#Programa que realiza un análisis Probit para estimar la DL50 y otros parámetros.
# José Gustavo Marín COntreras
#Cargar librerias
library(data.table)
#Cargar carpeta
setwd("C:/Users/gusta/Desktop/CursosR ladys/Pp/2019")
x <- read.csv("DatosProbit.csv")

#Cargar datos a analizar (Dosis vs Supervivinetes)
dosis <- c(0, 5, 10, 20, 40, 80, 160, 320, 640, 1280)
supervivientes <- c(10, 10, 9, 8, 6, 4, 2, 1, 0, 0)
#Número total de insectos por dosis, ojo, es importante ingresar los datos con la mortalidad ajustada en caso de que se requiera
Tt <- 10

#Tabla con datos para analiszar por probit
x <- data.frame(Dosis = dosis, Supervivientes = supervivientes)

Probit_JGMC <- function(x,Tt){
  #Función que realiza el analisis probit de los datos de dosis vs insectos supervivientes.
  #Es requerido conocer los datos de "Dosis" aplicada e insectos "Supervivientes" [x] y Total [Tt] de insectos por cada dosis
  #Se obtiene el valor Z de la distribución normal de acuerdo a las probabilidades observadas de mortalidad de cada dosis
  Tempdata <- cbind(x,Tp = qnorm((Tt-x$Supervivientes)/Tt))
  #Se elimina los datos que puedan generar errores
  rl.data <- Tempdata[Tempdata$Supervivientes<Tt & Tempdata$Supervivientes>0 ,c("Dosis","Tp")]
  #Modelo estádistico para graficar el resultado de los valores Z [Tp] vs Dosis [Dosis]
  rl.mod.datax <- lm(Tp~ Dosis, data = rl.data)
  #Modelo estadístico para obtener el DL50 de los datos ingresados
  rl.mod.datay <- lm(Dosis~ Tp, data = rl.data)
  
  #Graficación de los ressultados
  temp <-  predict(rl.mod.datax, interval="prediction")
  #my.pdf("8-SobXedad61-80B") # Como antes pero con el modelo y lim. conf
  plot(rl.data$Dosis, rl.data$Tp, col="red", pch=4, cex=1.5, xlab="Dosis aplicada (ug/mL)", ylab="Distribucion estandar acumulada (CDF)", main="Analisis Probit", sub="Estimación DL50")
  grid()
  # Vamos a incluir el modelo lineal
  abline(coef=rl.mod.datax$coefficients, lwd=2, col="blue")
  # Límite de confianza inferior:
  points(rl.data$Dosis,temp[,2] , type="l", lwd=2, lty=2, col="blue")
  # Límite de confianza superior:
  points(rl.data$Dosis,temp[,3], type="l", lwd=2, lty=2, col="blue")
  # Repintamos los datos como círculos:
  points(rl.data$Dosis, rl.data$Tp, col="red", pch=1, cex=1)
  #legend("topright", legend=c("Modelo: Y = 1.3255 -0.0062 x Edad", "Explica el 99.8% de la variación", "Líneas punteadas: LC 95%"))
  #dev.off()
  
  #Imprimirintervalos de confianza de la DL50
  DL50 <- data.frame(Tp = c(qnorm(.5),qnorm(.75),qnorm(.90),qnorm(.99)))
  temp <-  predict(rl.mod.datay, newdata = DL50 ,interval="confidence")
  cat("El valor DL50 es:",rl.mod.datay$coefficients[1],"ug/mL. \n Intervalos de confianza(95%): [",temp[1,2],",",temp[1,3],"] ug/mL. \n")
  print(summary(rl.mod.datax))
  return(rl.mod.datax)
  }
Res <- Probit_JGMC(x,Tt)



