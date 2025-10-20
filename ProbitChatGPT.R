# datos
dosis <- c(0,5,10,20,40,80,160,320,640,1280)
n     <- rep(10, length(dosis))
muertos <- c(0,0,1,2,4,6,8,9,10,10)

# preparar datos: evitar log(0) -> sustituimos 0 por 0.1 (valor pequeño)
dosis_for_log <- ifelse(dosis == 0, 0.1, dosis)
logd <- log(dosis_for_log)

df <- data.frame(dosis, logd, muertos, vivos = n - muertos)

# ajustar modelo probit con binomial (respuesta en cbind(muertos, vivos))
mod <- glm(cbind(muertos, vivos) ~ logd, family = binomial(link = "probit"), data = df)

summary(mod)

# Estimar DL50: para probit, DL50 ocurre cuando la función lineal = 0.
# Si el modelo es: eta = a + b*logd, para p=0.5 => eta = 0 => log(DL50) = -a/b
coef <- coef(mod)
a <- coef[1]; b <- coef[2]
logDL50 <- -a / b
DL50 <- exp(logDL50)
DL50

# Intervalo de confianza (aproximación delta en escala logarítmica)
vcovm <- vcov(mod)
var_a <- vcovm[1,1]
var_b <- vcovm[2,2]
cov_ab <- vcovm[1,2]

var_logDL50 <- var_a / (b^2) + (a^2 * var_b) / (b^4) - 2 * a * cov_ab / (b^3)
se_logDL50 <- sqrt(var_logDL50)

# IC 95% en escala de dosis
DL50_low  <- DL50 * exp(-1.96 * se_logDL50)
DL50_high <- DL50 * exp( 1.96 * se_logDL50)
c(DL50 = DL50, CI95 = c(DL50_low, DL50_high))


xp(logDL50)

# Graficar
plot(x$logd, x$Muertos/Tt, pch=19, col="blue", 
     xlab="Dosis (mg/L)", ylab="Mortalidad", log="x",
     ylim=c(0,1), main="Curva Probit - Estimación DL50")
lines(new_dosis, pred, col="red", lwd=2)
abline(h=0.5, col="darkgreen", lty=2)       # línea horizontal en 50%
abline(v=DL50, col="purple", lty=2)         # línea vertical en DL50
legend("bottomright", legend = paste("DL50 =", round(DL50,2), "mg/L"),
       col="purple", lty=2, bty="n")

