# =============================================================================
# Fase 1 - Construcción de datos funcionales (ítem 1 de la rúbrica)
# Proyecto final de Analítica de Datos: lluvia del Valle del Cauca, 1981-2003
#
# Qué hace este script:
#   1. Lee los datos limpios del preprocesamiento (estaciones diarias y satélite por péntadas).
#   2. Pone las dos fuentes sobre el MISMO dominio anual: t en días, de 0 a 365.
#   3. Construye la base B-spline A MANO (recursión de Cox-de Boor) y la matriz Phi.
#   4. Calcula los coeficientes por mínimos cuadrados: c = (Phi'Phi)^-1 Phi' y.
#   5. Verifica la matriz suavizadora S = Phi (Phi'Phi)^-1 Phi': simétrica, idempotente, tr(S) = K.
#   6. Elige K con ECM, VC y VCG, para cada fuente, sin usar los años de prueba.
#   7. Calcula sigma^2, la covarianza residual y las bandas de confianza PUNTO A PUNTO.
#   8. Guarda los coeficientes de todas las curvas para las fases 2, 3 y 4.
#
# Las gráficas se imprimen en pantalla (panel Plots de RStudio) y no se guardan como archivos.
#
# La librería fda se usa SOLO para verificar que lo programado a mano da lo mismo.
#
# Decisiones:
#   - Escala: raíz cuadrada de la lluvia en mm/día (sqrt(mm/día)), como pide el modelo de la fase 4.
#   - Años bisiestos: la lluvia del 29 de febrero se suma al 28 de febrero. Todos los años
#     quedan con 365 días y se conserva el total anual.
#   - Estaciones: cada día se ubica en su centro (día 1 -> t = 0,5; día 365 -> t = 364,5).
#   - Satélite: cada péntada se ubica en el centro de sus días (calendario no bisiesto).
#   - Base: B-splines cúbicos (orden 4) con nudos equiespaciados en [0, 365].
#   - Partición recomendada por el profesor Ferro (asesoría del 20 de noviembre): 1981-1996 para
#     entrenar y validar con validación cruzada temporal, y 1997-2003 para la prueba (incluye El Niño
#     1997, La Niña 1998-2000 y años neutros 2001-2003). K se elige con 1981-1996, sin la prueba.
#     Las etiquetas de salida siguen a la fase 4: entrenamiento 1981-1993, validación 1994-1996 y
#     prueba 1997-2003.
#
# Uso (desde la carpeta del repositorio):
#   source("fase1_construccion_funcional.R")       # en RStudio
#   Rscript fase1_construccion_funcional.R
#
# Paquetes: fda (solo para verificar)
# =============================================================================

rm(list = ls())
library(fda)

ENT <- "./salidas_preprocesamiento"    # datos limpios de la fase 0 (Hoyos)
OUT <- "./salidas_fase1"               # salidas de esta fase
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ENTRENA <- 1981:1996                   # años fuera de la prueba (1997-2003): se usan para elegir K
ORDEN   <- 4                           # B-splines cúbicos
RANGO   <- c(0, 365)                   # dominio anual común



#####################################################################################
# 1. Carga de datos
#####################################################################################

# ---- Estaciones: lluvia diaria limpia (incluye días imputados y desagregados)
d <- read.csv(file.path(ENT, "estaciones_diario_limpio.csv"))
d <- d[order(d$STATION_ID, d$YEAR, d$MONTH, d$DAY), ]

# 29 de febrero -> se suma al 28 de febrero de la misma estación y año
es29    <- d$MONTH == 2 & d$DAY == 29
clave   <- paste(d$STATION_ID, d$YEAR, d$MONTH, d$DAY)
pos28   <- match(paste(d$STATION_ID[es29], d$YEAR[es29], 2, 28), clave)
stopifnot(!anyNA(pos28))
d$VALUE[pos28] <- d$VALUE[pos28] + d$VALUE[es29]
d <- d[!es29, ]

# Claves de las curvas: una curva = una estación en un año
claves <- unique(d[, c("STATION_ID", "YEAR")])
rownames(claves) <- NULL
stopifnot(nrow(d) == 365 * nrow(claves))         # todos los años con 365 días

# Matriz Y: filas = días (365), columnas = curvas estación-año. Escala sqrt(mm/día)
Y <- matrix(sqrt(d$VALUE), nrow = 365)
N <- ncol(Y)
cat("Curvas estación-año:", N, "| estaciones:", length(unique(claves$STATION_ID)), "\n")

# ---- Satélite: lluvia por péntada en el píxel de cada estación
p <- read.csv(file.path(ENT, "estaciones_pentadas_limpio.csv"))
p <- p[order(p$STATION_ID, p$YEAR, p$PENTADA), ]
claves_sat <- unique(p[, c("STATION_ID", "YEAR")])
rownames(claves_sat) <- NULL
stopifnot(identical(claves$STATION_ID, claves_sat$STATION_ID),
          identical(claves$YEAR, claves_sat$YEAR))  # misma curva en la misma columna
stopifnot(nrow(p) == 72 * N, !anyNA(p$satelite_mm), all(p$dias %in% 3:6))

#Matriz X: filas = péntadas (72), columnas = las mismas curvas de Y
#El satélite trae el total de la péntada; se divide por sus días (3 a 6) para que quede en
#mm/día, la misma escala de las estaciones. Si no, la sexta péntada de cada mes (hasta 6 días)
#sale más alta solo por tener más días, y estación - satélite no se podría comparar
X <- matrix(sqrt(p$satelite_mm / p$dias), nrow = 72)

# Años fuera de la prueba
tr <- claves$YEAR %in% ENTRENA
cat("Curvas para elegir K (1981-1996):", sum(tr), "\n")


#####################################################################################
# 2. Dominio común: t en días [0, 365]
#####################################################################################

# Estaciones: centro de cada día
t_est <- (1:365) - 0.5

# Satélite: centro de cada péntada (días 1-5, 6-10, 11-15, 16-20, 21-25 y 26-fin de mes)
dias_mes <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
inicio   <- cumsum(c(0, dias_mes[-12]))
t_sat <- unlist(lapply(1:12, function(m) sapply(1:6, function(k) {
  a <- 5 * (k - 1) + 1
  b <- if (k < 6) 5 * k else dias_mes[m]
  inicio[m] + (a + b) / 2 - 0.5
})))
stopifnot(length(t_sat) == 72)

n_est <- length(t_est)   # 365 puntos por curva
n_sat <- length(t_sat)   # 72 puntos por curva

# Malla fina para dibujar las curvas
t_fino <- seq(RANGO[1], RANGO[2], length.out = 731)


#####################################################################################
# 3. Base B-spline programada a mano (recursión de Cox-de Boor)
#####################################################################################
# Nudos: los extremos se repiten ORDEN veces y los nudos interiores se reparten igual.
# Orden 1: "escalones" que valen 1 en su intervalo y 0 fuera.
# Orden m: cada función se arma combinando dos funciones de orden m-1:
#   B_{j,m}(t) = (t - u_j)/(u_{j+m-1} - u_j) B_{j,m-1}(t) + (u_{j+m} - t)/(u_{j+m} - u_{j+1}) B_{j+1,m-1}(t)

nudos_bspline <- function(K, orden = ORDEN, rango = RANGO) {
  internos <- seq(rango[1], rango[2], length.out = K - orden + 2)
  internos <- internos[-c(1, length(internos))]
  c(rep(rango[1], orden), internos, rep(rango[2], orden))
}

base_bspline <- function(t, K, orden = ORDEN, rango = RANGO) {
  u <- nudos_bspline(K, orden, rango)
  nt <- length(t)
  # Orden 1
  B <- matrix(0, nt, length(u) - 1)
  for (j in 1:(length(u) - 1)) B[, j] <- as.numeric(t >= u[j] & t < u[j + 1])
  B[t == rango[2], max(which(u < rango[2]))] <- 1      # el extremo derecho entra en el último tramo
  # Órdenes 2, ..., orden
  if (orden > 1) for (m in 2:orden) {
    Bn <- matrix(0, nt, length(u) - m)
    for (j in 1:(length(u) - m)) {
      d1 <- u[j + m - 1] - u[j]
      d2 <- u[j + m] - u[j + 1]
      izq <- if (d1 > 0) (t - u[j]) / d1 * B[, j] else 0
      der <- if (d2 > 0) (u[j + m] - t) / d2 * B[, j + 1] else 0
      Bn[, j] <- izq + der
    }
    B <- Bn
  }
  B                                                    # nt x K
}

# Ejemplo: base con K = 8 en la malla fina
Phi_ej <- base_bspline(t_fino, 8)
stopifnot(max(abs(rowSums(Phi_ej) - 1)) < 1e-12)       # las piezas suman 1 en todo t

matplot(t_fino, Phi_ej, type = "l", lty = 1, lwd = 2, col = rainbow(8),
        xlab = "Día del año", ylab = "Valor de la función base",
        main = "Base B-spline cúbica con K = 8 (programada con Cox-de Boor)")
abline(v = nudos_bspline(8), lty = 3, col = "grey50")
par(mfrow = c(1, 1))


#####################################################################################
# 4. Mínimos cuadrados paso a paso y matriz suavizadora S
#####################################################################################
# c_hat = (Phi'Phi)^-1 Phi' y      (una columna de coeficientes por curva)
# S     = Phi (Phi'Phi)^-1 Phi'    (y_hat = S y)

coeficientes <- function(Phi, Ymat) solve(crossprod(Phi), crossprod(Phi, Ymat))
matriz_S     <- function(Phi) Phi %*% solve(crossprod(Phi), t(Phi))

propiedades_S <- function(S, K) data.frame(
  K = K,
  simetria_max_dif     = max(abs(S - t(S))),          # S = S'
  idempotencia_max_dif = max(abs(S %*% S - S)),       # S S = S
  traza                = sum(diag(S)),                # tr(S) = K
  traza_menos_K        = sum(diag(S)) - K)


#####################################################################################
# 5. Selección de K con ECM, VC y VCG (solo 1981-1996, sin la prueba)
#####################################################################################
# ECM(K) = promedio de (y - y_hat)^2                    -> siempre baja con K; no sirve solo
# VC(K)  = promedio de ((y - y_hat) / (1 - s_ii))^2     -> atajo de dejar-uno-fuera con la diagonal de S
# VCG(K) = promedio sobre curvas de n * SSE / (n - K)^2 -> validación cruzada generalizada

criterios_K <- function(t, Ymat, Ks) {
  n <- length(t)
  res <- NULL
  for (K in Ks) {
    Phi  <- base_bspline(t, K)
    S    <- matriz_S(Phi)
    E    <- Ymat - S %*% Ymat
    s_ii <- diag(S)
    res <- rbind(res, data.frame(
      K   = K,
      ECM = mean(E^2),
      VC  = mean((E / (1 - s_ii))^2),
      VCG = mean(n * colSums(E^2) / (n - K)^2)))
  }
  res
}

# Estaciones: curvas de entrenamiento
Ks_est    <- seq(6, 60, by = 2)
tabla_est <- criterios_K(t_est, Y[, tr], Ks_est)

# Satélite: varias estaciones comparten píxel; para elegir K se usa cada píxel-año una sola vez
X_tr_unicas <- t(unique(t(X[, tr])))
cat("Curvas de satélite (píxel-año) de entrenamiento:", ncol(X_tr_unicas), "\n")
Ks_sat    <- seq(6, 40, by = 2)
tabla_sat <- criterios_K(t_sat, X_tr_unicas, Ks_sat)

K_est <- tabla_est$K[which.min(tabla_est$VCG)]
K_sat <- tabla_sat$K[which.min(tabla_sat$VCG)]
cat("\nK elegido por VCG -> estaciones:", K_est, "| satélite:", K_sat, "\n")
cat("K elegido por VC  -> estaciones:", tabla_est$K[which.min(tabla_est$VC)],
    "| satélite:", tabla_sat$K[which.min(tabla_sat$VC)], "\n")
print(round(tabla_est, 4)); print(round(tabla_sat, 4))

grafica_K <- function(tabla, K_elegido, fuente) {
  par(mfrow = c(1, 3))
  plot(tabla$K, tabla$ECM, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "ECM",
       main = paste("ECM de ajuste -", fuente))
  grid(5, 5, lwd = 2)
  plot(tabla$K, tabla$VC, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VC",
       main = paste("Validación cruzada -", fuente))
  abline(v = tabla$K[which.min(tabla$VC)], lty = 3, col = "red")
  grid(5, 5, lwd = 2)
  plot(tabla$K, tabla$VCG, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VCG",
       main = paste("VCG -", fuente))
  abline(v = K_elegido, lty = 3, col = "red")
  abline(h = min(tabla$VCG), lty = 2, col = "red")
  grid(5, 5, lwd = 2)
}
grafica_K(tabla_est, K_est, "estaciones")
par(mfrow = c(1, 1))
grafica_K(tabla_sat, K_sat, "satélite")
par(mfrow = c(1, 1))
# Nota: el ECM baja siempre al subir K (sobreajuste). La VC y la VCG penalizan la complejidad
# y tienen un mínimo interior. Con datos diarios la lluvia es muy intermitente: la VCG prefiere
# curvas suaves (K pequeño) y la curva conserva los dos picos del régimen bimodal (ver sección 8).


#####################################################################################
# 6. Ajuste final con el K elegido y verificación contra la librería fda
#####################################################################################

Phi_est <- base_bspline(t_est, K_est)              # 365 x K_est
Phi_sat <- base_bspline(t_sat, K_sat)              # 72  x K_sat
C_est   <- coeficientes(Phi_est, Y)                # K_est x N
C_sat   <- coeficientes(Phi_sat, X)                # K_sat x N
S_est   <- matriz_S(Phi_est)
S_sat   <- matriz_S(Phi_sat)

# Propiedades de S (deben dar ~0, ~0 y K)
verif_S <- rbind(cbind(fuente = "estaciones", propiedades_S(S_est, K_est)),
                 cbind(fuente = "satélite",   propiedades_S(S_sat, K_sat)))
print(verif_S)

# Verificación con fda: misma base y mismos coeficientes (diferencias del orden de 1e-12)
base_fda_est <- create.bspline.basis(RANGO, nbasis = K_est, norder = ORDEN)
base_fda_sat <- create.bspline.basis(RANGO, nbasis = K_sat, norder = ORDEN)
muestra <- 1:50
verif_fda <- data.frame(
  fuente = c("estaciones", "satélite"),
  max_dif_Phi = c(max(abs(Phi_est - eval.basis(t_est, base_fda_est))),
                  max(abs(Phi_sat - eval.basis(t_sat, base_fda_sat)))),
  max_dif_coeficientes = c(
    max(abs(C_est[, muestra] - smooth.basis(t_est, Y[, muestra], base_fda_est)$fd$coefs)),
    max(abs(C_sat[, muestra] - smooth.basis(t_sat, X[, muestra], base_fda_sat)$fd$coefs))))
print(verif_fda)

# Las dos fuentes evaluadas en la misma malla fina (mismo dominio anual)
Phi_fino_est <- base_bspline(t_fino, K_est)
Phi_fino_sat <- base_bspline(t_fino, K_sat)
curvas_est <- Phi_fino_est %*% C_est               # 731 x N
curvas_sat <- Phi_fino_sat %*% C_sat

par(mfrow = c(1, 2))
set.seed(2026)
m <- sample(which(tr), 200)
matplot(t_fino, curvas_est[, m], type = "l", lty = 1, col = adjustcolor("steelblue", 0.25),
        xlab = "Día del año", ylab = "sqrt(mm/día)",
        main = paste0("Estaciones: 200 curvas (K = ", K_est, ")"))
lines(t_fino, rowMeans(curvas_est[, tr]), lwd = 3, col = "navy")
matplot(t_fino, curvas_sat[, m], type = "l", lty = 1, col = adjustcolor("darkorange", 0.25),
        xlab = "Día del año", ylab = "sqrt(mm/día)",
        main = paste0("Satélite: las mismas 200 curvas (K = ", K_sat, ")"))
lines(t_fino, rowMeans(curvas_sat[, tr]), lwd = 3, col = "darkred")
par(mfrow = c(1, 1))


#####################################################################################
# 7. Varianza residual, covarianza residual y bandas de confianza punto a punto
#####################################################################################
# Bajo Var(y) = sigma^2 I:
#   sigma^2 (insesgado)     = SSE / (n - K)            (un valor por curva)
#   Var(c_hat)              = sigma^2 (Phi'Phi)^-1
#   Var(x_hat(t))           = sigma^2 phi(t)' (Phi'Phi)^-1 phi(t)
#   Banda punto a punto     = x_hat(t) +/- 2 EE(t)
# Es PUNTO A PUNTO: en cada t, por separado, cubre el valor de la curva con ~95 %.
# No es una banda simultánea para toda la curva a la vez.

sigma2 <- function(Ymat, S, K) colSums((Ymat - S %*% Ymat)^2) / (nrow(Ymat) - K)
sig2_est <- sigma2(Y, S_est, K_est)
sig2_sat <- sigma2(X, S_sat, K_sat)
cat("\nsigma^2 mediana -> estaciones:", round(median(sig2_est), 3),
    "| satélite:", round(median(sig2_sat), 3), "\n")

# h(t) = phi(t)' (Phi'Phi)^-1 phi(t) en la malla fina
h_fino <- function(Phi_fino, Phi) rowSums((Phi_fino %*% solve(crossprod(Phi))) * Phi_fino)
h_est  <- h_fino(Phi_fino_est, Phi_est)
h_sat  <- h_fino(Phi_fino_sat, Phi_sat)

# Covarianza residual: ¿se cumple Var(y) = sigma^2 I?  Si los residuos de días vecinos
# están poco correlacionados, el supuesto es razonable y las bandas son válidas.
E_est   <- (Y - S_est %*% Y)[, tr]
Cor_res <- cor(t(E_est))                                   # 365 x 365
cor_lag1 <- mean(diag(Cor_res[-1, -365]))
cat("Correlación media entre residuos de días consecutivos (estaciones):", round(cor_lag1, 3), "\n")

image(t_est, t_est, Cor_res, col = hcl.colors(50, "Blue-Red 3"), zlim = c(-1, 1),
      xlab = "Día s", ylab = "Día t", main = "Correlación residual entre días (estaciones, 1981-1996)")
par(mfrow = c(1, 1))

# Ejemplo de bandas: una estación-año de entrenamiento, ambas fuentes
j <- which(tr)[1]
grafica_banda <- function(t_obs, y_obs, curva, sig2, h, titulo, unidad) {
  ee <- sqrt(sig2 * h)
  plot(t_obs, y_obs, pch = 19, cex = 0.5, col = "grey40", xlab = "Día del año", ylab = unidad,
       ylim = range(y_obs, curva + 2 * ee, curva - 2 * ee), main = titulo)
  lines(t_fino, curva, lwd = 2, col = "blue")
  lines(t_fino, curva + 2 * ee, lty = 2, col = "blue")
  lines(t_fino, curva - 2 * ee, lty = 2, col = "blue")
  legend("topleft", c("Observado", "Curva estimada", "Banda +/- 2 EE (punto a punto)"),
         pch = c(19, NA, NA), lty = c(NA, 1, 2), col = c("grey40", "blue", "blue"), bty = "n", cex = 0.8)
}
par(mfrow = c(1, 2))
grafica_banda(t_est, Y[, j], curvas_est[, j], sig2_est[j], h_est,
              paste("Estación", claves$STATION_ID[j], "-", claves$YEAR[j]), "sqrt(mm/día)")
grafica_banda(t_sat, X[, j], curvas_sat[, j], sig2_sat[j], h_sat,
              paste("Satélite en su píxel -", claves$YEAR[j]), "sqrt(mm/día)")
par(mfrow = c(1, 1))


#####################################################################################
# 8. Control: ¿la curva media conserva el régimen bimodal?
#####################################################################################
media_est <- rowMeans(curvas_est[, tr])
picos <- which(diff(sign(diff(media_est))) == -2) + 1
cat("Picos de la curva media de estaciones (día del año):", round(t_fino[picos]), "\n")


#####################################################################################
# 9. Salidas para las fases 2, 3 y 4
#####################################################################################
coef_est <- data.frame(claves, t(C_est));  names(coef_est)[-(1:2)] <- paste0("c", 1:K_est)
coef_sat <- data.frame(claves, t(C_sat));  names(coef_sat)[-(1:2)] <- paste0("c", 1:K_sat)
coef_est$particion <- coef_sat$particion <- ifelse(claves$YEAR <= 1993, "entrenamiento",
                                             ifelse(claves$YEAR <= 1996, "validacion", "prueba"))
write.csv(coef_est, file.path(OUT, "coeficientes_estaciones.csv"), row.names = FALSE)
write.csv(coef_sat, file.path(OUT, "coeficientes_satelite.csv"),   row.names = FALSE)

write.csv(tabla_est, file.path(OUT, "seleccion_K_estaciones.csv"), row.names = FALSE)
write.csv(tabla_sat, file.path(OUT, "seleccion_K_satelite.csv"),   row.names = FALSE)
write.csv(verif_S,   file.path(OUT, "verificacion_matriz_S.csv"),  row.names = FALSE)
write.csv(verif_fda, file.path(OUT, "verificacion_vs_fda.csv"),    row.names = FALSE)

# Cómo reconstruir la base en las otras fases (mismos nudos, mismo orden)
base_info <- data.frame(fuente = c("estaciones", "satélite"), K = c(K_est, K_sat), orden = ORDEN,
                        rango = "0-365", escala = "sqrt(mm/día)",
                        nudos = c(paste(round(nudos_bspline(K_est), 4), collapse = " "),
                                  paste(round(nudos_bspline(K_sat), 4), collapse = " ")))
write.csv(base_info, file.path(OUT, "base_info.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_fase1.txt"))
cat("\nListo. Salidas en", OUT, "\n")
