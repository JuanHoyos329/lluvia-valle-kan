#####################################################################################
# Fase 1 - Complemento: ¿por qué B-spline, mínimos cuadrados ordinarios y K = 8?
#
# Tres revisiones que respaldan las decisiones de fase1_construccion_funcional.R:
#   A. Base de Fourier frente a B-spline (VCG y condición de Phi'Phi)
#   B. Varianza de los residuos por día del año: ¿hace falta mínimos cuadrados ponderados?
#   C. Suavizado penalizado (segunda derivada) con lambda elegido por VCG
#
# Usa solo las curvas de 1981-1996 (las mismas con que se eligió K). No modifica nada de la
# fase 1 ni guarda archivos: imprime los números en la consola y muestra las figuras.
# Correr desde la carpeta del repositorio:  source("fase1_complemento_bases.R")
#####################################################################################

rm(list = ls())
ENT     <- "./salidas_preprocesamiento"
ENTRENA <- 1981:1996
RANGO   <- c(0, 365)


#####################################################################################
# 0. Datos (igual que en la fase 1)
#####################################################################################

d <- read.csv(file.path(ENT, "estaciones_diario_limpio.csv"))
d <- d[order(d$STATION_ID, d$YEAR, d$MONTH, d$DAY), ]
es29  <- d$MONTH == 2 & d$DAY == 29
clave <- paste(d$STATION_ID, d$YEAR, d$MONTH, d$DAY)
pos28 <- match(paste(d$STATION_ID[es29], d$YEAR[es29], 2, 28), clave)
d$VALUE[pos28] <- d$VALUE[pos28] + d$VALUE[es29]
d <- d[!es29, ]
claves <- unique(d[, c("STATION_ID", "YEAR")])
Y <- matrix(sqrt(d$VALUE), nrow = 365)

p <- read.csv(file.path(ENT, "estaciones_pentadas_limpio.csv"))
p <- p[order(p$STATION_ID, p$YEAR, p$PENTADA), ]
X <- matrix(sqrt(p$satelite_mm / p$dias), nrow = 72)

tr  <- claves$YEAR %in% ENTRENA
Ytr <- Y[, tr]                       # estaciones 1981-1996
Xtr <- t(unique(t(X[, tr])))         # satélite: cada píxel-año una sola vez
cat("Curvas de estaciones:", ncol(Ytr), "| curvas píxel-año del satélite:", ncol(Xtr), "\n")

t_est <- (1:365) - 0.5
dias_mes <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
inicio   <- cumsum(c(0, dias_mes[-12]))
t_sat <- unlist(lapply(1:12, function(m) sapply(1:6, function(k) {
  a <- 5 * (k - 1) + 1
  b <- if (k < 6) 5 * k else dias_mes[m]
  inicio[m] + (a + b) / 2 - 0.5
})))


#####################################################################################
# 1. Bases programadas a mano
#####################################################################################

# B-spline cúbica (Cox-de Boor), la misma de la fase 1
nudos_bspline <- function(K, orden = 4, rango = RANGO) {
  internos <- seq(rango[1], rango[2], length.out = K - orden + 2)
  internos <- internos[-c(1, length(internos))]
  c(rep(rango[1], orden), internos, rep(rango[2], orden))
}
base_bspline <- function(t, K, orden = 4, rango = RANGO) {
  u <- nudos_bspline(K, orden, rango)
  nt <- length(t)
  B <- matrix(0, nt, length(u) - 1)
  for (j in 1:(length(u) - 1)) B[, j] <- as.numeric(t >= u[j] & t < u[j + 1])
  B[t == rango[2], max(which(u < rango[2]))] <- 1
  for (m in 2:orden) {
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
  B
}

# Fourier: 1, sen(wt), cos(wt), sen(2wt), cos(2wt), ...  con w = 2*pi/365.  K es impar.
base_fourier <- function(t, K, periodo = 365) {
  Phi <- matrix(1, length(t), 1)
  for (j in seq_len((K - 1) / 2)) {
    w <- 2 * pi * j / periodo
    Phi <- cbind(Phi, sin(w * t), cos(w * t))
  }
  Phi
}

# VCG(K) = promedio sobre curvas de n * SSE / (n - K)^2   (igual que en la fase 1)
vcg <- function(Phi, Ymat) {
  n <- nrow(Phi); K <- ncol(Phi)
  E <- Ymat - Phi %*% solve(crossprod(Phi), crossprod(Phi, Ymat))
  mean(n * colSums(E^2) / (n - K)^2)
}


#####################################################################################
# A. Fourier frente a B-spline
#####################################################################################

compara_bases <- function(t, Ymat, Ks) {
  res <- NULL
  for (K in Ks) {
    if (K %% 2 == 1) {                                   # Fourier solo con K impar
      Phi <- base_fourier(t, K)
      res <- rbind(res, data.frame(base = "Fourier", K = K, VCG = vcg(Phi, Ymat),
                                   condicion = kappa(crossprod(Phi), exact = TRUE)))
    }
    if (K >= 4) {
      Phi <- base_bspline(t, K)
      res <- rbind(res, data.frame(base = "B-spline", K = K, VCG = vcg(Phi, Ymat),
                                   condicion = kappa(crossprod(Phi), exact = TRUE)))
    }
  }
  res
}
bases_est <- compara_bases(t_est, Ytr, 3:31)
bases_sat <- compara_bases(t_sat, Xtr, 3:31)

mejor <- function(tab) do.call(rbind, lapply(split(tab, tab$base), function(s) s[which.min(s$VCG), ]))
cat("\n--- A. K de mínima VCG por base ---\nEstaciones:\n"); print(mejor(bases_est), digits = 5)
cat("Satélite:\n");                                       print(mejor(bases_sat), digits = 5)
cat("B-spline con K = 8 -> VCG estaciones:", round(vcg(base_bspline(t_est, 8), Ytr), 5),
    "| satélite:", round(vcg(base_bspline(t_sat, 8), Xtr), 5), "\n")

# Fourier obliga a que cada curva valga lo mismo el 1 de enero y el 31 de diciembre.
# ¿Pasa eso en los datos? Diferencia entre los primeros y los últimos 15 días de cada año.
dif_ext <- colMeans(Ytr[1:15, ]) - colMeans(Ytr[351:365, ])
cat("Primeros 15 días - últimos 15 días del mismo año: mediana |dif| =",
    round(median(abs(dif_ext)), 2), "frente a un nivel medio de", round(mean(Ytr), 2), "\n")

par(mfrow = c(1, 2))
for (fuente in c("estaciones", "satélite")) {
  tab <- if (fuente == "estaciones") bases_est else bases_sat
  f <- tab[tab$base == "Fourier", ]; b <- tab[tab$base == "B-spline", ]
  plot(b$K, b$VCG, type = "b", pch = 19, ylim = range(tab$VCG[tab$K >= 5]),
       xlab = "Número de bases (K)", ylab = "VCG",
       main = if (fuente == "estaciones") "(a) Estaciones" else "(b) Satélite")
  lines(f$K, f$VCG, type = "b", pch = 1, lty = 2)
  legend("topright", c("B-spline", "Fourier"), pch = c(19, 1), lty = c(1, 2), bty = "n")
}
par(mfrow = c(1, 1))


#####################################################################################
# B. Varianza de los residuos por día del año (¿mínimos cuadrados ponderados?)
#####################################################################################

Phi8 <- base_bspline(t_est, 8)
C8   <- solve(crossprod(Phi8), crossprod(Phi8, Ytr))
E8   <- Ytr - Phi8 %*% C8
v_dia <- rowMeans(E8^2)                                     # varianza residual de cada día
mes   <- rep(1:12, dias_mes)
v_mes <- tapply(v_dia, mes, mean)
cat("\n--- B. Varianza de los residuos por mes ---\n"); print(round(v_mes, 2))
cat("Mes más variable / menos variable:", round(max(v_mes) / min(v_mes), 2), "\n")

# Mínimos cuadrados ponderados (MCP) con pesos 1 / varianza del día (media móvil de 15 días)
v_suave <- stats::filter(c(tail(v_dia, 7), v_dia, head(v_dia, 7)), rep(1 / 15, 15))
v_suave <- as.numeric(v_suave[8:372])
w   <- 1 / v_suave
Cw  <- solve(crossprod(Phi8, w * Phi8), crossprod(Phi8, w * Ytr))
t_fino <- seq(0, 365, length.out = 731)
Pf  <- base_bspline(t_fino, 8)
m_mco <- rowMeans(Pf %*% C8); m_mcp <- rowMeans(Pf %*% Cw)
pico <- function(m, a, b) t_fino[t_fino > a & t_fino < b][which.max(m[t_fino > a & t_fino < b])]
cat("Máxima diferencia entre la curva media MCO y MCP:", round(max(abs(m_mco - m_mcp)), 3),
    "(nivel medio", round(mean(m_mco), 2), ")\n")
cat("Diferencia curva a curva (raíz del error cuadrático medio):",
    round(sqrt(mean((Pf %*% C8 - Pf %*% Cw)^2)), 3), "\n")
cat("Picos MCO:", round(pico(m_mco, 60, 200)), round(pico(m_mco, 200, 365)),
    "| picos MCP:", round(pico(m_mcp, 60, 200)), round(pico(m_mcp, 200, 365)), "\n")

par(mfrow = c(1, 2))
plot(t_est, v_dia, type = "l", col = "grey60", xlab = "Día del año",
     ylab = "Varianza de los residuos", main = "(a)")
lines(t_est, v_suave, lwd = 2.5, col = "black")
plot(t_fino, m_mco, type = "l", lwd = 2.5, xlab = "Día del año",
     ylab = "Precipitación media, √(mm/día)", main = "(b)")
lines(t_fino, m_mcp, lwd = 2, lty = 2, col = "grey40")
legend("topleft", c("Mínimos cuadrados ordinarios", "Mínimos cuadrados ponderados"),
       lty = c(1, 2), lwd = c(2.5, 2), col = c("black", "grey40"), bty = "n", cex = 0.8)
par(mfrow = c(1, 1))


#####################################################################################
# C. Suavizado penalizado: muchas bases y castigo a la curvatura
#####################################################################################
# c_lambda = (Phi'Phi + lambda R)^-1 Phi' y,   R = integral de phi''(t) phi''(t)'
# grados de libertad efectivos: df = tr(S_lambda).   VCG(lambda) = n SSE / (n - df)^2
# La segunda derivada se calcula con diferencias finitas sobre una malla fina.

matriz_R <- function(K, h = 0.01) {
  g  <- seq(h, 365 - h, by = 0.05)
  D2 <- (base_bspline(g + h, K) - 2 * base_bspline(g, K) + base_bspline(g - h, K)) / h^2
  crossprod(D2) * 0.05
}
vcg_penalizado <- function(t, Ymat, K, lambdas) {
  Phi <- base_bspline(t, K); R <- matriz_R(K); n <- length(t)
  t(sapply(lambdas, function(l) {
    S  <- Phi %*% solve(crossprod(Phi) + l * R, t(Phi))
    df <- sum(diag(S))
    E  <- Ymat - S %*% Ymat
    c(lambda = l, df = df, VCG = mean(n * colSums(E^2) / (n - df)^2))
  }))
}
lambdas <- 10^seq(0, 8, by = 0.5)
pen_est <- vcg_penalizado(t_est, Ytr, 40, lambdas)
pen_sat <- vcg_penalizado(t_sat, Xtr, 30, lambdas)
cat("\n--- C. Suavizado penalizado (mínimo de la VCG) ---\n")
print(round(pen_est[which.min(pen_est[, "VCG"]), ], 4))
print(round(pen_sat[which.min(pen_sat[, "VCG"]), ], 4))

par(mfrow = c(1, 2))
plot(pen_est[, "df"], pen_est[, "VCG"], type = "b", pch = 19, xlab = "Grados de libertad, tr(S)",
     ylab = "VCG", main = "(a) Estaciones")
points(8, vcg(base_bspline(t_est, 8), Ytr), pch = 4, cex = 1.8, lwd = 2)
plot(pen_sat[, "df"], pen_sat[, "VCG"], type = "b", pch = 19, xlab = "Grados de libertad, tr(S)",
     ylab = "VCG", main = "(b) Satélite")
points(8, vcg(base_bspline(t_sat, 8), Xtr), pch = 4, cex = 1.8, lwd = 2)
legend("topright", c("Penalizado (lambda variable)", "B-spline K = 8 sin penalizar"),
       pch = c(19, 4), bty = "n", cex = 0.8)
par(mfrow = c(1, 1))
