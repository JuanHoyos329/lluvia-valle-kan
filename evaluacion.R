# =============================================================================
# Evaluación
# =============================================================================

setwd("C:/Users/juana/OneDrive/Escritorio/lluvia-valle-kan")

library(writexl)

set.seed(2026)
B_BOOT  <- 2000                         # remuestreos bootstrap (la rúbrica pide >= 1000)
ORDEN   <- 4                            # B-splines cúbicos (igual que la fase 1)
RANGO   <- c(0, 365)
MODELOS <- c("satelite", "lineal", "kan", "mlp")
NOMBRE  <- c(satelite = "satélite sin corregir", lineal = "red - cabeza lineal",
             kan = "red - cabeza KAN", mlp = "red - cabeza MLP")
SEMILLAS <- 0:2

F1  <- "./salidas_fase1"
F4  <- "./salidas_fase4"
PRE <- "./salidas_preprocesamiento"
OUT <- "./salidas_evaluacion"
dir.create(OUT, showWarnings = FALSE)
guardar <- function(tabla, nombre) write_xlsx(tabla, file.path(OUT, paste0(nombre, ".xlsx")))


#####################################################################################
# 1. Carga de datos
#####################################################################################

# ---- Coeficientes de la fase 1 (escala sqrt(mm/día)): fila i = misma estación-año en los dos
coef_est <- read.csv(file.path(F1, "coeficientes_estaciones.csv"))
coef_sat <- read.csv(file.path(F1, "coeficientes_satelite.csv"))
stopifnot(identical(coef_est$STATION_ID, coef_sat$STATION_ID), identical(coef_est$YEAR, coef_sat$YEAR))
claves <- coef_est[, c("STATION_ID", "YEAR", "particion")]
N <- nrow(claves)
C_est <- as.matrix(coef_est[, grep("^c[0-9]+$", names(coef_est))])      # N x 8
C_sat <- as.matrix(coef_sat[, grep("^c[0-9]+$", names(coef_sat))])      # N x 12
K_est <- ncol(C_est); K_sat <- ncol(C_sat)

# ---- Coeficientes predichos por la red (3 semillas por cabeza), mismas filas que la fase 1
leer_pred <- function(tipo, s) {
  p <- read.csv(file.path(F4, sprintf("pred_coef_%s_semilla%d.csv", tipo, s)))
  stopifnot(identical(p$STATION_ID, claves$STATION_ID), identical(p$YEAR, claves$YEAR))
  as.matrix(p[, grep("^c[0-9]+$", names(p))])
}
C_pred <- lapply(setNames(MODELOS[-1], MODELOS[-1]), function(t) lapply(SEMILLAS, function(s) leer_pred(t, s)))

# ---- Curvas diarias en mm (365 días; el 29 de febrero va sumado al 28, como en la fase 1)
leer_mm <- function(archivo) {
  m <- read.csv(file.path(F4, archivo))
  stopifnot(identical(m$STATION_ID, claves$STATION_ID), identical(m$YEAR, claves$YEAR))
  m
}
obs <- leer_mm("mm_diario_observado.csv")
fase_enso <- obs$fase_enso
dcols <- paste0("d", 1:365)
Y_obs <- as.matrix(obs[, dcols])
Y_mm  <- lapply(setNames(MODELOS[-1], MODELOS[-1]),
                function(t) as.matrix(leer_mm(sprintf("mm_diario_pred_%s.csv", t))[, dcols]))

# ---- Satélite: total de cada péntada en el píxel de la estación
pent <- read.csv(file.path(PRE, "estaciones_pentadas_limpio.csv"))
pent <- pent[order(pent$STATION_ID, pent$YEAR, pent$PENTADA), ]
stopifnot(nrow(pent) == 72 * N,
          identical(pent$STATION_ID[seq(1, nrow(pent), 72)], claves$STATION_ID),
          identical(pent$YEAR[seq(1, nrow(pent), 72)], claves$YEAR))
SAT_mm   <- matrix(pent$satelite_mm, nrow = N, byrow = TRUE)               # N x 72, mm por péntada
SAT_dias <- matrix(pent$dias, nrow = N, byrow = TRUE)                      # días de cada péntada

# ---- Día del año (1..365) -> péntada (1..72): péntadas 1-5 de cada mes con 5 días, la 6.ª el resto
dias_mes <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
mes_dia  <- data.frame(mes = rep(1:12, dias_mes), dia = unlist(lapply(dias_mes, seq_len)))
pentada_de_dia <- (mes_dia$mes - 1) * 6 + pmin(ceiling(mes_dia$dia / 5), 6)
stopifnot(length(pentada_de_dia) == 365, all(table(pentada_de_dia) %in% 3:6))

# Satélite en mm/día para cada día: el total de su péntada repartido por igual entre sus días
SAT_dia <- (SAT_mm / SAT_dias)[, pentada_de_dia]                           # N x 365

# Índices de cada partición. El entrenamiento no se evalúa (ya se usó para ajustar).
i_va <- which(claves$particion == "validacion")
i_te <- which(claves$particion == "prueba")
cat("Curvas evaluadas: validación", length(i_va), "| prueba", length(i_te), "\n")


#####################################################################################
# 2. Bases B-spline y matrices de Gram (para el error L² en la escala de las curvas)
#####################################################################################
# Misma base de la fase 1. Las integrales se hacen por Gauss-Legendre de 4 puntos en cada tramo
# entre nudos; como el producto de dos cúbicas tiene grado 6 <= 7, la cuadratura es EXACTA.
nudos_bspline <- function(K, orden = ORDEN, rango = RANGO) {
  internos <- seq(rango[1], rango[2], length.out = K - orden + 2)
  internos <- internos[-c(1, length(internos))]
  c(rep(rango[1], orden), internos, rep(rango[2], orden))
}

base_bspline <- function(t, K, orden = ORDEN, rango = RANGO) {
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

# W = ∫ φ_1(t) φ_2(t)ᵀ dt, en la unión de los nudos de las dos bases
gram <- function(K1, K2) {
  u  <- sort(unique(c(nudos_bspline(K1), nudos_bspline(K2))))
  x  <- c(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
  w  <- c( 0.3478548451374538,  0.6521451548625461, 0.6521451548625461, 0.3478548451374538)
  a  <- head(u, -1); b <- u[-1]
  tq <- as.vector(outer(x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
  wq <- as.vector(outer(w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
  t(base_bspline(tq, K1)) %*% (wq * base_bspline(tq, K2))
}
W_est  <- gram(K_est, K_est)           # 8 x 8
W_sat  <- gram(K_sat, K_sat)           # 12 x 12
W_cruz <- gram(K_sat, K_est)           # 12 x 8: el satélite y la estación tienen bases distintas
stopifnot(abs(sum(W_est) - 365) < 1e-9, abs(sum(W_sat) - 365) < 1e-9)   # ∫ 1 dt = 365

forma <- function(A, W, B) rowSums((A %*% W) * B)                        # aᵢᵀ W bᵢ por fila


#####################################################################################
# 3. Las cuatro métricas por curva estación-año
#####################################################################################

# ---- 3.1 Error L² (escala sqrt(mm/día)):  e_i = ∫ (ŷ_i(t) - y_i(t))² dt
# Satélite: ∫ (x - y)² = xᵀW_sat x - 2 xᵀW_cruz y + yᵀW_est y
# Red: error de cada semilla y luego el promedio de las 3 semillas (igual que la fase 4).
# El L² RELATIVO de un grupo de curvas es 100 · sqrt(Σ e_i / Σ ∫ y_i²): se calcula al resumir.
norma <- forma(C_est, W_est, C_est)                                       # ∫ y_i² dt
err_L2 <- list(satelite = forma(C_sat, W_sat, C_sat) - 2 * forma(C_sat, W_cruz, C_est) + norma)
for (t in MODELOS[-1]) {
  err_L2[[t]] <- Reduce(`+`, lapply(C_pred[[t]], function(C) forma(C - C_est, W_est, C - C_est))) / length(SEMILLAS)
}

# ---- 3.2 KGE (Gupta et al., 2009) sobre las 72 péntadas del año, en mm/día
#   KGE = 1 - sqrt((r - 1)² + (α - 1)² + (β - 1)²)
#   r = correlación, α = sd(predicho)/sd(observado), β = media(predicho)/media(observado)
#   KGE = 1 es perfecto. Se usan péntadas (la resolución del satélite) y no días sueltos: la red
#   predice la lluvia ESPERADA, una curva suave, y la lluvia de un día cualquiera es casi ruido.
a_pentadas <- function(Y_dia, usar = NULL) {                              # N x 365 -> N x 72 (mm/día)
  if (is.null(usar)) usar <- matrix(TRUE, nrow(Y_dia), 365)
  Y_dia[!usar] <- 0
  suma <- t(rowsum(t(Y_dia), pentada_de_dia))
  n    <- t(rowsum(t(usar * 1), pentada_de_dia))
  ifelse(n > 0, suma / n, NA)
}
kge <- function(s, o) {
  ok <- !is.na(s) & !is.na(o)
  s <- s[ok]; o <- o[ok]
  if (length(o) < 3 || sd(o) == 0 || sd(s) == 0) return(NA_real_)
  1 - sqrt((cor(s, o) - 1)^2 + (sd(s) / sd(o) - 1)^2 + (mean(s) / mean(o) - 1)^2)
}
kge_curvas <- function(P_mod, P_obs) sapply(seq_len(nrow(P_obs)), function(i) kge(P_mod[i, ], P_obs[i, ]))

P_obs <- a_pentadas(Y_obs)
P_mod <- c(list(satelite = SAT_mm / SAT_dias), lapply(Y_mm, a_pentadas))
KGE   <- lapply(P_mod, kge_curvas, P_obs = P_obs)

# ---- 3.3 Acumulado anual (mm). El sesgo de un grupo es 100 · (Σ predicho - Σ observado) / Σ observado
tot_obs <- rowSums(Y_obs)
tot_mod <- c(list(satelite = rowSums(SAT_mm)), lapply(Y_mm, rowSums))   # satélite: suma de sus 72 péntadas

# ---- 3.4 Fechas de los dos picos de lluvia (días del año)
# Mismo criterio de la fase 2: máximo de la curva en [45, 200] (primera temporada, mar-may) y en
# [230, 350] (segunda, sep-nov). Las fechas no cambian entre sqrt(mm) y mm (la raíz es creciente),
# así que se usan las curvas de la fase 1. Para la red: el promedio de las 3 semillas.
# Un pico observado solo cuenta si es un máximo "real" (no queda pegado al borde de la ventana);
# si la estación no tuvo ese pico ese año, no hay fecha que comparar y la curva sale de esa métrica.
tt   <- seq(0, 365, by = 0.5)
P_tt_est <- base_bspline(tt, K_est)
P_tt_sat <- base_bspline(tt, K_sat)
VENTANAS <- list(pico1 = c(45, 200), pico2 = c(230, 350))
picos <- function(X) {                                                    # X: N x length(tt)
  sapply(VENTANAS, function(v) {
    w <- which(tt >= v[1] & tt <= v[2])
    i <- w[max.col(X[, w, drop = FALSE], ties.method = "first")]
    real <- X[cbind(seq_len(nrow(X)), i)] > X[cbind(seq_len(nrow(X)), i - 1)] &
      X[cbind(seq_len(nrow(X)), i)] > X[cbind(seq_len(nrow(X)), i + 1)]
    ifelse(real, tt[i], -tt[i])                                           # signo negativo = no es pico real
  })
}
pk_obs <- picos(C_est %*% t(P_tt_est))
pico_real_obs <- pk_obs > 0
pk_obs <- abs(pk_obs)
pk_mod <- list(satelite = abs(picos(C_sat %*% t(P_tt_sat))))
for (t in MODELOS[-1]) pk_mod[[t]] <- abs(picos((Reduce(`+`, C_pred[[t]]) / length(SEMILLAS)) %*% t(P_tt_est)))
err_pico <- lapply(pk_mod, function(p) ifelse(pico_real_obs, abs(p - pk_obs), NA))   # N x 2, en días

# Control de las ventanas con la curva media de ENTRENAMIENTO (no se ajusta nada con val/prueba)
media_tr <- colMeans(C_est[claves$particion == "entrenamiento", ]) %*% t(P_tt_est)
cat("Picos de la curva media de entrenamiento: día", abs(picos(media_tr)), "\n")
cat("Curvas de val+prueba con los dos picos observados bien definidos:",
    round(100 * mean(rowSums(pico_real_obs[c(i_va, i_te), ]) == 2), 1), "%\n")


#####################################################################################
# 4. Resumen por grupo: valor, IC bootstrap del 95 % y Wilcoxon pareado contra el satélite
#####################################################################################
# Las curvas de una misma estación no son independientes (comparten sitio, píxel y altitud), así que
# el bootstrap remuestrea ESTACIONES con reemplazo y se lleva todas sus curvas del grupo.
#
# Resumen de cada métrica en un conjunto de curvas idx:
resumen <- function(idx) {
  out <- list()
  for (m in MODELOS) {
    out[[m]] <- c(
      L2_relativo_pct   = 100 * sqrt(sum(err_L2[[m]][idx]) / sum(norma[idx])),
      KGE_mediana       = median(KGE[[m]][idx], na.rm = TRUE),
      sesgo_anual_pct   = 100 * (sum(tot_mod[[m]][idx]) - sum(tot_obs[idx])) / sum(tot_obs[idx]),
      error_pico1_dias  = mean(err_pico[[m]][idx, 1], na.rm = TRUE),
      error_pico2_dias  = mean(err_pico[[m]][idx, 2], na.rm = TRUE))
  }
  out
}

# Valor de cada métrica POR ESTACIÓN (para Wilcoxon): una cifra por estación y modelo
por_estacion <- function(idx, m) {
  g <- factor(claves$STATION_ID[idx])
  data.frame(
    L2_relativo_pct  = 100 * sqrt(tapply(err_L2[[m]][idx], g, sum) / tapply(norma[idx], g, sum)),
    KGE_mediana      = tapply(KGE[[m]][idx], g, median, na.rm = TRUE),
    sesgo_anual_pct  = abs(100 * (tapply(tot_mod[[m]][idx], g, sum) - tapply(tot_obs[idx], g, sum)) /
                             tapply(tot_obs[idx], g, sum)),                  # |sesgo|: cuánto se aleja
    error_pico1_dias = tapply(err_pico[[m]][idx, 1], g, mean, na.rm = TRUE),
    error_pico2_dias = tapply(err_pico[[m]][idx, 2], g, mean, na.rm = TRUE))
}

evaluar <- function(idx, particion, fase = "todas") {
  estaciones <- unique(claves$STATION_ID[idx])
  curvas_de  <- split(idx, factor(claves$STATION_ID[idx], levels = estaciones))
  obs_val    <- resumen(idx)
  # Bootstrap: B_BOOT réplicas, cada una con tantas estaciones como el grupo, sorteadas con reemplazo
  boot <- replicate(B_BOOT, {
    idx_b <- unlist(curvas_de[sample.int(length(estaciones), replace = TRUE)], use.names = FALSE)
    unlist(resumen(idx_b))
  })
  ic <- apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  # Wilcoxon pareado por estación: cada modelo contra el satélite sin corregir
  est_sat <- por_estacion(idx, "satelite")
  filas <- list()
  for (m in MODELOS) {
    est_m <- por_estacion(idx, m)
    for (met in names(obs_val[[m]])) {
      clave <- paste(m, met, sep = ".")
      dif <- est_m[[met]] - est_sat[[met]]
      ok  <- !is.na(dif)
      p   <- if (m == "satelite" || sum(ok) < 5) NA else
        suppressWarnings(wilcox.test(est_m[[met]][ok], est_sat[[met]][ok], paired = TRUE)$p.value)
      filas[[clave]] <- data.frame(
        particion = particion, fase_enso = fase, modelo = NOMBRE[[m]], metrica = met,
        valor = obs_val[[m]][[met]], IC95_inf = ic[1, clave], IC95_sup = ic[2, clave],
        n_curvas = length(idx), n_estaciones = length(estaciones),
        mediana_dif_estacion_vs_satelite = if (m == "satelite") NA else median(dif, na.rm = TRUE),
        wilcoxon_p_vs_satelite = p, row.names = NULL)
    }
  }
  do.call(rbind, filas)
}

# ---- 4.1 Validación y prueba completas
cat("\nBootstrap (", B_BOOT, " remuestreos de estaciones) en validación y prueba...\n", sep = "")
tabla_eval1 <- rbind(evaluar(i_va, "validacion"), evaluar(i_te, "prueba"))
tabla_eval1 <- tabla_eval1[order(tabla_eval1$particion != "validacion", tabla_eval1$metrica), ]
print(format(tabla_eval1[, c("particion", "metrica", "modelo", "valor", "IC95_inf", "IC95_sup",
                             "wilcoxon_p_vs_satelite")], digits = 3), row.names = FALSE)
guardar(tabla_eval1, "tabla_eval1_metricas")

# ---- 4.2 Por fase ENSO (dentro de cada partición; la prueba no tiene años Niña)
cat("\nBootstrap por fase ENSO...\n")
tabla_eval2 <- NULL
for (part in c("validacion", "prueba")) {
  idx_p <- if (part == "validacion") i_va else i_te
  for (f in c("Niña", "Neutro", "Niño")) {
    idx_f <- idx_p[fase_enso[idx_p] == f]
    if (length(idx_f) == 0) next
    tabla_eval2 <- rbind(tabla_eval2, evaluar(idx_f, part, f))
  }
}
tabla_eval2 <- tabla_eval2[order(tabla_eval2$particion != "validacion", tabla_eval2$metrica,
                                 tabla_eval2$fase_enso), ]
print(format(tabla_eval2[tabla_eval2$metrica == "L2_relativo_pct",
                         c("particion", "fase_enso", "modelo", "valor", "IC95_inf", "IC95_sup", "n_curvas")],
             digits = 3), row.names = FALSE)
guardar(tabla_eval2, "tabla_eval2_metricas_por_enso")


#####################################################################################
# 5. Sensibilidad: KGE y sesgo anual usando solo los días MEDIDOS
#####################################################################################
# Cerca del 1 % de los días de validación y prueba fueron imputados en el preprocesamiento (CCW).
# Si se quitan, ¿cambian las conclusiones? Se usa la misma máscara para observado, satélite y red.
diario <- read.csv(file.path(PRE, "estaciones_diario_limpio.csv"))[, c("STATION_ID", "YEAR", "MONTH", "DAY", "origen")]
diario <- diario[diario$YEAR >= 1996, ]
cum_mes <- cumsum(c(0, head(dias_mes, -1)))
diario$doy <- ifelse(diario$MONTH == 2 & diario$DAY == 29, 59, cum_mes[diario$MONTH] + diario$DAY)
diario$fila <- match(paste(diario$STATION_ID, diario$YEAR), paste(claves$STATION_ID, claves$YEAR))
stopifnot(!anyNA(diario$fila))
medido <- matrix(TRUE, N, 365)
no_med <- diario[diario$origen != "medido", ]
medido[cbind(no_med$fila, no_med$doy)] <- FALSE                 # 29 feb imputado -> el 28 tampoco cuenta
rm(diario, no_med)

Y_mod_dia <- c(list(satelite = SAT_dia), Y_mm)
P_obs_m   <- a_pentadas(Y_obs, medido)
tabla_eval3 <- NULL
for (part in c("validacion", "prueba")) {
  idx <- if (part == "validacion") i_va else i_te
  for (m in MODELOS) {
    k_todos <- median(KGE[[m]][idx], na.rm = TRUE)
    k_med   <- median(kge_curvas(a_pentadas(Y_mod_dia[[m]], medido)[idx, ], P_obs_m[idx, ]), na.rm = TRUE)
    o <- sum((Y_obs * medido)[idx, ]); s <- sum((Y_mod_dia[[m]] * medido)[idx, ])
    tabla_eval3 <- rbind(tabla_eval3, data.frame(
      particion = part, modelo = NOMBRE[[m]],
      pct_dias_imputados = 100 * mean(!medido[idx, ]),
      KGE_mediana_todos_los_dias = k_todos, KGE_mediana_solo_medidos = k_med,
      sesgo_anual_pct_todos_los_dias = 100 * (sum(tot_mod[[m]][idx]) - sum(tot_obs[idx])) / sum(tot_obs[idx]),
      sesgo_anual_pct_solo_medidos = 100 * (s - o) / o))
  }
}
print(format(tabla_eval3, digits = 3), row.names = FALSE)
guardar(tabla_eval3, "tabla_eval3_solo_dias_medidos")


#####################################################################################
# 6. Control: las cifras coinciden con las de la fase 4
#####################################################################################
# L² relativo del satélite y diferencias entre cabezas en prueba (ablacion_comparaciones_bootstrap.csv)
# y sesgo agregado (reconstruccion_mm.csv). Si algo no coincide, el script se detiene.
comp4 <- read.csv(file.path(F4, "ablacion_comparaciones_bootstrap.csv"))
rec4  <- read.csv(file.path(F4, "reconstruccion_mm.csv"))
L2te  <- sapply(MODELOS, function(m) 100 * sqrt(sum(err_L2[[m]][i_te]) / sum(norma[i_te])))
abl4  <- read.csv(file.path(F4, "ablacion_resultados.csv"))
stopifnot(abs(L2te[["satelite"]] - abl4$prueba_media[1]) < 1e-6,
          abs((L2te[["kan"]] - L2te[["lineal"]]) - comp4$dif_L2rel_pp[comp4$comparacion == "kan - lineal"]) < 1e-6,
          abs((L2te[["mlp"]] - L2te[["lineal"]]) - comp4$dif_L2rel_pp[comp4$comparacion == "mlp - lineal"]) < 1e-6)
for (m in MODELOS[-1]) for (part in c("validacion", "prueba")) {
  idx <- if (part == "validacion") i_va else i_te
  ref <- rec4$sesgo_anual_agregado[rec4$modelo == paste("red - cabeza", m) & rec4$particion == part]
  stopifnot(abs(100 * (sum(tot_mod[[m]][idx]) - sum(tot_obs[idx])) / sum(tot_obs[idx]) - ref) < 1e-2)
}
cat("\nControl OK: L² relativo y sesgo anual coinciden con la fase 4.\n")
cat("Tablas guardadas en", OUT, "\n")