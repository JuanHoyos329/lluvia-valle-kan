# =============================================================================
# Evaluación
# =============================================================================
paquetes <- c("writexl")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan) > 0) install.packages(faltan, repos = "https://cloud.r-project.org")
invisible(lapply(paquetes, library, character.only = TRUE))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
SEMILLA <- 2026
set.seed(SEMILLA)
cat("R:", R.version.string, "| writexl", as.character(packageVersion("writexl")), "\n")

B_BOOT  <- 2000                         # remuestreos bootstrap (la rúbrica pide >= 1000)
ORDEN   <- 4                            # B-splines cúbicos
RANGO   <- c(0, 365)
MODELOS <- c("satelite", "lineal", "kan", "mlp")
NOMBRE  <- c(satelite = "satélite sin corregir", lineal = "red - cabeza lineal",
             kan = "red - cabeza KAN", mlp = "red - cabeza MLP")
SEMILLAS <- 0:2
BLOQUES <- c("entrenamiento", "validacion", "prueba")

F1  <- "./salidas_fase1"
F4  <- "./salidas_fase4"
PRE <- "./salidas_preprocesamiento"
OUT <- "./salidas_evaluacion"
for (d in c(F1, F4, PRE)) if (!dir.exists(d))
  stop("No se encuentra la carpeta '", d, "' desde ", getwd(), ". Ponga el directorio de trabajo en la carpeta del proyecto.")
dir.create(OUT, showWarnings = FALSE)
guardar <- function(tabla, nombre) write_xlsx(tabla, file.path(OUT, paste0(nombre, ".xlsx")))
clave_de <- function(df) paste(df$STATION_ID, df$YEAR)


#####################################################################################
# 1. Carga de datos
#####################################################################################

# ---- Coeficientes predichos por la red (3 semillas por cabeza). De aquí sale la PARTICIÓN.
leer_pred <- function(tipo, s) read.csv(file.path(F4, sprintf("pred_coef_%s_semilla%d.csv", tipo, s)))
p0 <- leer_pred("kan", 0)
claves <- p0[, c("STATION_ID", "YEAR", "particion")]
stopifnot(all(claves$particion %in% BLOQUES))
N <- nrow(claves)
C_pred <- lapply(setNames(MODELOS[-1], MODELOS[-1]), function(t) lapply(SEMILLAS, function(s) {
  p <- leer_pred(t, s)
  if (!identical(clave_de(p), clave_de(claves)) || !identical(p$particion, claves$particion))
    stop("pred_coef_", t, "_semilla", s, ".csv no tiene las mismas curvas o partición que pred_coef_kan_semilla0.csv")
  as.matrix(p[, grep("^c[0-9]+$", names(p))])
}))

# ---- Coeficientes de la fase 1 (escala sqrt(mm/día)): solo para las fechas de los picos y el control
coef_est <- read.csv(file.path(F1, "coeficientes_estaciones.csv"))
coef_sat <- read.csv(file.path(F1, "coeficientes_satelite.csv"))
if (!identical(clave_de(coef_est), clave_de(claves)) || !identical(clave_de(coef_sat), clave_de(claves)))
  stop("Las salidas de la fase 1 y de la fase 4 no tienen las mismas curvas: corra de nuevo la fase 4.")
C_est <- as.matrix(coef_est[, grep("^c[0-9]+$", names(coef_est))])      # N x 8
C_sat <- as.matrix(coef_sat[, grep("^c[0-9]+$", names(coef_sat))])      # N x 12
K_est <- ncol(C_est); K_sat <- ncol(C_sat)

# ---- Curvas diarias en mm (365 días; el 29 de febrero va sumado al 28, como en la fase 1)
leer_mm <- function(archivo) {
  m <- read.csv(file.path(F4, archivo))
  if (!identical(clave_de(m), clave_de(claves))) stop(archivo, " no tiene las mismas curvas que pred_coef_*.csv")
  m
}
obs <- leer_mm("mm_diario_observado.csv")
fase_enso <- obs$fase_enso
dcols <- paste0("d", 1:365)
Y_obs <- as.matrix(obs[, dcols])
Y_mm  <- lapply(setNames(MODELOS[-1], MODELOS[-1]),
                function(t) as.matrix(leer_mm(sprintf("mm_diario_pred_%s.csv", t))[, dcols]))

# ---- Satélite: total de cada péntada (mm) en el píxel de la estación, unido por estación-año
pent <- read.csv(file.path(PRE, "estaciones_pentadas_limpio.csv"))
pent <- pent[order(pent$STATION_ID, pent$YEAR, pent$PENTADA), ]
fila_pent <- match(clave_de(claves), clave_de(pent[pent$PENTADA == 1, ]))
if (anyNA(fila_pent))
  stop(sum(is.na(fila_pent)), " de ", N, " curvas de la fase 4 no están en estaciones_pentadas_limpio.csv. ",
       "Las fases 1-4 se corrieron con un preprocesamiento anterior: córralas de nuevo antes de evaluar.")
SAT_mm   <- matrix(pent$satelite_mm, ncol = 72, byrow = TRUE)[fila_pent, ]   # N x 72, mm por péntada
SAT_dias <- matrix(pent$dias, ncol = 72, byrow = TRUE)[fila_pent, ]          # días de cada péntada
if (nrow(pent) != 72 * nrow(pent[pent$PENTADA == 1, ])) stop("estaciones_pentadas_limpio.csv no tiene 72 péntadas por curva")

# ---- Día del año (1..365) -> péntada (1..72): péntadas 1-5 de cada mes con 5 días, la 6.ª el resto
dias_mes <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
mes_dia  <- data.frame(mes = rep(1:12, dias_mes), dia = unlist(lapply(dias_mes, seq_len)))
pentada_de_dia <- (mes_dia$mes - 1) * 6 + pmin(ceiling(mes_dia$dia / 5), 6)
stopifnot(length(pentada_de_dia) == 365, all(table(pentada_de_dia) %in% 3:6))
n_pent <- as.vector(table(pentada_de_dia))                                 # días de cada péntada (año de 365)

# Satélite en mm/día para cada día: el total de su péntada repartido por igual entre sus días
SAT_dia <- (SAT_mm / SAT_dias)[, pentada_de_dia]                           # N x 365

# Índices de cada bloque de la partición
idx_de <- lapply(setNames(BLOQUES, BLOQUES), function(b) which(claves$particion == b))
cat("Curvas: entrenamiento", length(idx_de$entrenamiento), "| validación", length(idx_de$validacion),
    "| prueba", length(idx_de$prueba), "\n")


#####################################################################################
# 2. Bases B-spline (fechas de los picos) y matrices de Gram (control en escala raíz, sección 7)
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


#####################################################################################
# 3. Las cuatro métricas por curva estación-año, en mm
#####################################################################################
# Escala común: PÉNTADAS en mm/día (la resolución nativa del satélite). Observado y red se agregan de
# sus días; el satélite es el total de su péntada dividido entre sus días.
a_pentadas <- function(Y_dia, usar = NULL) {                              # N x 365 -> N x 72 (mm/día)
  if (is.null(usar)) usar <- matrix(TRUE, nrow(Y_dia), 365)
  Y_dia[!usar] <- 0
  suma <- t(rowsum(t(Y_dia), pentada_de_dia))
  n    <- t(rowsum(t(usar * 1), pentada_de_dia))
  ifelse(n > 0, suma / n, NA)
}
P_obs <- a_pentadas(Y_obs)
P_mod <- c(list(satelite = SAT_mm / SAT_dias), lapply(Y_mm, a_pentadas))

# ---- 3.1 Error L² relativo en mm:  e_i = Σ_p n_p (ŷ_ip - y_ip)²  ≈ ∫ (ŷ_i(t) - y_i(t))² dt
# (n_p = días de la péntada). El L² RELATIVO de un grupo es 100 · sqrt(Σ e_i / Σ ‖y_i‖²).
# Red: la curva en mm de la fase 4 (promedio de las 3 semillas + smearing).
pesos  <- matrix(n_pent, N, 72, byrow = TRUE)
norma  <- rowSums(pesos * P_obs^2)                                        # ‖y_i‖² en mm²/día
err_L2 <- lapply(P_mod, function(P) rowSums(pesos * (P - P_obs)^2))

# ---- 3.2 KGE (Gupta et al., 2009) sobre las 72 péntadas del año, en mm/día
#   KGE = 1 - sqrt((r - 1)² + (α - 1)² + (β - 1)²)
#   r = correlación, α = sd(predicho)/sd(observado), β = media(predicho)/media(observado)
# Escala: péntadas y no días. La red predice la lluvia ESPERADA de cada día (una curva suave), y la
# lluvia de un día cualquiera salta entre 0 y > 50 mm. Día a día, α de cualquier predicción suave
# sale muy por debajo de 1 y el KGE cae por variabilidad, no por error de la forma de la curva. En
# péntadas el ruido diario se promedia, y es la escala más fina en la que el satélite tiene dato.
# La sección 5 muestra el KGE en las tres escalas (día, péntada, mes) con sus componentes.
kge_comp <- function(s, o) {
  ok <- !is.na(s) & !is.na(o)
  s <- s[ok]; o <- o[ok]
  if (length(o) < 3 || sd(o) == 0 || sd(s) == 0) return(c(KGE = NA, r = NA, alfa = NA, beta = NA))
  r <- cor(s, o); a <- sd(s) / sd(o); b <- mean(s) / mean(o)
  c(KGE = 1 - sqrt((r - 1)^2 + (a - 1)^2 + (b - 1)^2), r = r, alfa = a, beta = b)
}
kge_curvas <- function(P_m, P_o) t(sapply(seq_len(nrow(P_o)), function(i) kge_comp(P_m[i, ], P_o[i, ])))
KGE <- lapply(P_mod, function(P) kge_curvas(P, P_obs)[, "KGE"])

# ---- 3.3 Acumulado anual (mm). El sesgo de un grupo es 100 · (Σ predicho - Σ observado) / Σ observado
tot_obs <- rowSums(Y_obs)
tot_mod <- c(list(satelite = rowSums(SAT_mm)), lapply(Y_mm, rowSums))   # satélite: suma de sus 72 péntadas

# ---- 3.4 Fechas de los dos picos de lluvia (días del año)
# Mismo criterio de la fase 2: máximo de la curva en [45, 200] (primera temporada, mar-may) y en
# [230, 350] (segunda, sep-nov). Las fechas de un máximo no cambian entre sqrt(mm) y mm (la raíz es
# creciente), así que se usan las curvas suaves de la fase 1. Para la red: el promedio de las 3 semillas.
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
media_tr <- colMeans(C_est[idx_de$entrenamiento, ]) %*% t(P_tt_est)
cat("Picos de la curva media de entrenamiento: día", abs(picos(media_tr)), "\n")
cat("Curvas con los dos picos observados bien definidos:",
    round(100 * mean(rowSums(pico_real_obs) == 2), 1), "%\n")


#####################################################################################
# 4. Resumen por bloque: valor e IC bootstrap del 95 %
#####################################################################################
# Las curvas de una misma estación no son independientes (comparten sitio, píxel y altitud), así que
# el bootstrap remuestrea ESTACIONES con reemplazo y se lleva todas sus curvas del grupo.
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

evaluar <- function(idx, bloque, fase = "todas") {
  estaciones <- unique(claves$STATION_ID[idx])
  curvas_de  <- split(idx, factor(claves$STATION_ID[idx], levels = estaciones))
  obs_val    <- resumen(idx)
  # Bootstrap: B_BOOT réplicas, cada una con tantas estaciones como el grupo, sorteadas con reemplazo
  boot <- replicate(B_BOOT, {
    idx_b <- unlist(curvas_de[sample.int(length(estaciones), replace = TRUE)], use.names = FALSE)
    unlist(resumen(idx_b))
  })
  ic <- apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  filas <- list()
  for (m in MODELOS) {
    for (met in names(obs_val[[m]])) {
      clave <- paste(m, met, sep = ".")
      filas[[clave]] <- data.frame(
        particion = bloque, fase_enso = fase, modelo = NOMBRE[[m]], metrica = met,
        valor = obs_val[[m]][[met]], IC95_inf = ic[1, clave], IC95_sup = ic[2, clave],
        n_curvas = length(idx), n_estaciones = length(estaciones), row.names = NULL)
    }
  }
  do.call(rbind, filas)
}
ordenar <- function(tb) tb[order(match(tb$particion, BLOQUES), tb$metrica, match(tb$fase_enso, c("todas", "Niña", "Neutro", "Niño")),
                                 match(tb$modelo, NOMBRE)), ]

# ---- 4.1 Los tres bloques completos
# Entrenamiento (1981-1993) = ajuste (dentro de muestra); validación (1994-1996) = elección de
# hiperparámetros; prueba (1997-2003) = desempeño. Solo la prueba sirve para concluir.
cat("\nBootstrap (", B_BOOT, " remuestreos de estaciones) en los tres bloques...\n", sep = "")
tabla_eval1 <- ordenar(do.call(rbind, lapply(BLOQUES, function(b) evaluar(idx_de[[b]], b))))
print(format(tabla_eval1[, c("particion", "metrica", "modelo", "valor", "IC95_inf", "IC95_sup")], digits = 3), row.names = FALSE)
guardar(tabla_eval1, "tabla_eval1_metricas")

# ---- 4.2 Por fase ENSO dentro de cada bloque
# La validación (1994-1996) solo tiene años Neutro. La prueba tiene las tres fases:
# Niño (1997, 2002), Niña (1998-2000) y Neutro (2001, 2003), según la regla de la fase 4.
cat("\nBootstrap por fase ENSO...\n")
tabla_eval2 <- NULL
for (b in BLOQUES) for (f in c("Niña", "Neutro", "Niño")) {
  idx_f <- idx_de[[b]][fase_enso[idx_de[[b]]] == f]
  if (length(idx_f) > 0) tabla_eval2 <- rbind(tabla_eval2, evaluar(idx_f, b, f))
}
tabla_eval2 <- ordenar(tabla_eval2)
print(format(tabla_eval2[tabla_eval2$metrica == "KGE_mediana",
                         c("particion", "fase_enso", "modelo", "valor", "IC95_inf", "IC95_sup", "n_curvas")],
             digits = 3), row.names = FALSE)
guardar(tabla_eval2, "tabla_eval2_metricas_por_enso")


#####################################################################################
# 5. Escala del KGE: día, péntada y mes
#####################################################################################
# Justifica la escala de 3.2. Para cada escala: mediana por curva del KGE y de sus tres componentes.
# Si el KGE diario de la red sale bajo con α << 1 mientras r y β están bien, la caída es por la
# variabilidad diaria que una predicción suave (valor esperado) no puede ni debe reproducir.
mes_de_dia <- mes_dia$mes
a_meses <- function(Y_dia) t(rowsum(t(Y_dia), mes_de_dia)) / matrix(dias_mes, nrow(Y_dia), 12, byrow = TRUE)
Y_mod_dia <- c(list(satelite = SAT_dia), Y_mm)
escalas <- list(
  `diaria (365 días)`    = list(obs = Y_obs, mod = Y_mod_dia),
  `péntadas (72)`        = list(obs = P_obs, mod = P_mod),
  `mensual (12 meses)`   = list(obs = a_meses(Y_obs), mod = lapply(Y_mod_dia, a_meses)))
tabla_eval4 <- NULL
for (e in names(escalas)) for (m in MODELOS) {
  comp <- kge_curvas(escalas[[e]]$mod[[m]], escalas[[e]]$obs)
  for (b in BLOQUES) {
    cb <- comp[idx_de[[b]], , drop = FALSE]
    tabla_eval4 <- rbind(tabla_eval4, data.frame(
      particion = b, escala = e, modelo = NOMBRE[[m]],
      KGE_mediana = median(cb[, "KGE"], na.rm = TRUE), r_mediana = median(cb[, "r"], na.rm = TRUE),
      alfa_mediana = median(cb[, "alfa"], na.rm = TRUE), beta_mediana = median(cb[, "beta"], na.rm = TRUE)))
  }
}
tabla_eval4 <- tabla_eval4[order(match(tabla_eval4$particion, BLOQUES), match(tabla_eval4$modelo, NOMBRE)), ]
print(format(tabla_eval4, digits = 3), row.names = FALSE)
guardar(tabla_eval4, "tabla_eval4_KGE_por_escala")


#####################################################################################
# 6. Sensibilidad: KGE y sesgo anual usando solo los días MEDIDOS
#####################################################################################
# Una parte de los días fue imputada en el preprocesamiento (IDW_altitud o, sin vecinas, satélite)
# o desagregada de acumulados (solo en entrenamiento).
# Si se quitan, ¿cambian las conclusiones? Se usa la misma máscara para observado, satélite y red.
diario <- read.csv(file.path(PRE, "estaciones_diario_limpio.csv"))[, c("STATION_ID", "YEAR", "MONTH", "DAY", "origen")]
cum_mes <- cumsum(c(0, head(dias_mes, -1)))
diario$doy <- ifelse(diario$MONTH == 2 & diario$DAY == 29, 59, cum_mes[diario$MONTH] + diario$DAY)
diario$fila <- match(clave_de(diario), clave_de(claves))
diario <- diario[!is.na(diario$fila), ]
if (length(unique(diario$fila)) != N)
  stop("estaciones_diario_limpio.csv no tiene todas las curvas de la fase 4: corra de nuevo las fases 1-4.")
medido <- matrix(TRUE, N, 365)
no_med <- diario[diario$origen != "medido", ]
medido[cbind(no_med$fila, no_med$doy)] <- FALSE                 # 29 feb imputado -> el 28 tampoco cuenta
rm(diario, no_med)

P_obs_m <- a_pentadas(Y_obs, medido)
tabla_eval3 <- NULL
for (b in BLOQUES) {
  idx <- idx_de[[b]]
  for (m in MODELOS) {
    k_med <- median(kge_curvas(a_pentadas(Y_mod_dia[[m]], medido)[idx, ], P_obs_m[idx, ])[, "KGE"], na.rm = TRUE)
    o <- sum((Y_obs * medido)[idx, ]); s <- sum((Y_mod_dia[[m]] * medido)[idx, ])
    tabla_eval3 <- rbind(tabla_eval3, data.frame(
      particion = b, modelo = NOMBRE[[m]],
      pct_dias_no_medidos = 100 * mean(!medido[idx, ]),
      KGE_mediana_todos_los_dias = median(KGE[[m]][idx], na.rm = TRUE), KGE_mediana_solo_medidos = k_med,
      sesgo_anual_pct_todos_los_dias = 100 * (sum(tot_mod[[m]][idx]) - sum(tot_obs[idx])) / sum(tot_obs[idx]),
      sesgo_anual_pct_solo_medidos = 100 * (s - o) / o))
  }
}
print(format(tabla_eval3, digits = 3), row.names = FALSE)
guardar(tabla_eval3, "tabla_eval3_solo_dias_medidos")


#####################################################################################
# 7. Control: las cifras coinciden con las de la fase 4
#####################################################################################
# La fase 4 reporta el L² en sqrt(mm/día) (la escala de entrenamiento de la red); aquí se recalcula en
# esa escala SOLO para comprobar que se leen las mismas curvas, y se imprime al lado del L² en mm para
# mostrar cuánto infla la escala raíz el error del satélite (Jensen, ver encabezado).
W_est  <- gram(K_est, K_est); W_sat <- gram(K_sat, K_sat); W_cruz <- gram(K_sat, K_est)
stopifnot(abs(sum(W_est) - 365) < 1e-9, abs(sum(W_sat) - 365) < 1e-9)   # ∫ 1 dt = 365
forma <- function(A, W, B) rowSums((A %*% W) * B)                        # aᵢᵀ W bᵢ por fila
norma_r <- forma(C_est, W_est, C_est)
err_r <- list(satelite = forma(C_sat, W_sat, C_sat) - 2 * forma(C_sat, W_cruz, C_est) + norma_r)
for (t in MODELOS[-1])
  err_r[[t]] <- Reduce(`+`, lapply(C_pred[[t]], function(C) forma(C - C_est, W_est, C - C_est))) / length(SEMILLAS)
L2_raiz <- function(m, idx) 100 * sqrt(sum(err_r[[m]][idx]) / sum(norma_r[idx]))
L2_mm   <- function(m, idx) 100 * sqrt(sum(err_L2[[m]][idx]) / sum(norma[idx]))
cat("\nL² relativo (%) en escala raíz (fase 4) y en mm (esta evaluación):\n")
print(do.call(rbind, lapply(BLOQUES, function(b) data.frame(
  particion = b, modelo = NOMBRE[MODELOS],
  L2_raiz_pct = sapply(MODELOS, L2_raiz, idx = idx_de[[b]]),
  L2_mm_pct   = sapply(MODELOS, L2_mm,   idx = idx_de[[b]]), row.names = NULL))), digits = 3, row.names = FALSE)

avisos <- character(0)
abl4  <- read.csv(file.path(F4, "ablacion_resultados.csv"))
comp4 <- read.csv(file.path(F4, "ablacion_comparaciones_bootstrap.csv"))
if ("subconjunto" %in% names(comp4)) comp4 <- comp4[comp4$subconjunto == "Todas", ]
i_te <- idx_de$prueba
if (abs(L2_raiz("satelite", i_te) - abl4$prueba_media[1]) > 1e-6) avisos <- c(avisos, "L² del satélite en prueba")
for (m in c("kan", "mlp")) {
  ref <- comp4$dif_L2rel_pp[comp4$comparacion == paste(m, "- lineal")]
  if (length(ref) != 1 || abs((L2_raiz(m, i_te) - L2_raiz("lineal", i_te)) - ref) > 1e-6)
    avisos <- c(avisos, paste("diferencia de L²", m, "- lineal en prueba"))
}
rec4 <- read.csv(file.path(F4, "reconstruccion_mm.csv"))
for (m in MODELOS[-1]) for (b in BLOQUES) {
  ref <- rec4$sesgo_anual_agregado[rec4$modelo == paste("red - cabeza", m) & rec4$particion == b]
  idx <- idx_de[[b]]
  if (length(ref) == 1 && abs(100 * (sum(tot_mod[[m]][idx]) - sum(tot_obs[idx])) / sum(tot_obs[idx]) - ref) > 1e-2)
    avisos <- c(avisos, paste("sesgo anual", m, b))
}
if (length(avisos) == 0) cat("\nControl OK: L² (escala raíz) y sesgo anual coinciden con la fase 4.\n") else
  warning("No coincide con la fase 4: ", paste(avisos, collapse = "; "),
          ". Revise que todas las salidas de la fase 4 sean de la misma corrida.")
cat("Tablas guardadas en", OUT, "\n")