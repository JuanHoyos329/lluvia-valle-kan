# =============================================================================
# Fase 3 - Comparación funcional por fuente (ítem 3) y regresión funcional de referencia
# Proyecto final de Analítica de Datos: lluvia del Valle del Cauca, 1981-2003
#
# Fragmentos (se irán agregando):
#   0. Curvas comparables: estación y satélite en la misma escala y con la MISMA base.
#   1. Emparejamiento estación-píxel documentado.
#   2. Curva de diferencia d(t) = estación - satélite con banda +- 2EE.
#   3. Correlación cruzada funcional y su diagonal.
#   4. Comparación de las covarianzas de cada fuente.
#   5. Sesgo por altitud / región y por fase ENSO.
#   6. Regresión funcional función-función (modelo de referencia).
#
# Por qué no se restan directamente las curvas de la fase 1:
#   la estación se suavizó con 365 valores diarios de sqrt(mm) y el satélite con 72 valores de
#   sqrt(promedio diario de la péntada). Con lluvia intermitente, sqrt(promedio) > promedio(sqrt),
#   así que esa resta mezclaría el sesgo del satélite con un efecto de la agregación.
#   Aquí las dos fuentes se agregan igual (total de la péntada / días) antes de sacar la raíz.
#
# Uso (desde la carpeta del repositorio):  source("fase3_comparacion_funcional.R")
# Paquetes: fda, readxl, terra (terra solo para verificar la extracción del píxel)
# =============================================================================

rm(list = ls())
library(fda)
library(readxl)
library(terra)

ENT0 <- "./salidas_preprocesamiento"     # datos limpios (Hoyos)
ENT1 <- "./salidas_fase1"                # coeficientes de la fase 1 (Daniela): orden de las curvas
DSAT <- "./Datos/ImagenesSatelitePentadal"
OUT  <- "./salidas_fase3"
FIG  <- file.path(OUT, "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)
ENTRENA <- 1981:1995                     # misma partición que las demás fases

guardar <- function(df, nombre) write.csv(df, file.path(OUT, paste0(nombre, ".csv")), row.names = FALSE)

# Dibuja en pantalla (si hay) y guarda el PNG en salidas_fase3/figuras
figura <- function(nombre, codigo, ancho = 10, alto = 4.5) {
  codigo <- substitute(codigo); entorno <- parent.frame()
  if (interactive()) { eval(codigo, entorno); par(mfrow = c(1, 1)) }
  png(file.path(FIG, paste0(nombre, ".png")), width = ancho, height = alto, units = "in", res = 150)
  on.exit(dev.off())
  eval(codigo, entorno)
}


#####################################################################################
# 0. Curvas comparables
#####################################################################################

# ---- 0.1 Curvas estación-año, en el mismo orden que la fase 1 (y que la fase 4)
ce     <- read.csv(file.path(ENT1, "coeficientes_estaciones.csv"))
claves <- ce[, c("STATION_ID", "YEAR", "particion")]
N      <- nrow(claves)

# ---- 0.2 Péntadas: lluvia de la estación y del satélite en su píxel (misma péntada)
p <- read.csv(file.path(ENT0, "estaciones_pentadas_limpio.csv"))
p <- p[order(p$STATION_ID, p$YEAR, p$PENTADA), ]
stopifnot(nrow(p) == 72 * N,
          all(p$STATION_ID[p$PENTADA == 1] == claves$STATION_ID),
          all(p$YEAR[p$PENTADA == 1] == claves$YEAR),
          !anyNA(p$satelite_mm))

# Las dos fuentes igual: total de la péntada / días de la péntada (3 a 6) -> mm/día -> raíz
Y <- matrix(sqrt(p$lluvia_mm   / p$dias), nrow = 72)   # estación  (72 x N)
X <- matrix(sqrt(p$satelite_mm / p$dias), nrow = 72)   # satélite  (72 x N)
colnames(Y) <- colnames(X) <- paste(claves$STATION_ID, claves$YEAR, sep = "_")

# Tabla de escala: en mm/día las fuentes casi coinciden; en raíz no, porque el satélite reparte
# la lluvia en más péntadas (menos péntadas secas y menos extremos)
escala <- data.frame(
  fuente             = c("Estación", "Satélite"),
  mm_dia_medio       = round(c(mean(Y^2), mean(X^2)), 3),
  raiz_media         = round(c(mean(Y), mean(X)), 3),
  pct_pentadas_secas = round(100 * c(mean(Y == 0), mean(X == 0)), 2),
  p99_mm_dia         = round(c(quantile(Y^2, 0.99), quantile(X^2, 0.99)), 2))
print(escala)
guardar(escala, "tabla_f3_01_escala_fuentes")

# ---- 0.3 Dominio común [0, 365]: cada péntada en el centro de sus días (año no bisiesto)
fin_mes <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
ini_mes <- c(0, fin_mes[-12])
t_pent  <- unlist(lapply(1:12, function(m) {
  ini <- ini_mes[m] + c(0, 5, 10, 15, 20, 25)          # días 1-5, 6-10, ..., 26-fin de mes
  fin <- c(ini[-1], fin_mes[m])
  (ini + fin) / 2
}))
stopifnot(length(t_pent) == 72)

# ---- 0.4 Una sola base B-spline cúbica para las dos fuentes; K por VCG solo con entrenamiento
tr <- claves$YEAR %in% ENTRENA
Ks <- c(4, 5, seq(6, 30, by = 2))          # desde 4 (mínimo para orden 4) para no cortar el mínimo
vcg <- as.data.frame(t(sapply(Ks, function(K) {
  b <- create.bspline.basis(c(0, 365), nbasis = K, norder = 4)
  c(K        = K,
    estacion = mean(smooth.basis(t_pent, Y[, tr], b)$gcv),
    satelite = mean(smooth.basis(t_pent, X[, tr], b)$gcv))
})))
# Las fuentes tienen escalas de error distintas: se suma la VCG de cada una relativa a su mínimo
vcg$rel_estacion <- vcg$estacion / min(vcg$estacion)
vcg$rel_satelite <- vcg$satelite / min(vcg$satelite)
vcg$rel_suma     <- vcg$rel_estacion + vcg$rel_satelite
K_com <- vcg$K[which.min(vcg$rel_suma)]
print(vcg)
cat("K por fuente: estación", vcg$K[which.min(vcg$estacion)], "| satélite", vcg$K[which.min(vcg$satelite)],
    "| K común:", K_com, "\n")
guardar(vcg, "tabla_f3_02_vcg_K_comun")

figura("fig_f3_01_vcg_K_comun", {
  matplot(vcg$K, vcg[, c("rel_estacion", "rel_satelite")], type = "b", pch = 19, lty = 1,
          col = c("slategray", "darkturquoise"), xlab = "K (número de funciones base)",
          ylab = "VCG / mínimo", main = "Selección de una base común (péntadas, 1981-1995)")
  abline(v = K_com, lty = 2, col = "red")
  legend("top", c("Estación", "Satélite", paste("K común =", K_com)),
         col = c("slategray", "darkturquoise", "red"), lty = c(1, 1, 2), pch = c(19, 19, NA), bty = "n")
}, ancho = 7, alto = 4.5)

# ---- 0.5 Datos funcionales de las dos fuentes con la misma base
base_com <- create.bspline.basis(c(0, 365), nbasis = K_com, norder = 4)
nombres  <- list("Día del año", "Estación-año", "sqrt(mm/día)")
Yfd <- smooth.basis(t_pent, Y, base_com, fdnames = nombres)$fd
Xfd <- smooth.basis(t_pent, X, base_com, fdnames = nombres)$fd

# Las mismas curvas en mm/día (escala en la que se lee el sesgo). En raíz el satélite queda por
# encima todo el año aunque en mm llueva un poco menos: como no tiene péntadas secas y tiene menos
# extremos, su media de sqrt es mayor (desigualdad de Jensen). Por eso d(t) se interpreta en
# mm/día y la escala sqrt (la del modelo de la fase 4) se reporta como complemento.
nombres_mm <- list("Día del año", "Estación-año", "mm/día")
Ymm_fd <- smooth.basis(t_pent, Y^2, base_com, fdnames = nombres_mm)$fd
Xmm_fd <- smooth.basis(t_pent, X^2, base_com, fdnames = nombres_mm)$fd

t_eval <- seq(0, 365, length.out = 366)
Yv <- eval.fd(t_eval, Yfd)                               # 366 x N
Xv <- eval.fd(t_eval, Xfd)
Ymm <- eval.fd(t_eval, Ymm_fd)
Xmm <- eval.fd(t_eval, Xmm_fd)

j_ej <- sample(N, 1)                                     # curva de ejemplo (fija)
figura("fig_f3_02_medias_fuentes", {
  par(mfrow = c(1, 3))
  plot(t_eval, rowMeans(Ymm), type = "l", lwd = 3, col = "slategray", ylim = range(0, rowMeans(Ymm), rowMeans(Xmm)),
       xlab = "Día del año", ylab = "mm/día", main = "Media funcional en mm/día")
  lines(t_eval, rowMeans(Xmm), lwd = 3, lty = 2, col = "darkturquoise")
  legend("top", c("Estación", "Satélite"), col = c("slategray", "darkturquoise"), lty = 1:2, lwd = 3, bty = "n")
  plot(t_eval, rowMeans(Yv), type = "l", lwd = 3, col = "slategray", ylim = range(rowMeans(Yv), rowMeans(Xv)),
       xlab = "Día del año", ylab = "sqrt(mm/día)", main = "Media funcional en sqrt(mm/día)")
  lines(t_eval, rowMeans(Xv), lwd = 3, lty = 2, col = "darkturquoise")
  legend("top", c("Estación", "Satélite"), col = c("slategray", "darkturquoise"), lty = 1:2, lwd = 3, bty = "n")
  j <- j_ej
  plot(t_pent, Y[, j], pch = 19, cex = 0.6, col = "slategray", ylim = range(Y[, j], X[, j]),
       xlab = "Día del año", ylab = "sqrt(mm/día)", main = paste("Ejemplo:", colnames(Y)[j]))
  points(t_pent, X[, j], pch = 17, cex = 0.6, col = "darkturquoise")
  lines(t_eval, Yv[, j], lwd = 2, col = "slategray"); lines(t_eval, Xv[, j], lwd = 2, lty = 2, col = "darkturquoise")
}, ancho = 14, alto = 4.5)


#####################################################################################
# 1. Emparejamiento estación-píxel
#####################################################################################

px  <- as.data.frame(read_excel(file.path(ENT0, "estacion_pixel.xlsx")))
rk  <- as.data.frame(read_excel(file.path(ENT0, "ranking_estaciones.xlsx")))
est <- unique(claves$STATION_ID)

emp <- px[match(est, px$STATION_ID), c("STATION_ID", "STATION_NAME", "Lon", "lat", "fuente_coord", "fila", "col")]
stopifnot(!anyNA(emp$STATION_ID))
i <- match(emp$STATION_ID, rk$STATION_ID)
emp$altitud    <- round(rk$altitud[i])
emp$fuente_alt <- rk$fuente_alt[i]
emp$confianza  <- rk$confianza[i]
emp$curvas     <- tabulate(match(claves$STATION_ID, emp$STATION_ID), nbins = nrow(emp))

# ---- 1.1 Regla de asignación: el píxel que CONTIENE la estación
# Malla de 0,05° (37 x 39), esquina superior izquierda en (-77,55; 5,00); fila 1 = norte
col_r  <- floor((emp$Lon + 77.55) / 0.05) + 1
fila_r <- floor((5.00 - emp$lat) / 0.05) + 1
stopifnot(all(col_r == emp$col), all(fila_r == emp$fila))
emp$lon_pixel      <- -77.55 + (emp$col - 0.5) * 0.05
emp$lat_pixel      <-   5.00 - (emp$fila - 0.5) * 0.05
emp$dist_centro_km <- round(111.32 * sqrt(((emp$Lon - emp$lon_pixel) * cos(emp$lat * pi / 180))^2 +
                                          (emp$lat - emp$lat_pixel)^2), 2)

# ---- 1.2 Estaciones que comparten píxel (reciben la MISMA curva satelital)
emp$pixel        <- paste0("F", emp$fila, "-C", emp$col)
emp$est_en_pixel <- ave(rep(1, nrow(emp)), emp$pixel, FUN = sum)
# Coordenadas no confiables (Hoyos, 17 nov): la triangulación por lluvia las ubicó mal (correlación
# con las estaciones de referencia <= 0,30). Lista PENDIENTE de confirmar con Hoyos; por ahora solo se
# marcan, no se excluyen.
NO_CONFIABLES <- c(5300800101, 5300900102, 5300600101, 5300000203, 2621500105, 5300000202,
                   5301300101, 5301000103, 5400000201, 2621500104, 5300900101)
emp$coord_no_confiable <- emp$STATION_ID %in% NO_CONFIABLES
guardar(emp, "tabla_f3_03_emparejamiento_estacion_pixel")

# ---- 1.3 Verificación: el satélite del archivo de péntadas es el del píxel asignado
verif_px <- do.call(rbind, lapply(1:5, function(k) {
  a <- p[sample(nrow(p), 1), ]
  r <- rast(file.path(DSAT, sprintf("ppt%d%02d.tif", a$YEAR, a$PENTADA)))
  q <- p[p$YEAR == a$YEAR & p$PENTADA == a$PENTADA, ]
  j <- match(q$STATION_ID, emp$STATION_ID)
  v <- r[cellFromRowCol(r, emp$fila[j], emp$col[j])][, 1]
  data.frame(anio = a$YEAR, pentada = a$PENTADA, estaciones = nrow(q), max_dif_mm = max(abs(v - q$satelite_mm)))
}))
print(verif_px)
guardar(verif_px, "tabla_f3_04_verificacion_extraccion_pixel")

# ---- 1.4 Resumen para el informe
resumen_emp <- data.frame(
  concepto = c("Estaciones", "Píxeles distintos", "Estaciones que comparten píxel",
               "Máximo de estaciones en un píxel", "Distancia mediana al centro del píxel (km)",
               "Coordenadas de catálogo", "Coordenadas estimadas por lluvia",
               "Coordenadas estimadas por lluvia + altitud", "Ubicadas solo por su código",
               "Altitud de la grilla (no de catálogo)"),
  valor = c(nrow(emp), length(unique(emp$pixel)), sum(emp$est_en_pixel > 1),
            max(emp$est_en_pixel), median(emp$dist_centro_km),
            sum(emp$fuente_coord == "catálogo"), sum(emp$fuente_coord == "lluvia"),
            sum(emp$fuente_coord == "lluvia + altitud"), sum(emp$fuente_coord == "ID"),
            sum(emp$fuente_alt != "catálogo")))
print(resumen_emp)
guardar(resumen_emp, "tabla_f3_05_resumen_emparejamiento")

figura("fig_f3_03_mapa_emparejamiento", {
  cols <- c("catálogo" = "black", "lluvia" = "orange", "lluvia + altitud" = "steelblue", "ID" = "red")
  par(mar = c(4.2, 4.2, 2.5, 1))
  plot(emp$Lon, emp$lat, type = "n", xlim = c(-77.55, -75.70), ylim = c(3.05, 5.00),
       asp = 1 / cos(4 * pi / 180), xlab = "Longitud", ylab = "Latitud",
       main = "Estaciones y malla del satélite (píxeles de 0,05°)")
  abline(v = -77.55 + 0.05 * (0:37), h = 5 - 0.05 * (0:39), col = "grey88")
  o <- order(emp$fuente_coord == "ID")                                  # las ubicadas por código, encima
  points(emp$Lon[o], emp$lat[o], pch = ifelse(emp$est_en_pixel[o] > 1, 17, 19),
         cex = ifelse(emp$fuente_coord[o] == "ID", 1.4, 0.9), col = cols[emp$fuente_coord[o]])
  legend("topleft", c(names(cols), "comparte píxel (triángulo)"), col = c(cols, "grey40"),
         pch = c(19, 19, 19, 19, 17), bg = "white", cex = 0.8)
}, ancho = 6.5, alto = 7.5)


#####################################################################################
# 2. Curva de diferencia d(t) = estación - satélite
#####################################################################################
# Signo: d(t) > 0 -> la estación midió más -> el satélite SUBESTIMA (y al revés).
# Escala principal: mm/día (la del sesgo). La escala sqrt se reporta como complemento.
#
# Unidad de análisis: la ESTACIÓN. Las curvas de una misma estación (hasta 23 años) comparten
# sitio y píxel y no son independientes; si se tratan como 3.829 curvas independientes la banda
# sale artificialmente angosta. Por eso primero se promedia d_i(t) dentro de cada estación y el
# error estándar se calcula entre estaciones: EE(t) = DE_entre_estaciones(t) / sqrt(n).
#
# Decisiones tomadas tras la primera corrida:
#   - La media de d(t) con todas las estaciones la dominan ~10 estaciones de la vertiente pacífica,
#     casi todas con coordenadas estimadas por lluvia, que quedaron en píxeles del lado del valle.
#     Por eso el análisis se separa por REGIÓN (valle del Cauca / vertiente pacífica), se agrega la
#     MEDIANA por estación y el % de estaciones con d > 0, y se repite solo con coordenadas de catálogo.
#   - Los primeros y últimos BORDE días del año se excluyen de los resúmenes (borde de los B-splines).

# ---- 2.1 Curvas de diferencia (resta de objetos fd, como en el script de clase)
Dmm_fd <- Ymm_fd - Xmm_fd
Dsq_fd <- Yfd - Xfd
Dmm <- eval.fd(t_eval, Dmm_fd)                           # 366 x N, mm/día
Dsq <- eval.fd(t_eval, Dsq_fd)                           # 366 x N, sqrt(mm/día)
stopifnot(max(abs(Dmm - (Ymm - Xmm))) < 1e-8)           # la resta de fd = resta de valores

# ---- 2.2 Promedio por estación, grupos y banda d_barra(t) +- 2 EE(t)
promedio_por_estacion <- function(M) {
  sapply(split(seq_len(N), claves$STATION_ID), function(j) rowMeans(M[, j, drop = FALSE]))
}
banda <- function(G) {
  m  <- rowMeans(G)
  ee <- apply(G, 1, sd) / sqrt(ncol(G))
  data.frame(t = t_eval, media = m, ee = ee, inf = m - 2 * ee, sup = m + 2 * ee,
             mediana = apply(G, 1, median),                    # resumen robusto
             pct_subestima = 100 * rowMeans(G > 0))            # % de estaciones con d > 0
}

G_dmm <- promedio_por_estacion(Dmm)                      # 366 x 244 (una columna por estación)
G_dsq <- promedio_por_estacion(Dsq)
G_ymm <- promedio_por_estacion(Ymm)
G_xmm <- promedio_por_estacion(Xmm)
n_est <- ncol(G_dmm)

ids <- colnames(G_dmm)
ie  <- match(ids, as.character(emp$STATION_ID)); stopifnot(!anyNA(ie))
info_est <- data.frame(STATION_ID = emp$STATION_ID[ie],
                       region = ifelse(substr(ids, 1, 2) == "26", "Cauca", "Pacífico"),   # zona hidrográfica
                       coord  = ifelse(emp$fuente_coord[ie] == "catálogo", "catálogo", "estimada"))

grupos <- list(
  "Todas"               = rep(TRUE, n_est),
  "Solo catálogo"       = info_est$coord == "catálogo",
  "Cauca"               = info_est$region == "Cauca",
  "Cauca - catálogo"    = info_est$region == "Cauca" & info_est$coord == "catálogo",
  "Pacífico"            = info_est$region == "Pacífico",
  "Pacífico - catálogo" = info_est$region == "Pacífico" & info_est$coord == "catálogo")
print(sapply(grupos, sum))

bandas <- lapply(grupos, function(g) banda(G_dmm[, g, drop = FALSE]))
b_sq   <- banda(G_dsq)

# ---- 2.3 Intervalos del año donde la banda no contiene el cero
meses <- c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")
dia_a_fecha <- function(t) {
  d <- as.Date(pmin(floor(t), 364), origin = "2019-01-01")
  paste(as.integer(format(d, "%d")), meses[as.integer(format(d, "%m"))])
}
BORDE <- 10   # días al inicio y al final del año afectados por el borde de los B-splines

intervalos <- function(b, analisis) {
  sentido <- ifelse(b$inf > 0, "satélite subestima (d > 0)",
             ifelse(b$sup < 0, "satélite sobreestima (d < 0)", ""))
  r <- rle(sentido)
  fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  k <- which(r$values != "")
  if (length(k) == 0) return(NULL)
  data.frame(analisis   = analisis,
             sentido    = r$values[k],
             desde      = dia_a_fecha(b$t[ini[k]]),
             hasta      = dia_a_fecha(b$t[fin[k]]),
             dias       = round(b$t[fin[k]] - b$t[ini[k]]),
             d_media    = sapply(k, function(q) round(mean(b$media[ini[q]:fin[q]]), 3)),
             d_mediana  = sapply(k, function(q) round(mean(b$mediana[ini[q]:fin[q]]), 3)),
             toca_borde = b$t[ini[k]] < BORDE | b$t[fin[k]] > 365 - BORDE)
}
tabla_int <- do.call(rbind, c(lapply(names(bandas), function(g) intervalos(bandas[[g]], paste("mm/día -", g))),
                             list(intervalos(b_sq, "sqrt(mm/día) - Todas"))))
print(tabla_int, row.names = FALSE)
guardar(tabla_int, "tabla_f3_06_intervalos_diferencia")

# ---- 2.4 Sesgo por temporada y grupo (sin los días de borde)
# Por estación: promedio de d(t) y de la lluvia observada en la temporada; luego media (con IC),
# mediana y % de estaciones que subestima. "Año completo" = sesgo del acumulado anual (ítem 4.4).
mes_t      <- as.integer(format(as.Date(t_eval[1:365], origin = "2019-01-01"), "%m"))
temporadas <- c("Dic-Feb (seca)", "Mar-May (lluvias)", "Jun-Ago (seca)", "Sep-Nov (lluvias)")
temp_t     <- temporadas[c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)[mes_t]]
interior   <- t_eval[1:365] >= BORDE & t_eval[1:365] <= 365 - BORDE

sesgo_temporada <- function(g, s) {
  f  <- interior & (if (s == "Año completo") TRUE else temp_t == s)
  dy <- colMeans(G_dmm[1:365, grupos[[g]], drop = FALSE][f, , drop = FALSE])   # una cifra por estación
  yy <- colMeans(G_ymm[1:365, grupos[[g]], drop = FALSE][f, , drop = FALSE])
  ee <- sd(dy) / sqrt(length(dy))
  data.frame(grupo = g, temporada = s, estaciones = length(dy),
             lluvia_estacion_mm_dia   = round(mean(yy), 2),
             d_media                  = round(mean(dy), 3),
             IC95_inf                 = round(mean(dy) - 2 * ee, 3),
             IC95_sup                 = round(mean(dy) + 2 * ee, 3),
             d_mediana                = round(median(dy), 3),
             sesgo_pct                = round(100 * mean(dy) / mean(yy), 1),
             pct_estaciones_subestima = round(100 * mean(dy > 0), 1))
}
tabla_temp <- do.call(rbind, lapply(names(grupos), function(g)
  do.call(rbind, lapply(c(temporadas, "Año completo"), function(s) sesgo_temporada(g, s)))))
print(tabla_temp, row.names = FALSE)
guardar(tabla_temp, "tabla_f3_07_sesgo_por_temporada_y_grupo")

# ---- 2.5 Estaciones influyentes: diferencia anual por estación
d_anual <- colMeans(G_dmm[1:365, , drop = FALSE][interior, , drop = FALSE])
y_anual <- colMeans(G_ymm[1:365, , drop = FALSE][interior, , drop = FALSE])
x_anual <- colMeans(G_xmm[1:365, , drop = FALSE][interior, , drop = FALSE])
influ <- data.frame(STATION_ID = info_est$STATION_ID, nombre = emp$STATION_NAME[ie], region = info_est$region,
                    fuente_coord = emp$fuente_coord[ie], confianza = emp$confianza[ie],
                    altitud = emp$altitud[ie], Lon = round(emp$Lon[ie], 3), lat = round(emp$lat[ie], 3),
                    pixel = emp$pixel[ie], curvas = emp$curvas[ie], coord_no_confiable = emp$coord_no_confiable[ie],
                    estacion_mm_dia = round(y_anual, 2), satelite_mm_dia = round(x_anual, 2),
                    d_mm_dia = round(d_anual, 2), razon_est_sat = round(y_anual / x_anual, 2))
influ <- influ[order(-abs(influ$d_mm_dia)), ]
influ$top10 <- seq_len(nrow(influ)) <= 10
print(head(influ, 15), row.names = FALSE)
cat("d anual media - todas:", round(mean(influ$d_mm_dia), 3),
    "| sin las 10 más influyentes:", round(mean(influ$d_mm_dia[!influ$top10]), 3),
    "| mediana:", round(median(influ$d_mm_dia), 3), "\n")
guardar(influ, "tabla_f3_08_diferencia_anual_por_estacion")

# ---- 2.6 Figuras
dibujar_banda <- function(b, ylab, main, col = "darkred", espacio = 0) {
  yl <- range(b$inf, b$sup, b$mediana, 0); yl[2] <- yl[2] + espacio * diff(yl)    # espacio para la leyenda
  plot(b$t, b$media, type = "n", ylim = yl, xlab = "Día del año",
       ylab = ylab, main = main, cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e6, c(BORDE, 365), 1e6, col = "grey93", border = NA)    # zona de borde
  polygon(c(b$t, rev(b$t)), c(b$inf, rev(b$sup)), col = adjustcolor(col, 0.25), border = NA)
  lines(b$t, b$media, lwd = 3, col = col)
  lines(b$t, b$mediana, lwd = 2, lty = 2, col = "black")
  abline(h = 0, lty = 3)
  y0 <- par("usr")[3] + 0.03 * diff(par("usr")[3:4])
  lines(b$t, ifelse(b$inf > 0, y0, NA), lwd = 6, col = "steelblue")                  # subestima
  lines(b$t, ifelse(b$sup < 0, y0, NA), lwd = 6, col = "orange")                     # sobreestima
}
leyenda_banda <- function(pos = "topright") {
  legend(pos, c("media", "mediana", "media ± 2 EE", "subestima", "sobreestima", "borde"),
         col = c("darkred", "black", adjustcolor("darkred", 0.25), "steelblue", "orange", "grey85"),
         lty = c(1, 2, 1, 1, 1, 1), lwd = c(3, 2, 8, 6, 6, 8), bty = "n", cex = 0.7)
}
comparar_bandas <- function(b1, b2, et1, et2, main) {
  plot(t_eval, b1$media, type = "n", ylim = range(b1$inf, b1$sup, b2$inf, b2$sup, 0),
       xlab = "Día del año", ylab = "mm/día", main = main, cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e6, c(BORDE, 365), 1e6, col = "grey93", border = NA)
  for (q in list(list(b1, "darkred"), list(b2, "darkgreen"))) {
    lines(t_eval, q[[1]]$media, lwd = 3, col = q[[2]])
    lines(t_eval, q[[1]]$inf, lty = 2, col = q[[2]]); lines(t_eval, q[[1]]$sup, lty = 2, col = q[[2]])
  }
  abline(h = 0, lty = 3)
  legend("topright", c(et1, et2), col = c("darkred", "darkgreen"), lwd = 3, bty = "n", cex = 0.75)
}
etiq <- function(g) paste0(g, " (n = ", sum(grupos[[g]]), ")")

figura("fig_f3_04_diferencia_mm", {
  par(mfrow = c(1, 2))
  col_reg <- ifelse(info_est$region == "Cauca", "grey75", adjustcolor("orange", 0.7))
  matplot(t_eval, G_dmm, type = "l", lty = 1, col = col_reg, xlab = "Día del año", ylab = "mm/día",
          main = "d(t) por estación (promedio de sus años)", cex.main = 0.95)
  lines(t_eval, bandas$Todas$media, lwd = 3, col = "darkred")
  lines(t_eval, bandas$Todas$mediana, lwd = 2, lty = 2)
  abline(h = 0, lty = 3)
  legend("topright", c("valle del Cauca", "vertiente pacífica", "media", "mediana"),
         col = c("grey60", "orange", "darkred", "black"), lty = c(1, 1, 1, 2), lwd = c(1, 1, 3, 2), bty = "n", cex = 0.75)
  dibujar_banda(bandas$Todas, "mm/día", paste("Todas las estaciones (n =", n_est, ")"), espacio = 0.6)
  leyenda_banda()
})

figura("fig_f3_05_diferencia_por_region", {
  par(mfrow = c(2, 2))
  dibujar_banda(bandas$Cauca, "mm/día", paste("Valle del Cauca:", etiq("Cauca")), espacio = 0.6)
  leyenda_banda()
  dibujar_banda(bandas$`Pacífico`, "mm/día", paste("Vertiente pacífica:", etiq("Pacífico")))
  comparar_bandas(bandas$Cauca, bandas$`Cauca - catálogo`, etiq("Cauca"), etiq("Cauca - catálogo"),
                  "Cauca: todas vs. solo catálogo")
  comparar_bandas(bandas$`Pacífico`, bandas$`Pacífico - catálogo`, etiq("Pacífico"), etiq("Pacífico - catálogo"),
                  "Pacífico: todas vs. solo catálogo")
}, ancho = 12, alto = 8.5)

figura("fig_f3_06_diferencia_sqrt_y_pct", {
  par(mfrow = c(1, 2))
  dibujar_banda(b_sq, "sqrt(mm/día)", "Todas, en sqrt(mm/día)")
  plot(t_eval, bandas$Cauca$pct_subestima, type = "l", lwd = 3, col = "grey40", ylim = c(0, 100),
       xlab = "Día del año", ylab = "% de estaciones", main = "Estaciones donde el satélite subestima", cex.main = 0.95)
  lines(t_eval, bandas$`Pacífico`$pct_subestima, lwd = 3, col = "orange")
  lines(t_eval, bandas$`Solo catálogo`$pct_subestima, lwd = 2, lty = 2, col = "darkgreen")
  abline(h = 50, lty = 3)
  legend("topright", c(etiq("Cauca"), etiq("Pacífico"), etiq("Solo catálogo")),
         col = c("grey40", "orange", "darkgreen"), lty = c(1, 1, 2), lwd = c(3, 3, 2), bty = "n", cex = 0.75)
})


#####################################################################################
# 3. Correlación cruzada funcional
#####################################################################################
# r_XY(s, t) = correlación, entre curvas, del satélite en el día s con la estación en el día t.
# La diagonal r_XY(t, t) dice si en la misma época del año las dos fuentes suben y bajan juntas.
# Dos versiones (escala mm/día):
#   - Total (como en el script de clase): incluye las diferencias ENTRE sitios. Una estación muy
#     lluviosa con un píxel muy lluvioso sube la correlación aunque el satélite no siga los años.
#   - Anomalías: a cada curva se le resta la curva media de su estación (y al satélite la de su
#     píxel-estación). Mide si el satélite sigue la variación de UN AÑO A OTRO en el mismo sitio,
#     que es lo que la red de la fase 4 debe corregir. Solo estaciones con >= 3 años.

# ---- 3.1 Anomalías y grupos de curvas
est_curva    <- match(as.character(claves$STATION_ID), colnames(G_dmm))     # columna de su estación
region_curva <- info_est$region[est_curva]
n_curvas_est <- emp$curvas[match(claves$STATION_ID, emp$STATION_ID)]
Xmm_an <- Xmm - G_xmm[, est_curva]
Ymm_an <- Ymm - G_ymm[, est_curva]
con_anom <- n_curvas_est >= 3

curvas_grupo <- list("Todas" = rep(TRUE, N), "Cauca" = region_curva == "Cauca",
                     "Pacífico" = region_curva == "Pacífico")

# ---- 3.2 Correlación cruzada paso a paso y verificación con cor.fd (todas, total)
cor_cruzada <- function(Xm, Ym) cor(t(Xm), t(Ym))       # filas = s (satélite), columnas = t (estación)
R_fd <- cor.fd(t_eval, Xmm_fd, t_eval, Ymm_fd)
cat("Correlación cruzada: máx. diferencia paso a paso vs cor.fd =",
    signif(max(abs(cor_cruzada(Xmm, Ymm) - R_fd)), 3), "\n")

analisis_cc <- list(); n_cc <- c()
for (g in names(curvas_grupo)) {
  k  <- curvas_grupo[[g]]
  ka <- k & con_anom
  analisis_cc[[paste(g, "- total")]]     <- cor_cruzada(Xmm[, k], Ymm[, k])
  analisis_cc[[paste(g, "- anomalías")]] <- cor_cruzada(Xmm_an[, ka], Ymm_an[, ka])
  n_cc[paste(g, "- total")] <- sum(k); n_cc[paste(g, "- anomalías")] <- sum(ka)
}
diagonales <- sapply(analisis_cc, diag)                  # 366 x 6

# Desfase: para cada día t de la estación, el día s (a +-30 días) donde el satélite correlaciona más.
# desfase > 0: el satélite "llega tarde" (se parece más a la estación de días anteriores).
desfase <- function(R, ventana = 30) {
  sapply(seq_along(t_eval), function(j) {
    i <- which(abs(t_eval - t_eval[j]) <= ventana)
    t_eval[i[which.max(R[i, j])]] - t_eval[j]
  })
}

# ---- 3.3 Resumen de la diagonal por temporada (sin días de borde)
resumen_cc <- do.call(rbind, lapply(names(analisis_cc), function(a) {
  dg <- diagonales[1:365, a]
  ds <- desfase(analisis_cc[[a]])[1:365]
  do.call(rbind, lapply(c(temporadas, "Año completo"), function(s) {
    f <- interior & (if (s == "Año completo") TRUE else temp_t == s)
    data.frame(analisis = a, temporada = s, curvas = n_cc[a],
               r_diag_media = round(mean(dg[f]), 3), r_diag_min = round(min(dg[f]), 3),
               r_diag_max = round(max(dg[f]), 3), desfase_medio_dias = round(mean(ds[f]), 1))
  }))
}))
print(resumen_cc, row.names = FALSE)
guardar(resumen_cc, "tabla_f3_09_correlacion_cruzada_diagonal")

# ---- 3.4 Figura
figura("fig_f3_07_correlacion_cruzada", {
  par(mfrow = c(1, 3))
  for (a in c("Cauca - total", "Cauca - anomalías")) {
    image(t_eval, t_eval, analisis_cc[[a]], zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
          xlab = "s: satélite (día del año)", ylab = "t: estación (día del año)",
          main = paste0("r(s,t) ", a, " (", n_cc[a], " curvas)"), cex.main = 0.95)
    contour(t_eval, t_eval, analisis_cc[[a]], levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
    abline(0, 1, lty = 3, col = "black")
  }
  cols_g <- c("Todas" = "black", "Cauca" = "grey50", "Pacífico" = "orange")
  matplot(t_eval, diagonales, type = "l", ylim = c(-0.2, 1), xlab = "Día del año", ylab = "r(t, t)",
          lty = rep(c(1, 2), 3), lwd = 2, col = rep(cols_g, each = 2),
          main = "Diagonal de la correlación cruzada", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -2, c(BORDE, 365), 2, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("bottomright", colnames(diagonales), col = rep(cols_g, each = 2), lty = rep(c(1, 2), 3),
         lwd = 2, bty = "n", cex = 0.7)
}, ancho = 15, alto = 5)


#####################################################################################
# 4. Comparación de las covarianzas de cada fuente
#####################################################################################
# ¿El satélite reproduce la variabilidad de la lluvia real o la suaviza? Se compara, en mm/día:
#   - la desviación estándar funcional DE(t) de cada fuente (diagonal de la covarianza),
#   - la varianza total integrada y cuánta concentra la primera componente (forma de la covarianza),
#   - la correlación de cada fuente consigo misma a 30 días (qué tan "suave" es en el tiempo).
# Total y anomalías, igual que en el fragmento 3.

# ---- 4.1 Covarianza con var.fd y verificación paso a paso (como en clase)
Cov_Y_fd <- eval.bifd(t_eval, t_eval, var.fd(Ymm_fd))
Cov_X_fd <- eval.bifd(t_eval, t_eval, var.fd(Xmm_fd))
cat("Covarianza: máx. diferencia paso a paso vs var.fd: estación =",
    signif(max(abs(cov(t(Ymm)) - Cov_Y_fd)), 3), "| satélite =", signif(max(abs(cov(t(Xmm)) - Cov_X_fd)), 3), "\n")

# ---- 4.2 Resumen por grupo
w <- rep(1, length(t_eval)); w[c(1, length(t_eval))] <- 0.5     # pesos del trapecio (paso de 1 día)
pct_componentes <- function(C) {
  v <- pmax(eigen(sqrt(w) * t(sqrt(w) * C), symmetric = TRUE, only.values = TRUE)$values, 0)
  100 * v / sum(v)
}
autocor30 <- function(C) {
  R <- cov2cor(C); i <- which(t_eval >= BORDE & t_eval + 30 <= 365 - BORDE)
  mean(R[cbind(i, i + 30)])
}

cov_grupo <- list()
for (g in names(curvas_grupo)) {
  k <- curvas_grupo[[g]]; ka <- k & con_anom
  cov_grupo[[paste(g, "- total")]]     <- list(Y = cov(t(Ymm[, k])),     X = cov(t(Xmm[, k])))
  cov_grupo[[paste(g, "- anomalías")]] <- list(Y = cov(t(Ymm_an[, ka])), X = cov(t(Xmm_an[, ka])))
}

tabla_cov <- do.call(rbind, lapply(names(cov_grupo), function(a) {
  CY <- cov_grupo[[a]]$Y; CX <- cov_grupo[[a]]$X
  deY <- sqrt(diag(CY))[1:365]; deX <- sqrt(diag(CX))[1:365]
  pY <- pct_componentes(CY); pX <- pct_componentes(CX)
  raz_temp <- sapply(temporadas, function(s) round(mean(deX[interior & temp_t == s] / deY[interior & temp_t == s]), 2))
  data.frame(analisis = a, curvas = n_cc[a],
             var_total_estacion = round(sum(w * diag(CY)), 0), var_total_satelite = round(sum(w * diag(CX)), 0),
             razon_var_total_sat_est = round(sum(w * diag(CX)) / sum(w * diag(CY)), 3),
             razon_DE_DicFeb = raz_temp[1], razon_DE_MarMay = raz_temp[2],
             razon_DE_JunAgo = raz_temp[3], razon_DE_SepNov = raz_temp[4],
             pct_comp1_estacion = round(pY[1], 1), pct_comp1_satelite = round(pX[1], 1),
             pct_3comp_estacion = round(sum(pY[1:3]), 1), pct_3comp_satelite = round(sum(pX[1:3]), 1),
             autocor30_estacion = round(autocor30(CY), 3), autocor30_satelite = round(autocor30(CX), 3),
             row.names = NULL)
}))
print(tabla_cov, row.names = FALSE)
guardar(tabla_cov, "tabla_f3_10_comparacion_covarianzas")

# ---- 4.3 Figura
figura("fig_f3_08_covarianzas", {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1))
  CY <- cov_grupo[["Cauca - total"]]$Y; CX <- cov_grupo[["Cauca - total"]]$X
  niveles <- pretty(range(CY, CX), 12)
  for (q in list(list(CY, "Covarianza estación - Cauca, total (mm²/día²)"),
                 list(CX, "Covarianza satélite - Cauca, total (misma escala)"))) {
    image(t_eval, t_eval, q[[1]], zlim = range(niveles), col = hcl.colors(40, "YlOrRd", rev = TRUE),
          xlab = "s (día)", ylab = "t (día)", main = q[[2]], cex.main = 0.95)
    contour(t_eval, t_eval, q[[1]], levels = niveles, labcex = 0.6, add = TRUE)
  }
  de <- sapply(c("Cauca - total", "Cauca - anomalías"), function(a)
    cbind(sqrt(diag(cov_grupo[[a]]$Y)), sqrt(diag(cov_grupo[[a]]$X))))
  de <- matrix(de, ncol = 4)
  matplot(t_eval, de, type = "l", lty = c(1, 2, 1, 2), lwd = c(3, 3, 1.5, 1.5),
          col = c("slategray", "darkturquoise", "slategray", "darkturquoise"), ylim = c(0, max(de)),
          xlab = "Día del año", ylab = "mm/día", main = "Desviación estándar funcional - Cauca", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = adjustcolor("grey90", 0.6), border = NA)
  legend("topright", c("estación total", "satélite total", "estación anomalías", "satélite anomalías"),
         col = c("slategray", "darkturquoise", "slategray", "darkturquoise"), lty = c(1, 2, 1, 2),
         lwd = c(3, 3, 1.5, 1.5), bty = "n", cex = 0.75)
  pc <- rbind(pct_componentes(cov_grupo[["Cauca - anomalías"]]$Y)[1:6],
              pct_componentes(cov_grupo[["Cauca - anomalías"]]$X)[1:6])
  barplot(pc, beside = TRUE, names.arg = 1:6, col = c("slategray", "darkturquoise"),
          xlab = "Componente", ylab = "% de varianza", main = "Componentes de la covarianza - Cauca, anomalías",
          cex.main = 0.95, legend.text = c("Estación", "Satélite"), args.legend = list(x = "topright", bty = "n"))
}, ancho = 12, alto = 9)


#####################################################################################
# 5. Sesgo por altitud y por fase ENSO
#####################################################################################
# 5A. Altitud: se analiza DENTRO del valle del Cauca, porque en todo el conjunto la altitud se
#     confunde con la región (casi todas las estaciones bajas son del Pacífico). Franjas iguales a la
#     fase 2 (Annie). Prueba: ANOVA funcional por permutación sobre las curvas promedio por estación
#     (script de clase), con F_max para controlar el error en todo el año.
# 5B. ENSO: fase de cada año según el ONI (NOAA CPC; tabla de Fabián, fase 4). Como una estación
#     aporta años de las tres fases, la comparación es PAREADA dentro de cada estación:
#     Delta_i(t) = d_i(t) en años Niño (o Niña) - d_i(t) en años neutros. Así se quita el efecto del sitio.

# ---- 5.0 Resumen de sesgo para cualquier grupo de estaciones (columnas = estaciones)
resumen_sesgo <- function(Gd, Gy, etiqueta) {
  do.call(rbind, lapply(c(temporadas, "Año completo"), function(s) {
    f  <- interior & (if (s == "Año completo") TRUE else temp_t == s)
    dy <- colMeans(Gd[1:365, , drop = FALSE][f, , drop = FALSE])
    yy <- colMeans(Gy[1:365, , drop = FALSE][f, , drop = FALSE])
    ee <- if (length(dy) > 1) sd(dy) / sqrt(length(dy)) else NA
    data.frame(grupo = etiqueta, temporada = s, estaciones = length(dy),
               lluvia_estacion_mm_dia   = round(mean(yy), 2),
               d_media                  = round(mean(dy), 3),
               IC95_inf                 = round(mean(dy) - 2 * ee, 3),
               IC95_sup                 = round(mean(dy) + 2 * ee, 3),
               d_mediana                = round(median(dy), 3),
               sesgo_pct                = round(100 * mean(dy) / mean(yy), 1),
               pct_estaciones_subestima = round(100 * mean(dy > 0), 1))
  }))
}
cauca <- info_est$region == "Cauca"

# ---- 5A.1 Franjas de altitud
cortes_alt     <- c(-Inf, 800, 1100, 1500, 2000, Inf)
nombres_franja <- c("< 800 m", "800-1100 m", "1100-1500 m", "1500-2000 m", "> 2000 m")
info_est$altitud <- emp$altitud[ie]
info_est$franja  <- cut(info_est$altitud, cortes_alt, labels = nombres_franja, right = FALSE)
print(table(info_est$region, info_est$franja))

franjas_cauca <- nombres_franja[sapply(nombres_franja, function(fr) sum(cauca & info_est$franja == fr) >= 3)]
tabla_alt <- do.call(rbind, c(
  lapply(franjas_cauca, function(fr) {
    k <- cauca & info_est$franja == fr
    resumen_sesgo(G_dmm[, k, drop = FALSE], G_ymm[, k, drop = FALSE], paste("Cauca", fr))
  }),
  list(resumen_sesgo(G_dmm[, !cauca], G_ymm[, !cauca], "Pacífico (todas las altitudes)"))))
print(tabla_alt[tabla_alt$temporada %in% c("Jun-Ago (seca)", "Año completo"), ], row.names = FALSE)
guardar(tabla_alt, "tabla_f3_11_sesgo_por_altitud")

# ---- 5A.2 Sesgo anual por estación vs. altitud y vs. lluvia de la estación
info_est$d_anual    <- d_anual
info_est$sesgo_pct  <- 100 * d_anual / y_anual
info_est$lluvia_est <- y_anual
info_est$razon      <- y_anual / x_anual
info_est$no_conf    <- info_est$STATION_ID %in% NO_CONFIABLES
cor_alt <- cor.test(info_est$altitud[cauca], info_est$sesgo_pct[cauca], method = "spearman", exact = FALSE)
# Compresión espacial: ¿el satélite reproduce las diferencias de lluvia entre sitios?
# Se ajusta log(satélite) ~ log(estación) con la lluvia media de cada estación (sin las no confiables).
# Pendiente = 1: el satélite reproduce las diferencias; pendiente < 1: las comprime (sobreestima los
# sitios secos y subestima los lluviosos). No se correlaciona la razón est/sat con la lluvia de la
# estación porque esa correlación es en parte espuria (la estación está en los dos lados).
ok_conf  <- !info_est$no_conf
fit_comp <- lm(log(x_anual[ok_conf]) ~ log(y_anual[ok_conf]))
ic_comp  <- confint(fit_comp)[2, ]
tabla_cor_alt <- data.frame(
  relacion   = c("Cauca: altitud vs. sesgo anual (%) - rho de Spearman",
                 "Todas sin no confiables: pendiente de log(satélite) ~ log(estación)"),
  estaciones = c(sum(cauca), sum(ok_conf)),
  valor      = round(c(cor_alt$estimate, coef(fit_comp)[2]), 3),
  IC95_inf   = c(NA, round(ic_comp[1], 3)),
  IC95_sup   = c(NA, round(ic_comp[2], 3)),
  valor_p    = signif(c(cor_alt$p.value, summary(fit_comp)$coefficients[2, 4]), 3))
print(tabla_cor_alt, row.names = FALSE)
guardar(tabla_cor_alt, "tabla_f3_12_correlacion_sesgo_altitud_lluvia")

# ---- 5A.3 ANOVA funcional por franja (Cauca), por permutación (script de clase ANOVA funcional)
anova_funcional <- function(X, grupo, B = 2000) {
  grupo <- droplevels(factor(grupo)); k <- nlevels(grupo); n <- ncol(X)
  Z  <- model.matrix(~ grupo - 1)
  Mg <- rowMeans(X); SST <- rowSums((X - Mg)^2)
  ssr <- function(Zb) { Mb <- X %*% Zb %*% diag(1 / colSums(Zb)); as.numeric((Mb - Mg)^2 %*% colSums(Zb)) }
  Fdet <- function(S) (S / (k - 1)) / ((SST - S) / (n - k))
  F_t <- Fdet(ssr(Z))
  Fmax_nulos <- replicate(B, max(Fdet(ssr(Z[sample(n), ]))))
  list(F_t = F_t, c_max = unname(quantile(Fmax_nulos, 0.95)),
       p_global = (1 + sum(Fmax_nulos >= max(F_t))) / (B + 1), k = k, n = n)
}
# Con B permutaciones el menor valor p posible es 1/(B+1); si se alcanza se reporta como "< 0,001"
texto_p <- function(p, B = 2000) if (p <= 1 / (B + 1) + 1e-12) "< 0,001" else as.character(signif(p, 3))
t_int <- t_eval[1:365][interior]
X_alt <- G_dmm[1:365, cauca, drop = FALSE][interior, , drop = FALSE]
g_alt <- info_est$franja[cauca]
mantener <- g_alt %in% franjas_cauca
an_alt <- anova_funcional(X_alt[, mantener], g_alt[mantener])
sig <- an_alt$F_t > an_alt$c_max
r <- rle(sig); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1; kk <- which(r$values)
tabla_anova_alt <- data.frame(prueba = "ANOVA funcional d(t) ~ franja de altitud (Cauca)",
                              grupos = an_alt$k, estaciones = an_alt$n, F_max_obs = round(max(an_alt$F_t), 2),
                              F_critico_max_5pct = round(an_alt$c_max, 2), p_global = texto_p(an_alt$p_global),
                              dias_significativos = if (length(kk) == 0) "ninguno" else
                                paste(paste(dia_a_fecha(t_int[ini[kk]]), "-", dia_a_fecha(t_int[fin[kk]])), collapse = "; "))
print(tabla_anova_alt, row.names = FALSE)
guardar(tabla_anova_alt, "tabla_f3_13_anova_altitud")

# ---- 5B.1 Fase ENSO de cada curva
enso <- read.csv(file.path(".", "salidas_fase4", "fase_enso_por_anio.csv"), fileEncoding = "UTF-8")
fase_curva <- enso$fase_enso[match(claves$YEAR, enso$YR)]
stopifnot(!anyNA(fase_curva))
fases <- c("Niño", "Neutro", "Niña")
print(tapply(enso$YR, enso$fase_enso, paste, collapse = ", "))

promedio_sub <- function(M, sel) {
  sapply(split(which(sel), claves$STATION_ID[sel]), function(j) rowMeans(M[, j, drop = FALSE]))
}
G_fase_d <- lapply(setNames(fases, fases), function(f) promedio_sub(Dmm, fase_curva == f))
G_fase_y <- lapply(setNames(fases, fases), function(f) promedio_sub(Ymm, fase_curva == f))
es_cauca <- function(G) substr(colnames(G), 1, 2) == "26"

tabla_enso <- do.call(rbind, c(
  lapply(fases, function(f) resumen_sesgo(G_fase_d[[f]][, es_cauca(G_fase_d[[f]]), drop = FALSE],
                                          G_fase_y[[f]][, es_cauca(G_fase_y[[f]]), drop = FALSE], paste("Cauca -", f))),
  lapply(fases, function(f) resumen_sesgo(G_fase_d[[f]][, !es_cauca(G_fase_d[[f]]), drop = FALSE],
                                          G_fase_y[[f]][, !es_cauca(G_fase_y[[f]]), drop = FALSE], paste("Pacífico -", f)))))
print(tabla_enso[tabla_enso$temporada %in% c("Jun-Ago (seca)", "Año completo"), ], row.names = FALSE)
guardar(tabla_enso, "tabla_f3_14_sesgo_por_fase_enso")

# ---- 5B.2 Comparación pareada dentro de cada estación (Cauca): Niño - Neutro y Niña - Neutro
pareada <- function(f) {
  comunes <- intersect(colnames(G_fase_d[[f]]), colnames(G_fase_d[["Neutro"]]))
  comunes <- comunes[substr(comunes, 1, 2) == "26"]
  list(D = G_fase_d[[f]][, comunes] - G_fase_d[["Neutro"]][, comunes], Y = G_fase_y[["Neutro"]][, comunes])
}
par_nino <- pareada("Niño"); par_nina <- pareada("Niña")
b_nino <- banda(par_nino$D); b_nina <- banda(par_nina$D)
tabla_par <- rbind(resumen_sesgo(par_nino$D, par_nino$Y, "Cauca: cambio del sesgo Niño - Neutro"),
                   resumen_sesgo(par_nina$D, par_nina$Y, "Cauca: cambio del sesgo Niña - Neutro"))
names(tabla_par)[names(tabla_par) == "d_media"]   <- "cambio_d_media"
names(tabla_par)[names(tabla_par) == "sesgo_pct"] <- "cambio_pct_de_lluvia_neutra"
print(tabla_par, row.names = FALSE)
guardar(tabla_par, "tabla_f3_15_cambio_sesgo_enso_pareado")

# ¿Cuál fuente responde al ENSO? Cambio de la lluvia de cada fuente respecto a los años neutros,
# pareado por estación (Cauca, año completo sin bordes). Satélite = estación - d.
respuesta_enso <- function(f) {
  com <- intersect(colnames(G_fase_y[[f]]), colnames(G_fase_y[["Neutro"]])); com <- com[substr(com, 1, 2) == "26"]
  cy  <- colMeans((G_fase_y[[f]][, com] - G_fase_y[["Neutro"]][, com])[1:365, ][interior, ])
  cx  <- colMeans(((G_fase_y[[f]][, com] - G_fase_d[[f]][, com]) -
                   (G_fase_y[["Neutro"]][, com] - G_fase_d[["Neutro"]][, com]))[1:365, ][interior, ])
  ic  <- function(v) paste0("[", round(mean(v) - 2 * sd(v) / sqrt(length(v)), 3), "; ",
                            round(mean(v) + 2 * sd(v) / sqrt(length(v)), 3), "]")
  data.frame(comparacion = paste(f, "- Neutro"), estaciones = length(com),
             cambio_estacion_mm_dia = round(mean(cy), 3), IC95_estacion = ic(cy),
             cambio_satelite_mm_dia = round(mean(cx), 3), IC95_satelite = ic(cx))
}
tabla_resp_enso <- rbind(respuesta_enso("Niño"), respuesta_enso("Niña"))
print(tabla_resp_enso, row.names = FALSE)
guardar(tabla_resp_enso, "tabla_f3_17_respuesta_enso_por_fuente")
tabla_int_enso <- rbind(intervalos(b_nino, "Cauca: Niño - Neutro"), intervalos(b_nina, "Cauca: Niña - Neutro"))
print(tabla_int_enso, row.names = FALSE)
guardar(tabla_int_enso, "tabla_f3_16_intervalos_enso")

# ---- 5.4 Figuras
cols_fr <- setNames(hcl.colors(length(franjas_cauca), "Viridis"), franjas_cauca)
figura("fig_f3_09_sesgo_altitud", {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1))
  med_fr <- sapply(franjas_cauca, function(fr) rowMeans(G_dmm[, cauca & info_est$franja == fr, drop = FALSE]))
  yl <- range(med_fr); yl[2] <- yl[2] + 0.6 * diff(yl)
  matplot(t_eval, med_fr, type = "l", lty = 1, lwd = 2.5, col = cols_fr, xlab = "Día del año", ylab = "mm/día",
          ylim = yl, main = "Cauca: d(t) medio por franja de altitud", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("topright", paste0(franjas_cauca, " (n = ", sapply(franjas_cauca, function(fr) sum(cauca & info_est$franja == fr)), ")"),
         col = cols_fr, lwd = 2.5, bty = "n", cex = 0.75, ncol = 2)
  plot(t_int, an_alt$F_t, type = "l", lwd = 2, xlab = "Día del año", ylab = "F(t)",
       main = paste0("ANOVA funcional por franja (p global ", texto_p(an_alt$p_global), ")"), cex.main = 0.95,
       ylim = c(0, max(an_alt$F_t, an_alt$c_max) * 1.45))
  abline(h = an_alt$c_max, lty = 2, col = "red")
  legend("topright", c("F(t) observado", "valor crítico máximo (5 %)"), col = c("black", "red"), lty = 1:2, bty = "n", cex = 0.75)
  plot(info_est$altitud, info_est$sesgo_pct, pch = ifelse(info_est$no_conf, 4, 19), cex = 0.8,
       col = ifelse(cauca, "grey40", "orange"), xlab = "Altitud (m)", ylab = "Sesgo anual (%)",
       main = paste0("Sesgo anual por estación (Cauca: rho = ", round(cor_alt$estimate, 2), ")"), cex.main = 0.95)
  abline(h = 0, lty = 3)
  legend("bottomright", c("valle del Cauca", "vertiente pacífica", "coordenadas no confiables"),
         col = c("grey40", "orange", "black"), pch = c(19, 19, 4), bg = "white", cex = 0.75)
  plot(y_anual, x_anual, log = "xy", pch = ifelse(info_est$no_conf, 4, 19), cex = 0.8,
       col = ifelse(cauca, "grey40", "orange"), xlab = "Lluvia media de la estación (mm/día, log)",
       ylab = "Lluvia media del satélite en su píxel (mm/día, log)",
       main = paste0("Compresión espacial: pendiente = ", round(coef(fit_comp)[2], 2), " (1 = sin compresión)"),
       cex.main = 0.95)
  abline(0, 1, lty = 3)
  xs <- range(y_anual[ok_conf]); lines(xs, exp(coef(fit_comp)[1]) * xs^coef(fit_comp)[2], lwd = 2, col = "darkred")
  legend("topleft", c("1:1 (sin sesgo)", "ajuste sin no confiables"), lty = c(3, 1), lwd = c(1, 2),
         col = c("black", "darkred"), bty = "n", cex = 0.75)
}, ancho = 12, alto = 9)

cols_fase <- c("Niño" = "firebrick", "Neutro" = "grey40", "Niña" = "steelblue")
figura("fig_f3_10_sesgo_enso", {
  par(mfrow = c(1, 3))
  med_fase <- sapply(fases, function(f) rowMeans(G_fase_d[[f]][, es_cauca(G_fase_d[[f]]), drop = FALSE]))
  matplot(t_eval, med_fase, type = "l", lty = 1, lwd = 2.5, col = cols_fase, xlab = "Día del año", ylab = "mm/día",
          main = "Cauca: d(t) medio por fase ENSO", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("topright", fases, col = cols_fase, lwd = 2.5, bty = "n", cex = 0.8)
  dibujar_banda(b_nino, "mm/día", paste0("Cauca: Niño - Neutro, pareado (n = ", ncol(par_nino$D), ")"), col = "firebrick")
  dibujar_banda(b_nina, "mm/día", paste0("Cauca: Niña - Neutro, pareado (n = ", ncol(par_nina$D), ")"), col = "steelblue")
}, ancho = 15, alto = 5)


#####################################################################################
# 6. Regresión funcional función-función (modelo de referencia de la fase 4)
#####################################################################################
# Y_i(t) = beta0(t) + ∫ beta1(s, t) X_i(s) ds + e_i(t)       (script de clase: linmod)
#
# Para que sea comparable con las redes de la fase 4 se usan EXACTAMENTE sus entradas y su métrica:
#   - X = curva satelital de la fase 1 (K = 12, sqrt(mm/día)); Y = curva de la estación de la fase 1
#     (K = 8, sqrt(mm), diaria). Mismo orden de curvas y misma partición.
#   - Error L² relativo = 100 * sqrt( suma ||y_hat - y||² / suma ||y||² ), con la matriz de Gram de la
#     base de las estaciones (igual que error_L2 y L2_relativo en fase4_red_funcional.py).
# beta1 en base B-spline (12 en s x 8 en t), penalización de la segunda derivada en s y en t con un
# lambda común elegido con VALIDACIÓN (1996-2000). La prueba (2001-2003) solo se usa para reportar.
# Diferencia con la cabeza lineal de Fabián: linmod no admite covariables escalares (altitud, ONI).

# ---- 6.1 Curvas de la fase 1
cs1   <- read.csv(file.path(ENT1, "coeficientes_satelite.csv"))
info1 <- read.csv(file.path(ENT1, "base_info.csv"), fileEncoding = "UTF-8")
stopifnot(all(cs1$STATION_ID == ce$STATION_ID), all(cs1$YEAR == ce$YEAR))
K_e1 <- info1$K[info1$fuente == "estaciones"]; K_s1 <- info1$K[info1$fuente != "estaciones"]
base_e1 <- create.bspline.basis(c(0, 365), nbasis = K_e1, norder = 4)
base_s1 <- create.bspline.basis(c(0, 365), nbasis = K_s1, norder = 4)
C_e1 <- t(as.matrix(ce[,  paste0("c", 1:K_e1)]))        # K_e1 x N
C_s1 <- t(as.matrix(cs1[, paste0("c", 1:K_s1)]))        # K_s1 x N
Yfd1 <- fd(C_e1, base_e1)
Xfd1 <- fd(C_s1, base_s1)
i_tr <- which(claves$particion == "entrenamiento")
i_va <- which(claves$particion == "validacion")
i_te <- which(claves$particion == "prueba")
cat("Curvas: entrenamiento", length(i_tr), "| validación", length(i_va), "| prueba", length(i_te), "\n")

# ---- 6.2 Matrices de Gram y error L² (misma métrica de la fase 4)
W_e  <- inprod(base_e1, base_e1)
W_s  <- inprod(base_s1, base_s1)
W_se <- inprod(base_s1, base_e1)                         # K_s1 x K_e1 (Gram cruzada)
err_curvas <- function(Chat, idx) { E <- Chat - C_e1[, idx, drop = FALSE]; colSums(E * (W_e %*% E)) }
norma_y    <- function(idx) colSums(C_e1[, idx, drop = FALSE] * (W_e %*% C_e1[, idx, drop = FALSE]))
L2rel      <- function(err, idx) 100 * sqrt(sum(err) / sum(norma_y(idx)))
err_sat <- function(idx) {                               # satélite SIN corregir (bases distintas)
  cs_ <- C_s1[, idx, drop = FALSE]; ce_ <- C_e1[, idx, drop = FALSE]
  colSums(cs_ * (W_s %*% cs_)) - 2 * colSums(cs_ * (W_se %*% ce_)) + colSums(ce_ * (W_e %*% ce_))
}
# Control de consistencia con la fase 4 (Fabián reporta 65,66 % en validación y 84,48 % en prueba)
cat("Satélite sin corregir: validación", round(L2rel(err_sat(i_va), i_va), 2),
    "% | prueba", round(L2rel(err_sat(i_te), i_te), 2), "%\n")

# ---- 6.3 Ajuste con linmod y predicción paso a paso
# Con beta1(s,t) = phi_s(s)' B phi_t(t) y X(s) = phi_s(s)' c:  ∫ X(s) beta1(s,t) ds = c' W_s B phi_t(t),
# así que los coeficientes de y_hat en la base de las estaciones son  b0 + B' W_s c.
ajustar_ff <- function(lam, idx = i_tr) {
  b0Par <- fdPar(base_e1, int2Lfd(2), 0)
  b1Par <- bifdPar(bifd(matrix(0, K_s1, K_e1), base_s1, base_e1), int2Lfd(2), int2Lfd(2), lam, lam)
  linmod(Xfd1[idx], Yfd1[idx], list(b0Par, b1Par))
}
predecir_ff <- function(fit, idx) {
  B  <- fit$beta1estbifd$coefs                         # K_s1 x K_e1
  b0 <- as.vector(fit$beta0estfd$coefs)                # K_e1
  b0 + t(B) %*% W_s %*% C_s1[, idx, drop = FALSE]      # K_e1 x n
}

# ---- 6.4 Elección de lambda con validación (y lambda = 0 como referencia sin penalizar)
loglam <- seq(-4, 14, by = 1)       # hasta 1e14: en la primera corrida el óptimo cayó en el borde (1e8)
tabla_lam <- do.call(rbind, lapply(c(NA, loglam), function(l) {
  fit <- ajustar_ff(if (is.na(l)) 0 else 10^l)
  data.frame(log10_lambda = l,
             L2rel_entrenamiento = round(L2rel(err_curvas(predecir_ff(fit, i_tr), i_tr), i_tr), 3),
             L2rel_validacion    = round(L2rel(err_curvas(predecir_ff(fit, i_va), i_va), i_va), 3))
}))
tabla_lam$log10_lambda[is.na(tabla_lam$log10_lambda)] <- -Inf      # lambda = 0
print(tabla_lam, row.names = FALSE)
guardar(tabla_lam, "tabla_f3_18_regresion_eleccion_lambda")
fila_opt <- which.min(tabla_lam$L2rel_validacion)
lam_opt  <- if (is.infinite(tabla_lam$log10_lambda[fila_opt])) 0 else 10^tabla_lam$log10_lambda[fila_opt]
cat("lambda elegido con validación:", lam_opt, "\n")

fit_ff <- ajustar_ff(lam_opt)
C_hat_ff <- predecir_ff(fit_ff, seq_len(N))              # todas las curvas (para exportar)

# Verificación: la predicción paso a paso coincide con el ajuste de linmod (entrenamiento)
dif_linmod <- max(abs(eval.fd(t_eval, fit_ff$yhatfdobj) - eval.fd(t_eval, fd(C_hat_ff[, i_tr], base_e1))))
cat("Predicción paso a paso vs linmod (entrenamiento): máx. diferencia =", signif(dif_linmod, 3), "\n")

# ---- 6.5 Comparación con la fase 4 (mismas curvas, misma métrica)
leer_pred <- function(tipo, s) {
  P <- read.csv(file.path(".", "salidas_fase4", sprintf("pred_coef_%s_semilla%d.csv", tipo, s)))
  stopifnot(all(P$STATION_ID == claves$STATION_ID), all(P$YEAR == claves$YEAR))
  t(as.matrix(P[, paste0("c", 1:K_e1)]))
}
err_todos <- list(satelite = err_sat(seq_len(N)), regresion = err_curvas(C_hat_ff, seq_len(N)))
for (tipo in c("lineal", "kan", "mlp")) {
  e_sem <- sapply(0:2, function(s) err_curvas(leer_pred(tipo, s), seq_len(N)))      # N x 3
  err_todos[[tipo]] <- rowMeans(e_sem)                    # promedio de las 3 semillas (como Fabián)
  attr(err_todos[[tipo]], "sd_val")  <- sd(apply(e_sem, 2, function(e) L2rel(e[i_va], i_va)))
  attr(err_todos[[tipo]], "sd_test") <- sd(apply(e_sem, 2, function(e) L2rel(e[i_te], i_te)))
}
nombres_mod <- c(satelite = "satélite sin corregir", regresion = "regresión funcional (linmod)",
                 lineal = "red - cabeza lineal", kan = "red - cabeza KAN", mlp = "red - cabeza MLP")
tabla_mod <- data.frame(modelo = nombres_mod[names(err_todos)],
  L2rel_validacion = sapply(err_todos, function(e) round(L2rel(e[i_va], i_va), 2)),
  L2rel_prueba     = sapply(err_todos, function(e) round(L2rel(e[i_te], i_te), 2)),
  sd_semillas_prueba = sapply(err_todos, function(e) if (is.null(attr(e, "sd_test"))) NA else round(attr(e, "sd_test"), 2)),
  row.names = NULL)
# % de valores diarios negativos en prueba (escala sqrt): se recortan a 0 al pasar a mm (rubrica 4.2)
P_dia <- eval.basis(t_eval[1:365] + 0.5, base_e1)                    # 365 x K_e1 (centro de cada día)
pct_neg <- function(Chat) round(100 * mean(P_dia %*% Chat[, i_te, drop = FALSE] < 0), 2)
tabla_mod$pct_negativos_prueba <- c(NA, pct_neg(C_hat_ff), pct_neg(leer_pred("lineal", 0)),
                                    pct_neg(leer_pred("kan", 0)), pct_neg(leer_pred("mlp", 0)))
print(tabla_mod, row.names = FALSE)
guardar(tabla_mod, "tabla_f3_19_regresion_vs_redes")

# IC bootstrap por estación (prueba) y Wilcoxon pareado, con el mismo esquema de la fase 4
est_te <- claves$STATION_ID[i_te]
sum_est <- function(e) tapply(e[i_te], est_te, sum)
E_est  <- lapply(err_todos, sum_est)
nor_te <- tapply(norma_y(i_te), est_te, sum)
cont   <- t(replicate(2000, tabulate(sample(length(nor_te), replace = TRUE), nbins = length(nor_te))))
l2b    <- function(E) 100 * sqrt(as.vector(cont %*% E) / as.vector(cont %*% nor_te))
# Además del agregado (bootstrap) se reporta lo que pasa en la estación TÍPICA (mediana de la diferencia
# de error relativo por estación y % de estaciones donde gana el modelo a), y se repite por subconjunto:
# así se ve si una ganancia viene de toda la región o de unas pocas estaciones (p. ej. las no confiables).
ids_te <- as.numeric(names(nor_te))
subconj <- list("Todas"            = rep(TRUE, length(nor_te)),
                "Sin no confiables" = !(ids_te %in% NO_CONFIABLES),
                "Cauca"            = substr(names(nor_te), 1, 2) == "26",
                "Pacífico"         = substr(names(nor_te), 1, 2) != "26")
comparar <- function(a, b, sub = "Todas") {
  k  <- subconj[[sub]]
  Ea <- E_est[[a]][k]; Eb <- E_est[[b]][k]; nt <- nor_te[k]; cb <- cont[, k, drop = FALSE]
  dif_b <- 100 * (sqrt(as.vector(cb %*% Ea) / as.vector(cb %*% nt)) - sqrt(as.vector(cb %*% Eb) / as.vector(cb %*% nt)))
  ic <- quantile(dif_b, c(0.025, 0.975))
  ra <- sqrt(Ea / nt); rb <- sqrt(Eb / nt)
  data.frame(comparacion = paste(nombres_mod[a], "-", nombres_mod[b]), estaciones = sub, n = sum(k),
             dif_L2rel_pp = round(100 * (sqrt(sum(Ea) / sum(nt)) - sqrt(sum(Eb) / sum(nt))), 2),
             IC95_inf = round(ic[1], 2), IC95_sup = round(ic[2], 2), IC_excluye_0 = ic[1] > 0 | ic[2] < 0,
             mediana_dif_por_estacion_pp = round(100 * median(ra - rb), 2),
             pct_estaciones_gana_a = round(100 * mean(ra < rb), 1),
             wilcoxon_p = signif(wilcox.test(ra, rb, paired = TRUE)$p.value, 3), row.names = NULL)
}
tabla_comp_ff <- do.call(rbind, lapply(names(subconj), function(sb)
  rbind(comparar("regresion", "satelite", sb), comparar("lineal", "regresion", sb),
        comparar("kan", "regresion", sb), comparar("mlp", "regresion", sb))))
print(tabla_comp_ff, row.names = FALSE)
guardar(tabla_comp_ff, "tabla_f3_20_regresion_comparaciones_bootstrap")

# ---- 6.6 R²(t) en prueba (como en el script de clase)
y_te  <- eval.fd(t_eval, Yfd1[i_te])
R2_t  <- function(Chat) { yh <- eval.fd(t_eval, fd(Chat, base_e1)); 1 - rowSums((y_te - yh)^2) / rowSums((y_te - rowMeans(y_te))^2) }
R2_ff  <- R2_t(C_hat_ff[, i_te])
R2_lin <- R2_t(leer_pred("lineal", 0)[, i_te])
R2_kan <- R2_t(leer_pred("kan", 0)[, i_te])
tabla_R2 <- data.frame(modelo = c("regresión funcional", "red - cabeza lineal (semilla 0)", "red - cabeza KAN (semilla 0)"),
                       R2_medio = round(c(mean(R2_ff[interior]), mean(R2_lin[interior]), mean(R2_kan[interior])), 3),
                       R2_min = round(c(min(R2_ff[interior]), min(R2_lin[interior]), min(R2_kan[interior])), 3),
                       R2_max = round(c(max(R2_ff[interior]), max(R2_lin[interior]), max(R2_kan[interior])), 3))
print(tabla_R2, row.names = FALSE)
guardar(tabla_R2, "tabla_f3_21_regresion_R2_prueba")

# ---- 6.7 Predicciones exportadas con el formato de la fase 4 (para la evaluación 4.4 de Hoyos)
pred_ff <- data.frame(claves[, c("STATION_ID", "YEAR", "particion")], t(C_hat_ff))
names(pred_ff)[-(1:3)] <- paste0("c", 1:K_e1)
write.csv(pred_ff, file.path(OUT, "pred_coef_regresion_funcional.csv"), row.names = FALSE)

# ---- 6.8 Figuras
semanas <- seq(0, 365, length.out = 53)
B1mat   <- eval.bifd(semanas, semanas, fit_ff$beta1estbifd)
j_med   <- i_te[which.min(abs(err_todos$regresion[i_te] - median(err_todos$regresion[i_te])))]   # curva típica
figura("fig_f3_11_regresion_funcional", {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1))
  persp(semanas, semanas, B1mat, theta = 30, phi = 30, col = "lightblue", border = "grey40",
        xlab = "s (satélite)", ylab = "t (estación)", zlab = "beta1(s,t)", main = "Superficie beta1(s, t)",
        ticktype = "simple", cex.main = 0.95)
  image(semanas, semanas, B1mat, col = hcl.colors(40, "Blue-Red 3"),
        zlim = c(-1, 1) * max(abs(B1mat)), xlab = "s: día del satélite", ylab = "t: día de la estación",
        main = paste0("beta1(s, t), lambda = ", signif(lam_opt, 2)), cex.main = 0.95)
  contour(semanas, semanas, B1mat, nlevels = 12, labcex = 0.6, add = TRUE)
  abline(0, 1, lty = 2)
  matplot(t_eval, cbind(R2_ff, R2_lin, R2_kan), type = "l", lty = c(1, 2, 2), lwd = c(3, 2, 2),
          col = c("darkred", "grey40", "darkgreen"), ylim = c(min(0, R2_ff, R2_lin, R2_kan), 1),
          xlab = "Día del año", ylab = expression(R^2*(t)), main = "R²(t) en prueba (2001-2003)", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -10, c(BORDE, 365), 10, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("bottomright", c("regresión funcional", "red lineal", "red KAN"), col = c("darkred", "grey40", "darkgreen"),
         lty = c(1, 2, 2), lwd = c(3, 2, 2), bty = "n", cex = 0.8)
  yv <- eval.fd(t_eval, Yfd1[j_med]); xv <- eval.fd(t_eval, Xfd1[j_med])
  hv <- eval.fd(t_eval, fd(C_hat_ff[, j_med, drop = FALSE], base_e1))
  yl <- range(yv, xv, hv); yl[2] <- yl[2] + 0.25 * diff(yl)
  matplot(t_eval, cbind(yv, xv, hv), type = "l", lty = c(1, 2, 1), lwd = c(3, 2, 2), ylim = yl,
          col = c("slategray", "darkturquoise", "darkred"), xlab = "Día del año", ylab = "sqrt(mm/día)",
          main = paste("Curva típica de prueba:", claves$STATION_ID[j_med], claves$YEAR[j_med]), cex.main = 0.95)
  legend("top", c("estación", "satélite sin corregir", "regresión funcional"), ncol = 3,
         col = c("slategray", "darkturquoise", "darkred"), lty = c(1, 2, 1), lwd = c(3, 2, 2), bty = "n", cex = 0.8)
}, ancho = 12, alto = 9)
