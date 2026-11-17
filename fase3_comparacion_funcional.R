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
#   5. Sesgo por altitud / región y por fase ENSO.                             (pendiente)
#   6. Regresión funcional función-función (modelo de referencia).             (pendiente)
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
  contour(t_eval, t_eval, CY, levels = niveles, labcex = 0.6, xlab = "s (día)", ylab = "t (día)",
          main = "Covarianza estación - Cauca, total (mm²/día²)", cex.main = 0.95)
  contour(t_eval, t_eval, CX, levels = niveles, labcex = 0.6, xlab = "s (día)", ylab = "t (día)",
          main = "Covarianza satélite - Cauca, total (mismos niveles)", cex.main = 0.95)
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
          cex.main = 0.95, legend.text = c("Estación", "Satélite"), args.legend = list(bty = "n"))
}, ancho = 12, alto = 9)
