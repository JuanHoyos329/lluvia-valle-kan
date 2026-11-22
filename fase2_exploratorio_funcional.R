
# Análisis exploratorio de datos funcionales 
# Proyecto final de Analítica de Datos: lluvia del Valle del Cauca, 1981-2003
#
# Nos guiamos con la estructura del Script_AnalisisExploratorio.R del profesor Johann Ospina Galindez, es una adaptación 
# de lo visto en clase a los requerimientos del presente proyecto. 

# En este script vamos a manejar la siguiente estructura:
#   0. Carga: coeficientes de la fase 1 -> objetos fd (estaciones y satélite)
#   1. Representatividad: se usan TODAS las curvas de 1981-2003 (tope de la rúbrica)
#   2. Media funcional: paso a paso 
#   3. Varianza, desviación estándar y coeficiente de variación funcional
#   4. Covarianza y correlación funcional de cada fuente (paso a paso y con var.fd / cor.fd)
#   5. Derivadas
#   6. Profundidad (Fraiman-Muniz, modal y MBD), mediana, media recortada y boxplot funcional
#      - curvas crudas      -> atípicos explicados por la ubicación (orografía)
#      - anomalías por sitio -> atípicos explicados por el año (ENSO)
#   7. Orografía (franjas de altitud) y ENSO (fase del año)
#   8. Tabla resumen
#
# Entradas (rutas relativas a la carpeta del repositorio):
#   salidas_fase1/coeficientes_estaciones.csv, coeficientes_satelite.csv, base_info.csv 
#   salidas_preprocesamiento/ranking_estaciones.xlsx (coordenadas y altitud), estacion_pixel.xlsx,
#                            estaciones_diario_limpio.csv, estaciones_pentadas_limpio.csv 
#   Datos/oni_cpc.csv (ONI de la NOAA CPC). La fase ENSO de cada año se calcula aquí con la misma
#     regla de la fase 4, para que el script de R no dependa de una salida de Python; si existe
#     salidas_fase4/fase_enso_por_anio.csv, se comprueba que las dos clasificaciones coinciden.

#Limpiar área de trabajo
rm(list=ls())

# Instalar los paquetes solo si faltan (fda.usc necesita dependencias como SuppDists y kSamples)
for (pq in c("fda", "fda.usc", "readxl")) {
  if (!requireNamespace(pq, quietly = TRUE)) install.packages(pq, dependencies = TRUE)
}
library(fda)
library(fda.usc)
library(readxl)

ENT0 <- "./salidas_preprocesamiento"  
ENT1 <- "./salidas_fase1"              
ENT4 <- "./salidas_fase4"              
DATOS <- "./Datos"                     # ONI de la NOAA CPC

#Carpeta para guardar las figuras y tablas 
OUT  <- "./salidas_fase2"
FIG  <- file.path(OUT, "figuras")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

GUARDAR_FIGURAS <- TRUE                # FALSE: las figuras solo se imprimen en pantalla
set.seed(2026)

# Dibuja una figura en pantalla (si la sesión es interactiva) y la guarda en salidas_fase2/figuras
# Si la sesión no tiene pantalla gráfica, no se detiene:
figura <- function(nombre, codigo, ancho = 10, alto = 4.5) {
  codigo <- substitute(codigo); entorno <- parent.frame()
  if (interactive()) tryCatch({
    if (dev.cur() == 1) dev.new()
    eval(codigo, entorno); par(mfrow = c(1, 1))
  }, error = function(e) message("Sin pantalla gráfica (", conditionMessage(e), "): ", nombre,
                                 " se guarda solo en ", FIG))
  if (GUARDAR_FIGURAS) {
    png(file.path(FIG, paste0(nombre, ".png")), width = ancho, height = alto, units = "in", res = 150)
    on.exit(dev.off())
    eval(codigo, entorno)
  }
}
guardar <- function(df, nombre) write.csv(df, file.path(OUT, paste0(nombre, ".csv")), row.names = FALSE)

# Día del año (t en [0, 365), día k centrado en k - 0.5) 
meses <- c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")
fecha <- function(t) {
  d <- as.Date(pmin(floor(t), 364), origin = "2019-01-01")
  ifelse(is.na(t), NA, paste(as.integer(format(d, "%d")), meses[as.integer(format(d, "%m"))]))
}


#Cargar Datos 

ce   <- read.csv(file.path(ENT1, "coeficientes_estaciones.csv"))
cs   <- read.csv(file.path(ENT1, "coeficientes_satelite.csv"))
info <- read.csv(file.path(ENT1, "base_info.csv"))
K_est <- info$K[info$fuente == "estaciones"]
K_sat <- info$K[info$fuente != "estaciones"]
stopifnot(identical(ce$STATION_ID, cs$STATION_ID), identical(ce$YEAR, cs$YEAR))

C_est <- t(as.matrix(ce[, paste0("c", 1:K_est)]))     # K_est x N
C_sat <- t(as.matrix(cs[, paste0("c", 1:K_sat)]))     # K_sat x N
N <- ncol(C_est)

# Información de cada estacióny fase ENSO de cada año
ranking <- as.data.frame(read_excel(file.path(ENT0, "ranking_estaciones.xlsx")))   # altitud y coordenadas
pixeles <- as.data.frame(read_excel(file.path(ENT0, "estacion_pixel.xlsx")))       # fila y columna del píxel

# Fase ENSO de cada año con la regla de episodios de la NOAA CPC (la misma de la fase 4):
# un episodio son >= 5 temporadas trimestrales seguidas con ONI >= 0,5 (Niño) o <= -0,5 (Niña);
# el año toma la fase si el episodio cubre >= 6 de sus 12 temporadas; si no, es Neutro.
oni <- read.csv(file.path(DATOS, "oni_cpc.csv"))
en_episodio <- function(cumple, minimo = 5) {
  r <- rle(cumple); r$values <- r$values & r$lengths >= minimo; inverse.rle(r)
}
oni$nino <- en_episodio(oni$ANOM >= 0.5); oni$nina <- en_episodio(oni$ANOM <= -0.5)
enso <- data.frame(YR = sort(unique(oni$YR)))
enso$ONI_medio       <- as.vector(tapply(oni$ANOM, oni$YR, mean))
enso$temporadas_nino <- as.vector(tapply(oni$nino, oni$YR, sum))
enso$temporadas_nina <- as.vector(tapply(oni$nina, oni$YR, sum))
enso$fase_enso <- ifelse(enso$temporadas_nino >= 6, "Niño", ifelse(enso$temporadas_nina >= 6, "Niña", "Neutro"))
enso <- enso[enso$YR >= 1981 & enso$YR <= 2003, ]
print(table(enso$fase_enso))
guardar(enso, "tabla00_fase_enso_por_anio")
if (file.exists(file.path(ENT4, "fase_enso_por_anio.csv"))) {
  enso4 <- read.csv(file.path(ENT4, "fase_enso_por_anio.csv"))
  stopifnot("La fase ENSO no coincide con la de la fase 4" =
              identical(enso$fase_enso, enso4$fase_enso[match(enso$YR, enso4$YR)]))
  cat("Fase ENSO: coincide con la de la fase 4 en los 23 años\n")
}

meta <- data.frame(STATION_ID = ce$STATION_ID, YEAR = ce$YEAR, particion = ce$particion)
im <- match(meta$STATION_ID, ranking$STATION_ID)
stopifnot("Hay estaciones sin ubicación en ranking_estaciones.xlsx" = !anyNA(im))
meta$nombre     <- ifelse(is.na(ranking$STATION_NAME[im]), "(sin nombre)", ranking$STATION_NAME[im])
meta$altitud    <- ranking$altitud[im]
meta$fuente_alt <- ranking$fuente_alt[im]
meta$Lon        <- ranking$Lon[im]
meta$lat        <- ranking$lat[im]

# Nombres de las zonas acordados por el grupo (zona hidrográfica IDEAM por los 2 primeros dígitos del código)
REG_CAUCA  <- "Cuenca del río Cauca"                              # 26
REG_CALIMA <- "Vertiente del Pacífico del Valle: Calima-Garrapatas" # 54
REG_DAGUA  <- "Vertiente del Pacífico del Valle: Dagua-Anchicayá"   # 53
ETQ_REGION <- c("Cuenca del río Cauca", "Pacífico: Calima-Garrapatas", "Pacífico: Dagua-Anchicayá")  # etiquetas cortas para leyendas
zona_a_region <- function(id) {
  z <- substr(as.character(id), 1, 2)
  factor(ifelse(z == "26", REG_CAUCA, ifelse(z == "53", REG_DAGUA, REG_CALIMA)),
         levels = c(REG_CAUCA, REG_CALIMA, REG_DAGUA))
}
meta$region     <- zona_a_region(meta$STATION_ID)
cortes_alt <- c(-Inf, 800, 1100, 1500, 2000, Inf)
nombres_franja <- c("< 800 m", "800-1100 m", "1100-1500 m", "1500-2000 m", "> 2000 m")
cols_region <- c("grey45", "darkorange", "steelblue")      # Cauca, Calima-Garrapatas, Dagua-Anchicayá
meta$franja    <- cut(meta$altitud, cortes_alt, labels = nombres_franja, right = FALSE)
meta$fase_enso <- factor(enso$fase_enso[match(meta$YEAR, enso$YR)], levels = c("Niño", "Neutro", "Niña"))
ip <- match(meta$STATION_ID, pixeles$STATION_ID); stopifnot(!anyNA(ip))
meta$pixel <- paste(pixeles$fila[ip], pixeles$col[ip], sep = "-")
stopifnot(!anyNA(meta$altitud), !anyNA(meta$fase_enso))

# Bases (las mismas de la fase 1) y objetos fd
rangval  <- c(0, 365)
base_est <- create.bspline.basis(rangval, nbasis = K_est, norder = 4)
base_sat <- create.bspline.basis(rangval, nbasis = K_sat, norder = 4)
Y_fd <- fd(C_est, base_est, list("Día del año", "Estación-año", "sqrt(mm/día)"))

# Satélite: una curva por píxel-año (las estaciones que comparten píxel traen la misma curva)
clave_px <- paste(meta$pixel, meta$YEAR)
unico    <- !duplicated(clave_px)
dup_dif  <- tapply(seq_len(N), clave_px, function(j) max(abs(C_sat[, j] - C_sat[, j[1]])))
stopifnot(max(dup_dif) < 1e-10)                       # las curvas repetidas son idénticas
X_fd <- fd(C_sat[, unico], base_sat, list("Día del año", "Píxel-año", "sqrt(mm/día)"))
meta_sat <- meta[unico, c("pixel", "YEAR", "particion", "fase_enso")]
alt_px   <- tapply(meta$altitud, meta$pixel, mean)    # altitud media de las estaciones del píxel
nest_px  <- tapply(meta$STATION_ID, meta$pixel, function(s) length(unique(s)))
meta_sat$altitud <- as.numeric(alt_px[meta_sat$pixel])
meta_sat$n_estaciones_pixel <- as.integer(nest_px[meta_sat$pixel])
reg_px <- tapply(as.character(meta$region), meta$pixel, function(r) names(which.max(table(r))))   # región mayoritaria del píxel
meta_sat$region <- factor(reg_px[meta_sat$pixel], levels = levels(meta$region))
meta_sat$franja <- cut(meta_sat$altitud, cortes_alt, labels = nombres_franja, right = FALSE)
N_sat <- ncol(X_fd$coefs)
cat("Curvas: estaciones =", N, "| satélite (píxel-año únicos) =", N_sat, "\n")

# Verificación de coherencia con la fase 1: smooth.basis() 
d <- read.csv(file.path(ENT0, "estaciones_diario_limpio.csv"))
d <- d[order(d$STATION_ID, d$YEAR, d$MONTH, d$DAY), ]
es29  <- d$MONTH == 2 & d$DAY == 29                   # misma convención de la fase 1
clave <- paste(d$STATION_ID, d$YEAR, d$MONTH, d$DAY)
pos28 <- match(paste(d$STATION_ID[es29], d$YEAR[es29], 2, 28), clave)
d$VALUE[pos28] <- d$VALUE[pos28] + d$VALUE[es29]
d <- d[!es29, ]


# Los coeficientes de la fase 1 deben venir de los datos del preprocesamiento: mismas estaciones-año y en el mismo orden.

pares_d <- unique(d[, c("STATION_ID", "YEAR")])
stopifnot("Los coeficientes de la fase 1 no corresponden a los datos del preprocesamiento: vuelva a correr la fase 1" =
            nrow(pares_d) == N && all(pares_d$STATION_ID == ce$STATION_ID) && all(pares_d$YEAR == ce$YEAR))
stopifnot(nrow(d) == 365 * N)
t_dia <- (1:365) - 0.5
Y_obs <- matrix(sqrt(d$VALUE), nrow = 365)           # 365 x N, mismo orden que ce
muestra <- sample(N, 100)
coef_fda <- smooth.basis(t_dia, Y_obs[, muestra], base_est)$fd$coefs
verif_fase1 <- max(abs(coef_fda - C_est[, muestra]))
cat("Coherencia con la fase 1 (100 curvas, max |dif| de coeficientes):", signif(verif_fase1, 3), "\n")

# Mallas de evaluación
t_fino <- seq(rangval[1], rangval[2], length.out = 731)   # para medias, varianzas y figuras
t_prof <- seq(2.5, 362.5, by = 5)                         # 73 puntos (cada 5 días) para profundidades
Y_eval <- eval.fd(t_fino, Y_fd)                            # 731 x N
X_eval <- eval.fd(t_fino, X_fd)                            # 731 x N_sat


# Representatividad: se analizan todas las curvas de 1981-2003
# Aquí entran todas las curvas que dejó el preprocesamiento, de los 23 años y de todas las altitudes.

repres <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  curvas = c(N, N_sat),
  unidades = c(length(unique(meta$STATION_ID)), length(unique(meta_sat$pixel))),
  anios = c(length(unique(meta$YEAR)), length(unique(meta_sat$YEAR))),
  min_curvas_por_anio = c(min(table(meta$YEAR)), min(table(meta_sat$YEAR))),
  max_curvas_por_anio = c(max(table(meta$YEAR)), max(table(meta_sat$YEAR))))
print(repres)
guardar(repres, "tabla01_representatividad")
tabla_franja <- as.data.frame.matrix(table(meta$franja, meta$particion))
tabla_franja <- cbind(franja = rownames(tabla_franja), tabla_franja,
                      estaciones = as.integer(tapply(meta$STATION_ID, meta$franja, function(s) length(unique(s)))))
print(tabla_franja)
guardar(tabla_franja, "tabla02_curvas_por_franja")
tabla_region <- data.frame(region = levels(meta$region),
  estaciones = as.integer(tapply(meta$STATION_ID, meta$region, function(s) length(unique(s)))),
  curvas = as.vector(table(meta$region)),
  altitud_mediana = round(as.vector(tapply(meta$altitud[!duplicated(meta$STATION_ID)], meta$region[!duplicated(meta$STATION_ID)], median))),
  altitud_min = round(as.vector(tapply(meta$altitud, meta$region, min))),
  altitud_max = round(as.vector(tapply(meta$altitud, meta$region, max))))
print(tabla_region)
guardar(tabla_region, "tabla02b_curvas_por_region")

figura("fase2_fig01_representatividad", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  tab_anio <- t(table(meta$YEAR, meta$region))
  barplot(tab_anio, col = cols_region, border = NA, las = 2, cex.names = 0.7,
          ylim = c(0, max(colSums(tab_anio)) * 1.45),            
          ylab = "Curvas estación-año", main = "Curvas por año (estaciones)")
  legend("topleft", ETQ_REGION, fill = cols_region, border = NA, bty = "n", cex = 0.7)
  box()
  est1 <- !duplicated(meta$STATION_ID)
  plot(meta$Lon[est1], meta$altitud[est1], pch = 19, col = cols_region[as.integer(meta$region[est1])],
       xlab = "Longitud", ylab = "Altitud (m)", main = "Estaciones: longitud y altitud")
  legend("topleft", ETQ_REGION, col = cols_region, pch = 19, bty = "n", cex = 0.7)
})


#Media funcional: x_barra(t) = c_barra' phi(t)

# Paso a paso: media de los coeficientes por la matriz de bases
bvals_est <- eval.basis(t_fino, base_est)                
bvals_sat <- eval.basis(t_fino, base_sat)
media_est_paso <- bvals_est %*% rowMeans(C_est)
media_sat_paso <- bvals_sat %*% rowMeans(X_fd$coefs)

# Ahora Media funcional con la función mean.fd()
media_est_fd <- mean.fd(Y_fd); media_sat_fd <- mean.fd(X_fd)
media_est <- as.vector(eval.fd(t_fino, media_est_fd))
media_sat <- as.vector(eval.fd(t_fino, media_sat_fd))

#Comprobación: la diferencia debe ser cero
verif_media <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  dif_paso_vs_mean.fd = c(max(abs(media_est_paso - media_est)), max(abs(media_sat_paso - media_sat))),
  dif_paso_vs_promedio_puntual = c(max(abs(media_est_paso - rowMeans(Y_eval))),
                                   max(abs(media_sat_paso - rowMeans(X_eval)))))
print(verif_media)
guardar(verif_media, "tabla03_verificacion_media")


# El suavizado por mínimos cuadrados no impone positividad: en meses secos algunas curvas bajan de 0
cat("Valores negativos de las curvas suavizadas (% de la malla) -> estaciones:", round(100 * mean(Y_eval < 0), 2),
    "| satélite:", round(100 * mean(X_eval < 0), 2), "\n")

figura("fase2_fig02_media_funcional", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  matplot(t_fino, Y_eval, type = "l", lty = 1, col = adjustcolor("grey50", 0.04), ylim = c(-0.7, 6.5),
          xlab = "Día del año", ylab = "sqrt(mm/día)", main = paste0("Estaciones (", N, " curvas)"))
  lines(t_fino, media_est, lwd = 3, col = "navy"); abline(h = 0, lty = 3)
  matplot(t_fino, X_eval, type = "l", lty = 1, col = adjustcolor("grey50", 0.04), ylim = c(-0.7, 6.5),
          xlab = "Día del año", ylab = "sqrt(mm/día)", main = paste0("Satélite (", N_sat, " curvas píxel-año)"))
  lines(t_fino, media_sat, lwd = 3, col = "darkred"); abline(h = 0, lty = 3)
})

# Climatología en mm/día con los DATOS OBSERVADOS (no con las curvas al cuadrado)
mes_dia <- as.integer(format(as.Date(0:364, origin = "2019-01-01"), "%m"))     # mes de cada día
p <- read.csv(file.path(ENT0, "estaciones_pentadas_limpio.csv"))
p <- p[order(p$STATION_ID, p$YEAR, p$PENTADA), ]
p$pixel <- meta$pixel[match(p$STATION_ID, meta$STATION_ID)]
p <- p[!duplicated(paste(p$pixel, p$YEAR, p$PENTADA)), ]                       # píxel-año únicos
p$mes <- ceiling(p$PENTADA / 6)
clim_mm <- data.frame(
  mes = meses,
  estaciones_mm_dia = as.vector(tapply(as.vector(Y_obs^2), rep(mes_dia, N), mean)),
  estaciones_cauca = as.vector(tapply(as.vector(Y_obs[, meta$region == REG_CAUCA]^2),
                                      rep(mes_dia, sum(meta$region == REG_CAUCA)), mean)),
  estaciones_calima_garrapatas = as.vector(tapply(as.vector(Y_obs[, meta$region == REG_CALIMA]^2),
                                                  rep(mes_dia, sum(meta$region == REG_CALIMA)), mean)),
  estaciones_dagua_anchicaya = as.vector(tapply(as.vector(Y_obs[, meta$region == REG_DAGUA]^2),
                                                rep(mes_dia, sum(meta$region == REG_DAGUA)), mean)),
  satelite_mm_dia = as.vector(tapply(p$satelite_mm, p$mes, sum) / tapply(p$dias, p$mes, sum)))
clim_mm[, -1] <- round(clim_mm[, -1], 2)
print(clim_mm)
guardar(clim_mm, "tabla04_climatologia_mensual_mm")
cat("Media anual en mm/día -> estaciones:", round(mean(Y_obs^2), 2),
    "| satélite:", round(sum(p$satelite_mm) / sum(p$dias), 2),
    "| (media de sqrt)^2 estaciones:", round(mean(Y_obs)^2, 2), "(sesgo por elevar al cuadrado)\n")

figura("fase2_fig03_climatologia_mm", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  barplot(t(as.matrix(clim_mm[, c("estaciones_mm_dia", "satelite_mm_dia")])), beside = TRUE,
          names.arg = meses, col = c("navy", "darkred"), ylab = "mm/día", cex.names = 0.8, las = 2,
          main = "Climatología mensual (datos observados)",
          legend.text = c("Estaciones", "Satélite"), args.legend = list(x = "topleft", bty = "n", cex = 0.8))
  barplot(t(as.matrix(clim_mm[, c("estaciones_cauca", "estaciones_calima_garrapatas", "estaciones_dagua_anchicaya")])),
          beside = TRUE, names.arg = meses, col = cols_region, ylab = "mm/día", cex.names = 0.8, las = 2,
          ylim = c(0, max(clim_mm[, 3:5]) * 1.35),
          main = "Estaciones por región (zona hidrográfica)", legend.text = ETQ_REGION,
          args.legend = list(x = "topleft", bty = "n", cex = 0.8))
})





#Varianza y desviación estándar funcional
#Var(t) = (N-1)^{-1} sum_i [x_i(t) - x_barra(t)]^2

#coeficiente de variación funcional
#Var(t) = (N-1)^{-1} sum_i [x_i(t) - x_barra(t)]^2 = phi(t)' Sigma_c phi(t)


Sigma_c_est <- cov(t(C_est))                               # K x K, denominador N - 1
Sigma_c_sat <- cov(t(X_fd$coefs))

var_est_paso <- rowSums((Y_eval - media_est)^2) / (N - 1)                   # punto a punto
var_est_coef <- rowSums((bvals_est %*% Sigma_c_est) * bvals_est)            # phi' Sigma_c phi
var_sat_paso <- rowSums((X_eval - media_sat)^2) / (N_sat - 1)
var_sat_coef <- rowSums((bvals_sat %*% Sigma_c_sat) * bvals_sat)
de_est <- sqrt(var_est_paso); de_sat <- sqrt(var_sat_paso)
cv_est <- de_est / media_est;  cv_sat <- de_sat / media_sat


interior <- t_fino >= 15 & t_fino <= 350
tabla_de <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  DE_max_interior = round(c(max(de_est[interior]), max(de_sat[interior])), 3),
  fecha_DE_max = c(fecha(t_fino[interior][which.max(de_est[interior])]), fecha(t_fino[interior][which.max(de_sat[interior])])),
  DE_min_interior = round(c(min(de_est[interior]), min(de_sat[interior])), 3),
  fecha_DE_min = c(fecha(t_fino[interior][which.min(de_est[interior])]), fecha(t_fino[interior][which.min(de_sat[interior])])),
  DE_en_bordes = round(c(max(de_est[!interior]), max(de_sat[!interior])), 3),
  CV_max_interior = round(c(max(cv_est[interior]), max(cv_sat[interior])), 3),
  fecha_CV_max = c(fecha(t_fino[interior][which.max(cv_est[interior])]), fecha(t_fino[interior][which.max(cv_sat[interior])])),
  CV_min_interior = round(c(min(cv_est[interior]), min(cv_sat[interior])), 3),
  fecha_CV_min = c(fecha(t_fino[interior][which.min(cv_est[interior])]), fecha(t_fino[interior][which.min(cv_sat[interior])])))
print(tabla_de)
guardar(tabla_de, "tabla05b_DE_CV_extremos")

#Con la función std.fd() se obtiene una versión suavizada de la desviación estándar
de_est_stdfd <- as.vector(eval.fd(t_fino, std.fd(Y_fd)))
de_sat_stdfd <- as.vector(eval.fd(t_fino, std.fd(X_fd)))
verif_var <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  dif_var_puntual_vs_coeficientes = c(max(abs(var_est_paso - var_est_coef)), max(abs(var_sat_paso - var_sat_coef))),
  dif_media_DE_vs_std.fd = c(mean(abs(de_est - de_est_stdfd)), mean(abs(de_sat - de_sat_stdfd))))
print(verif_var)
guardar(verif_var, "tabla05_verificacion_varianza")

figura("fase2_fig04_media_DE", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  plot(t_fino, media_est, type = "l", lwd = 3, col = "navy", ylim = c(0, 5), xlab = "Día del año",
       ylab = "sqrt(mm/día)", main = "Estaciones: media y media ± 2 DE")
  lines(t_fino, media_est + 2 * de_est, lty = 2, lwd = 2, col = "red")
  lines(t_fino, pmax(media_est - 2 * de_est, 0), lty = 2, lwd = 2, col = "red")
  lines(t_fino, de_est, lwd = 2, col = "darkgreen")
  legend("topleft", c("Media", "Media ± 2 DE (truncada en 0)", "DE"), col = c("navy", "red", "darkgreen"),
         lty = c(1, 2, 1), lwd = c(3, 2, 2), bty = "n", cex = 0.8)
  plot(t_fino, media_sat, type = "l", lwd = 3, col = "darkred", ylim = c(0, 5), xlab = "Día del año",
       ylab = "sqrt(mm/día)", main = "Satélite: media y media ± 2 DE")
  lines(t_fino, media_sat + 2 * de_sat, lty = 2, lwd = 2, col = "red")
  lines(t_fino, pmax(media_sat - 2 * de_sat, 0), lty = 2, lwd = 2, col = "red")
  lines(t_fino, de_sat, lwd = 2, col = "darkgreen")
})

figura("fase2_fig05_DE_CV", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  matplot(t_fino, cbind(de_est, de_sat), type = "l", lty = c(1, 2), lwd = 3, col = c("navy", "darkred"),
          ylim = c(0, max(de_est, de_sat) * 1.15), xlab = "Día del año", ylab = "sqrt(mm/día)",
          main = "Desviación estándar funcional")
  legend("topleft", c("Estaciones", "Satélite"), col = c("navy", "darkred"), lty = c(1, 2), lwd = 3, bty = "n")
  abline(v = c(15, 350), lty = 3, col = "grey50")          # fuera de estas líneas: efecto de borde
  matplot(t_fino, cbind(cv_est, cv_sat), type = "l", lty = c(1, 2), lwd = 3, col = c("navy", "darkred"),
          ylim = c(0, max(cv_est, cv_sat) * 1.15), xlab = "Día del año", ylab = "CV(t)",
          main = "Coeficiente de variación funcional")
  abline(v = c(15, 350), lty = 3, col = "grey50")
})


# Covarianza funcional paso a paso: Cov(s,t) = phi(s)' Sigma_c phi(t)
# correlación funcional de cada fuente: r(s,t) = v(s,t) / sqrt(v(s,s) v(t,t))

t_cov <- seq(0, 365, length.out = 147)                     # malla para las superficies
b_cov_est <- eval.basis(t_cov, base_est); b_cov_sat <- eval.basis(t_cov, base_sat)
cov_est_paso <- b_cov_est %*% Sigma_c_est %*% t(b_cov_est)
cov_sat_paso <- b_cov_sat %*% Sigma_c_sat %*% t(b_cov_sat)
cov_est <- eval.bifd(t_cov, t_cov, var.fd(Y_fd))
cov_sat <- eval.bifd(t_cov, t_cov, var.fd(X_fd))
cor_est <- cor.fd(t_cov, Y_fd)
cor_sat <- cor.fd(t_cov, X_fd)
cor_est_paso <- cov_est_paso / sqrt(diag(cov_est_paso) %o% diag(cov_est_paso))

#Covarianza funcional con la función var.fd()
verif_cov <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  dif_cov_paso_vs_var.fd = c(max(abs(cov_est_paso - cov_est)), max(abs(cov_sat_paso - cov_sat))),
  dif_diagonal_vs_varianza = c(max(abs(diag(cov_est_paso) - rowSums((b_cov_est %*% Sigma_c_est) * b_cov_est))),
                               max(abs(diag(cov_sat_paso) - rowSums((b_cov_sat %*% Sigma_c_sat) * b_cov_sat)))),
  dif_cor_paso_vs_cor.fd = c(max(abs(cor_est_paso - cor_est)), NA),
  simetria_cov = c(max(abs(cov_est - t(cov_est))), max(abs(cov_sat - t(cov_sat)))))
print(verif_cov)
guardar(verif_cov, "tabla06_verificacion_covarianza")

# Correlaciones entre fechas clave (calculadas exactamente con las bases)
cor_entre <- function(fdobj, s, t) {
  b <- eval.basis(c(s, t), fdobj$basis); S <- cov(t(fdobj$coefs)); v <- b %*% S %*% t(b)
  v[1, 2] / sqrt(v[1, 1] * v[2, 2])
}
fechas_clave <- c(ene = 15.5, abr = 115.5, jul = 195.5, oct = 297.5)
pares <- combn(names(fechas_clave), 2)
cor_clave <- data.frame(
  par = apply(pares, 2, paste, collapse = "-"),
  estaciones = apply(pares, 2, function(q) cor_entre(Y_fd, fechas_clave[q[1]], fechas_clave[q[2]])),
  satelite   = apply(pares, 2, function(q) cor_entre(X_fd, fechas_clave[q[1]], fechas_clave[q[2]])))
cor_clave[, 2:3] <- round(cor_clave[, 2:3], 3)
print(cor_clave)
guardar(cor_clave, "tabla07_correlacion_fechas_clave")

# "Memoria" temporal: distancia media (días) a la que la correlación cae por debajo de 0.5
memoria <- function(R, tt) mean(sapply(seq_along(tt), function(i) {
  dist <- abs(tt - tt[i]); max(dist[R[i, ] >= 0.5]) }))
cat("Memoria (r >= 0.5) -> estaciones:", round(memoria(cor_est, t_cov)), "días | satélite:",
    round(memoria(cor_sat, t_cov)), "días\n")

# Descomposición espectral de la covarianza 
pca_est <- pca.fd(Y_fd, nharm = 4); pca_sat <- pca.fd(X_fd, nharm = 4)
var_modos <- data.frame(modo = 1:4, estaciones_pct = round(100 * pca_est$varprop, 1),
                        satelite_pct = round(100 * pca_sat$varprop, 1))
print(var_modos)
guardar(var_modos, "tabla08_modos_covarianza")

# Correlación del puntaje del modo 1 con el nivel medio de cada curva y con la región
nivel_est <- colMeans(Y_eval)
cat("Estaciones: cor(puntaje modo 1, nivel medio de la curva) =",
    round(cor(pca_est$scores[, 1], nivel_est), 3), "\n")

figura("fase2_fig06_cov_cor_estaciones", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  persp(t_cov, t_cov, cov_est, theta = -45, phi = 25, r = 3, expand = 0.6, col = "grey85", border = NA,
        shade = 0.5, xlab = "s (día)", ylab = "t (día)", zlab = "Covarianza", main = "Estaciones: covarianza")
  image(t_cov, t_cov, cor_est, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
        xlab = "s (día)", ylab = "t (día)", main = "Estaciones: correlación")
  contour(t_cov, t_cov, cor_est, levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
}, ancho = 10, alto = 4.8)

figura("fase2_fig07_cov_cor_satelite", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  persp(t_cov, t_cov, cov_sat, theta = -45, phi = 25, r = 3, expand = 0.6, col = "grey85", border = NA,
        shade = 0.5, xlab = "s (día)", ylab = "t (día)", zlab = "Covarianza", main = "Satélite: covarianza")
  image(t_cov, t_cov, cor_sat, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
        xlab = "s (día)", ylab = "t (día)", main = "Satélite: correlación")
  contour(t_cov, t_cov, cor_sat, levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
}, ancho = 10, alto = 4.8)

# Derivadas de los datos funcionales: D x(t) = c' D phi(t)
# inicio, pico y fin de las temporadas de lluvia


D_media_est <- as.vector(eval.fd(t_fino, media_est_fd, 1))
D_media_sat <- as.vector(eval.fd(t_fino, media_sat_fd, 1))
DY_fino <- eval.fd(t_fino, Y_fd, 1); DX_fino <- eval.fd(t_fino, X_fd, 1)
verif_deriv <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  media_derivadas_vs_derivada_media = c(max(abs(rowMeans(DY_fino) - D_media_est)),
                                        max(abs(rowMeans(DX_fino) - D_media_sat))),
  
  # Lo mismo con deriv.fd()
  igualdad_con_deriv.fd = c(max(abs(eval.fd(t_fino, mean.fd(deriv.fd(Y_fd, 1))) - eval.fd(t_fino, deriv.fd(media_est_fd, 1)))),
                            max(abs(eval.fd(t_fino, mean.fd(deriv.fd(X_fd, 1))) - eval.fd(t_fino, deriv.fd(media_sat_fd, 1))))),
  # pero se aleja de la derivada exacta
  deriv.fd_vs_derivada_exacta = c(max(abs(eval.fd(t_fino, deriv.fd(media_est_fd, 1)) - D_media_est)),
                                  max(abs(eval.fd(t_fino, deriv.fd(media_sat_fd, 1)) - D_media_sat))),
  derivada_media_max_abs = c(max(abs(D_media_est)), max(abs(D_media_sat))))
print(verif_deriv)
guardar(verif_deriv, "tabla09a_verificacion_derivadas")

# Fechas de una curva x (con su derivada D) evaluadas en la malla tt.

fechas_temporadas <- function(x, D, tt) {
  en <- function(a, b) which(tt >= a & tt <= b)
  w1 <- en(45, 200);  i_p1 <- w1[which.max(x[w1])]
  w2 <- en(230, 350); i_p2 <- w2[which.max(x[w2])]
  wv <- en(tt[i_p1], tt[i_p2]); i_v <- wv[which.min(x[wv])]
  wa <- en(10, tt[i_p1]);       i_i1 <- wa[which.max(D[wa])]
  wb <- en(tt[i_p1], tt[i_v]);  i_f1 <- wb[which.min(D[wb])]
  wc <- en(tt[i_v], tt[i_p2]);  i_i2 <- wc[which.max(D[wc])]
  wd <- en(tt[i_p2], 355);      i_f2 <- wd[which.min(D[wd])]
  
  # un pico es "real" si la derivada cambia de signo (no queda pegado al borde de la ventana)
  real <- function(i) i > 1 && i < length(tt) && D[i - 1] > 0 && D[i + 1] < 0
  
  # si el extremo de la derivada cae en el límite de la búsqueda (10 o 355), la fecha no queda
  # determinada dentro del año (la lluvia sigue cayendo hacia el borde): se deja como NA
  paso <- tt[2] - tt[1]
  c(inicio1 = if (tt[i_i1] <= 10 + paso) NA else tt[i_i1], pico1 = tt[i_p1], fin1 = tt[i_f1], valle = tt[i_v],
    inicio2 = tt[i_i2], pico2 = tt[i_p2], fin2 = if (tt[i_f2] >= 355 - paso) NA else tt[i_f2],
    pico1_real = real(i_p1), pico2_real = real(i_p2))
}

# Definición complementaria (robustez): una temporada es el tramo en que la curva está por encima
# de su propio nivel medio anual; inicio y fin son los cruces con ese nivel.
umbral_temporadas <- function(x, tt, pico1, pico2) {
  arriba <- x > mean(x); r <- rle(arriba); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  tramo <- function(p) { k <- which(r$values & tt[ini] <= p & tt[fin] >= p)
    if (length(k) == 0) return(c(NA, NA)); c(if (ini[k] == 1) NA else tt[ini[k]], if (fin[k] == length(tt)) NA else tt[fin[k]]) }
  s1 <- tramo(pico1); s2 <- tramo(pico2)
  c(inicio1_umbral = s1[1], fin1_umbral = s1[2], inicio2_umbral = s2[1], fin2_umbral = s2[2])
}

f_media_est <- fechas_temporadas(media_est, D_media_est, t_fino)
f_media_sat <- fechas_temporadas(media_sat, D_media_sat, t_fino)
nombres_f <- c("inicio1", "pico1", "fin1", "valle", "inicio2", "pico2", "fin2")
tabla_fechas <- data.frame(evento = nombres_f,
                           estaciones_dia = round(f_media_est[nombres_f]), estaciones_fecha = fecha(f_media_est[nombres_f]),
                           satelite_dia = round(f_media_sat[nombres_f]), satelite_fecha = fecha(f_media_sat[nombres_f]))
u_est <- umbral_temporadas(media_est, t_fino, f_media_est["pico1"], f_media_est["pico2"])
u_sat <- umbral_temporadas(media_sat, t_fino, f_media_sat["pico1"], f_media_sat["pico2"])
tabla_umbral <- data.frame(evento = c("inicio1", "fin1", "inicio2", "fin2"),
                           estaciones_umbral = fecha(u_est), satelite_umbral = fecha(u_sat))
tabla_fechas <- merge(tabla_fechas, tabla_umbral, by = "evento", all.x = TRUE, sort = FALSE)
tabla_fechas <- tabla_fechas[match(nombres_f, tabla_fechas$evento), ]
print(tabla_fechas)
guardar(tabla_fechas, "tabla09_fechas_temporadas_media")

# Fechas de cada curva (malla diaria)
DY_dia <- eval.fd(t_dia, Y_fd, 1); Y_dia <- eval.fd(t_dia, Y_fd)
F_est <- t(sapply(1:N, function(j) fechas_temporadas(Y_dia[, j], DY_dia[, j], t_dia)))
F_est <- data.frame(meta[, c("STATION_ID", "YEAR", "region", "franja", "fase_enso")], F_est)
cat("Curvas con los dos picos bien definidos:", round(100 * mean(F_est$pico1_real & F_est$pico2_real), 1), "%\n")
cat("Fechas no determinadas dentro del año (borde) -> inicio1:", round(100 * mean(is.na(F_est$inicio1)), 1),
    "% | fin2:", round(100 * mean(is.na(F_est$fin2)), 1), "%\n")
guardar(F_est, "fechas_temporadas_por_curva")

# Resumen: mediana de cada fecha por región y por fase ENSO (curvas con los dos picos definidos)
ok <- F_est$pico1_real == 1 & F_est$pico2_real == 1
resumen_fechas <- function(grupo) {
  r <- aggregate(F_est[ok, nombres_f], list(grupo = grupo[ok]), function(v) median(v, na.rm = TRUE))
  r$n_curvas <- as.vector(table(grupo[ok])[as.character(r$grupo)]); r
}
fechas_grupo <- rbind(cbind(tipo = "region", resumen_fechas(F_est$region)),
                      cbind(tipo = "enso",   resumen_fechas(F_est$fase_enso)))
fechas_grupo_txt <- fechas_grupo
fechas_grupo_txt[, nombres_f] <- lapply(fechas_grupo[, nombres_f], fecha)
print(fechas_grupo_txt)
guardar(fechas_grupo_txt, "tabla10_fechas_por_region_y_enso")

figura("fase2_fig08_derivadas_temporadas", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 4.2), oma = c(2, 0, 0, 0))
  for (fuente in c("est", "sat")) {
    m <- if (fuente == "est") media_est else media_sat
    Dm <- if (fuente == "est") D_media_est else D_media_sat
    ff <- if (fuente == "est") f_media_est else f_media_sat
    plot(t_fino, m, type = "l", lwd = 3, col = if (fuente == "est") "navy" else "darkred",
         xlab = "Día del año", ylab = "Media sqrt(mm/día)", ylim = range(m) + c(-0.1, 0.1),
         main = if (fuente == "est") "Estaciones: media y su derivada" else "Satélite: media y su derivada")
    abline(v = ff[c("inicio1", "inicio2")], col = "darkgreen", lty = 2, lwd = 2)
    abline(v = ff[c("fin1", "fin2")], col = "orange", lty = 2, lwd = 2)
    abline(v = ff[c("pico1", "pico2")], col = "grey40", lty = 3)
    par(new = TRUE)
    plot(t_fino, Dm, type = "l", lwd = 2, col = "purple", axes = FALSE, xlab = "", ylab = "")
    abline(h = 0, lty = 3, col = "purple"); axis(4, col.axis = "purple"); mtext("D x(t)", side = 4, line = 2.5, col = "purple")
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("bottom", c("Inicio (máx. de D)", "Fin (mín. de D)", "Pico (D = 0)", "Derivada D x(t)"),
         col = c("darkgreen", "orange", "grey40", "purple"), lty = c(2, 2, 3, 1), lwd = 2, bty = "n", cex = 0.85, horiz = TRUE)
})

figura("fase2_fig09_fechas_por_enso", {
  par(mfrow = c(1, 4), mar = c(4.2, 4.2, 2.5, 1))
  for (ev in c("inicio1", "fin1", "inicio2", "fin2")) {
    boxplot(F_est[ok, ev] ~ F_est$fase_enso[ok], col = c("tomato", "grey80", "steelblue"),
            xlab = "", ylab = "Día del año", main = c(inicio1 = "Inicio 1.ª temporada", fin1 = "Fin 1.ª temporada",
                                                       inicio2 = "Inicio 2.ª temporada", fin2 = "Fin 2.ª temporada")[ev],
            outline = FALSE)
  }
}, ancho = 11, alto = 4)

# Profundidad funcional, mediana, media recortada y boxplot funcional
# Fraiman-Muniz 
# fda.usc::depth.FM devuelve la versión reescalada a [0, 1]: 2 * D_FM - 1 (se verifica abajo).
# MBD: proporción promedio del dominio en que la curva queda dentro de la banda de cada par.

Y_prof <- eval.fd(t_prof, Y_fd)                            # 73 x N
X_prof <- eval.fd(t_prof, X_fd)                            # 73 x N_sat

# Fraiman-Muniz paso a paso y con fda.usc

FM_paso <- function(M) {                                   # M: puntos x curvas
  Fn <- t(apply(M, 1, function(v) rank(v, ties.method = "max") / length(v)))
  colMeans(1 - abs(0.5 - Fn))                              # regla del rectángulo en malla regular
}
fm_est <- depth.FM(fdata(t(Y_prof), t_prof, rangval), trim = 0.25)
fm_sat <- depth.FM(fdata(t(X_prof), t_prof, rangval), trim = 0.25)
verif_FM <- c(estaciones = max(abs(fm_est$dep - (2 * FM_paso(Y_prof) - 1))),
              satelite   = max(abs(fm_sat$dep - (2 * FM_paso(X_prof) - 1))))
cat("FM: depth.FM vs 2*D_FM(paso a paso) - 1 (max |dif|):\n"); print(verif_FM)

#MBD por definición (pares de curvas) en una submuestra, contra fda:::fMBD
MBD_definicion <- function(M) {
  n <- ncol(M); dep <- numeric(n); np <- 0
  for (a in 1:(n - 1)) for (b in (a + 1):n) {
    lo <- pmin(M[, a], M[, b]); hi <- pmax(M[, a], M[, b])
    dep <- dep + colMeans(M >= lo & M <= hi); np <- np + 1
  }
  dep / np
}
sub <- sample(N, 80)
verif_MBD <- max(abs(MBD_definicion(Y_prof[, sub]) - fbplot(Y_prof[, sub], t_prof, plot = FALSE)$depth))
cat("MBD por definición vs fbplot (80 curvas, max |dif|):", signif(verif_MBD, 3), "\n")

# Profundidad modal
mode_est <- depth.mode(fdata(t(Y_prof), t_prof, rangval), trim = 0.25)
mode_sat <- depth.mode(fdata(t(X_prof), t_prof, rangval), trim = 0.25)

# ¿Ordenan igual las tres profundidades? 
mbd_est <- fbplot(Y_prof, t_prof, method = "MBD", plot = FALSE)$depth
mbd_sat <- fbplot(X_prof, t_prof, method = "MBD", plot = FALSE)$depth
concordancia <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  FM_vs_MBD = c(cor(fm_est$dep, mbd_est, method = "spearman"), cor(fm_sat$dep, mbd_sat, method = "spearman")),
  FM_vs_modal = c(cor(fm_est$dep, mode_est$dep, method = "spearman"), cor(fm_sat$dep, mode_sat$dep, method = "spearman")),
  MBD_vs_modal = c(cor(mbd_est, mode_est$dep, method = "spearman"), cor(mbd_sat, mode_sat$dep, method = "spearman")),
  mediana_FM = c(paste(meta$nombre[fm_est$lmed], meta$YEAR[fm_est$lmed]), paste("píxel", meta_sat$pixel[fm_sat$lmed], meta_sat$YEAR[fm_sat$lmed])),
  mediana_MBD = c(paste(meta$nombre[which.max(mbd_est)], meta$YEAR[which.max(mbd_est)]), paste("píxel", meta_sat$pixel[which.max(mbd_sat)], meta_sat$YEAR[which.max(mbd_sat)])),
  mediana_modal = c(paste(meta$nombre[mode_est$lmed], meta$YEAR[mode_est$lmed]), paste("píxel", meta_sat$pixel[mode_sat$lmed], meta_sat$YEAR[mode_sat$lmed])))
concordancia[, 2:4] <- round(concordancia[, 2:4], 3)
print(concordancia)
guardar(concordancia, "tabla11b_concordancia_profundidades")

#Medidas de centralidad: media, mediana (curva más profunda FM) y media recortada al 25 %
centralidad <- function(fm, M, tt, media_vec, tf) data.frame(
  medida = c("Media", "Mediana FM", "Media recortada 25 %"),
  maximo = c(max(media_vec), max(fm$median$data), max(fm$mtrim$data)),
  dia_maximo = c(round(tf[which.max(media_vec)]), round(tt[which.max(fm$median$data)]), round(tt[which.max(fm$mtrim$data)])))
cent <- rbind(cbind(fuente = "Estaciones", centralidad(fm_est, Y_prof, t_prof, media_est, t_fino)),
              cbind(fuente = "Satélite",   centralidad(fm_sat, X_prof, t_prof, media_sat, t_fino)))
cent$maximo <- round(cent$maximo, 3)
print(cent)
cat("Curva mediana FM de estaciones:", meta$STATION_ID[fm_est$lmed], meta$nombre[fm_est$lmed],
    meta$YEAR[fm_est$lmed], "\n")
guardar(cent, "tabla11_centralidad")

# Boxplot funcional con MBD sobre las curvas crudas
# fbplot recibe una matriz puntos x curvas

figura("fase2_fig10_boxplot_crudo", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  bx_est <- fbplot(Y_prof, t_prof, method = "MBD", xlim = rangval, ylim = c(0, 7),
                    xlab = "Día del año", ylab = "sqrt(mm/día)")
  title("Estaciones (curvas crudas)")
  bx_sat <- fbplot(X_prof, t_prof, method = "MBD", xlim = rangval, ylim = c(0, 7),
                    xlab = "Día del año", ylab = "sqrt(mm/día)")
  title("Satélite (curvas crudas)")
})
if (!exists("bx_est")) {                      
  bx_est <- fbplot(Y_prof, t_prof, method = "MBD", plot = FALSE)
  bx_sat <- fbplot(X_prof, t_prof, method = "MBD", plot = FALSE)
}
cat("Spearman FM vs MBD -> estaciones:", round(cor(fm_est$dep, bx_est$depth, method = "spearman"), 3),
    "| satélite:", round(cor(fm_sat$dep, bx_sat$depth, method = "spearman"), 3), "\n")

# Nivel de cada curva respecto a la media (para saber si un atípico es húmedo o seco)
nivel_rel <- function(M) colMeans(M) - mean(M)

atipicos_crudo <- data.frame(meta[bx_est$outpoint, c("STATION_ID", "nombre", "YEAR", "altitud", "fuente_alt",
                                                     "region", "franja", "fase_enso")],
                             prof_FM = round(fm_est$dep[bx_est$outpoint], 3),
                             prof_MBD = round(bx_est$depth[bx_est$outpoint], 4),
                             tipo = ifelse(nivel_rel(Y_prof)[bx_est$outpoint] > 0, "más lluviosa", "más seca"))
atipicos_crudo <- atipicos_crudo[order(atipicos_crudo$prof_MBD), ]
guardar(atipicos_crudo, "tabla12_atipicos_estaciones_crudo")
cat("\nAtípicos (crudo, estaciones):", nrow(atipicos_crudo), "de", N, "curvas\n")
resumen_crudo <- data.frame(
  franja = levels(meta$franja),
  curvas = as.vector(table(meta$franja)),
  atipicas = as.vector(table(factor(atipicos_crudo$franja, levels = levels(meta$franja)))))
resumen_crudo$pct_atipicas <- round(100 * resumen_crudo$atipicas / resumen_crudo$curvas, 1)
print(resumen_crudo)
guardar(resumen_crudo, "tabla13_atipicos_crudo_por_franja")
resumen_crudo_reg <- data.frame(region = levels(meta$region), curvas = as.vector(table(meta$region)),
  atipicas = as.vector(table(factor(atipicos_crudo$region, levels = levels(meta$region)))))
resumen_crudo_reg$pct_atipicas <- round(100 * resumen_crudo_reg$atipicas / resumen_crudo_reg$curvas, 1)
print(resumen_crudo_reg)
guardar(resumen_crudo_reg, "tabla13b_atipicos_crudo_por_region")
cat("Estaciones distintas entre las curvas atípicas (crudo):", length(unique(atipicos_crudo$STATION_ID)), "\n")

atipicos_sat_crudo <- data.frame(meta_sat[bx_sat$outpoint, ], prof_MBD = round(bx_sat$depth[bx_sat$outpoint], 4),
                                 tipo = ifelse(nivel_rel(X_prof)[bx_sat$outpoint] > 0, "más lluviosa", "más seca"))
guardar(atipicos_sat_crudo, "tabla14_atipicos_satelite_crudo")
cat("Atípicos (crudo, satélite):", nrow(atipicos_sat_crudo), "de", N_sat, "\n")

# Anomalías: curva menos la curva media de su propio sitio
# Quita la diferencia de nivel entre sitios (orografía) y deja la variación de un año a otro.
# Se usan los sitios con al menos 3 años (con 1 o 2 años la anomalía es casi nula por construcción).
anomalias <- function(C, sitio, min_anios = 3) {
  n_sitio <- table(sitio); usar <- sitio %in% names(n_sitio)[n_sitio >= min_anios]
  medias <- sapply(split(seq_along(sitio), sitio), function(j) rowMeans(C[, j, drop = FALSE]))
  list(C = C[, usar] - medias[, as.character(sitio[usar])], usar = usar)
}
an_est <- anomalias(C_est, meta$STATION_ID)
an_sat <- anomalias(X_fd$coefs, meta_sat$pixel)
A_fd  <- fd(an_est$C, base_est); AX_fd <- fd(an_sat$C, base_sat)
meta_an <- meta[an_est$usar, ]; meta_an_sat <- meta_sat[an_sat$usar, ]
cat("Anomalías: estaciones =", ncol(an_est$C), "curvas (", sum(!an_est$usar), "excluidas) | satélite =",
    ncol(an_sat$C), "(", sum(!an_sat$usar), "excluidas)\n")
A_prof <- eval.fd(t_prof, A_fd); AX_prof <- eval.fd(t_prof, AX_fd)
fm_an <- depth.FM(fdata(t(A_prof), t_prof, rangval), trim = 0.25)

# Boxplot funcional (MBD, factor 1.5) de las anomalías de cada estación-año y de cada píxel-año
bxa_est <- fbplot(A_prof, t_prof, method = "MBD", plot = FALSE)
bxa_sat <- fbplot(AX_prof, t_prof, method = "MBD", plot = FALSE)
cat("Atípicos por la regla del boxplot en anomalías -> estaciones:", length(bxa_est$outpoint), "de", ncol(A_prof),
    "| satélite:", length(bxa_sat$outpoint), "de", ncol(AX_prof), "\n")

atip_an_box <- data.frame(meta_an[bxa_est$outpoint, c("STATION_ID", "nombre", "YEAR", "altitud", "region", "fase_enso")],
                          prof_MBD = round(bxa_est$depth[bxa_est$outpoint], 4),
                          anomalia_media = round(colMeans(A_prof)[bxa_est$outpoint], 3))
print(atip_an_box)
guardar(atip_an_box, "tabla15b_atipicos_boxplot_anomalias")


# (a) Una curva por AÑO: anomalía media de las estaciones de la cuenca del Cauca (y de sus píxeles) en ese año.
#     Con 23 curvas, el boxplot funcional señala los años atípicos de toda la región.

anom_anual <- function(A_eval, anios) sapply(split(seq_along(anios), anios), function(j) rowMeans(A_eval[, j, drop = FALSE]))
and_an   <- meta_an$region == REG_CAUCA
anual_and_prof <- anom_anual(A_prof[, and_an], meta_an$YEAR[and_an])        # 73 x 23
and_an_sat <- meta_an_sat$region == REG_CAUCA                                # mismos criterios que las estaciones
anual_sat_prof <- anom_anual(AX_prof[, and_an_sat], meta_an_sat$YEAR[and_an_sat])
bx_anio_est <- fbplot(anual_and_prof, t_prof, method = "MBD", plot = FALSE)
bx_anio_sat <- fbplot(anual_sat_prof, t_prof, method = "MBD", plot = FALSE)
fm_anio_est <- depth.FM(fdata(t(anual_and_prof), t_prof, rangval), trim = 0.25)
anios_tab <- data.frame(YEAR = as.integer(colnames(anual_and_prof)))
anios_tab$fase_enso <- enso$fase_enso[match(anios_tab$YEAR, enso$YR)]
anios_tab$ONI_medio <- round(enso$ONI_medio[match(anios_tab$YEAR, enso$YR)], 2)
anios_tab$anomalia_media_estaciones <- round(colMeans(anual_and_prof), 3)
anios_tab$prof_MBD_estaciones <- round(bx_anio_est$depth, 3)
anios_tab$prof_FM_estaciones  <- round(fm_anio_est$dep, 3)
anios_tab$atipico_estaciones  <- seq_len(ncol(anual_and_prof)) %in% bx_anio_est$outpoint
anios_tab$anomalia_media_satelite <- round(colMeans(anual_sat_prof), 3)
anios_tab$prof_MBD_satelite <- round(bx_anio_sat$depth, 3)
anios_tab$atipico_satelite  <- seq_len(ncol(anual_sat_prof)) %in% bx_anio_sat$outpoint
anios_tab$rango_profundidad <- rank(anios_tab$prof_MBD_estaciones)          # 1 = año menos profundo
print(anios_tab[order(anios_tab$prof_MBD_estaciones), ])
guardar(anios_tab, "tabla15_anios_profundidad_y_atipicos")
cat("Correlación entre la anomalía media anual de estaciones y de satélite (23 años):",
    round(cor(anios_tab$anomalia_media_estaciones, anios_tab$anomalia_media_satelite), 3), "\n")
cat("Correlación de la anomalía media anual de estaciones con el ONI medio del año:",
    round(cor(anios_tab$anomalia_media_estaciones, anios_tab$ONI_medio), 3), "\n")

# Correlación funcional de las ANOMALÍAS: sin la diferencia de nivel entre sitios, ¿cuánta memoria
# temporal queda dentro del año?

cor_an_est <- cor.fd(t_cov, A_fd); cor_an_sat <- cor.fd(t_cov, AX_fd)
cat("Memoria (r >= 0.5) de las anomalías -> estaciones:", round(memoria(cor_an_est, t_cov)), "días | satélite:",
    round(memoria(cor_an_sat, t_cov)), "días\n")
figura("fase2_fig07b_correlacion_anomalias", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  image(t_cov, t_cov, cor_an_est, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
        xlab = "s (día)", ylab = "t (día)", main = "Estaciones: correlación de las anomalías")
  contour(t_cov, t_cov, cor_an_est, levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
  image(t_cov, t_cov, cor_an_sat, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
        xlab = "s (día)", ylab = "t (día)", main = "Satélite: correlación de las anomalías")
  contour(t_cov, t_cov, cor_an_sat, levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
}, ancho = 10, alto = 4.8)

# (b) Estación-año: el 5 % de curvas de anomalía MENOS profundas (MBD). No es la regla del boxplot,
#     sino la cola de la profundidad; sirve para ver en qué años y fases se concentran los casos extremos.

signo_an <- colMeans(A_prof)                               # >0: año más lluvioso que lo normal en esa estación
cola <- bxa_est$depth <= quantile(bxa_est$depth, 0.05)
extremas_an <- data.frame(meta_an[cola, c("STATION_ID", "nombre", "YEAR", "altitud", "region", "fase_enso")],
                          prof_FM = round(fm_an$dep[cola], 3), prof_MBD = round(bxa_est$depth[cola], 4),
                          anomalia_media = round(signo_an[cola], 3),
                          tipo = ifelse(signo_an[cola] > 0, "año húmedo", "año seco"))
extremas_an <- extremas_an[order(extremas_an$prof_MBD), ]
guardar(extremas_an, "tabla16_anomalias_extremas_5pct_estacion_anio")
por_fase <- data.frame(fase_enso = levels(meta_an$fase_enso),
                       curvas = as.vector(table(meta_an$fase_enso)),
                       extremas_humedas = as.vector(table(factor(extremas_an$fase_enso, levels(meta_an$fase_enso))[extremas_an$tipo == "año húmedo"])),
                       extremas_secas = as.vector(table(factor(extremas_an$fase_enso, levels(meta_an$fase_enso))[extremas_an$tipo == "año seco"])))
por_fase$pct_humedas <- round(100 * por_fase$extremas_humedas / por_fase$curvas, 1)
por_fase$pct_secas   <- round(100 * por_fase$extremas_secas / por_fase$curvas, 1)
print(por_fase)                                            # si ENSO no importara, cada fase tendría ~5 % en total
guardar(por_fase, "tabla17_anomalias_extremas_por_enso")
top_anios <- sort(table(extremas_an$YEAR), decreasing = TRUE)
cat("Años con más estaciones en la cola del 5 %:\n"); print(head(top_anios, 8))

figura("fase2_fig11_boxplot_anomalias", {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.6, 1))
  atipicos_txt <- function(at) paste("Años atípicos:", if (length(at)) paste(at, collapse = ", ") else "ninguno")
  fbplot(A_prof, t_prof, method = "MBD", xlim = rangval, ylim = c(-3, 3), xlab = "Día del año",
         ylab = "Anomalía sqrt(mm/día)")
  title(paste0("Estaciones: anomalía de cada estación-año (", ncol(A_prof), ")"), line = 1.8); abline(h = 0, lty = 3)
  fbplot(anual_and_prof, t_prof, method = "MBD", xlim = rangval, ylim = c(-0.8, 0.8), xlab = "Día del año",
         ylab = "Anomalía sqrt(mm/día)")
  title("Estaciones: anomalía media anual, cuenca del río Cauca", line = 1.8); abline(h = 0, lty = 3)
  fbplot(AX_prof, t_prof, method = "MBD", xlim = rangval, ylim = c(-3, 3), xlab = "Día del año",
         ylab = "Anomalía sqrt(mm/día)")
  title(paste0("Satélite: anomalía de cada píxel-año (", ncol(AX_prof), ")"), line = 1.8); abline(h = 0, lty = 3)
  fbplot(anual_sat_prof, t_prof, method = "MBD", xlim = rangval, ylim = c(-0.8, 0.8), xlab = "Día del año",
         ylab = "Anomalía sqrt(mm/día)")
  title("Satélite: anomalía media anual, cuenca del río Cauca", line = 1.8); abline(h = 0, lty = 3)
 
}, ancho = 11, alto = 8)

figura("fase2_fig12_profundidad_centralidad", {
  par(mfrow = c(1, 2), mar = c(8, 4.2, 2.5, 2))             # margen inferior para la leyenda
  orden <- order(fm_est$dep, decreasing = TRUE)
  grises <- gray(seq(0.85, 0.15, length.out = N))         # claras = más profundas (se dibujan encima)
  matplot(t_prof, Y_prof[, rev(orden)], type = "l", lty = 1, col = adjustcolor(rev(grises), 0.25),
          xlab = "Día del año", ylab = "sqrt(mm/día)", ylim = c(0, 7), main = "Estaciones ordenadas por profundidad FM")
  legend("bottom", c("Más profunda (centro)", "Menos profunda (extremos)"), col = c("grey85", "grey15"),
         lwd = 3, bty = "n", cex = 0.8, ncol = 2, xpd = NA, inset = c(0, -0.42))
  matplot(t_prof, Y_prof, type = "l", lty = 1, col = adjustcolor("grey70", 0.05), ylim = c(0, 4),
          xlab = "Día del año", ylab = "sqrt(mm/día)", main = "Medidas de centralidad (estaciones)")
  lines(t_fino, media_est, lwd = 3, col = "blue")
  lines(t_prof, fm_est$median$data, lwd = 3, col = "red")
  lines(t_prof, fm_est$mtrim$data, lwd = 3, lty = 2, col = "darkgreen")
  legend("bottom", c("Media", "Mediana (FM)", "Media recortada 25 %"), col = c("blue", "red", "darkgreen"),
         lty = c(1, 1, 2), lwd = 3, bty = "n", cex = 0.75, ncol = 2, xpd = NA, inset = c(0, -0.42))
}, alto = 5.2)



#Orografía: media y DE por franja de altitud y por región, y lluvia anual observada en mm

media_franja <- sapply(levels(meta$franja), function(f) rowMeans(Y_eval[, meta$franja == f, drop = FALSE]))
de_franja    <- sapply(levels(meta$franja), function(f) apply(Y_eval[, meta$franja == f, drop = FALSE], 1, sd))
anual_obs <- colSums(Y_obs^2)                                # mm/año observados de cada curva (con 29 feb)
orografia <- data.frame(
  franja = levels(meta$franja),
  estaciones = as.integer(tapply(meta$STATION_ID, meta$franja, function(s) length(unique(s)))),
  curvas = as.vector(table(meta$franja)),
  lluvia_anual_mediana_mm = round(as.vector(tapply(anual_obs, meta$franja, median))),
  pico1 = fecha(apply(media_franja, 2, function(m) fechas_temporadas(m, c(diff(m), 0), t_fino)["pico1"])),
  pico2 = fecha(apply(media_franja, 2, function(m) fechas_temporadas(m, c(diff(m), 0), t_fino)["pico2"])),
  razon_pico_valle = round(apply(media_franja, 2, function(m) {
    f <- fechas_temporadas(m, c(diff(m), 0), t_fino)
    mean(m[t_fino %in% f[c("pico1", "pico2")]]) / m[t_fino == f["valle"]] }), 2))
print(orografia)
guardar(orografia, "tabla19_orografia")
media_region <- sapply(levels(meta$region), function(g) rowMeans(Y_eval[, meta$region == g, drop = FALSE]))
region_tab <- data.frame(
  region = levels(meta$region),
  estaciones = as.integer(tapply(meta$STATION_ID, meta$region, function(s) length(unique(s)))),
  lluvia_anual_mediana_mm = round(as.vector(tapply(anual_obs, meta$region, median))),
  pico1 = fecha(apply(media_region, 2, function(m) fechas_temporadas(m, c(diff(m), 0), t_fino)["pico1"])),
  pico2 = fecha(apply(media_region, 2, function(m) fechas_temporadas(m, c(diff(m), 0), t_fino)["pico2"])),
  razon_pico_valle = round(apply(media_region, 2, function(m) {
    f <- fechas_temporadas(m, c(diff(m), 0), t_fino)
    mean(m[t_fino %in% f[c("pico1", "pico2")]]) / m[t_fino == f["valle"]] }), 2))

# dentro de la cuenca del Cauca: ¿crece la lluvia con la altitud?

cauca <- meta$region == REG_CAUCA
alt_cauca <- data.frame(franja = levels(meta$franja),
  estaciones = as.integer(tapply(meta$STATION_ID[cauca], meta$franja[cauca], function(s) length(unique(s)))),
  lluvia_anual_mediana_mm = round(as.vector(tapply(anual_obs[cauca], meta$franja[cauca], median))))
print(region_tab); print(alt_cauca)
guardar(region_tab, "tabla19b_region")
guardar(alt_cauca, "tabla19c_altitud_dentro_del_cauca")

# Robustez: solo estaciones con altitud de catálogo (las demás vienen de la grilla corregida)

es_cat <- grepl("^cat", meta$fuente_alt)
orografia_cat <- tapply(anual_obs[es_cat], meta$franja[es_cat], median)
cat("Lluvia anual mediana por franja, solo altitud de catálogo:\n"); print(round(orografia_cat))

cols_franja <- c("steelblue", "black", "darkgreen", "orange", "purple")
figura("fase2_fig13_orografia", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  matplot(t_fino, media_franja, type = "l", lty = 1, lwd = 3, col = cols_franja, ylim = c(0, max(media_franja) * 1.2),
          xlab = "Día del año", ylab = "sqrt(mm/día)", main = "Media funcional por franja de altitud")
  legend("topleft", levels(meta$franja), col = cols_franja, lwd = 3, bty = "n", cex = 0.7, ncol = 2)
  matplot(t_fino, de_franja, type = "l", lty = 1, lwd = 3, col = cols_franja, ylim = c(0, max(de_franja) * 1.1),
          xlab = "Día del año", ylab = "sqrt(mm/día)", main = "Desviación estándar por franja de altitud")
})

# ENSO: curva de anomalía media de cada año (unidad de réplica = año) y media por fase
# Promediar primero por año evita contar 160 veces el mismo año: las bandas usan los 23 años.

A_eval  <- eval.fd(t_fino, A_fd)
AX_eval <- eval.fd(t_fino, AX_fd)
anual_and <- anom_anual(A_eval[, meta_an$region == REG_CAUCA], meta_an$YEAR[meta_an$region == REG_CAUCA])
anual_sat <- anom_anual(AX_eval[, meta_an_sat$region == REG_CAUCA], meta_an_sat$YEAR[meta_an_sat$region == REG_CAUCA])
fase_anio <- enso$fase_enso[match(as.integer(colnames(anual_and)), enso$YR)]

banda_fase <- function(M, fase) lapply(c("Niño", "Neutro", "Niña"), function(f) {
  X <- M[, fase == f, drop = FALSE]
  list(media = rowMeans(X), ee = apply(X, 1, sd) / sqrt(ncol(X)), n = ncol(X)) })
bf_and <- banda_fase(anual_and, fase_anio)
bf_sat <- banda_fase(anual_sat, enso$fase_enso[match(as.integer(colnames(anual_sat)), enso$YR)])

cols_enso <- c("tomato", "grey40", "steelblue")
figura("fase2_fig14_enso", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  for (k in 1:2) {
    bf <- if (k == 1) bf_and else bf_sat
    plot(NA, xlim = rangval, ylim = c(-0.6, 0.6), xlab = "Día del año", ylab = "Anomalía sqrt(mm/día)",
         main = if (k == 1) "Estaciones (cuenca del río Cauca) por fase" else "Satélite (píxeles de la cuenca del río Cauca)")
    for (f in 1:3) {
      polygon(c(t_fino, rev(t_fino)), c(bf[[f]]$media + 2 * bf[[f]]$ee, rev(bf[[f]]$media - 2 * bf[[f]]$ee)),
              col = adjustcolor(cols_enso[f], 0.15), border = NA)
      lines(t_fino, bf[[f]]$media, lwd = 3, col = cols_enso[f])
    }
    abline(h = 0, lty = 3)
    legend("topright", paste0(c("Niño", "Neutro", "Niña"), " (", sapply(bf, `[[`, "n"), " años)"),
           col = cols_enso, lwd = 3, bty = "n", cex = 0.8)
  }
})

# ENSO en mm con los datos observados: anomalía porcentual de la lluvia por trimestre
# (cada estación contra su propio promedio; se promedia por año y luego por fase)

trim_dia <- c(DEF = 1, MAM = 2, JJA = 3, SON = 4)[c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)][mes_dia]
lluvia_trim <- sapply(1:4, function(q) colSums(Y_obs[trim_dia == q, ]^2))      # N x 4 (mm)
colnames(lluvia_trim) <- c("DEF", "MAM", "JJA", "SON")
lluvia_trim <- cbind(lluvia_trim, anual = rowSums(lluvia_trim))
clim_est <- apply(lluvia_trim, 2, function(v) tapply(v, meta$STATION_ID, mean))
anom_pct <- 100 * (lluvia_trim / clim_est[as.character(meta$STATION_ID), ] - 1)
anom_pct <- anom_pct[an_est$usar & meta$region == REG_CAUCA, ]
anios_and <- meta$YEAR[an_est$usar & meta$region == REG_CAUCA]
anom_pct_anio <- apply(anom_pct, 2, function(v) tapply(v, anios_and, mean))
enso_mm <- aggregate(anom_pct_anio, list(fase_enso = enso$fase_enso[match(as.integer(rownames(anom_pct_anio)), enso$YR)]), mean)
enso_mm[, -1] <- round(enso_mm[, -1], 1)
print(enso_mm)
guardar(enso_mm, "tabla20_enso_anomalia_pct_por_trimestre")

# Lo mismo para el satélite en los píxeles de la cuenca del Cauca (datos observados por péntada)

p$region <- meta_sat$region[match(paste(p$pixel, p$YEAR), paste(meta_sat$pixel, meta_sat$YEAR))]
p$trim <- c("DEF", "DEF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DEF")[p$mes]
pc <- p[p$region == REG_CAUCA, ]
tot_sat <- tapply(pc$satelite_mm, list(paste(pc$pixel, pc$YEAR), pc$trim), sum)[, c("DEF", "MAM", "JJA", "SON")]
tot_sat <- cbind(tot_sat, anual = rowSums(tot_sat))
px_sat <- sub(" .*", "", rownames(tot_sat)); anio_sat <- as.integer(sub(".* ", "", rownames(tot_sat)))
clim_px <- apply(tot_sat, 2, function(v) tapply(v, px_sat, mean))
an_sat_pct <- 100 * (tot_sat / clim_px[px_sat, ] - 1)
an_sat_anio <- apply(an_sat_pct, 2, function(v) tapply(v, anio_sat, mean))
enso_mm_sat <- aggregate(an_sat_anio, list(fase_enso = enso$fase_enso[match(as.integer(rownames(an_sat_anio)), enso$YR)]), mean)
enso_mm_sat[, -1] <- round(enso_mm_sat[, -1], 1)
cat("ENSO en mm: estaciones (arriba) vs satélite (abajo), cuenca del río Cauca:\n"); print(enso_mm); print(enso_mm_sat)
guardar(enso_mm_sat, "tabla20b_enso_anomalia_pct_satelite")
cat("Correlación año a año de la anomalía anual en mm, estaciones vs satélite (cuenca del río Cauca):",
    round(cor(anom_pct_anio[, "anual"], an_sat_anio[rownames(anom_pct_anio), "anual"]), 3), "\n")
anom_anio_txt <- data.frame(YEAR = as.integer(rownames(anom_pct_anio)), round(anom_pct_anio, 1))
anom_anio_txt$satelite_anual <- round(an_sat_anio[as.character(anom_anio_txt$YEAR), "anual"], 1)
anom_anio_txt$fase_enso <- enso$fase_enso[match(anom_anio_txt$YEAR, enso$YR)]
guardar(anom_anio_txt, "tabla21_anomalia_pct_por_anio_cauca")


#Resumen 

resumen <- data.frame(
  fuente = c("Estaciones", "Satélite"),
  curvas = c(N, N_sat),
  K = c(K_est, K_sat),
  lluvia_media_mm_dia_observada = round(c(mean(Y_obs^2), sum(p$satelite_mm) / sum(p$dias)), 2),
  media_max_sqrt = round(c(max(media_est), max(media_sat)), 3),
  pico1 = c(fecha(f_media_est["pico1"]), fecha(f_media_sat["pico1"])),
  pico2 = c(fecha(f_media_est["pico2"]), fecha(f_media_sat["pico2"])),
  DE_max_interior = tabla_de$DE_max_interior,
  fecha_DE_max = tabla_de$fecha_DE_max,
  varianza_modo1_pct = c(var_modos$estaciones_pct[1], var_modos$satelite_pct[1]),
  atipicas_crudo = c(length(bx_est$outpoint), length(bx_sat$outpoint)),
  atipicas_anomalias = c(length(bxa_est$outpoint), length(bxa_sat$outpoint)),
  anios_atipicos = c(paste(anios_tab$YEAR[anios_tab$atipico_estaciones], collapse = " "),
                     paste(anios_tab$YEAR[anios_tab$atipico_satelite], collapse = " ")))
print(resumen)
guardar(resumen, "tabla22_resumen")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_fase2.txt"))
cat("\nListo. Tablas en", OUT, "y figuras en", FIG, "\n")
