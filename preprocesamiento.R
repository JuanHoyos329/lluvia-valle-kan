# =============================================================================
#  Analítica de Datos 2026-2S, proyecto final
# Juan Hoyos, Daniela Marin, Miguel Arroyo, Annie Correa, Fabian Gomezcasseres
# =============================================================================
install.packages("writexl")
install.packges("terra")

library("writexl")
library("terra")

setwd("C:/Users/juana/OneDrive/Escritorio/lluvia-valle-kan")

# ---- 0.1 Paquetes ------------------------------------------------------------
# Solo se usan tres paquetes y únicamente para LEER y ESCRIBIR archivos. Todo
# el cálculo (incluidos la FPCA y el VAE) se hace con R base.
for (paquete in c("readxl", "terra", "writexl")) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop("Falta el paquete '", paquete, "'. Instálelo con install.packages('", paquete, "')")
  }
}

# ---- 0.2 Rutas -----------------------------------------------------------------
RUTA_DATOS   <- "Datos"
RUTA_SALIDAS <- "salidas"
CARPETAS <- c(pre_eda = "pre_eda", prep = "preprocesamiento", limpieza = "limpieza",
              imputacion = "imputacion")
for (carpeta in CARPETAS) dir.create(file.path(RUTA_SALIDAS, carpeta), recursive = TRUE, showWarnings = FALSE)
options(scipen = 999)        # mostrar números sin notación científica

# ---- 0.3 Parámetros --------------------------------------------------------
# Todos los umbrales están aquí arriba para poder cambiarlos y justificarlos.

# Periodo
ANIO_INICIO     <- 1970      # antes de 1970 hay muy pocas estaciones (ver Parte 3.1)
ANIO_FIN        <- 2003      # último año de estaciones (indicado por el profesor)
PERIODO_COMUN   <- 1981:2003 # años en que hay estaciones Y satélite
ANIO_REFERENCIA <- ANIO_FIN  # último año de ENTRENAMIENTO: lo "normal" y los modelos de
# imputación se aprenden solo hasta aquí. Todavía no hay
# partición oficial, así que se usan todos los años. Cuando
# el profesor la defina, basta con cambiar este número para
# que ninguna decisión use datos de prueba.

# Valores imposibles y exagerados (Parte 4.5)
LIMITE_FISICO     <- 1825    # mm en un día: récord mundial de lluvia en 24 h (OMM, La
# Reunión, 1966). Más que esto, o negativo, es imposible.
FACTOR_EXAGERADO  <- 3       # un valor es sospechoso si supera 3 veces el percentil 99
# de los días con lluvia de su estación...
RESPALDO_MIN      <- .25     # ...y ninguna vecina (día anterior, mismo día o siguiente)
# registró al menos el 25 % de ese valor
FACTOR_SIN_VECINAS <- 5      # si ese día no hay vecinas con dato: exagerado si supera 5 veces el p99

# Rachas de ceros
RACHA_CEROS     <- 30        # días seguidos en cero para considerar una racha
LLUVIA_RED      <- 1         # mm: desde aquí un día cuenta como "día con lluvia"
RADIO_KM        <- 20        # se buscan estaciones vecinas a esta distancia...
RADIO_MAX_KM    <- 30        # ...y si hay menos de MIN_VECINAS se amplía hasta aquí
MIN_VECINAS     <- 3
DIF_ALTURA      <- 500       # m: vecinas con más diferencia de altura no se usan

# Acumulados y estaciones mensuales
HUECO_MIN       <- 2         # días sin dato antes de un valor para sospechar acumulado
PERCENTIL_ACUM  <- .99       # y el valor debe superar este percentil de la estación
MAX_DIAS_MES    <- 5         # mediana de días con dato por mes <= esto: lectura mensual

# Años aptos para construir curvas
COB_MIN         <- 90        # % mínimo de días con dato en el año
HUECO_MAX       <- 30        # hueco más largo permitido dentro del año (días)
COB_TEMPORADA   <- 80        # % mínimo en cada temporada de lluvias (mar-may y sep-nov)

# Vecinas por correlación para estaciones sin coordenadas (Parte 4.3)
MIN_DIAS_COMUNES <- 365      # días con dato en ambas estaciones para calcular su correlación
N_VECINAS_CORR   <- 5        # vecinas por correlación: las 5 estaciones más parecidas
# (validado en 4.3: con 5 quedan a ~12 km reales en mediana)
R_MIN_VECINA     <- .3       # correlación mínima para aceptar una vecina

# Imputación (Parte 5)
MIN_PENTADAS_IMPUTAR <- 36   # una curva se imputa si tiene al menos la mitad de sus 72 péntadas
SEMILLA          <- 2026     # para que los sorteos (validación, VAE) den siempre lo mismo
VAE_OCULTA       <- 32       # neuronas de la capa oculta del VAE
VAE_LATENTE      <- 4        # dimensión del espacio latente z
VAE_TASA         <- 1e-3     # tasa de aprendizaje (Adam)
VAE_LOTE         <- 64       # curvas por lote
VAE_EPOCAS_MAX   <- 300
REPETICIONES_VALIDACION <- 20 # sorteos de péntadas escondidas para comparar los métodos

# ---- 0.4 Colores de las gráficas ---------------------------------------------
COL_PRINCIPAL <- "#1f6fb2"   # azul: dato normal
COL_CERO      <- "#c77c12"   # ocre: ceros / sequía
COL_ROJO      <- "#bf3636"   # rojo: crítico / sospechoso
COL_GRIS      <- "#b8c4cc"   # gris: revisar / sin dato
COLORES_COB   <- c("#bf3636", "#e8a33d", "#9cc3e6", COL_PRINCIPAL)
CORTES_COB    <- c(-1, 50, 80, 90, 100.1)
NIVELES_COB   <- c("< 50 %", "50-80 %", "80-90 %", ">= 90 %")

# Clases de las rachas de ceros (ver Parte 3.5)
CLASES <- c("Probable sequía: la red también estuvo seca", "Revisar", "Sospechoso: la red llovió normal")
CLASES_CORTAS <- c("Probable sequía", "Revisar", "Sospechoso")

# ---- 0.5 Funciones de apoyo --------------------------------------------------

# Muestra una gráfica en pantalla. Recibe una función que dibuja.
graficar <- function(dibujar) invisible(dibujar())

# Guarda una tabla en CSV con punto y coma y coma decimal (abre bien en Excel en español).
guardar_tabla <- function(x, carpeta, nombre) {
  write.csv2(x, file.path(RUTA_SALIDAS, CARPETAS[[carpeta]], paste0(nombre, ".csv")), row.names = FALSE)
}

# Número de días de un año (366 si es bisiesto).
dias_del_anio <- function(a) ifelse((a %% 4 == 0 & a %% 100 != 0) | a %% 400 == 0, 366, 365)

# Días que dura una péntada. El satélite trae 6 péntadas por mes: las cinco
# primeras duran 5 días y la sexta completa el mes (de 3 a 6 días).
dias_pentada <- function(anio, pentada) {
  mes <- ceiling(pentada / 6)
  posicion <- pentada - (mes - 1) * 6
  siguiente <- as.Date(sprintf("%d-%02d-01", ifelse(mes == 12, anio + 1, anio),
                               ifelse(mes == 12, 1, mes + 1)), format = "%Y-%m-%d")
  dias_mes <- as.integer(format(siguiente - 1, "%d"))
  ifelse(posicion < 6, 5, dias_mes - 25)
}

# Moda: el valor que más se repite.
moda <- function(x) {
  frecuencias <- table(x)
  as.numeric(names(frecuencias)[which.max(frecuencias)])
}

# Estadísticos básicos de un vector. "esperados" es cuántos datos debería haber,
# para poder contar los nulos (los días que faltan no vienen como fila).
describir <- function(x, esperados = NA) {
  datos <- x[!is.na(x)]
  if (!length(datos)) {                     # sin ningún dato: solo se cuentan los nulos
    return(c(esperados = esperados, con_dato = 0, nulos = esperados, pct_nulos = 100,
             ceros = NA, pct_ceros = NA, media = NA, mediana = NA, moda = NA, desv_estandar = NA,
             minimo = NA, q1 = NA, q3 = NA, p99 = NA, maximo = NA))
  }
  c(esperados     = esperados,
    con_dato      = length(datos),
    nulos         = esperados - length(datos),
    pct_nulos     = (esperados - length(datos)) / esperados * 100,
    ceros         = sum(datos == 0),
    pct_ceros     = mean(datos == 0) * 100,
    media         = mean(datos),
    mediana       = median(datos),
    moda          = moda(round(datos, 2)),
    desv_estandar = sd(datos),
    minimo        = min(datos),
    q1            = quantile(datos, .25, names = FALSE),
    q3            = quantile(datos, .75, names = FALSE),
    p99           = quantile(datos, .99, names = FALSE),
    maximo        = max(datos))
}

# Une años seguidos en rangos: c(1991, 1992, 1993, 2001) -> "1991-1993, 2001".
rangos <- function(a) {
  if (!length(a)) return("")
  grupos <- split(a, cumsum(c(1, diff(a) != 1)))
  paste(sapply(grupos, function(x) if (length(x) > 1) paste0(min(x), "-", max(x)) else x),
        collapse = ", ")
}

# Distancia en km entre dos puntos (longitud, latitud) sobre la Tierra (haversine).
dist_km <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(sqrt(a))
}

# Largo del tramo más largo de TRUE seguidos (x = TRUE cuando falta el dato).
racha_mas_larga <- function(x) {
  r <- rle(x)
  if (any(r$values)) max(r$lengths[r$values]) else 0
}

# Regla de valor imposible: negativo o por encima del récord mundial en 24 h.
es_imposible <- function(x) !is.na(x) & (x < 0 | x > LIMITE_FISICO)

# Trazabilidad: cada paso que cambia los datos deja una fila con cuántas
# estaciones, pares estación-año y datos quedan. Así se ve qué se perdió y
# dónde. "unidad" dice si se cuentan días (datos diarios) o péntadas (curvas).
trazabilidad <- data.frame()
registrar <- function(paso, id, anio, valor, unidad = "días") {
  con_dato <- !is.na(valor)
  clave <- paste(id, anio)
  fila <- data.frame(
    paso                = paso,
    unidad              = unidad,
    estaciones          = length(unique(id[con_dato])),
    pares_estacion_anio = length(unique(clave[con_dato])),
    con_dato            = sum(con_dato),
    sin_dato            = sum(!con_dato),
    pct_sin_dato        = round(mean(!con_dato) * 100, 2))
  trazabilidad <<- rbind(trazabilidad, fila)
  print(fila, row.names = FALSE)
}


# =============================================================================
# PARTE 1. LECTURA DE LOS DATOS
# =============================================================================
cat("\n######## PARTE 1. LECTURA ########\n")

# ---- 1.1 Estaciones ----------------------------------------------------------
# Vienen en tres archivos: un Excel 1930-1979 y dos CSV 1980-1990 y 1991-2003.
# El código de estación (10 dígitos) se lee como texto para no perder ceros
# ni que R lo escriba en notación científica.
leer_csv <- function(archivo) {
  read.csv(file.path(RUTA_DATOS, archivo), colClasses = c(STATION_ID = "character"))
}
historico <- as.data.frame(readxl::read_excel(
  file.path(RUTA_DATOS, "P_1930_1979.xlsx"), sheet = "Hoja1",
  col_types = c("numeric", "text", "numeric", "numeric", "numeric", "numeric", "text")))
historico$STATION_ID <- sprintf("%.0f", historico$STATION_ID)

columnas <- c("STATION_ID", "YEAR", "MONTH", "DAY", "VALUE")
est <- rbind(historico[, columnas],
             leer_csv("P_1980_1990.csv")[, columnas],
             leer_csv("P_1991_2003.csv")[, columnas])
# Solo el Excel trae unidad y sensor; se guardan para revisar unidades (Parte 3.2)
est$UNIDAD <- c(historico$UNITS, rep(NA_character_, nrow(est) - nrow(historico)))
est$SENSOR <- c(historico$SENSORNAME, rep(NA_character_, nrow(est) - nrow(historico)))
rm(historico)
est$fecha <- as.Date(sprintf("%04d-%02d-%02d", est$YEAR, est$MONTH, est$DAY), format = "%Y-%m-%d")
cat(sprintf("Estaciones: %d filas, %d estaciones, %d-%d\n",
            nrow(est), length(unique(est$STATION_ID)), min(est$YEAR), max(est$YEAR)))

# ---- 1.2 Satélite --------------------------------------------------------------
# Una imagen GeoTIFF por péntada (archivo pptAAAAPP.tif: año y péntada 01-72).
# Si la carpeta no existe, se descomprime el zip (solo los .tif reales, sin
# los archivos basura "._" que deja macOS).
carpeta_sat <- file.path(RUTA_DATOS, "ImagenesSatelitePentadal")
if (!dir.exists(carpeta_sat)) {
  zip_sat <- paste0(carpeta_sat, ".zip")
  en_zip <- unzip(zip_sat, list = TRUE)$Name
  unzip(zip_sat, files = grep("/ppt[0-9]{6}\\.tif$", en_zip, value = TRUE), exdir = RUTA_DATOS)
}
archivos_sat <- sort(list.files(carpeta_sat, pattern = "^ppt[0-9]{6}\\.tif$", full.names = TRUE))
anio_sat <- as.integer(substr(basename(archivos_sat), 4, 7))
pent_sat <- as.integer(substr(basename(archivos_sat), 8, 9))
mes_sat  <- ceiling(pent_sat / 6)
dias_sat <- dias_pentada(anio_sat, pent_sat)

# M: matriz píxel x péntada (filas = 1443 píxeles, columnas = imágenes), en mm/péntada.
plantilla_sat <- terra::rast(archivos_sat[1])          # para coordenadas de los píxeles
M <- unname(sapply(archivos_sat, function(f) terra::values(terra::rast(f))[, 1]))
M[M == -9999] <- NA                                    # -9999 = sin dato
M_dia <- sweep(M, 2, dias_sat, "/")                    # mm/día: se divide por los días de la péntada
xy_pixel <- terra::xyFromCell(plantilla_sat, seq_len(nrow(M)))
cat(sprintf("Satélite: %d imágenes, %d-%d, %d píxeles\n",
            ncol(M), min(anio_sat), max(anio_sat), nrow(M)))


# =============================================================================
# PARTE 2. PRE-EDA: CÓMO LLEGAN LOS DATOS (no se cambia nada)
# =============================================================================
# Objetivo: describir los datos tal como vienen, desde 1930, para decidir
# después qué hay que corregir. Esta parte también justifica el corte de años.
cat("\n######## PARTE 2. PRE-EDA ########\n")

# ---- 2.1 Días esperados ------------------------------------------------------
# Un día faltante NO viene como NA: simplemente no hay fila. Por eso los días
# esperados se cuentan desde la primera hasta la última fecha de cada estación.
valido <- !is.na(est$VALUE) & !es_imposible(est$VALUE)
rango_est <- tapply(as.integer(est$fecha), est$STATION_ID, function(x) max(x) - min(x) + 1)
esperados_est <- sum(rango_est)
anios_sat <- sort(unique(anio_sat))
esperados_sat <- nrow(M) * 72 * length(anios_sat)      # 72 péntadas por año

# ---- 2.2 Tabla descriptiva ---------------------------------------------------
v_est <- est$VALUE[valido]
v_sat <- as.vector(M_dia)
descriptiva <- rbind(
  "Estaciones (mm/día)"                = describir(est$VALUE, esperados_est),
  "Estaciones sin valor imposible"     = describir(v_est, esperados_est),
  "Estaciones, solo días con lluvia"   = describir(v_est[v_est > 0]),
  "Satélite (mm/péntada)"              = describir(as.vector(M), esperados_sat),
  "Satélite (mm/día)"                  = describir(v_sat, esperados_sat),
  "Satélite, solo péntadas con lluvia" = describir(v_sat[!is.na(v_sat) & v_sat > 0]))
cat("\n== 2.2 TABLA DESCRIPTIVA ==\n")
print(round(t(descriptiva), 2))
guardar_tabla(data.frame(fuente = rownames(descriptiva), round(descriptiva, 3)), "pre_eda", "tabla_descriptiva")

# ---- 2.3 Cobertura de estaciones por año (1930-2003) -------------------------
# Cobertura = % de los días del año que tienen dato.
# m_dias_pre: matriz estación x año con el número de días con dato.
anios_pre <- seq(min(est$YEAR), max(est$YEAR))
con_dato <- !is.na(est$VALUE)
m_dias_pre <- unclass(table(est$STATION_ID[con_dato], factor(est$YEAR[con_dato], levels = anios_pre)))
m_cob_pre  <- sweep(m_dias_pre, 2, dias_del_anio(anios_pre), "/") * 100
en_comun   <- anios_pre %in% PERIODO_COMUN

# Estado de cada estación según su cobertura en el periodo común 1981-2003
resumen_estaciones <- data.frame(
  STATION_ID        = rownames(m_dias_pre),
  primer_anio       = anios_pre[apply(m_dias_pre > 0, 1, function(x) min(which(x)))],
  ultimo_anio       = anios_pre[apply(m_dias_pre > 0, 1, function(x) max(which(x)))],
  anios_con_dato    = rowSums(m_dias_pre > 0),
  cobertura_periodo = rowSums(m_dias_pre[, en_comun]) / sum(dias_del_anio(PERIODO_COMUN)) * 100,
  anios_periodo_80  = rowSums(m_cob_pre[, en_comun] >= 80),
  pct_ceros         = as.vector(tapply(est$VALUE == 0, est$STATION_ID, mean, na.rm = TRUE)[rownames(m_dias_pre)]) * 100,
  stringsAsFactors  = FALSE)
resumen_estaciones$estado <- as.character(cut(
  resumen_estaciones$cobertura_periodo, breaks = c(-1, 0, 50, 80, 90, 100),
  labels = c("Sin datos en el periodo", "Crítica", "Regular", "Buena", "Completa")))
resumen_estaciones <- resumen_estaciones[order(-resumen_estaciones$cobertura_periodo), ]
rownames(resumen_estaciones) <- NULL

guardar_tabla(data.frame(STATION_ID = rownames(m_cob_pre), round(m_cob_pre, 1), check.names = FALSE),
              "pre_eda", "cobertura_estacion_anio_1930_2003")
guardar_tabla(transform(resumen_estaciones, cobertura_periodo = round(cobertura_periodo, 1),
                        pct_ceros = round(pct_ceros, 1)), "pre_eda", "resumen_estaciones")

con_periodo <- resumen_estaciones[resumen_estaciones$cobertura_periodo > 0, ]
cat(sprintf("\n== 2.3 COBERTURA DE ESTACIONES %d-%d ==\n", min(PERIODO_COMUN), max(PERIODO_COMUN)))
print(table(factor(resumen_estaciones$estado,
                   levels = c("Completa", "Buena", "Regular", "Crítica", "Sin datos en el periodo"))))
cat("\nLas 15 más completas:\n")
print(head(con_periodo[, c("STATION_ID", "cobertura_periodo", "anios_periodo_80")], 15), row.names = FALSE)
cat("\nLas 15 más críticas (con algún dato en el periodo):\n")
print(tail(con_periodo[, c("STATION_ID", "cobertura_periodo", "anios_periodo_80")], 15), row.names = FALSE)

# Mapa de calor: cada fila es una estación, cada columna un año, el color es la cobertura
mapa_cobertura <- function(z, anios, titulo, eje_y) {
  z[z == 0] <- NA
  paleta <- hcl.colors(10, "Blues 3", rev = TRUE)
  image(x = anios, y = seq_len(nrow(z)), z = t(z), col = paleta, breaks = seq(0, 100, 10),
        xlab = "Año", ylab = eje_y, main = titulo)
  legend("topleft", fill = paleta[c(2, 5, 10)], legend = c("20 %", "50 %", "100 %"),
         bg = "white", box.col = NA, horiz = TRUE)
}
graficar(function() {   # Gráfica: cobertura por estación y año desde 1930
  op <- par(mar = c(4, 4, 3, 1)); on.exit(par(op))
  mapa_cobertura(m_cob_pre[rev(resumen_estaciones$STATION_ID), ], anios_pre,
                 "Cobertura por estación y año (arriba, las más completas en 1981-2003)", "Estación")
  abline(v = min(PERIODO_COMUN) - .5, lty = 2, lwd = 2)
})
graficar(function() {   # Gráfica: las 15 estaciones más completas y las 15 más críticas
  op <- par(mfrow = c(1, 2), mar = c(4, 7, 3, 1)); on.exit(par(op))
  arriba <- head(con_periodo, 15)
  abajo  <- tail(con_periodo, 15)
  barplot(rev(arriba$cobertura_periodo), names.arg = rev(arriba$STATION_ID), horiz = TRUE, las = 1,
          cex.names = .65, col = COL_PRINCIPAL, border = NA, xlim = c(0, 100),
          main = "15 más completas", xlab = "% de días con dato, 1981-2003")
  barplot(abajo$cobertura_periodo, names.arg = abajo$STATION_ID, horiz = TRUE, las = 1,
          cex.names = .65, col = COL_ROJO, border = NA, xlim = c(0, 100),
          main = "15 más críticas", xlab = "% de días con dato, 1981-2003")
})

# ---- 2.4 Cobertura del satélite por píxel y año ------------------------------
m_sat <- sapply(anios_sat, function(a) rowSums(!is.na(M[, anio_sat == a, drop = FALSE])))
colnames(m_sat) <- anios_sat
resumen_pixeles <- data.frame(
  pixel     = sprintf("P%04d", seq_len(nrow(M))), lon = xy_pixel[, 1], lat = xy_pixel[, 2],
  cobertura = rowSums(m_sat) / (72 * length(anios_sat)) * 100,
  pct_ceros = rowMeans(M == 0, na.rm = TRUE) * 100)
guardar_tabla(transform(resumen_pixeles, cobertura = round(cobertura, 1), pct_ceros = round(pct_ceros, 2)),
              "pre_eda", "resumen_pixeles")
cat("\n== 2.4 COBERTURA DEL SATÉLITE: péntadas por año ==\n")
print(table(anio_sat))
cat(sprintf("Años incompletos: %s\n", paste(names(which(table(anio_sat) < 72)), collapse = ", ")))
graficar(function() {   # Gráfica: cobertura por píxel y año
  op <- par(mar = c(4, 4, 3, 1)); on.exit(par(op))
  mapa_cobertura(m_sat / 72 * 100, anios_sat, "Cobertura del satélite por píxel y año", "Píxel")
})

# ---- 2.5 Comportamiento de la lluvia contando los ceros ----------------------
# Un cero no es un dato faltante: es un día sin lluvia. Aquí se ve cuánto pesan.
TRAMOS <- c("0", "0-1", "1-5", "5-10", "10-20", "20-50", "50-100", ">100")
por_tramo <- function(x) {
  x <- x[!is.na(x)]
  prop.table(table(cut(x, c(-Inf, 0, 1, 5, 10, 20, 50, 100, Inf), labels = TRAMOS))) * 100
}
ciclo <- function(x, mes) {        # lluvia media por mes, con y sin ceros
  list(con_ceros = tapply(x, mes, mean, na.rm = TRUE),
       sin_ceros = tapply(x[x > 0 & !is.na(x)], mes[x > 0 & !is.na(x)], mean))
}
ciclo_est <- ciclo(v_est, est$MONTH[valido])
ciclo_sat <- ciclo(v_sat, rep(mes_sat, each = nrow(M)))
tope_ciclo <- max(unlist(ciclo_est), unlist(ciclo_sat)) * 1.1
graficar(function() {   # Gráfica: reparto de los días por tramos y ciclo anual
  op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  panel_tramos <- function(x, titulo, eje) {
    barplot(por_tramo(x), col = c(COL_CERO, rep(COL_PRINCIPAL, 7)), border = NA, ylim = c(0, 100),
            main = titulo, xlab = "mm/día", ylab = eje, cex.names = .85)
  }
  panel_ciclo <- function(cc, titulo) {
    plot(1:12, cc$sin_ceros, type = "o", pch = 16, col = COL_PRINCIPAL, lwd = 2, ylim = c(0, tope_ciclo),
         xaxt = "n", xlab = "", ylab = "mm/día", main = titulo)
    lines(1:12, cc$con_ceros, type = "o", pch = 16, col = COL_CERO, lwd = 2)
    axis(1, at = 1:12, labels = c("E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
    legend("topleft", col = c(COL_PRINCIPAL, COL_CERO), lwd = 2, pch = 16, bty = "n",
           legend = c("Solo con lluvia", "Contando los ceros"))
  }
  panel_tramos(v_est, "Estaciones: reparto de los días", "% de días")
  panel_ciclo(ciclo_est, "Estaciones: lluvia media por mes")
  panel_tramos(v_sat, "Satélite: reparto de las péntadas", "% de péntadas")
  panel_ciclo(ciclo_sat, "Satélite: lluvia media por mes")
})

# ---- 2.6 Porcentaje de ceros por estación y por píxel ------------------------
cat("\n== 2.6 PORCENTAJE DE CEROS ==\n")
print(summary(resumen_estaciones$pct_ceros))
print(summary(resumen_pixeles$pct_ceros))
graficar(function() {   # Gráfica: % de ceros, de menor a mayor
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  panel <- function(x, titulo, eje) {
    x <- sort(x)
    plot(seq_along(x), x, type = "h", col = COL_CERO, ylim = c(0, 100), main = titulo,
         xlab = eje, ylab = "% de registros en cero")
    abline(h = mean(x), lty = 2, lwd = 2)
    legend("topleft", lty = 2, lwd = 2, bty = "n", legend = sprintf("Media: %.1f %%", mean(x)))
  }
  panel(resumen_estaciones$pct_ceros, "Días en cero por estación", "Estaciones, de menor a mayor")
  panel(resumen_pixeles$pct_ceros, "Péntadas en cero por píxel", "Píxeles, de menor a mayor")
})
rm(v_est, v_sat, valido, con_dato)


# =============================================================================
# PARTE 3. PREPROCESAMIENTO
# =============================================================================
cat("\n######## PARTE 3. PREPROCESAMIENTO ########\n")
registrar("1. Datos originales (1930-2003)", est$STATION_ID, est$YEAR, est$VALUE)

# ---- 3.1 Corte de años: se quitan 1930-1969 ----------------------------------
# Justificación: se cuenta cuántas estaciones reportan cada año desde 1930 y
# cuántas tienen 80 % de los días o más (sale de la matriz del Pre-EDA).
red_por_anio <- data.frame(
  anio                = anios_pre,
  estaciones_con_dato = colSums(m_dias_pre > 0),
  estaciones_80       = colSums(m_cob_pre >= 80),
  dias_con_dato       = colSums(m_dias_pre))
guardar_tabla(red_por_anio, "prep", "red_por_anio_1930_2003")
decada <- paste0(anios_pre %/% 10 * 10, "s")
cat("\n== 3.1 CORTE DE AÑOS: promedio de estaciones por década ==\n")
print(data.frame(
  estaciones_con_dato = round(tapply(red_por_anio$estaciones_con_dato, decada, mean)),
  estaciones_80       = round(tapply(red_por_anio$estaciones_80, decada, mean))))

graficar(function() {   # Gráfica: estaciones que reportan cada año, con la línea del corte
  op <- par(mar = c(4, 5, 4, 1)); on.exit(par(op))
  bp <- barplot(red_por_anio$estaciones_con_dato, names.arg = anios_pre, col = "#c9d9ea",
                border = NA, las = 2, cex.names = .65, ylab = "estaciones",
                main = "Estaciones que reportan cada año, 1930-2003")
  barplot(red_por_anio$estaciones_80, col = COL_PRINCIPAL, border = NA, add = TRUE, axes = FALSE)
  corte <- mean(bp[anios_pre %in% c(ANIO_INICIO - 1, ANIO_INICIO)])
  abline(v = corte, lty = 2, lwd = 2)
  text(corte, max(red_por_anio$estaciones_con_dato) * .6, sprintf("corte en %d", ANIO_INICIO), pos = 2, cex = .9)
  legend("topleft", fill = c("#c9d9ea", COL_PRINCIPAL), bty = "n",
         legend = c("con algún dato", "con 80 % o más de los días"))
})

# El corte: desde aquí "est" solo tiene 1970-2003
est <- est[est$YEAR >= ANIO_INICIO & est$YEAR <= ANIO_FIN, ]
registrar(sprintf("2. Quitar años antes de %d", ANIO_INICIO), est$STATION_ID, est$YEAR, est$VALUE)

# ---- 3.2 Unidades de las estaciones ------------------------------------------
# Solo el Excel dice la unidad. Para los CSV se revisa la magnitud: la mediana
# de un día lluvioso de cada estación contra la de su zona. Si una estación
# estuviera en cm o en décimas de mm, quedaría unas 10 veces por debajo o encima.
cat("\n== 3.2 UNIDADES DE LAS ESTACIONES ==\n")
print(table(UNIDAD = est$UNIDAD, SENSOR = est$SENSOR, useNA = "ifany"))

# Zona según los dos primeros dígitos del código (26 = una zona; 52, 53, 54 = otra)
est$zona <- ifelse(substr(est$STATION_ID, 1, 2) == "26", "Código 26", "Códigos 52-54")
lluvioso <- !is.na(est$VALUE) & est$VALUE >= LLUVIA_RED & !es_imposible(est$VALUE)
mediana_lluvia <- tapply(est$VALUE[lluvioso], est$STATION_ID[lluvioso], median)
zona_est <- tapply(est$zona, est$STATION_ID, `[`, 1)        # zona de cada estación
mediana_zona <- tapply(mediana_lluvia, zona_est[names(mediana_lluvia)], median)

unidades <- data.frame(
  STATION_ID           = names(mediana_lluvia),
  zona                 = zona_est[names(mediana_lluvia)],
  mediana_dia_lluvioso = as.vector(mediana_lluvia),
  pct_decimales        = as.vector(tapply(est$VALUE != floor(est$VALUE), est$STATION_ID, mean,
                                          na.rm = TRUE)[names(mediana_lluvia)]) * 100,
  stringsAsFactors     = FALSE)
unidades$razon_con_su_zona <- unidades$mediana_dia_lluvioso / mediana_zona[unidades$zona]
unidades$revisar <- unidades$razon_con_su_zona > 4 | unidades$razon_con_su_zona < .25
guardar_tabla(transform(unidades, razon_con_su_zona = round(razon_con_su_zona, 2),
                        pct_decimales = round(pct_decimales, 1)), "prep", "unidades_estaciones")
cat(sprintf("Mediana de un día lluvioso por zona: %s\n",
            paste(names(mediana_zona), round(mediana_zona, 1), "mm", collapse = " | ")))
cat(sprintf("Estaciones con magnitud 4 veces mayor o menor que su zona: %d\n", sum(unidades$revisar)))
if (any(unidades$revisar)) print(unidades[unidades$revisar, ], row.names = FALSE)

# ---- 3.3 Calendario completo -------------------------------------------------
# Cada estación pasa a tener TODOS los días entre su primer y su último registro.
# Los días que no venían en el archivo quedan con VALUE = NA. Sin esto no se
# pueden contar huecos ni rachas.
cat("\n== 3.3 CALENDARIO COMPLETO ==\n")
est <- est[order(est$STATION_ID, est$fecha), ]
por_estacion <- lapply(split(seq_len(nrow(est)), est$STATION_ID), function(i) {
  dias  <- as.integer(est$fecha[i])
  todos <- seq(min(dias), max(dias))
  valor <- rep(NA_real_, length(todos))
  valor[dias - todos[1] + 1L] <- est$VALUE[i]
  list(dia = todos, valor = valor)
})
cal <- data.frame(
  STATION_ID = rep(names(por_estacion), sapply(por_estacion, function(x) length(x$dia))),
  fecha      = as.Date(unlist(lapply(por_estacion, `[[`, "dia"), use.names = FALSE)),
  VALUE      = unlist(lapply(por_estacion, `[[`, "valor"), use.names = FALSE),
  stringsAsFactors = FALSE)
rm(por_estacion)
cal$anio <- as.integer(format(cal$fecha, "%Y"))
cal$mes  <- as.integer(format(cal$fecha, "%m"))
cal$dia  <- as.integer(format(cal$fecha, "%d"))
cal$zona <- zona_est[cal$STATION_ID]
cal$sin_dato <- is.na(cal$VALUE)
registrar("3. Calendario completo (días faltantes como NA)", cal$STATION_ID, cal$anio, cal$VALUE)

sin_dato_mes <- tapply(cal$sin_dato, list(cal$anio, cal$mes), mean) * 100
para_grafica <- !cal$sin_dato & !es_imposible(cal$VALUE)
lluvia_mes <- tapply(cal$VALUE[para_grafica], list(cal$anio[para_grafica], cal$mes[para_grafica]), mean)
tiempo <- as.vector(t(outer(as.integer(rownames(sin_dato_mes)), (1:12 - .5) / 12, "+")))

graficar(function() {   # Gráfica: % de días sin dato, mes a mes y según el mes del año
  op <- par(mfrow = c(2, 1), mar = c(4, 5, 3, 1)); on.exit(par(op))
  plot(tiempo, as.vector(t(sin_dato_mes)), type = "l", col = COL_PRINCIPAL, lwd = 1.5,
       xlab = "", ylab = "% de días sin dato", main = "Días sin dato en cada mes, toda la red")
  barplot(colMeans(sin_dato_mes, na.rm = TRUE), names.arg = month.abb, col = COL_PRINCIPAL,
          border = NA, ylab = "% de días sin dato", main = "Días sin dato según el mes del año")
})
graficar(function() {   # Gráfica: lluvia media de la red por mes y por año
  op <- par(mfrow = c(2, 1), mar = c(4, 5, 3, 1)); on.exit(par(op))
  plot(tiempo, as.vector(t(lluvia_mes)), type = "l", col = COL_PRINCIPAL, lwd = 1.2,
       xlab = "", ylab = "mm/día", main = "Lluvia media de la red en cada mes")
  anual <- tapply(cal$VALUE[para_grafica], cal$anio[para_grafica], mean)
  barplot(anual, col = COL_PRINCIPAL, border = NA, las = 2, cex.names = .7, ylab = "mm/día",
          main = "Lluvia media de la red en cada año")
})

# ---- 3.4 Matriz estación x día -----------------------------------------------
# Para las rachas y la cobertura conviene tener los datos como una matriz:
# filas = estaciones, columnas = días de 1970 a 2003. Así se pueden mirar
# muchas estaciones en los mismos días de un solo golpe.
ids <- sort(unique(cal$STATION_ID))
primer_dia <- min(as.integer(cal$fecha))
n_dias <- max(as.integer(cal$fecha)) - primer_dia + 1
fechas <- as.Date(primer_dia + seq_len(n_dias) - 1)     # la fecha de cada columna
doy    <- as.integer(format(fechas, "%j"))              # día del año (1 a 366)
anio_col <- as.integer(format(fechas, "%Y"))
mes_col  <- as.integer(format(fechas, "%m"))
de_referencia <- anio_col <= ANIO_REFERENCIA

fila_est <- match(cal$STATION_ID, ids)                  # fila de la matriz de cada día del calendario
col_dia  <- as.integer(cal$fecha) - primer_dia + 1      # columna de la matriz
primera_fila <- match(ids, cal$STATION_ID)              # primera fila de cal de cada estación
primera_col  <- col_dia[primera_fila]
ultima_col   <- as.vector(tapply(col_dia, fila_est, max))
fila_de <- function(k, cols) primera_fila[k] + cols - primera_col[k]   # (estación, columna) -> fila de cal

V <- matrix(NA_real_, length(ids), n_dias, dimnames = list(ids, NULL))
V[cbind(fila_est, col_dia)] <- ifelse(es_imposible(cal$VALUE), NA, cal$VALUE)

# ---- 3.5 Rachas largas de ceros (primera versión: comparando con la zona) ----
# Un cero puede ser real (no llovió) o falso (el equipo falló o nadie anotó).
# Para cada racha de RACHA_CEROS días o más en cero se pregunta: en esos mismos
# días, ¿llovió en las demás estaciones de su zona?
#   razón = % de días con lluvia en la zona durante la racha
#           ÷ % normal de días con lluvia en la zona en esas fechas (±7 días, 1970-1997)
#   razón < 0,5      -> Probable sequía (la red también estuvo seca)
#   razón 0,5 a 0,8  -> Revisar
#   razón >= 0,8     -> Sospechoso (la red llovió normal)
# En la Parte 4 se repite con las estaciones vecinas, que es más preciso.
cat("\n== 3.5 RACHAS LARGAS DE CEROS (por zona) ==\n")
zona_fila <- zona_est[ids]
red <- lapply(split(seq_along(ids), zona_fila), function(filas) {
  list(reportan = colSums(!is.na(V[filas, , drop = FALSE])),
       llueve   = colSums(V[filas, , drop = FALSE] >= LLUVIA_RED, na.rm = TRUE))
})
# Fracción normal de estaciones con lluvia en cada día del año, suavizada con ±7 días
normal <- lapply(red, function(z) {
  f <- tapply(z$llueve[de_referencia], doy[de_referencia], sum) /
    tapply(z$reportan[de_referencia], doy[de_referencia], sum)
  f <- f[as.character(1:366)]
  sapply(1:366, function(k) mean(f[((k - 8):(k + 6)) %% 366 + 1], na.rm = TRUE))
})

# Buscar las rachas: rle() parte la serie en tramos de valores iguales
rachas <- do.call(rbind, lapply(seq_along(ids), function(k) {
  cols <- primera_col[k]:ultima_col[k]
  v <- V[k, cols]
  r <- rle(!is.na(v) & v == 0)
  fin <- cumsum(r$lengths)
  ini <- fin - r$lengths + 1
  j <- which(r$values & r$lengths >= RACHA_CEROS)
  if (!length(j)) return(NULL)
  data.frame(STATION_ID = ids[k], k = k, zona = zona_fila[[k]],
             c_ini = cols[ini[j]], c_fin = cols[fin[j]], dias = r$lengths[j],
             stringsAsFactors = FALSE)
}))
rachas$inicio <- fechas[rachas$c_ini]
rachas$fin    <- fechas[rachas$c_fin]
rachas$mes_central <- mes_col[round((rachas$c_ini + rachas$c_fin) / 2)]
rachas$temporada <- ifelse(rachas$mes_central %in% c(12, 1, 2, 6, 7, 8), "Seca", "Lluviosa")
rachas$lluvia_zona <- NA_real_
rachas$lluvia_zona_normal <- NA_real_
for (j in seq_len(nrow(rachas))) {
  d <- rachas$c_ini[j]:rachas$c_fin[j]
  z <- red[[rachas$zona[j]]]
  otras <- z$reportan[d] - 1                    # sin contar la propia estación
  rachas$lluvia_zona[j] <- if (sum(otras) > 0) sum(z$llueve[d]) / sum(otras) else NA
  rachas$lluvia_zona_normal[j] <- mean(normal[[rachas$zona[j]]][doy[d]])
}
rachas$razon_zona <- rachas$lluvia_zona / rachas$lluvia_zona_normal
rachas$clasificacion_zona <- ifelse(is.na(rachas$razon_zona), NA,
                                    ifelse(rachas$razon_zona < .5, CLASES[1], ifelse(rachas$razon_zona >= .8, CLASES[3], CLASES[2])))
cat(sprintf("Rachas de %d o más días en cero: %d\n", RACHA_CEROS, nrow(rachas)))
print(table(rachas$clasificacion_zona, rachas$temporada))

graficar(function() {   # Gráfica: rachas según su mes central y su clase (por zona)
  tabla <- table(factor(rachas$clasificacion_zona, levels = CLASES), factor(rachas$mes_central, levels = 1:12))
  colores <- c(COL_CERO, COL_GRIS, COL_ROJO)
  op <- par(mar = c(4, 5, 3, 1)); on.exit(par(op))
  barplot(tabla, names.arg = month.abb, col = colores, border = NA, ylab = "rachas",
          main = sprintf("Rachas de %d o más ceros según su mes central (comparando con la zona)", RACHA_CEROS))
  legend("topright", fill = colores, legend = CLASES_CORTAS, bty = "n")
})

# ---- 3.6 Cobertura por estación y año (1970-2003) ----------------------------
# Para cada estación y año: % de días con dato. Los años en que la estación
# no operaba quedan como NA (no cuentan como problema).
cat("\n== 3.6 COBERTURA POR ESTACIÓN Y AÑO ==\n")
anios <- ANIO_INICIO:ANIO_FIN
m_dias <- unclass(table(cal$STATION_ID[!cal$sin_dato], factor(cal$anio[!cal$sin_dato], levels = anios)))
m_cob  <- sweep(m_dias, 2, dias_del_anio(anios), "/") * 100
primero <- tapply(cal$anio, cal$STATION_ID, min)[rownames(m_cob)]
ultimo  <- tapply(cal$anio, cal$STATION_ID, max)[rownames(m_cob)]
for (k in seq_len(nrow(m_cob))) m_cob[k, anios < primero[k] | anios > ultimo[k]] <- NA

en_operacion <- tapply(!cal$sin_dato, cal$STATION_ID, mean)[rownames(m_cob)] * 100
cobertura <- data.frame(
  STATION_ID             = rownames(m_cob),
  zona                   = zona_est[rownames(m_cob)],
  primer_anio            = primero, ultimo_anio = ultimo,
  anios_operando         = ultimo - primero + 1,
  cobertura_en_operacion = round(en_operacion, 1),
  cobertura_1970_2003    = round(rowSums(m_dias) / sum(dias_del_anio(anios)) * 100, 1),
  anios_bajo_80          = rowSums(m_cob < 80, na.rm = TRUE),
  anios_con_problema     = apply(m_cob, 1, function(x) rangos(anios[!is.na(x) & x < 80])),
  stringsAsFactors       = FALSE)
# ¿Los años malos están juntos (un periodo) o repartidos?
bloques <- sapply(strsplit(cobertura$anios_con_problema, ", "), function(x) sum(nzchar(x)))
cobertura$patron <- ifelse(cobertura$anios_bajo_80 == 0, "Sin años problemáticos",
                           ifelse(bloques == 1, "Concentrado en un periodo", "Repartido en varios periodos"))
cobertura <- cobertura[order(cobertura$cobertura_en_operacion), ]
rownames(cobertura) <- NULL
guardar_tabla(cbind(cobertura, round(m_cob[cobertura$STATION_ID, ], 1)), "prep", "cobertura_por_estacion_y_anio")
print(table(cobertura$patron))
cat("\nLas 20 estaciones con peor cobertura mientras operaban:\n")
print(head(cobertura[, c("STATION_ID", "primer_anio", "ultimo_anio", "cobertura_en_operacion",
                         "anios_bajo_80", "anios_con_problema")], 20), row.names = FALSE)

# Vista 1, año por año: cuántas estaciones caen en cada nivel de cobertura
por_anio_nivel <- sapply(colnames(m_cob), function(a) table(cut(m_cob[, a], CORTES_COB, labels = NIVELES_COB)))
guardar_tabla(data.frame(anio = anios, t(por_anio_nivel), check.names = FALSE), "prep", "estaciones_por_nivel_y_anio")
graficar(function() {   # Gráfica: una grafiquita por año
  op <- par(mfrow = c(5, 7), mar = c(1, 2.5, 1.8, .3), oma = c(0, 0, 3, 0)); on.exit(par(op))
  tope <- max(por_anio_nivel)
  for (a in seq_along(anios)) {
    bp <- barplot(por_anio_nivel[, a], col = COLORES_COB, border = NA, ylim = c(0, tope * 1.2),
                  names.arg = FALSE, las = 1, cex.axis = .7, main = anios[a], cex.main = .9)
    text(bp, por_anio_nivel[, a], por_anio_nivel[, a], pos = 3, cex = .6)
  }
  plot.new()
  legend("center", fill = COLORES_COB, legend = NIVELES_COB, bty = "n", title = "Cobertura")
  mtext("Número de estaciones según su cobertura, año por año", outer = TRUE, line = 1, font = 2)
})

# Vista 2, por tandas: de la peor a la mejor estación, 40 por gráfica, con su código
TAMANO_TANDA <- 40
tandas <- split(cobertura$STATION_ID, ceiling(seq_len(nrow(cobertura)) / TAMANO_TANDA))
for (k in seq_along(tandas)) graficar(function() {   # Gráfica: tanda k
  ids_tanda <- rev(tandas[[k]])                       # la peor de la tanda queda arriba
  z <- m_cob[ids_tanda, , drop = FALSE]
  fila <- seq_along(ids_tanda)
  etiquetas <- sprintf("%s (%.0f %%)", ids_tanda,
                       cobertura$cobertura_en_operacion[match(ids_tanda, cobertura$STATION_ID)])
  op <- par(mar = c(4, 9.5, 5, 1)); on.exit(par(op))
  image(anios, fila, t(ifelse(is.na(z), 1, NA)), col = "#eeeeee", zlim = c(0, 1),
        xlab = "Año", ylab = "", yaxt = "n")
  image(anios, fila, t(z), col = COLORES_COB, breaks = CORTES_COB, add = TRUE)
  abline(v = anios + .5, h = fila + .5, col = "white", lwd = .6)
  axis(2, at = fila, labels = etiquetas, las = 1, cex.axis = .6)
  inicio <- (k - 1) * TAMANO_TANDA + 1
  title(main = sprintf("Tanda %d de %d: estaciones %d a %d de %d (de peor a mejor cobertura)",
                       k, length(tandas), inicio, inicio + length(ids_tanda) - 1, nrow(cobertura)), line = 3)
  legend("top", inset = c(0, -.07), xpd = NA, fill = c(COLORES_COB, "#eeeeee"), bty = "n", cex = .8,
         legend = c(NIVELES_COB, "No operaba"), horiz = TRUE)
})

# ---- 3.7 Satélite: 1981-2003 y unidades --------------------------------------
# Se deja solo el periodo del estudio. Luego se comprueba que los valores son
# mm ACUMULADOS en la péntada: la última péntada del mes (3 a 6 días) se compara
# con las péntadas 4 y 5 del mismo mes (5 días). Si es acumulado, la razón cruda
# sigue a los días (6/5 = 1,2; 3/5 = 0,6) y al pasar a mm/día queda cerca de 1.
cat("\n== 3.7 SATÉLITE ==\n")
en_estudio <- anio_sat <= ANIO_FIN
M_est <- M[, en_estudio]
anio_e <- anio_sat[en_estudio]; pent_e <- pent_sat[en_estudio]
mes_e  <- mes_sat[en_estudio];  dias_e <- dias_sat[en_estudio]
posicion <- pent_e - (mes_e - 1) * 6                  # 1 a 6 dentro del mes

media_pent <- colMeans(M_est, na.rm = TRUE)
clave_mes <- paste(anio_e, mes_e)
referencia_mes <- tapply(media_pent[posicion %in% 4:5], clave_mes[posicion %in% 4:5], mean)
ultima <- posicion == 6
razon_cruda <- media_pent[ultima] / referencia_mes[clave_mes[ultima]]
unidades_sat <- data.frame(
  dias_ultima_pentada      = sort(unique(dias_e[ultima])),
  esperado_si_es_acumulado = round(sort(unique(dias_e[ultima])) / 5, 2),
  razon_cruda              = round(tapply(razon_cruda, dias_e[ultima], median), 2),
  razon_en_mm_dia          = round(tapply(razon_cruda * 5 / dias_e[ultima], dias_e[ultima], median), 2))
cat("Última péntada del mes contra las péntadas 4 y 5 del mismo mes (medianas):\n")
print(unidades_sat, row.names = FALSE)
guardar_tabla(unidades_sat, "prep", "unidades_satelite")
cat(sprintf("Satélite en estudio: %d imágenes, %d-%d | celdas en el tope de 200 mm/péntada: %d\n",
            ncol(M_est), min(anio_e), max(anio_e), sum(M_est == 200, na.rm = TRUE)))

# Formato largo: una fila por píxel y péntada
satelite <- data.frame(
  pixel      = rep(sprintf("P%04d", seq_len(nrow(M_est))), times = ncol(M_est)),
  lon        = rep(xy_pixel[, 1], times = ncol(M_est)),
  lat        = rep(xy_pixel[, 2], times = ncol(M_est)),
  anio       = rep(anio_e, each = nrow(M_est)),
  pentada    = rep(pent_e, each = nrow(M_est)),
  dias       = rep(dias_e, each = nrow(M_est)),
  mm_pentada = as.vector(M_est),
  mm_dia     = round(as.vector(M_est) / rep(dias_e, each = nrow(M_est)), 4))


# =============================================================================
# PARTE 4. LIMPIEZA
# =============================================================================
cat("\n######## PARTE 4. LIMPIEZA ########\n")

# ---- 4.1 Ubicación de las estaciones -----------------------------------------
# Coordenadas: archivo de la CVC. Altitud: catálogo de estaciones.
coord  <- read.csv(file.path(RUTA_DATOS, "CVC_Wheather_Stations_Long_Lat.csv"),
                   colClasses = c(CODIGO = "character"))
altura <- as.data.frame(readxl::read_excel(file.path(RUTA_DATOS, "CatalogoEstacionesElevacion.xlsx"),
                                           sheet = "Hoja1", col_types = c("text", "text", "text", "text", "numeric")))
i_coord <- match(ids, coord$CODIGO)
i_alt   <- match(ids, altura$codigo)
info <- data.frame(
  STATION_ID = ids,
  zona       = zona_est[ids],
  nombre     = ifelse(is.na(i_coord), altura$nombre[i_alt], coord$STATION_NAME[i_coord]),
  municipio  = ifelse(is.na(i_coord), altura$municipio[i_alt], coord$Municipality[i_coord]),
  lon        = coord$Lon[i_coord],
  lat        = coord$lat[i_coord],
  altitud    = altura$alt[i_alt],
  stringsAsFactors = FALSE)
info$con_coordenadas <- !is.na(info$lon)
cat("\n== 4.1 UBICACIÓN ==\n")
cat(sprintf("Estaciones: %d | con coordenadas: %d | con altitud: %d | con ambas: %d\n",
            nrow(info), sum(info$con_coordenadas), sum(!is.na(info$altitud)),
            sum(info$con_coordenadas & !is.na(info$altitud))))
print(table(zona = info$zona, con_coordenadas = info$con_coordenadas))

# Lista de estaciones sin coordenadas, para pedírselas al profesor
sin_coord <- merge(info[!info$con_coordenadas, c("STATION_ID", "zona", "nombre", "municipio", "altitud")],
                   cobertura[, c("STATION_ID", "primer_anio", "ultimo_anio", "cobertura_en_operacion")],
                   by = "STATION_ID")
sin_coord$cuenca <- altura$cuenca[match(sin_coord$STATION_ID, altura$codigo)]
sin_coord <- sin_coord[order(-sin_coord$cobertura_en_operacion), ]
guardar_tabla(sin_coord, "limpieza", "estaciones_sin_coordenadas")

# ---- 4.2 Lecturas mensuales en vez de diarias --------------------------------
# Algunas estaciones (o algunos años de una estación) no anotan cada día:
# escriben el total del mes en uno o pocos días. Se revisa AÑO POR AÑO:
#   un año es "mensual" si tiene al menos 3 meses con dato y, en la mitad de
#   esos meses, MAX_DIAS_MES días con dato o menos.
# Decisión:
#   - si la mitad o más de los años de la estación son mensuales, se elimina
#     la estación completa (no sirve para curvas diarias ni como vecina);
#   - si solo son algunos años, se quitan solo esos años.
cat("\n== 4.2 LECTURAS MENSUALES ==\n")
clave_sy <- paste(cal$STATION_ID, cal$anio)                       # estación-año de cada día
dias_mes <- tapply(!cal$sin_dato, list(clave_sy, cal$mes), sum)   # filas: estación-año, columnas: mes
meses_con_dato <- rowSums(dias_mes > 0, na.rm = TRUE)
mediana_mes <- apply(dias_mes, 1, function(x) { x <- x[!is.na(x) & x > 0]; if (length(x)) median(x) else NA })
anios_revisados <- data.frame(
  STATION_ID           = substr(rownames(dias_mes), 1, 10),
  anio                 = as.integer(substr(rownames(dias_mes), 12, 15)),
  meses_con_dato       = meses_con_dato,
  mediana_dias_por_mes = mediana_mes,
  stringsAsFactors     = FALSE)
anios_revisados$mensual <- !is.na(mediana_mes) & meses_con_dato >= 3 & mediana_mes <= MAX_DIAS_MES

lecturas_mensuales <- do.call(rbind, lapply(split(anios_revisados, anios_revisados$STATION_ID), function(d) {
  data.frame(STATION_ID = d$STATION_ID[1], anios = nrow(d), anios_mensuales = sum(d$mensual),
             cuales = rangos(d$anio[d$mensual]), stringsAsFactors = FALSE)
}))
lecturas_mensuales <- lecturas_mensuales[lecturas_mensuales$anios_mensuales > 0, ]
lecturas_mensuales$decision <- ifelse(lecturas_mensuales$anios_mensuales / lecturas_mensuales$anios >= .5,
                                      "Eliminar la estación", "Quitar solo esos años")
lecturas_mensuales <- lecturas_mensuales[order(-lecturas_mensuales$anios_mensuales), ]
rownames(lecturas_mensuales) <- NULL
guardar_tabla(lecturas_mensuales, "limpieza", "lecturas_mensuales")
cat(sprintf("Estaciones con algún año mensual: %d (años mensuales en total: %d de %d)\n",
            nrow(lecturas_mensuales), sum(lecturas_mensuales$anios_mensuales), nrow(anios_revisados)))
print(lecturas_mensuales, row.names = FALSE)

info$lecturas_mensuales <- ids %in% lecturas_mensuales$STATION_ID[lecturas_mensuales$decision == "Eliminar la estación"]
info$anios_mensuales <- as.vector(tapply(anios_revisados$mensual, anios_revisados$STATION_ID, sum)[ids])
sy_mensual <- with(anios_revisados[anios_revisados$mensual, ], paste(STATION_ID, anio))
dia_mensual <- clave_sy %in% sy_mensual | cal$STATION_ID %in% ids[info$lecturas_mensuales]
# Esos días tampoco se usan como "vecina" ni para buscar acumulados
V[cbind(fila_est[dia_mensual], col_dia[dia_mensual])] <- NA
rm(dias_mes, clave_sy)

# ---- 4.3 Vecinas de cada estación: por distancia o por correlación -----------
# Para revisar los ceros falsos (4.4) y para imputar (Parte 5), cada estación
# necesita sus "vecinas": estaciones que viven la misma lluvia.
#   - Estaciones CON coordenadas: vecinas por DISTANCIA. Las que están a
#     RADIO_KM o menos (hasta RADIO_MAX_KM si hay menos de MIN_VECINAS), con
#     menos de DIF_ALTURA m de diferencia de altura.
#   - Estaciones SIN coordenadas: NO se les inventa una ubicación. Sus vecinas
#     se escogen por CORRELACIÓN de la lluvia diaria (triangulación): dos
#     estaciones cercanas ven los mismos aguaceros el mismo día, así que la
#     correlación funciona como una medida de cercanía. Son vecinas las
#     N_VECINAS_CORR estaciones más correlacionadas (con al menos R_MIN_VECINA).
# El número N no se escoge a ojo: a las estaciones que SÍ tienen coordenadas se
# les buscan vecinas por correlación como si no las tuvieran, y se mide a qué
# distancia real quedaron, para N = 3, 5, 8, 10 y 15. Con pocas vecinas quedan
# más cerca pero el promedio es más ruidoso; con muchas, se cuelan estaciones
# lejanas. N = 5 equilibra las dos cosas.
# (Se probó también un umbral fijo de correlación, r >= 0,38, que es la típica
# a 20 km: daba hasta 177 vecinas por estación y solo el 62 % a 30 km o menos,
# porque estaciones lejanas también pueden superar ese umbral.)
# Las estaciones sin coordenadas no se emparejan con el satélite (ver 4.7):
# para eso sí hace falta saber en qué píxel están.
cat("\n== 4.3 VECINAS DE CADA ESTACIÓN ==\n")
activas  <- which(!info$lecturas_mensuales)
ubicadas <- which(info$con_coordenadas & !info$lecturas_mensuales)

# a) Correlación de la lluvia diaria entre todos los pares de estaciones.
#    Raíz cuadrada para que los aguaceros no dominen; solo días con dato en ambas.
RAIZ  <- sqrt(V)
TIENE <- !is.na(RAIZ)
CORR <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
for (e1 in activas) for (e2 in activas) {
  if (e2 <= e1) next
  comun <- TIENE[e1, ] & TIENE[e2, ]
  if (sum(comun) < MIN_DIAS_COMUNES) next
  CORR[e1, e2] <- CORR[e2, e1] <- cor(RAIZ[e1, comun], RAIZ[e2, comun])
}
rm(RAIZ, TIENE)

# b) Calibración: correlación típica según la distancia (estaciones con coordenadas)
pares_u <- t(combn(ubicadas, 2))
calib <- data.frame(
  d = dist_km(info$lon[pares_u[, 1]], info$lat[pares_u[, 1]], info$lon[pares_u[, 2]], info$lat[pares_u[, 2]]),
  r = CORR[pares_u])
calib <- calib[!is.na(calib$r), ]
BANDAS_KM <- c(0, 5, 10, 20, 30, 50, 100, Inf)
banda <- cut(calib$d, BANDAS_KM)
calibracion <- data.frame(distancia_km = levels(banda), pares = as.vector(table(banda)),
                          correlacion_mediana = round(as.vector(tapply(calib$r, banda, median)), 3))
guardar_tabla(calibracion, "limpieza", "calibracion_correlacion_distancia")
cat("Correlación mediana de la lluvia diaria según la distancia entre estaciones:\n")
print(calibracion, row.names = FALSE)

graficar(function() {   # Gráfica: correlación contra distancia
  op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  boxplot(r ~ banda, data = data.frame(r = calib$r, banda = banda), col = COL_GRIS, border = "#56636f",
          outline = FALSE, xlab = "Distancia entre estaciones (km)", ylab = "Correlación de la lluvia diaria",
          main = "Cuanto más cerca, más se parece la lluvia (estaciones con coordenadas)")
})

# c) Las dos reglas para escoger vecinas
vecinas_por_distancia <- function(k) {
  d <- dist_km(info$lon[k], info$lat[k], info$lon, info$lat)
  d[k] <- NA
  d[info$lecturas_mensuales] <- NA
  altura_ok <- is.na(info$altitud) | is.na(info$altitud[k]) |
    abs(info$altitud - info$altitud[k]) <= DIF_ALTURA
  v <- which(!is.na(d) & altura_ok & d <= RADIO_KM)
  if (length(v) < MIN_VECINAS) v <- which(!is.na(d) & altura_ok & d <= RADIO_MAX_KM)
  v
}
vecinas_por_correlacion <- function(k, candidatas, n = N_VECINAS_CORR) {
  r <- CORR[k, candidatas]
  ok <- which(!is.na(r) & r >= R_MIN_VECINA)
  candidatas[ok[order(-r[ok])][seq_len(min(n, length(ok)))]]
}

# d) Validación: a las estaciones con coordenadas se les buscan vecinas por
#    correlación (como si no tuvieran coordenadas) y se mide la distancia real,
#    para varios valores de N
validacion_vecinas <- do.call(rbind, lapply(c(3, 5, 8, 10, 15), function(n) {
  medida <- sapply(ubicadas, function(k) {
    v <- vecinas_por_correlacion(k, setdiff(ubicadas, k), n)
    if (!length(v)) return(c(NA, NA, NA))
    d <- dist_km(info$lon[k], info$lat[k], info$lon[v], info$lat[v])
    c(median(d), mean(d <= RADIO_MAX_KM), mean(v %in% vecinas_por_distancia(k)))
  })
  data.frame(N = n, distancia_mediana_km = round(median(medida[1, ], na.rm = TRUE), 1),
             pct_a_30_km_o_menos = round(mean(medida[2, ], na.rm = TRUE) * 100, 1),
             pct_tambien_vecina_por_distancia = round(mean(medida[3, ], na.rm = TRUE) * 100, 1))
}))
guardar_tabla(validacion_vecinas, "limpieza", "validacion_vecinas_por_correlacion")
cat("Validación (estaciones con coordenadas, vecinas buscadas solo por correlación):\n")
print(validacion_vecinas, row.names = FALSE)
cat(sprintf("Se usa N = %d\n", N_VECINAS_CORR))

# e) Vecinas definitivas de cada estación
vecinas <- lapply(seq_along(ids), function(k) {
  if (info$lecturas_mensuales[k]) return(integer(0))
  if (info$con_coordenadas[k]) vecinas_por_distancia(k) else vecinas_por_correlacion(k, setdiff(activas, k))
})
info$tipo_vecinas <- ifelse(info$lecturas_mensuales, "Ninguna (lecturas mensuales)",
                            ifelse(info$con_coordenadas, "Por distancia", "Por correlación"))
info$r_max_vecina <- sapply(seq_along(ids), function(k) {
  r <- CORR[k, vecinas[[k]]]
  if (length(r) && any(!is.na(r))) round(max(r, na.rm = TRUE), 3) else NA
})

# f) Localidad de las estaciones sin coordenadas: el código indica la cuenca
#    (primeros 7 dígitos). El nombre de la cuenca sale del catálogo; si la
#    estación no está en él, se toma el de otra estación con el mismo código.
info$codigo_cuenca <- substr(ids, 1, 7)
info$cuenca <- altura$cuenca[i_alt]
nombre_cuenca <- tapply(info$cuenca, info$codigo_cuenca, function(x) {
  x <- x[!is.na(x)]
  if (length(x)) names(sort(table(x), decreasing = TRUE))[1] else NA
})
info$cuenca[is.na(info$cuenca)] <- nombre_cuenca[info$codigo_cuenca[is.na(info$cuenca)]]

vecinas_sin_coordenadas <- data.frame(
  info[!info$con_coordenadas & !info$lecturas_mensuales, c("STATION_ID", "zona", "nombre", "municipio",
                                                           "codigo_cuenca", "cuenca", "r_max_vecina")],
  n_vecinas = lengths(vecinas)[!info$con_coordenadas & !info$lecturas_mensuales],
  vecinas = sapply(vecinas, function(v) paste(ids[v], collapse = " "))[!info$con_coordenadas & !info$lecturas_mensuales],
  stringsAsFactors = FALSE)
guardar_tabla(vecinas_sin_coordenadas, "limpieza", "vecinas_estaciones_sin_coordenadas")
cat(sprintf("Estaciones sin coordenadas: %d | con vecinas por correlación: %d | sin ninguna vecina: %d\n",
            nrow(vecinas_sin_coordenadas), sum(vecinas_sin_coordenadas$n_vecinas > 0),
            sum(vecinas_sin_coordenadas$n_vecinas == 0)))
cat(sprintf("Estaciones sin coordenadas con nombre de cuenca: %d\n", sum(!is.na(vecinas_sin_coordenadas$cuenca))))

# ---- 4.4 Rachas de ceros con estaciones vecinas ------------------------------
# Misma pregunta que en 3.5 (¿llovió alrededor mientras esta estación marcaba
# cero?), pero comparando con las vecinas de 4.3 en lugar de toda la zona.
cat("\n== 4.4 RACHAS DE CEROS CON VECINAS ==\n")
info$n_vecinas <- lengths(vecinas)
info$vecinas <- sapply(vecinas, function(v) paste(ids[v], collapse = " "))
cat("Vecinas por estación ubicada:\n")
print(tapply(info$n_vecinas, info$tipo_vecinas, summary))

REPORTA <- !is.na(V)                      # TRUE si la estación tiene dato ese día
LLUEVE  <- !is.na(V) & V >= LLUVIA_RED    # TRUE si ese día llovió 1 mm o más

# Para una racha (estación k, columnas c1 a c2):
#   observado = % de registros de las vecinas con lluvia en esos días
#   normal    = el mismo % en las mismas fechas (±7 días) de 1970-1997
comparar_con_vecinas <- function(k, c1, c2) {
  N <- vecinas[[k]]
  if (!length(N)) return(c(NA, NA, 0))
  d <- c1:c2
  reportan <- sum(REPORTA[N, d])
  observado <- if (reportan > 0) sum(LLUEVE[N, d]) / reportan else NA
  ventana <- unique((rep(doy[d], each = 15) + rep(-7:7, length(d)) - 1) %% 366 + 1)
  ref <- which(de_referencia & doy %in% ventana)
  normal_vec <- sum(LLUEVE[N, ref]) / sum(REPORTA[N, ref])
  c(observado, normal_vec, length(N))
}
resultado <- t(mapply(comparar_con_vecinas, rachas$k, rachas$c_ini, rachas$c_fin))
rachas$n_vecinas             <- resultado[, 3]
rachas$lluvia_vecinas        <- resultado[, 1]
rachas$lluvia_vecinas_normal <- resultado[, 2]
rachas$razon_vecinas         <- rachas$lluvia_vecinas / rachas$lluvia_vecinas_normal
rachas$clasificacion_vecinas <- ifelse(is.na(rachas$razon_vecinas), NA,
                                       ifelse(rachas$razon_vecinas < .5, CLASES[1], ifelse(rachas$razon_vecinas >= .8, CLASES[3], CLASES[2])))

# Clasificación final: con vecinas si se pudo; si no, la de la zona
rachas$metodo <- ifelse(!is.na(rachas$clasificacion_vecinas),
                        ifelse(info$con_coordenadas[rachas$k], "Vecinas por distancia", "Vecinas por correlación"),
                        "Zona (sin vecinas con dato)")
rachas$clasificacion <- ifelse(is.na(rachas$clasificacion_vecinas), rachas$clasificacion_zona,
                               rachas$clasificacion_vecinas)

cat("Método usado:\n"); print(table(rachas$metodo))
con_vecinas <- grepl("^Vecinas", rachas$metodo)
cambio <- table(zona    = factor(rachas$clasificacion_zona[con_vecinas], CLASES, CLASES_CORTAS),
                vecinas = factor(rachas$clasificacion_vecinas[con_vecinas], CLASES, CLASES_CORTAS))
cat("\nClasificación por zona (filas) contra vecinas (columnas):\n"); print(cambio)
cat(sprintf("Rachas que cambian de clase: %d de %d\n", sum(cambio) - sum(diag(cambio)), sum(cambio)))
cat("\nClasificación final:\n"); print(table(rachas$clasificacion, rachas$temporada))

guardar_tabla(transform(rachas[, c("STATION_ID", "zona", "inicio", "fin", "dias", "mes_central", "temporada",
                                   "lluvia_zona", "lluvia_zona_normal", "razon_zona", "clasificacion_zona",
                                   "n_vecinas", "lluvia_vecinas", "lluvia_vecinas_normal", "razon_vecinas",
                                   "clasificacion_vecinas", "metodo", "clasificacion")],
                        lluvia_zona = round(lluvia_zona * 100, 1), lluvia_zona_normal = round(lluvia_zona_normal * 100, 1),
                        razon_zona = round(razon_zona, 2),
                        lluvia_vecinas = round(lluvia_vecinas * 100, 1),
                        lluvia_vecinas_normal = round(lluvia_vecinas_normal * 100, 1),
                        razon_vecinas = round(razon_vecinas, 2)), "limpieza", "rachas_de_ceros")
guardar_tabla(as.data.frame.matrix(cambio), "limpieza", "rachas_zona_contra_vecinas")

# Se marca cada día del calendario con la clase de su racha
cal$racha_ceros <- NA_character_
for (j in seq_len(nrow(rachas))) {
  cal$racha_ceros[fila_de(rachas$k[j], rachas$c_ini[j]:rachas$c_fin[j])] <- rachas$clasificacion[j]
}
marcados <- !is.na(cal$racha_ceros)
meses_marcados <- aggregate(
  list(dias_en_racha = rep(1, sum(marcados))),
  by = list(STATION_ID = cal$STATION_ID[marcados], anio = cal$anio[marcados],
            mes = cal$mes[marcados], clasificacion = cal$racha_ceros[marcados]),
  FUN = sum)
meses_marcados <- meses_marcados[order(meses_marcados$STATION_ID, meses_marcados$anio, meses_marcados$mes), ]
guardar_tabla(meses_marcados, "limpieza", "meses_marcados_por_ceros")

graficar(function() {   # Gráfica: cuántas rachas hay en cada clase, por zona y por vecinas
  op <- par(mar = c(4, 5, 4, 1)); on.exit(par(op))
  m <- rbind("Comparando con la zona" = rowSums(cambio), "Comparando con vecinas" = colSums(cambio))
  colnames(m) <- CLASES_CORTAS
  bp <- barplot(m, beside = TRUE, col = c(COL_GRIS, COL_PRINCIPAL), border = NA, ylim = c(0, max(m) * 1.15),
                ylab = "rachas", main = "Rachas con estaciones ubicadas: zona contra vecinas cercanas")
  text(bp, m, m, pos = 3, cex = .85)
  legend("topright", fill = c(COL_GRIS, COL_PRINCIPAL), legend = rownames(m), bty = "n")
})
graficar(function() {   # Gráfica: ejemplo de la racha sospechosa más larga contra sus vecinas
  s <- rachas[rachas$metodo == "Vecinas por distancia" & rachas$clasificacion == CLASES[3], ]
  s <- s[which.max(s$dias), ]
  d <- max(1, s$c_ini - 30):min(n_dias, s$c_fin + 30)
  N <- vecinas[[s$k]]
  media_vecinas <- colMeans(V[N, d, drop = FALSE], na.rm = TRUE)
  op <- par(mar = c(4, 5, 4, 1)); on.exit(par(op))
  plot(fechas[d], media_vecinas, type = "h", col = COL_PRINCIPAL, lwd = 2,
       ylim = c(0, max(media_vecinas, V[s$k, d], na.rm = TRUE)), xlab = "", ylab = "mm/día",
       main = sprintf("Racha sospechosa: estación %s, %d días en cero", s$STATION_ID, s$dias))
  rect(fechas[s$c_ini], -1, fechas[s$c_fin], 1e4, col = adjustcolor(COL_ROJO, .12), border = NA)
  lines(fechas[d], V[s$k, d], type = "h", col = COL_ROJO, lwd = 2)
  legend("topleft", bty = "n", col = c(COL_PRINCIPAL, COL_ROJO), lwd = 2,
         legend = c(sprintf("Promedio de sus %d vecinas", length(N)), "La estación"))
})

# ---- 4.5 Pasar a NA lo que no es lluvia real ---------------------------------
# VALUE se queda igual. VALUE_LIMPIO es la copia limpia y "motivo" dice por
# qué un día quedó sin dato. El orden importa: un día se marca una sola vez.
cat("\n== 4.5 VALORES QUE PASAN A NA ==\n")
cal$VALUE_LIMPIO <- cal$VALUE
cal$motivo <- ifelse(cal$sin_dato, "Sin dato original", NA_character_)
marcar <- function(filas, texto) {
  filas <- filas[!is.na(cal$VALUE_LIMPIO[filas])]
  cal$VALUE_LIMPIO[filas] <<- NA
  cal$motivo[filas] <<- texto
  cat(sprintf("%-40s %8d días\n", texto, length(filas)))
}

# Cada tipo de dato malo se quita en su propio paso y queda en la trazabilidad.

# a) Lecturas mensuales (ver 4.2): estaciones eliminadas completas y años sueltos.
#    Se eliminan: son pocas (ver tabla de 4.2) y repartir un total mensual en
#    días inventaría el momento de cada aguacero.
marcar(which(cal$STATION_ID %in% ids[info$lecturas_mensuales]), "Estación con lecturas mensuales")
marcar(which(dia_mensual), "Año con lecturas mensuales")
registrar("4a. Quitar lecturas mensuales", cal$STATION_ID, cal$anio, cal$VALUE_LIMPIO)

# b) Valores imposibles: negativos o por encima de LIMITE_FISICO (récord mundial)
marcar(which(es_imposible(cal$VALUE)), "Valor imposible (negativo o > récord mundial)")
registrar("4b. Quitar valores imposibles", cal$STATION_ID, cal$anio, cal$VALUE_LIMPIO)

# c) Acumulados: tras HUECO_MIN días o más sin dato llega un valor mayor que el
#    percentil 99 de la estación. Es la lluvia de varios días anotada en uno.
#    Se revisan antes que los exagerados porque, si el valor viene después de
#    un hueco, la explicación más probable es que sea un acumulado.
acumulados <- do.call(rbind, lapply(which(!info$lecturas_mensuales), function(k) {
  cols <- primera_col[k]:ultima_col[k]
  v <- V[k, cols]
  lluvia <- v[!is.na(v) & v >= LLUVIA_RED]
  if (length(lluvia) < 30) return(NULL)
  umbral <- quantile(lluvia, PERCENTIL_ACUM, names = FALSE)
  r <- rle(is.na(v))
  fin <- cumsum(r$lengths)
  j <- which(r$values & r$lengths >= HUECO_MIN & fin < length(v))   # hueco seguido de un dato
  if (!length(j)) return(NULL)
  despues <- fin[j] + 1
  alto <- v[despues] > umbral
  if (!any(alto)) return(NULL)
  data.frame(STATION_ID = ids[k], k = k, col = cols[despues[alto]], fecha = fechas[cols[despues[alto]]],
             valor = v[despues[alto]], dias_hueco = r$lengths[j][alto], umbral_p99 = round(umbral, 1),
             stringsAsFactors = FALSE)
}))
acumulados <- acumulados[order(-acumulados$valor), ]
guardar_tabla(acumulados[, setdiff(names(acumulados), c("k", "col"))], "limpieza", "acumulados_tras_hueco")
cat(sprintf("Acumulados tras hueco: %d en %d estaciones (mediana del hueco: %g días)\n",
            nrow(acumulados), length(unique(acumulados$STATION_ID)), median(acumulados$dias_hueco)))
marcar(fila_de(acumulados$k, acumulados$col), "Acumulado tras un hueco")
registrar("4c. Quitar acumulados tras un hueco", cal$STATION_ID, cal$anio, cal$VALUE_LIMPIO)
V[cbind(acumulados$k, acumulados$col)] <- NA

# d) Valores exagerados (posibles errores de digitación, como un cero de más):
#    un valor es exagerado si
#      1. supera FACTOR_EXAGERADO veces el percentil 99 de los días con lluvia
#         de su estación, Y
#      2. ninguna vecina registró al menos RESPALDO_MIN de ese valor el día
#         anterior, el mismo día o el siguiente (un aguacero real de ese tamaño
#         casi siempre deja rastro alrededor).
#    Si ese día ninguna vecina tiene dato, se exige más: FACTOR_SIN_VECINAS
#    veces el percentil 99.
#    Por qué 3 veces: entre los valores de 150 mm o más, los que pasan de 3
#    veces el p99 casi nunca tienen respaldo de las vecinas; por debajo de 3
#    veces, lo normal es que sí lo tengan (ver tabla valores_exagerados).
p99_est <- apply(V, 1, function(v) {
  v <- v[!is.na(v) & v >= LLUVIA_RED]
  if (length(v) < 30) NA else quantile(v, .99, names = FALSE)
})
candidatos <- which(V >= FACTOR_EXAGERADO * p99_est, arr.ind = TRUE)   # p99_est se recicla por fila
exagerados <- do.call(rbind, lapply(seq_len(nrow(candidatos)), function(j) {
  k <- candidatos[j, 1]; c <- candidatos[j, 2]
  N <- vecinas[[k]]
  alrededor <- if (length(N)) V[N, max(1, c - 1):min(n_dias, c + 1)] else NA
  respaldo <- if (all(is.na(alrededor))) NA else max(alrededor, na.rm = TRUE)
  valor <- V[k, c]
  data.frame(STATION_ID = ids[k], k = k, col = c, fecha = fechas[c], valor = valor,
             p99_estacion = round(p99_est[k], 1), veces_p99 = round(valor / p99_est[k], 1),
             max_vecinas_alrededor = respaldo,
             exagerado = if (is.na(respaldo)) valor >= FACTOR_SIN_VECINAS * p99_est[k] else respaldo < RESPALDO_MIN * valor,
             stringsAsFactors = FALSE)
}))
exagerados <- exagerados[order(-exagerados$valor), ]
guardar_tabla(exagerados[, setdiff(names(exagerados), c("k", "col"))], "limpieza", "valores_exagerados")
cat(sprintf("Valores de %d veces el p99 o más: %d | exagerados (sin respaldo de vecinas): %d | se conservan: %d\n",
            FACTOR_EXAGERADO, nrow(exagerados), sum(exagerados$exagerado), sum(!exagerados$exagerado)))
quitar <- exagerados[exagerados$exagerado, ]
marcar(fila_de(quitar$k, quitar$col), "Valor exagerado sin respaldo de vecinas")
registrar("4d. Quitar valores exagerados", cal$STATION_ID, cal$anio, cal$VALUE_LIMPIO)
V[cbind(quitar$k, quitar$col)] <- NA

# e) Ceros falsos: los días de las rachas sospechosas (4.4). Se borran aquí y
#    después se imputan en la Parte 5, como cualquier otro hueco.
marcar(which(!is.na(cal$racha_ceros) & cal$racha_ceros == CLASES[3]), "Cero falso (racha sospechosa)")
registrar("4e. Quitar ceros falsos", cal$STATION_ID, cal$anio, cal$VALUE_LIMPIO)

motivos <- table(cal$motivo)
guardar_tabla(data.frame(motivo = names(motivos), dias = as.vector(motivos)), "limpieza", "motivos_de_na")

graficar(function() {   # Gráfica: ejemplo del acumulado más grande
  a <- acumulados[1, ]
  d <- max(1, a$col - 40):min(n_dias, a$col + 20)
  op <- par(mar = c(4, 5, 4, 1)); on.exit(par(op))
  plot(fechas[d], V[a$k, d], type = "h", lwd = 2, col = COL_PRINCIPAL, xlab = "", ylab = "mm/día",
       main = sprintf("Acumulado tras un hueco: estación %s", a$STATION_ID))
  rect(fechas[a$col - a$dias_hueco], -10, fechas[a$col - 1], 1e5, col = adjustcolor(COL_GRIS, .5), border = NA)
  points(fechas[a$col], a$valor, pch = 19, col = COL_ROJO, cex = 1.4)
  abline(h = a$umbral_p99, lty = 2)
  legend("topleft", bty = "n", pch = c(15, 19, NA), lty = c(NA, NA, 2), col = c(COL_GRIS, COL_ROJO, 1),
         legend = c(sprintf("%d días sin dato", a$dias_hueco), sprintf("%g mm de golpe", a$valor),
                    sprintf("Percentil 99 de la estación (%g mm)", a$umbral_p99)))
})
graficar(function() {   # Gráfica: días que pasan a NA por cada motivo
  m <- sort(motivos[names(motivos) != "Sin dato original"])
  op <- par(mar = c(4, 17, 4, 2)); on.exit(par(op))
  bp <- barplot(m, horiz = TRUE, las = 1, col = COL_ROJO, border = NA, xlim = c(0, max(m) * 1.2),
                xlab = "días que pasan a NA", main = "Datos que se quitan y por qué")
  text(m, bp, format(m, big.mark = ".", decimal.mark = ","), pos = 4, cex = .85)
})

# ---- 4.6 Pares estación-año aptos para curvas --------------------------------
# Una curva anual necesita el año casi completo. Un par estación-año es apto si:
#   - tiene COB_MIN % o más de días con dato (ya limpio),
#   - su hueco más largo no pasa de HUECO_MAX días,
#   - cada temporada de lluvias (mar-may y sep-nov) tiene COB_TEMPORADA % o más.
# No se borra nada: la columna "apto" servirá de filtro más adelante.
cat("\n== 4.6 PARES ESTACIÓN-AÑO APTOS ==\n")
OK <- matrix(FALSE, length(ids), n_dias)
OK[cbind(fila_est, col_dia)] <- !is.na(cal$VALUE_LIMPIO)
pares <- unique(data.frame(STATION_ID = cal$STATION_ID, anio = cal$anio, stringsAsFactors = FALSE))
pares <- do.call(rbind, lapply(split(pares, pares$anio), function(p) {
  cols <- which(anio_col == p$anio[1])
  O <- OK[match(p$STATION_ID, ids), cols, drop = FALSE]
  mam <- mes_col[cols] %in% 3:5
  son <- mes_col[cols] %in% 9:11
  data.frame(p,
             cobertura       = round(rowMeans(O) * 100, 1),
             hueco_mas_largo = apply(!O, 1, racha_mas_larga),
             cobertura_MAM   = round(rowMeans(O[, mam, drop = FALSE]) * 100, 1),
             cobertura_SON   = round(rowMeans(O[, son, drop = FALSE]) * 100, 1))
}))
pares$apto <- pares$cobertura >= COB_MIN & pares$hueco_mas_largo <= HUECO_MAX &
  pares$cobertura_MAM >= COB_TEMPORADA & pares$cobertura_SON >= COB_TEMPORADA
pares$falla <- ifelse(pares$apto, "",
                      paste0(ifelse(pares$cobertura < COB_MIN, "cobertura ", ""),
                             ifelse(pares$hueco_mas_largo > HUECO_MAX, "hueco ", ""),
                             ifelse(pares$cobertura_MAM < COB_TEMPORADA | pares$cobertura_SON < COB_TEMPORADA, "temporada", "")))
pares <- pares[order(pares$STATION_ID, pares$anio), ]
rownames(pares) <- NULL
guardar_tabla(pares, "limpieza", "pares_estacion_anio")
cat(sprintf("Pares aptos: %d de %d (%.1f %%) en %d estaciones\n", sum(pares$apto), nrow(pares),
            mean(pares$apto) * 100, length(unique(pares$STATION_ID[pares$apto]))))
cat("Por qué fallan los demás (puede ser más de una razón):\n")
print(c(cobertura = sum(grepl("cobertura", pares$falla)), hueco = sum(grepl("hueco", pares$falla)),
        temporada = sum(grepl("temporada", pares$falla))))
info$anios_aptos <- as.vector(tapply(pares$apto, pares$STATION_ID, sum)[ids])

cal$apto <- pares$apto[match(paste(cal$STATION_ID, cal$anio), paste(pares$STATION_ID, pares$anio))]
registrar("5. Solo pares estación-año aptos (referencia, no se borran)",
          cal$STATION_ID, cal$anio, ifelse(cal$apto, cal$VALUE_LIMPIO, NA))

graficar(function() {   # Gráfica: pares aptos y no aptos en cada año
  m <- table(factor(pares$apto, c(TRUE, FALSE), c("Aptos", "No aptos")), pares$anio)
  op <- par(mar = c(4, 5, 4, 1)); on.exit(par(op))
  barplot(m, col = c(COL_PRINCIPAL, COL_GRIS), border = NA, las = 2, cex.names = .75, ylab = "estaciones",
          main = sprintf("Pares estación-año aptos para curvas (>= %d %%, hueco <= %d días, MAM y SON >= %d %%)",
                         COB_MIN, HUECO_MAX, COB_TEMPORADA), cex.main = .95)
  legend("topleft", fill = c(COL_PRINCIPAL, COL_GRIS), legend = rownames(m), bty = "n")
})

# ---- 4.7 Píxel del satélite de cada estación ---------------------------------
# Regla de asignación: cada estación CON coordenadas se empareja con el píxel
# de 0,05° que contiene su punto (terra::cellFromXY). Estas parejas son las
# que se comparan en la sección 3 del proyecto y con las que la red KAN
# aprende a corregir el satélite. Las estaciones sin coordenadas no se
# emparejan: no se sabe en qué píxel están.
cat("\n== 4.7 ESTACIÓN Y PÍXEL DEL SATÉLITE ==\n")
celda  <- terra::cellFromXY(plantilla_sat, cbind(info$lon, info$lat))
centro <- terra::xyFromCell(plantilla_sat, ifelse(is.na(celda), 1, celda))
info$pixel <- ifelse(is.na(celda), NA, sprintf("P%04d", celda))
info$km_al_centro_pixel <- ifelse(is.na(celda), NA,
                                  round(dist_km(info$lon, info$lat, centro[, 1], centro[, 2]), 2))
info$ubicacion <- ifelse(!info$con_coordenadas, "Sin coordenadas: no se empareja con el satélite",
                         ifelse(is.na(celda), "Fuera de la imagen", "Dentro de la imagen"))
info$comparte_pixel <- !is.na(info$pixel) & info$pixel %in% info$pixel[duplicated(info$pixel) & !is.na(info$pixel)]
print(table(info$ubicacion))
cat(sprintf("Píxeles distintos con estación: %d | píxeles con 2 o más estaciones: %d (%d estaciones)\n",
            length(unique(na.omit(info$pixel))), sum(table(info$pixel) >= 2), sum(info$comparte_pixel)))
cat(sprintf("Estaciones dentro de la imagen y con algún año apto en %d-%d: %d\n",
            min(PERIODO_COMUN), max(PERIODO_COMUN),
            length(unique(pares$STATION_ID[pares$apto & pares$anio %in% PERIODO_COMUN &
                                             pares$STATION_ID %in% ids[!is.na(info$pixel)]]))))
emparejamiento <- info[!is.na(info$pixel), c("STATION_ID", "nombre", "lon", "lat", "altitud", "pixel",
                                             "km_al_centro_pixel", "comparte_pixel")]
emparejamiento <- emparejamiento[order(emparejamiento$pixel), ]
guardar_tabla(emparejamiento, "limpieza", "emparejamiento_estacion_pixel")

graficar(function() {   # Gráfica: mapa de estaciones con coordenadas y área del satélite
  e <- terra::ext(plantilla_sat)
  u <- info[info$con_coordenadas, ]
  op <- par(mar = c(4, 4, 4, 1)); on.exit(par(op))
  plot(u$lon, u$lat, type = "n", asp = 1, xlab = "Longitud", ylab = "Latitud",
       xlim = range(u$lon, e[1:2]), ylim = range(u$lat, e[3:4]),
       main = sprintf("Estaciones con coordenadas (%d de %d) y área del satélite", nrow(u), nrow(info)))
  rect(e[1], e[3], e[2], e[4], border = COL_CERO, lwd = 2, lty = 2)
  apto <- !is.na(u$anios_aptos) & u$anios_aptos > 0
  points(u$lon, u$lat, pch = 19, cex = .9,
         col = ifelse(u$lecturas_mensuales, COL_ROJO, ifelse(apto, COL_PRINCIPAL, COL_GRIS)))
  legend("bottomleft", bty = "n", pch = c(19, 19, 19, NA), lty = c(NA, NA, NA, 2), lwd = c(NA, NA, NA, 2),
         col = c(COL_PRINCIPAL, COL_GRIS, COL_ROJO, COL_CERO), cex = .85,
         legend = c("Con algún año apto", "Sin años aptos", "Lecturas mensuales", "Imagen del satélite"))
})


# =============================================================================
# PARTE 5. IMPUTACIÓN: FPCA CONDICIONAL Y VAE
# =============================================================================
# Escala: PÉNTADAS, igual que el satélite (6 por mes, 72 por año), en mm/día.
# Una "curva" es un año de una estación: 72 valores.
#
# Modelo que usan los dos métodos (en raíz cuadrada, para que la lluvia sea
# menos asimétrica):
#     raíz(lluvia) = climatología de la estación
#                  + a · anomalía de sus vecinas ese mismo año
#                  + residuo
#   - climatología: la curva media de la estación (años de entrenamiento).
#   - vecinas: lo que llovió de más o de menos en las estaciones cercanas en
#     esas mismas péntadas. Es la pista más fuerte para rellenar un hueco.
#   - residuo: lo que las vecinas no explican. ESTO es lo que imputan la FPCA
#     y el VAE, usando las péntadas observadas de la propia curva.
#
# Imputación oficial: FPCA condicional (ganó en la validación y es más simple
# e interpretable). El VAE se mantiene solo como COMPARACIÓN.
#
# Reglas contra la fuga de datos: la climatología, el coeficiente a, las
# componentes de la FPCA y los pesos del VAE se aprenden SOLO con las curvas
# de ajuste (años hasta ANIO_REFERENCIA, sin las curvas apartadas para validar).
cat("\n######## PARTE 5. IMPUTACIÓN ########\n")
set.seed(SEMILLA)

# ---- 5.1 De días a péntadas --------------------------------------------------
# Péntada del día: días 1-5 -> 1, 6-10 -> 2, ..., 26 al fin de mes -> 6.
# Una péntada cuenta como observada si tiene al menos el 60 % de sus días con
# dato limpio; su valor es el promedio de esos días (mm/día).
cal$pentada <- (cal$mes - 1) * 6 + pmin(ceiling(cal$dia / 5), 6)
usar <- !cal$STATION_ID %in% ids[info$lecturas_mensuales]
curvas <- unique(data.frame(STATION_ID = cal$STATION_ID[usar], anio = cal$anio[usar], stringsAsFactors = FALSE))
curvas <- curvas[order(curvas$STATION_ID, curvas$anio), ]
rownames(curvas) <- NULL
clave_curva <- paste(curvas$STATION_ID, curvas$anio)
fila_curva <- match(paste(cal$STATION_ID, cal$anio), clave_curva)
con_valor <- !is.na(cal$VALUE_LIMPIO) & !is.na(fila_curva)
casilla <- (fila_curva[con_valor] - 1) * 72 + cal$pentada[con_valor]     # posición en la matriz curva x péntada
n_casillas <- nrow(curvas) * 72
dias_con_dato <- tabulate(casilla, n_casillas)
suma <- rowsum(cal$VALUE_LIMPIO[con_valor], casilla)
total <- numeric(n_casillas)
total[as.integer(rownames(suma))] <- suma[, 1]
Y <- matrix(total / dias_con_dato, ncol = 72, byrow = TRUE)               # mm/día; filas = curvas
largo_pent <- sapply(anios, function(a) dias_pentada(a, 1:72))            # 72 x años
dias_largo <- t(largo_pent[, match(curvas$anio, anios)])
Y[matrix(dias_con_dato, ncol = 72, byrow = TRUE) < ceiling(.6 * dias_largo)] <- NA
rm(suma, total, casilla)

curvas$n_obs <- rowSums(!is.na(Y))
curvas$imputable <- curvas$n_obs >= MIN_PENTADAS_IMPUTAR
curvas$entrenamiento <- curvas$anio <= ANIO_REFERENCIA
cat(sprintf("Curvas (estación-año): %d | completas: %d | imputables (>= %d péntadas): %d | no imputables: %d\n",
            nrow(curvas), sum(curvas$n_obs == 72), MIN_PENTADAS_IMPUTAR, sum(curvas$imputable),
            sum(!curvas$imputable)))
registrar("6. Pasar a péntadas (antes de imputar)", rep(curvas$STATION_ID, 72), rep(curvas$anio, 72),
          as.vector(Y), unidad = "péntadas")

# ---- 5.2 Conjuntos para aprender y para validar ------------------------------
# Para saber qué tan bien imputa cada método se ESCONDEN péntadas que sí se
# conocen (un bloque de 6 seguidas, como un mes perdido, y 6 sueltas) y se
# mide qué tan bien se recuperan. Así se compara con datos reales.
#   - prueba:  curvas completas que ningún método ve al aprender. Mientras no
#              haya partición oficial, son el 15 % de las curvas completas,
#              sorteadas. Si ANIO_REFERENCIA < ANIO_FIN, son las curvas
#              completas de los años posteriores.
#   - interna: otro 15 % de las curvas completas, para escoger K (FPCA) y
#              saber cuándo parar el VAE.
#   - ajuste:  todas las demás curvas imputables, con las que se aprende.
esconder <- function(observado) {
  H <- observado & FALSE
  for (i in seq_len(nrow(observado))) {
    inicio <- sample(1:67, 1)
    H[i, inicio:(inicio + 5)] <- TRUE
    H[i, sample(72, 6)] <- TRUE
  }
  H & observado
}
entrena <- which(curvas$entrenamiento & curvas$imputable)
completas <- entrena[curvas$n_obs[entrena] == 72]
if (any(!curvas$entrenamiento)) {
  prueba <- which(!curvas$entrenamiento & curvas$n_obs == 72)
} else {
  prueba <- sample(completas, round(.15 * length(completas)))
}
interna <- sample(setdiff(completas, prueba), round(.15 * length(completas)))
ajuste  <- setdiff(entrena, c(prueba, interna))
cat(sprintf("Curvas para aprender: %d | para ajustar parámetros: %d | de prueba: %d\n",
            length(ajuste), length(interna), length(prueba)))

# ---- 5.3 Climatología de cada estación y anomalías ---------------------------
# La climatología se calcula SOLO con las curvas de ajuste (si una estación
# tiene menos de 3, con todas sus curvas imputables que no sean de prueba).
RAIZ_Y <- sqrt(Y)
clim <- t(sapply(split(seq_len(nrow(curvas)), curvas$STATION_ID), function(i) {
  j <- intersect(i, ajuste)
  if (length(j) < 3) j <- setdiff(i[curvas$imputable[i]], prueba)
  if (length(j) == 0) j <- setdiff(i, prueba)
  if (length(j) == 0) j <- i
  colMeans(RAIZ_Y[j, , drop = FALSE], na.rm = TRUE)
}))
clim_red <- colMeans(clim, na.rm = TRUE)                  # si a una estación le falta una péntada, se usa la red
for (j in 1:72) clim[is.na(clim[, j]), j] <- clim_red[j]
CLIM <- clim[curvas$STATION_ID, ]                         # climatología de cada curva
X <- RAIZ_Y - CLIM                                        # anomalía: cuánto se aparta de lo normal

# ---- 5.4 Lo que dicen las vecinas --------------------------------------------
# Vecinas de imputación: las mismas de 4.3 (por distancia si la estación tiene
# coordenadas, por correlación si no). Si una estación se quedó sin ninguna,
# se usan sus 5 estaciones más correlacionadas.
vecinas_imp <- lapply(seq_along(ids), function(k) {
  if (length(vecinas[[k]])) return(ids[vecinas[[k]]])
  r <- CORR[k, ]
  o <- order(-r)
  ids[head(o[!is.na(r[o])], 5)]
})
names(vecinas_imp) <- ids
# Curva de las vecinas: promedio de sus anomalías en el mismo año, péntada a péntada
VEC <- t(sapply(seq_len(nrow(curvas)), function(i) {
  f <- match(paste(vecinas_imp[[curvas$STATION_ID[i]]], curvas$anio[i]), clave_curva)
  f <- f[!is.na(f)]
  if (!length(f)) return(rep(NA_real_, 72))
  m <- colMeans(X[f, , drop = FALSE], na.rm = TRUE)
  m[is.nan(m)] <- NA
  m
}))
VEC0 <- ifelse(is.na(VEC), 0, VEC)                        # sin información de vecinas = 0 (no mueve nada)
cat(sprintf("Curvas con información de vecinas: %d de %d\n", sum(rowSums(!is.na(VEC)) > 0), nrow(curvas)))
H_int <- esconder(!is.na(X[interna, ]))
H_pru <- esconder(!is.na(X[prueba, ]))

# ---- 5.5 Coeficiente a: cuánto se parece la estación a sus vecinas -----------
# Mínimos cuadrados por estación: a = suma(x·v) / suma(v²). Si la estación
# tiene pocos datos de entrenamiento, se usa el a de toda la red.
filas_aj <- split(ajuste, curvas$STATION_ID[ajuste])
num <- sapply(filas_aj, function(i) sum(X[i, ] * VEC0[i, ], na.rm = TRUE))
den <- sapply(filas_aj, function(i) sum(VEC0[i, ]^2 * !is.na(X[i, ])))
a_red <- sum(num) / sum(den)
a_est <- ifelse(den >= 50, num / den, a_red)
a_curva <- ifelse(is.na(a_est[curvas$STATION_ID]), a_red, a_est[curvas$STATION_ID])
cat(sprintf("Coeficiente a de la red: %.2f (mediana por estación: %.2f)\n", a_red, median(a_est)))
ESPACIAL <- as.vector(a_curva) * VEC0                     # lo que explican las vecinas
RES <- X - ESPACIAL                                       # residuo que imputan FPCA y VAE

# Devuelve la lluvia (mm/día) a partir del residuo: se deshace la raíz y se
# suma "corr", la varianza que se pierde al elevar al cuadrado (si no, la
# lluvia imputada sale sistemáticamente baja).
a_lluvia <- function(idx, residuo, corr) pmax(CLIM[idx, , drop = FALSE] + ESPACIAL[idx, , drop = FALSE] + residuo, 0)^2 + corr
varianza_perdida <- function(idx, H, residuo) {
  mean(((CLIM[idx, ] + ESPACIAL[idx, ] + residuo)[H] - RAIZ_Y[idx, ][H])^2)
}
error_rmse <- function(idx, H, y_est) sqrt(mean((y_est[H] - Y[idx, ][H])^2))

# ---- 5.6 FPCA condicional ----------------------------------------------------
# 1. Covarianza del residuo entre péntadas, C(s, t), con las curvas de ajuste.
# 2. Ruido (como en PACE): la diagonal de C trae además el ruido de cada
#    péntada. Se estima comparando la diagonal con sus vecinos inmediatos,
#    sigma² = promedio de [C(t, t) - (C(t, t-1) + C(t, t+1)) / 2].
# 3. Autofunciones y autovalores de la covarianza sin ruido (malla regular de
#    72 péntadas, pesos iguales).
# 4. Para una curva con péntadas observadas O, los puntajes se estiman con la
#    esperanza condicional:
#       puntajes = (Psi_O' Psi_O + sigma² Lambda^-1)^-1 Psi_O' (r_O - mu_O)
#    y las péntadas faltantes son mu + Psi · puntajes.
fpca_ajustar <- function(R) {
  mu <- colMeans(R, na.rm = TRUE)
  C <- cov(R, use = "pairwise.complete.obs")
  p <- ncol(C)
  sin_ruido <- sapply(1:p, function(t) mean(C[t, c(t - 1, t + 1)[c(t > 1, t < p)]]))
  sigma2 <- max(mean(diag(C) - sin_ruido), 1e-6)
  diag(C) <- sin_ruido
  e <- eigen(C, symmetric = TRUE)
  positivos <- e$values > 0
  list(mu = mu, lambda = e$values[positivos], psi = e$vectors[, positivos, drop = FALSE], sigma2 = sigma2)
}
fpca_imputar <- function(modelo, R, K) {
  psi <- modelo$psi[, 1:K, drop = FALSE]
  t(sapply(seq_len(nrow(R)), function(i) {
    r <- R[i, ]
    O <- !is.na(r)
    if (all(O)) return(r)
    psi_O <- psi[O, , drop = FALSE]
    puntajes <- solve(crossprod(psi_O) + modelo$sigma2 * diag(1 / modelo$lambda[1:K], K),
                      crossprod(psi_O, r[O] - modelo$mu[O]))
    r[!O] <- (modelo$mu + psi %*% puntajes)[!O]
    r
  }))
}
con_huecos <- function(idx, H) { R <- RES[idx, ]; R[H] <- NA; R }

fpca <- fpca_ajustar(RES[ajuste, ])
cat(sprintf("\nFPCA: sigma² del ruido = %.3f | componentes con varianza positiva: %d\n",
            fpca$sigma2, length(fpca$lambda)))
cat("Varianza explicada acumulada (%) por 1, 2, 3, 5 y 10 componentes:",
    round(cumsum(fpca$lambda)[c(1, 2, 3, 5, 10)] / sum(fpca$lambda) * 100, 1), "\n")
K_PROBAR <- c(1:8, 10, 15, 20)
error_K <- sapply(K_PROBAR, function(K) {
  r_est <- fpca_imputar(fpca, con_huecos(interna, H_int), K)
  corr <- varianza_perdida(interna, H_int, r_est)
  error_rmse(interna, H_int, a_lluvia(interna, r_est, corr))
})
K_FPCA <- K_PROBAR[which.min(error_K)]
corr_fpca <- varianza_perdida(interna, H_int, fpca_imputar(fpca, con_huecos(interna, H_int), K_FPCA))
cat(sprintf("K escogido: %d componentes (error en el conjunto interno: %.3f mm/día)\n", K_FPCA, min(error_K)))

graficar(function() {   # Gráfica: primeras autofunciones del residuo
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  hasta <- min(length(fpca$lambda), max(15, K_FPCA + 3))
  plot(1:hasta, cumsum(fpca$lambda)[1:hasta] / sum(fpca$lambda) * 100, type = "b", pch = 19, col = COL_PRINCIPAL,
       xlab = "Número de componentes", ylab = "% de varianza explicada", main = "Varianza explicada (FPCA)")
  abline(v = K_FPCA, lty = 2)
  legend("bottomright", lty = 2, bty = "n", legend = sprintf("K escogido = %d (%.0f %%)", K_FPCA,
                                                             sum(fpca$lambda[1:K_FPCA]) / sum(fpca$lambda) * 100))
  matplot(1:72, fpca$psi[, 1:3], type = "l", lty = 1, lwd = 2, col = c(COL_PRINCIPAL, COL_CERO, COL_ROJO),
          xlab = "Péntada del año", ylab = "Autofunción", main = "Tres primeras autofunciones")
  abline(h = 0, col = COL_GRIS)
  legend("topright", lty = 1, lwd = 2, col = c(COL_PRINCIPAL, COL_CERO, COL_ROJO), bty = "n",
         legend = paste("Componente", 1:3))
})

# ---- 5.7 VAE (autocodificador variacional) -----------------------------------
# Red neuronal escrita a mano en R base (sin paquetes):
#   codificador: [residuo con huecos en 0, máscara de observado] -> capa oculta
#                (tanh) -> media y log-varianza de z (espacio latente)
#   z = media + exp(logvar/2) · ruido normal      (truco de reparametrización)
#   decodificador: z -> capa oculta (tanh) -> curva de 72 péntadas
# Pérdida (ELBO negativo): error cuadrático SOLO en las péntadas observadas
#   + divergencia KL de z contra una normal estándar.
# Para que aprenda a rellenar, en cada lote se le esconden péntadas que sí se
# conocen (como en la validación) y se le pide reconstruirlas.
# Los gradientes se calculan con retropropagación y los pesos se actualizan con
# Adam. Se para cuando el error del conjunto interno deja de mejorar.
vae_iniciar <- function(p, h, q) {
  w <- function(a, b) matrix(rnorm(a * b, sd = sqrt(2 / (a + b))), a, b)   # inicialización de Xavier
  list(W1 = w(2 * p, h), b1 = rep(0, h),    # codificador
       Wm = w(h, q),     bm = rep(0, q),    # media de z
       Wv = w(h, q),     bv = rep(0, q),    # log-varianza de z
       W2 = w(q, h),     b2 = rep(0, h),    # decodificador
       W3 = w(h, p),     b3 = rep(0, p))    # salida
}
vae_adelante <- function(red, Z, M, ruido = NULL) {
  entrada <- cbind(Z * M, M)
  h1 <- tanh(sweep(entrada %*% red$W1, 2, red$b1, "+"))
  mu <- sweep(h1 %*% red$Wm, 2, red$bm, "+")
  lv <- sweep(h1 %*% red$Wv, 2, red$bv, "+")
  z  <- if (is.null(ruido)) mu else mu + exp(lv / 2) * ruido
  h2 <- tanh(sweep(z %*% red$W2, 2, red$b2, "+"))
  salida <- sweep(h2 %*% red$W3, 2, red$b3, "+")
  list(entrada = entrada, h1 = h1, mu = mu, lv = lv, z = z, h2 = h2, salida = salida)
}
# M_entra: lo que la red ve; M_eval: lo observado de verdad (donde se mide el error)
vae_gradientes <- function(red, Z, M_entra, M_eval, ruido) {
  n <- nrow(Z)
  a <- vae_adelante(red, Z, M_entra, ruido)
  err <- M_eval * (a$salida - Z)
  d_sal <- err / n
  g <- list()
  g$W3 <- t(a$h2) %*% d_sal;  g$b3 <- colSums(d_sal)
  d_a2 <- (d_sal %*% t(red$W3)) * (1 - a$h2^2)
  g$W2 <- t(a$z) %*% d_a2;    g$b2 <- colSums(d_a2)
  d_z  <- d_a2 %*% t(red$W2)
  d_mu <- d_z + a$mu / n
  d_lv <- d_z * ruido * 0.5 * exp(a$lv / 2) + 0.5 * (exp(a$lv) - 1) / n
  g$Wm <- t(a$h1) %*% d_mu;   g$bm <- colSums(d_mu)
  g$Wv <- t(a$h1) %*% d_lv;   g$bv <- colSums(d_lv)
  d_a1 <- (d_mu %*% t(red$Wm) + d_lv %*% t(red$Wv)) * (1 - a$h1^2)
  g$W1 <- t(a$entrada) %*% d_a1; g$b1 <- colSums(d_a1)
  g[names(red)]
}
# Escala del VAE: residuo centrado y dividido por su desviación estándar
centro_vae <- fpca$mu
escala_vae <- sd(RES[entrena, ], na.rm = TRUE)
Z_vae <- sweep(RES, 2, centro_vae) / escala_vae
OBS <- !is.na(Z_vae) * 1
Z_vae[is.na(Z_vae)] <- 0
vae_imputar <- function(red, idx, H = NULL) {
  M <- OBS[idx, , drop = FALSE]; Z <- Z_vae[idx, , drop = FALSE]
  if (!is.null(H)) { M[H] <- 0; Z[H] <- 0 }
  rec <- vae_adelante(red, Z, M)$salida * escala_vae + rep(centro_vae, each = length(idx))
  r <- RES[idx, , drop = FALSE]
  if (!is.null(H)) r[H] <- NA
  r[is.na(r)] <- rec[is.na(r)]
  r
}

red <- vae_iniciar(72, VAE_OCULTA, VAE_LATENTE)
m1 <- lapply(red, function(x) x * 0); m2 <- m1; paso <- 0      # memorias de Adam
mejor_error <- Inf; mejor_red <- red; sin_mejorar <- 0; historia <- c()
for (epoca in 1:VAE_EPOCAS_MAX) {
  orden <- sample(ajuste)
  for (lote in split(orden, ceiling(seq_along(orden) / VAE_LOTE))) {
    M_eval <- OBS[lote, ]
    M_entra <- M_eval * !esconder(M_eval == 1)              # se esconden péntadas conocidas
    ruido <- matrix(rnorm(length(lote) * VAE_LATENTE), length(lote))
    g <- vae_gradientes(red, Z_vae[lote, ], M_entra, M_eval, ruido)
    paso <- paso + 1
    for (w in names(red)) {                                  # Adam
      m1[[w]] <- .9 * m1[[w]] + .1 * g[[w]]
      m2[[w]] <- .999 * m2[[w]] + .001 * g[[w]]^2
      red[[w]] <- red[[w]] - VAE_TASA * (m1[[w]] / (1 - .9^paso)) / (sqrt(m2[[w]] / (1 - .999^paso)) + 1e-8)
    }
  }
  if (epoca %% 10 == 0) {                                    # cada 10 épocas: ¿mejoró?
    r_est <- vae_imputar(red, interna, H_int)
    e <- error_rmse(interna, H_int, a_lluvia(interna, r_est, varianza_perdida(interna, H_int, r_est)))
    historia <- c(historia, e)
    if (e < mejor_error - 1e-4) { mejor_error <- e; mejor_red <- red; sin_mejorar <- 0 }
    else { sin_mejorar <- sin_mejorar + 1; if (sin_mejorar >= 4) break }
  }
}
red <- mejor_red
corr_vae <- varianza_perdida(interna, H_int, vae_imputar(red, interna, H_int))
cat(sprintf("VAE: %d épocas | error en el conjunto interno: %.3f mm/día\n", epoca, mejor_error))

graficar(function() {   # Gráfica: aprendizaje del VAE
  op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  plot(seq_along(historia) * 10, historia, type = "b", pch = 19, col = COL_PRINCIPAL,
       xlab = "Época", ylab = "Error en el conjunto interno (mm/día)", main = "Entrenamiento del VAE")
})

# ---- 5.8 Validación: péntadas escondidas de las curvas de prueba -------------
# Se comparan cuatro formas de rellenar las péntadas escondidas de prueba:
#   climatología sola, climatología + vecinas, + FPCA, + VAE.
cat("\n== 5.8 VALIDACIÓN (péntadas escondidas de las curvas de prueba) ==\n")
corr_vec <- varianza_perdida(interna, H_int, matrix(0, length(interna), 72))
imputar_prueba <- function(H) list(
  "Climatología de la estación"  = CLIM[prueba, ]^2,
  "Climatología + vecinas"       = a_lluvia(prueba, matrix(0, length(prueba), 72), corr_vec),
  "FPCA condicional (+ vecinas)" = a_lluvia(prueba, fpca_imputar(fpca, con_huecos(prueba, H), K_FPCA), corr_fpca),
  "VAE (+ vecinas)"              = a_lluvia(prueba, vae_imputar(red, prueba, H), corr_vae))
y_pru <- imputar_prueba(H_pru)
verdad <- Y[prueba, ][H_pru]
validacion <- data.frame(
  metodo      = names(y_pru),
  RMSE        = sapply(y_pru, function(y) sqrt(mean((y[H_pru] - verdad)^2))),
  MAE         = sapply(y_pru, function(y) mean(abs(y[H_pru] - verdad))),
  sesgo       = sapply(y_pru, function(y) mean(y[H_pru] - verdad)),
  correlacion = sapply(y_pru, function(y) cor(y[H_pru], verdad)),
  row.names = NULL)
validacion[, -1] <- round(validacion[, -1], 3)
print(validacion, row.names = FALSE)
guardar_tabla(validacion, "imputacion", "validacion_imputacion")

# Validación repetida: se esconden otras péntadas REPETICIONES_VALIDACION veces
# (los modelos ya entrenados no cambian) para ver si la diferencia entre
# métodos se mantiene o fue suerte del sorteo. Se reporta el intervalo que
# contiene el 95 % de las diferencias de RMSE.
rmse_rep <- t(sapply(seq_len(REPETICIONES_VALIDACION), function(r) {
  H <- esconder(!is.na(X[prueba, ]))
  y <- imputar_prueba(H)
  sapply(y, function(m) sqrt(mean((m[H] - Y[prueba, ][H])^2)))
}))
resumir_dif <- function(d) c(media = mean(d), ic_inferior = quantile(d, .025, names = FALSE),
                             ic_superior = quantile(d, .975, names = FALSE), sorteos_que_gana_el_primero = sum(d < 0))
diferencias <- data.frame(
  comparacion = c("FPCA - VAE", "FPCA - Climatología + vecinas", "VAE - Climatología + vecinas"),
  round(rbind(resumir_dif(rmse_rep[, 3] - rmse_rep[, 4]),
              resumir_dif(rmse_rep[, 3] - rmse_rep[, 2]),
              resumir_dif(rmse_rep[, 4] - rmse_rep[, 2])), 4),
  row.names = NULL)
cat(sprintf("\nValidación repetida (%d sorteos). Diferencia de RMSE (negativa = gana el primero):\n",
            REPETICIONES_VALIDACION))
print(diferencias, row.names = FALSE)
guardar_tabla(diferencias, "imputacion", "validacion_repetida")

graficar(function() {   # Gráfica: error de cada método en la validación
  op <- par(mar = c(4.5, 15, 3, 2)); on.exit(par(op))
  bp <- barplot(rev(validacion$RMSE), names.arg = rev(validacion$metodo), horiz = TRUE, las = 1,
                col = rev(c(COL_GRIS, COL_GRIS, COL_PRINCIPAL, COL_CERO)), border = NA,
                xlim = c(0, max(validacion$RMSE) * 1.15), xlab = "RMSE (mm/día) en péntadas escondidas",
                main = "¿Qué tan bien se recuperan datos reales escondidos?")
  text(rev(validacion$RMSE), bp, format(rev(validacion$RMSE), decimal.mark = ","), pos = 4, cex = .85)
})
graficar(function() {   # Gráfica: RMSE de cada sorteo de la validación repetida
  op <- par(mar = c(4.5, 15, 3, 2)); on.exit(par(op))
  boxplot(rmse_rep[, 4:1], horizontal = TRUE, las = 1, col = c(COL_CERO, COL_PRINCIPAL, COL_GRIS, COL_GRIS),
          names = rev(colnames(rmse_rep)), xlab = "RMSE (mm/día)",
          main = sprintf("Validación repetida (%d sorteos de péntadas escondidas)", REPETICIONES_VALIDACION))
})
graficar(function() {   # Gráfica: ejemplo de una curva con un mes escondido
  i <- which.max(rowSums(H_pru))
  op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  plot(1:72, Y[prueba[i], ], type = "h", lwd = 3, col = COL_GRIS, xlab = "Péntada del año", ylab = "mm/día",
       ylim = c(0, max(Y[prueba[i], ], y_pru[[3]][i, ], y_pru[[4]][i, ], na.rm = TRUE)),
       main = sprintf("Estación %s, año %d: péntadas escondidas y su imputación",
                      curvas$STATION_ID[prueba[i]], curvas$anio[prueba[i]]))
  escondidas <- which(H_pru[i, ])
  points(escondidas, Y[prueba[i], escondidas], pch = 19, col = 1)
  points(escondidas, y_pru[[3]][i, escondidas], pch = 17, col = COL_PRINCIPAL, cex = 1.2)
  points(escondidas, y_pru[[4]][i, escondidas], pch = 15, col = COL_CERO, cex = 1.1)
  legend("topright", bty = "n", pch = c(NA, 19, 17, 15), lty = c(1, NA, NA, NA), lwd = c(3, NA, NA, NA),
         col = c(COL_GRIS, 1, COL_PRINCIPAL, COL_CERO),
         legend = c("Observado", "Valor real escondido", "FPCA", "VAE"))
})

# ---- 5.9 Imputación final con FPCA -------------------------------------------
# Se rellenan con FPCA las péntadas faltantes de TODAS las curvas imputables.
# Lo observado no se toca. Las curvas no imputables (menos de la mitad de las
# péntadas) se quedan con sus huecos. El VAE también se aplica, pero solo
# para comparar sus estadísticas con las de la FPCA (5.10).
cat("\n== 5.9 IMPUTACIÓN FINAL (FPCA) ==\n")
imp <- which(curvas$imputable & curvas$n_obs < 72)
Y_FPCA <- Y; Y_VAE <- Y
falta <- is.na(Y[imp, , drop = FALSE])
Y_FPCA[imp, ][falta] <- a_lluvia(imp, fpca_imputar(fpca, RES[imp, , drop = FALSE], K_FPCA), corr_fpca)[falta]
Y_VAE[imp, ][falta]  <- a_lluvia(imp, vae_imputar(red, imp), corr_vae)[falta]
IMPUTADA <- is.na(Y) & !is.na(Y_FPCA)

# Por qué faltaba cada péntada: el motivo más frecuente entre sus días sin dato
sin_valor <- is.na(cal$VALUE_LIMPIO) & !is.na(fila_curva)
casilla_na <- (fila_curva[sin_valor] - 1) * 72 + cal$pentada[sin_valor]
motivo_casilla <- tapply(cal$motivo[sin_valor], casilla_na, function(x) names(which.max(table(x))))
pos <- as.integer(names(motivo_casilla))
MOTIVO <- matrix(NA_character_, nrow(curvas), 72)
MOTIVO[cbind((pos - 1) %/% 72 + 1, (pos - 1) %% 72 + 1)] <- motivo_casilla
# Péntadas sin ningún día en el calendario: caen antes del primer o después del
# último registro de la estación (en su primer o último año de operación)
MOTIVO[IMPUTADA & is.na(MOTIVO)] <- "Fuera del periodo de operación de la estación"
imputado_por_motivo <- as.data.frame(table(motivo = MOTIVO[IMPUTADA]), stringsAsFactors = FALSE)
names(imputado_por_motivo)[2] <- "pentadas_imputadas"
imputado_por_motivo <- imputado_por_motivo[order(-imputado_por_motivo$pentadas_imputadas), ]
cat("Péntadas imputadas según por qué faltaban:\n")
print(imputado_por_motivo, row.names = FALSE)
guardar_tabla(imputado_por_motivo, "imputacion", "imputado_por_motivo")

resumen_imputacion <- data.frame(
  concepto = c("Curvas (estación-año)", "Curvas completas sin imputar", "Curvas imputadas",
               "Curvas no imputables (menos de la mitad de las péntadas)",
               "Péntadas observadas", "Péntadas imputadas con FPCA", "Péntadas que siguen sin dato"),
  valor = c(nrow(curvas), sum(curvas$n_obs == 72), length(imp), sum(!curvas$imputable),
            sum(!is.na(Y)), sum(IMPUTADA), sum(is.na(Y_FPCA))))
print(resumen_imputacion, row.names = FALSE)
guardar_tabla(resumen_imputacion, "imputacion", "resumen_imputacion")

# Trazabilidad de la imputación (en péntadas)
registrar("7. Imputar con FPCA (curvas con al menos la mitad de sus péntadas)",
          rep(curvas$STATION_ID, 72), rep(curvas$anio, 72), as.vector(Y_FPCA), unidad = "péntadas")
completa_final <- rowSums(is.na(Y_FPCA)) == 0
registrar("8. Curvas completas listas para construir datos funcionales",
          rep(curvas$STATION_ID, 72), rep(curvas$anio, 72),
          ifelse(rep(completa_final, 72), as.vector(Y_FPCA), NA), unidad = "péntadas")

# ---- 5.10 Tablas descriptivas en Excel: antes, FPCA y VAE --------------------
# Un libro de Excel con:
#   Resumen      estadísticos de toda la red: antes (diario y péntadas), FPCA y VAE
#   Antes / FPCA / VAE   los mismos estadísticos por estación (péntadas, mm/día)
#   Por_mes      lluvia media de cada mes en las tres versiones
#   Validacion   error de cada método con datos escondidos (y repetida)
#   Imputacion   cuántas péntadas se observaron e imputaron, y por qué faltaban
#   Trazabilidad todos los pasos, de los datos originales a las curvas completas
# La imputación oficial es la FPCA; el VAE aparece solo para comparar.
cat("\n== 5.10 TABLAS DESCRIPTIVAS (Excel) ==\n")
NOMBRES_ESTAD <- c("esperados", "con_dato", "nulos", "pct_nulos", "ceros", "pct_ceros", "media", "mediana",
                   "moda", "desv_estandar", "minimo", "q1", "q3", "p99", "maximo")
por_estacion <- function(M) {
  t(sapply(split(seq_len(nrow(curvas)), curvas$STATION_ID), function(i) describir(as.vector(M[i, ]), length(i) * 72)))
}
tabla_estacion <- function(M) {
  d <- por_estacion(M)
  data.frame(STATION_ID = rownames(d), round(d, 3), row.names = NULL, check.names = FALSE)
}
resumen_excel <- data.frame(
  version = c("Antes de imputar (diario)", "Antes de imputar (péntadas)", "Imputado con FPCA (péntadas)",
              "Imputado con VAE (péntadas)"),
  round(rbind(describir(cal$VALUE_LIMPIO[usar], sum(usar)),
              describir(as.vector(Y), length(Y)),
              describir(as.vector(Y_FPCA), length(Y)),
              describir(as.vector(Y_VAE), length(Y))), 3),
  check.names = FALSE)
mes_de_pentada <- rep(ceiling(1:72 / 6), each = nrow(curvas))
por_mes <- data.frame(
  mes          = month.abb,
  antes        = round(tapply(as.vector(Y), mes_de_pentada, mean, na.rm = TRUE), 3),
  FPCA         = round(tapply(as.vector(Y_FPCA), mes_de_pentada, mean, na.rm = TRUE), 3),
  VAE          = round(tapply(as.vector(Y_VAE), mes_de_pentada, mean, na.rm = TRUE), 3),
  pct_imputado = round(tapply(as.vector(is.na(Y) & !is.na(Y_FPCA)), mes_de_pentada, mean) * 100, 2))
print(resumen_excel, row.names = FALSE)
writexl::write_xlsx(
  list(Resumen = resumen_excel, Antes = tabla_estacion(Y), FPCA = tabla_estacion(Y_FPCA),
       VAE = tabla_estacion(Y_VAE), Por_mes = por_mes, Validacion = validacion,
       Validacion_repetida = diferencias, Imputacion = resumen_imputacion,
       Imputado_por_motivo = imputado_por_motivo, Trazabilidad = trazabilidad),
  file.path(RUTA_SALIDAS, CARPETAS[["imputacion"]], "estadisticas_descriptivas.xlsx"))
# Segundo libro: revisión de estaciones (para el informe y para el profesor)
writexl::write_xlsx(
  list(Lecturas_mensuales = lecturas_mensuales,
       Valores_exagerados = exagerados[, setdiff(names(exagerados), c("k", "col"))],
       Calibracion_correlacion = calibracion,
       Validacion_vecinas = validacion_vecinas, Vecinas_sin_coordenadas = vecinas_sin_coordenadas,
       Emparejamiento_pixel = emparejamiento, Sin_coordenadas = sin_coord),
  file.path(RUTA_SALIDAS, CARPETAS[["limpieza"]], "revision_estaciones.xlsx"))

graficar(function() {   # Gráfica: lluvia media por mes antes y después de imputar
  op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
  m <- t(as.matrix(por_mes[, c("antes", "FPCA", "VAE")]))
  barplot(m, beside = TRUE, names.arg = por_mes$mes, col = c(COL_GRIS, COL_PRINCIPAL, COL_CERO), border = NA,
          ylab = "mm/día", main = "Lluvia media por mes: antes y después de imputar")
  legend("topleft", fill = c(COL_GRIS, COL_PRINCIPAL, COL_CERO), legend = c("Antes", "FPCA", "VAE"), bty = "n")
})


# =============================================================================
# PARTE 6. GUARDAR
# =============================================================================
# Estos archivos son la entrada de los pasos siguientes (curvas, comparación, KAN).
cat("\n######## PARTE 6. GUARDAR ########\n")
guardar_tabla(info[, c("STATION_ID", "zona", "nombre", "municipio", "codigo_cuenca", "cuenca", "lon", "lat",
                       "altitud", "ubicacion", "pixel", "km_al_centro_pixel", "comparte_pixel",
                       "tipo_vecinas", "n_vecinas", "r_max_vecina", "vecinas", "lecturas_mensuales",
                       "anios_mensuales", "anios_aptos")], "limpieza", "estaciones_info")

# Datos diarios limpios (sin imputar)
estaciones_limpias <- cal[, c("STATION_ID", "zona", "fecha", "anio", "mes", "dia", "pentada", "VALUE",
                              "VALUE_LIMPIO", "motivo", "racha_ceros", "apto")]
carpeta <- file.path(RUTA_SALIDAS, CARPETAS[["limpieza"]])
saveRDS(estaciones_limpias, file.path(carpeta, "estaciones_limpias_1970_2003.rds"))
write.csv(estaciones_limpias, gzfile(file.path(carpeta, "estaciones_limpias_1970_2003.csv.gz")), row.names = FALSE)
saveRDS(satelite, file.path(carpeta, "satelite_1981_2003.rds"))
write.csv(satelite, gzfile(file.path(carpeta, "satelite_1981_2003.csv.gz")), row.names = FALSE)

# Péntadas: el dato observado y el imputado con FPCA, en columnas separadas
pentadas <- data.frame(
  STATION_ID      = rep(curvas$STATION_ID, times = 72),
  anio            = rep(curvas$anio, times = 72),
  pentada         = rep(1:72, each = nrow(curvas)),
  mm_dia_obs      = round(as.vector(Y), 4),
  mm_dia_imputado = round(as.vector(Y_FPCA), 4),       # observado + FPCA donde faltaba
  imputado        = as.vector(IMPUTADA),
  motivo_hueco    = as.vector(MOTIVO),
  curva_completa  = rep(completa_final, times = 72),
  curva_apta      = rep(pares$apto[match(clave_curva, paste(pares$STATION_ID, pares$anio))], times = 72),
  stringsAsFactors = FALSE)
pentadas$motivo_hueco[!is.na(pentadas$mm_dia_obs)] <- NA
pentadas <- pentadas[order(pentadas$STATION_ID, pentadas$anio, pentadas$pentada), ]
carpeta <- file.path(RUTA_SALIDAS, CARPETAS[["imputacion"]])
saveRDS(pentadas, file.path(carpeta, "estaciones_pentadas_imputadas.rds"))
write.csv(pentadas, gzfile(file.path(carpeta, "estaciones_pentadas_imputadas.csv.gz")), row.names = FALSE)
guardar_tabla(trazabilidad, "limpieza", "trazabilidad")

cat("\nTrazabilidad:\n")
print(trazabilidad, row.names = FALSE)
cat("\nListo. Resultados en:", normalizePath(RUTA_SALIDAS), "\n")
