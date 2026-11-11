# =============================================================================
#  LLUVIA EN EL VALLE DEL CAUCA: SCRIPT ÚNICO DE R DEL PROYECTO
#  Analítica de Datos 2026-2S, proyecto final
# =============================================================================
#
#  Qué hay en este archivo (se corre completo, de arriba a abajo):
#
#    PARTE 0. Configuración: parámetros, paquetes y funciones de apoyo
#    PARTE 1. Lectura de los datos (estaciones y satélite), una sola vez
#    PARTE 2. Pre-EDA: describe los datos tal como llegan, sin cambiar nada
#    PARTE 3. Preprocesamiento: corte de años, unidades, calendario,
#             rachas de ceros, cobertura y satélite 1981-2003
#    PARTE 4. Limpieza: coordenadas, rachas con vecinas, valores a NA,
#             años aptos para curvas y píxel de cada estación
#    PARTE 5. Guardar los datos limpios y la trazabilidad
#
#  Reglas que sigue todo el script:
#    - Los datos originales NUNCA se modifican. Lo que se quita queda como NA
#      en una columna aparte (VALUE_LIMPIO) y con su motivo.
#    - Todavía NO se imputa. Eso se hará después, con FPCA condicional.
#    - Las gráficas solo se muestran en pantalla (panel Plots de RStudio).
#      Las tablas sí se guardan como CSV en la carpeta "salidas".
#
#  Cómo correrlo: abrir el proyecto en RStudio (la carpeta del repositorio
#  debe ser el directorio de trabajo) y ejecutar todo (Ctrl + Shift + Enter).
#  Tarda unos 5 minutos. Paquetes: readxl (leer Excel) y terra (leer GeoTIFF).
# =============================================================================


# =============================================================================
# PARTE 0. CONFIGURACIÓN
# =============================================================================

# ---- 0.1 Paquetes ------------------------------------------------------------
# Solo se usan dos paquetes y únicamente para LEER archivos. Todo el cálculo
# se hace con R base.
for (paquete in c("readxl", "terra")) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop("Falta el paquete '", paquete, "'. Instálelo con install.packages('", paquete, "')")
  }
}

# ---- 0.2 Rutas -----------------------------------------------------------------
RUTA_DATOS   <- "Datos"
RUTA_SALIDAS <- "salidas"
CARPETAS <- c(pre_eda = "pre_eda", prep = "preprocesamiento", limpieza = "limpieza")
for (carpeta in CARPETAS) dir.create(file.path(RUTA_SALIDAS, carpeta), recursive = TRUE, showWarnings = FALSE)
options(scipen = 999)        # mostrar números sin notación científica

# ---- 0.3 Parámetros --------------------------------------------------------
# Todos los umbrales están aquí arriba para poder cambiarlos y justificarlos.

# Periodo
ANIO_INICIO     <- 1970      # antes de 1970 hay muy pocas estaciones (ver Parte 3.1)
ANIO_FIN        <- 2003      # último año de estaciones (indicado por el profesor)
PERIODO_COMUN   <- 1981:2003 # años en que hay estaciones Y satélite
ANIO_REFERENCIA <- 1997      # lo "normal" se calcula solo hasta aquí, para no usar
# años que podrían quedar de prueba (evita fuga de datos)

# Valores raros
VALOR_IMPOSIBLE <- 859375    # mm en un día: imposible, es un error de digitación

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

# Trazabilidad: cada paso que cambia los datos deja una fila con cuántas
# estaciones, pares estación-año y días quedan. Así se ve qué se perdió y dónde.
trazabilidad <- data.frame()
registrar <- function(paso, id, anio, valor) {
  con_dato <- !is.na(valor)
  fila <- data.frame(
    paso                = paso,
    estaciones          = length(unique(id[con_dato])),
    pares_estacion_anio = nrow(unique(data.frame(id, anio)[con_dato, ])),
    dias_con_dato       = sum(con_dato),
    dias_sin_dato       = sum(!con_dato))
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
valido <- !is.na(est$VALUE) & est$VALUE != VALOR_IMPOSIBLE
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
lluvioso <- !is.na(est$VALUE) & est$VALUE >= LLUVIA_RED & est$VALUE != VALOR_IMPOSIBLE
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
para_grafica <- !cal$sin_dato & cal$VALUE != VALOR_IMPOSIBLE
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
V[cbind(fila_est, col_dia)] <- ifelse(cal$VALUE == VALOR_IMPOSIBLE, NA, cal$VALUE)

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

# ---- 4.2 Estaciones con lecturas mensuales -----------------------------------
# Algunas estaciones no anotan cada día: escriben el total del mes en un solo
# día. Se detectan porque la mitad de sus meses tiene MAX_DIAS_MES días con
# dato o menos. No sirven para curvas diarias ni como vecinas.
dias_por_mes <- tapply(!cal$sin_dato, list(cal$STATION_ID, cal$anio * 100 + cal$mes), sum)
info$mediana_dias_por_mes <- apply(dias_por_mes[ids, , drop = FALSE], 1, function(x) median(x[!is.na(x) & x > 0]))
info$lecturas_mensuales <- info$mediana_dias_por_mes <= MAX_DIAS_MES
cat("Estaciones con lecturas mensuales:", paste(ids[info$lecturas_mensuales], collapse = ", "), "\n")
rm(dias_por_mes)

# ---- 4.3 Rachas de ceros con estaciones vecinas ------------------------------
# Misma pregunta que en 3.5, pero comparando con las estaciones cercanas en
# lugar de toda la zona. Vecinas de una estación: las que están a RADIO_KM o
# menos (hasta RADIO_MAX_KM si hay menos de MIN_VECINAS), con menos de
# DIF_ALTURA m de diferencia de altura y que no sean de lecturas mensuales.
cat("\n== 4.3 RACHAS DE CEROS CON VECINAS ==\n")
vecinas <- lapply(seq_along(ids), function(k) {
  if (!info$con_coordenadas[k]) return(integer(0))
  d <- dist_km(info$lon[k], info$lat[k], info$lon, info$lat)
  d[k] <- NA
  d[info$lecturas_mensuales] <- NA
  altura_ok <- is.na(info$altitud) | is.na(info$altitud[k]) |
    abs(info$altitud - info$altitud[k]) <= DIF_ALTURA
  v <- which(!is.na(d) & altura_ok & d <= RADIO_KM)
  if (length(v) < MIN_VECINAS) v <- which(!is.na(d) & altura_ok & d <= RADIO_MAX_KM)
  v
})
info$n_vecinas <- lengths(vecinas)
info$vecinas <- sapply(vecinas, function(v) paste(ids[v], collapse = " "))
cat("Vecinas por estación ubicada:\n")
print(summary(info$n_vecinas[info$con_coordenadas]))

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
rachas$metodo <- ifelse(!is.na(rachas$clasificacion_vecinas), "Vecinas",
                        ifelse(info$con_coordenadas[rachas$k], "Zona (sin vecinas con dato)", "Zona (sin coordenadas)"))
rachas$clasificacion <- ifelse(is.na(rachas$clasificacion_vecinas), rachas$clasificacion_zona,
                               rachas$clasificacion_vecinas)

cat("Método usado:\n"); print(table(rachas$metodo))
con_vecinas <- rachas$metodo == "Vecinas"
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
  s <- rachas[rachas$metodo == "Vecinas" & rachas$clasificacion == CLASES[3], ]
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

# ---- 4.4 Pasar a NA lo que no es lluvia real ---------------------------------
# VALUE se queda igual. VALUE_LIMPIO es la copia limpia y "motivo" dice por
# qué un día quedó sin dato. El orden importa: un día se marca una sola vez.
cat("\n== 4.4 VALORES QUE PASAN A NA ==\n")
cal$VALUE_LIMPIO <- cal$VALUE
cal$motivo <- ifelse(cal$sin_dato, "Sin dato original", NA_character_)
marcar <- function(filas, texto) {
  filas <- filas[!is.na(cal$VALUE_LIMPIO[filas])]
  cal$VALUE_LIMPIO[filas] <<- NA
  cal$motivo[filas] <<- texto
  cat(sprintf("%-40s %8d días\n", texto, length(filas)))
}

# a) Estaciones con lecturas mensuales (ver 4.2)
marcar(which(cal$STATION_ID %in% ids[info$lecturas_mensuales]), "Estación con lecturas mensuales")

# b) Valor imposible
marcar(which(cal$VALUE == VALOR_IMPOSIBLE), "Valor imposible (859375)")

# c) Acumulados: tras HUECO_MIN días o más sin dato llega un valor mayor que el
#    percentil 99 de la estación. Es la lluvia de varios días anotada en uno.
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

# d) Ceros falsos: los días de las rachas sospechosas
marcar(which(!is.na(cal$racha_ceros) & cal$racha_ceros == CLASES[3]), "Cero falso (racha sospechosa)")
registrar("4. Pasar a NA lo que no es lluvia real", cal$STATION_ID, cal$anio, cal$VALUE_LIMPIO)

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

# ---- 4.5 Pares estación-año aptos para curvas --------------------------------
# Una curva anual necesita el año casi completo. Un par estación-año es apto si:
#   - tiene COB_MIN % o más de días con dato (ya limpio),
#   - su hueco más largo no pasa de HUECO_MAX días,
#   - cada temporada de lluvias (mar-may y sep-nov) tiene COB_TEMPORADA % o más.
# No se borra nada: la columna "apto" servirá de filtro más adelante.
cat("\n== 4.5 PARES ESTACIÓN-AÑO APTOS ==\n")
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

# ---- 4.6 Píxel del satélite de cada estación ---------------------------------
# Cada estación ubicada se empareja con el píxel de 0,05° que la contiene.
# Estas parejas son las que se comparan en la sección 3 del proyecto y con las
# que la red KAN aprende a corregir el satélite.
cat("\n== 4.6 ESTACIÓN Y PÍXEL DEL SATÉLITE ==\n")
celda  <- terra::cellFromXY(plantilla_sat, cbind(info$lon, info$lat))
centro <- terra::xyFromCell(plantilla_sat, ifelse(is.na(celda), 1, celda))
info$pixel <- ifelse(is.na(celda), NA, sprintf("P%04d", celda))
info$km_al_centro_pixel <- ifelse(is.na(celda), NA,
                                  round(dist_km(info$lon, info$lat, centro[, 1], centro[, 2]), 2))
info$ubicacion <- ifelse(!info$con_coordenadas, "Sin coordenadas",
                         ifelse(is.na(celda), "Fuera de la imagen", "Dentro de la imagen"))
print(table(info$ubicacion))
cat(sprintf("Píxeles distintos con estación: %d | píxeles con 2 o más estaciones: %d\n",
            length(unique(na.omit(info$pixel))), sum(table(info$pixel) >= 2)))
cat(sprintf("Estaciones dentro de la imagen y con algún año apto en %d-%d: %d\n",
            min(PERIODO_COMUN), max(PERIODO_COMUN),
            length(unique(pares$STATION_ID[pares$apto & pares$anio %in% PERIODO_COMUN &
                                             pares$STATION_ID %in% ids[!is.na(info$pixel)]]))))

graficar(function() {   # Gráfica: mapa de estaciones y área del satélite
  e <- terra::ext(plantilla_sat)
  u <- info[info$con_coordenadas, ]
  op <- par(mar = c(4, 4, 4, 1)); on.exit(par(op))
  plot(u$lon, u$lat, type = "n", asp = 1, xlab = "Longitud", ylab = "Latitud",
       xlim = range(u$lon, e[1:2]), ylim = range(u$lat, e[3:4]),
       main = sprintf("Estaciones ubicadas (%d de %d) y área del satélite", nrow(u), nrow(info)))
  rect(e[1], e[3], e[2], e[4], border = COL_CERO, lwd = 2, lty = 2)
  apto <- !is.na(u$anios_aptos) & u$anios_aptos > 0
  points(u$lon, u$lat, pch = 19, cex = .9,
         col = ifelse(u$lecturas_mensuales, COL_ROJO, ifelse(apto, COL_PRINCIPAL, COL_GRIS)))
  legend("bottomleft", bty = "n", pch = c(19, 19, 19, NA), lty = c(NA, NA, NA, 2), lwd = c(NA, NA, NA, 2),
         col = c(COL_PRINCIPAL, COL_GRIS, COL_ROJO, COL_CERO), cex = .85,
         legend = c("Con algún año apto", "Sin años aptos", "Lecturas mensuales", "Imagen del satélite"))
})


# =============================================================================
# PARTE 5. GUARDAR
# =============================================================================
# Estos archivos son la entrada de los pasos siguientes (imputación, curvas).
cat("\n######## PARTE 5. GUARDAR ########\n")
guardar_tabla(info[, c("STATION_ID", "zona", "nombre", "municipio", "lon", "lat", "altitud", "ubicacion",
                       "pixel", "km_al_centro_pixel", "n_vecinas", "vecinas", "mediana_dias_por_mes",
                       "lecturas_mensuales", "anios_aptos")], "limpieza", "estaciones_info")

estaciones_limpias <- cal[, c("STATION_ID", "zona", "fecha", "anio", "mes", "dia", "VALUE", "VALUE_LIMPIO",
                              "motivo", "racha_ceros", "apto")]
carpeta <- file.path(RUTA_SALIDAS, CARPETAS[["limpieza"]])
saveRDS(estaciones_limpias, file.path(carpeta, "estaciones_limpias_1970_2003.rds"))
write.csv(estaciones_limpias, gzfile(file.path(carpeta, "estaciones_limpias_1970_2003.csv.gz")), row.names = FALSE)
saveRDS(satelite, file.path(carpeta, "satelite_1981_2003.rds"))
write.csv(satelite, gzfile(file.path(carpeta, "satelite_1981_2003.csv.gz")), row.names = FALSE)
guardar_tabla(trazabilidad, "limpieza", "trazabilidad")

cat("\nTrazabilidad:\n")
print(trazabilidad, row.names = FALSE)
cat("\nListo. Resultados en:", normalizePath(RUTA_SALIDAS), "\n")
