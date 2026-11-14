# =============================================================================
# Preprocesamiento de la lluvia del Valle del Cauca (estaciones CVC diarias + satélite por péntadas)
# Proyecto final de Analítica: los datos limpios alimentan los modelos KAN y de redes funcionales.
#
# Partición del proyecto: entrenamiento 1981-1995, validación 1996-2000, prueba 2001-2003.
# Todo lo que se CALIBRA (umbrales, correlaciones, vecinas, variograma, métricas de imputación) usa
# solo 1981-1995, para no filtrar información de validación ni de prueba.
#
# Pasos (cada uno es una sección del script) y por qué se hacen:
#   1. Carga: se unen los tres archivos de estaciones (1930-2003), el catálogo y la máscara del Valle.
#   2. Ubicación y altitud: hay estaciones sin coordenadas; se ubican con la correlación de su lluvia
#      con estaciones conocidas (error típico ~5 km, un píxel del satélite). Sin fuentes externas.
#   3. Solo el Valle: el proyecto es el Valle del Cauca; las de fuera solo sirven para ubicar.
#   4. Verificación de formato: confirmar que las estaciones están en mm/día y el satélite en
#      mm/péntada antes de comparar nada.
#   5. Completitud 1930-2003: ver qué años tienen datos y justificar el corte (G1, G2).
#   6. Corte 1981-2003: es el periodo del satélite y de la partición.
#   7. Satélite: solo píxeles del Valle; cada estación se empareja con su píxel.
#   8. Limpieza 1 - registros acumulados: algunas estaciones reportaron cada semana o cada mes; ese
#      valor es real y se reparte en días con la forma de la vecina más correlacionada.
#   9. Limpieza 2 - días de más de 300 mm y totales mensuales en un día: se quitan si ninguna vecina
#      los confirma (en la costa Pacífica pueden ser reales, por eso no hay tope fijo).
#  10. Limpieza 3 - ceros falsos: rachas de 15+ días en 0 que las vecinas (o el satélite) desmienten.
#      Las sequías reales se dejan.
#  11. Limpieza 4 - meses inconsistentes: errores sostenidos que no tienen ningún día imposible
#      pero sí un total mensual imposible frente a las vecinas. Además se quita un año sin vecinas
#      (5301100101 en 2001) que el satélite y la propia estación muestran como error.
#  12. Completitud 1981-2003 con los datos ya limpios.
#  13. Ranking de estaciones por cobertura.
#  14. Sensibilidad del criterio de cobertura: cuántos datos quedan con otros umbrales.
#  15. Criterio de cobertura: un año entra si tiene >= 90 % de días y ningún hueco > 31 días
#      (la curva anual necesita el año completo y sus dos picos de lluvia).
#  16. Imputación: se comparan IDW, Kriging ordinario, SBE y CCW escondiendo días conocidos;
#      se usa el de mayor KGE (métrica del proyecto).
#  17. Datos limpios en CSV: diario y por péntadas (misma escala que el satélite) y trazabilidad.
#  18. Estadísticas descriptivas.
#  19. Gráficas (solo se imprimen).
#
# Salidas: tablas en Excel (.xlsx); los datos limpios en CSV porque no caben en Excel.
#          Las gráficas se imprimen en pantalla y no se guardan.
#
# Uso (desde la carpeta del proyecto, la que contiene ./Datos):
#   source("preprocesamiento.R")                         # en RStudio
#   Rscript preprocesamiento.R [RUTA_DATOS] [RUTA_SALIDA]
#
# Paquetes: dplyr, tidyr, lubridate, ggplot2, readxl, writexl, terra
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate); library(ggplot2)
  library(readxl); library(writexl); library(terra)
})

args <- commandArgs(trailingOnly = TRUE)
D   <- if (length(args) >= 1) args[1] else "./Datos"
OUT <- if (length(args) >= 2) args[2] else "./salidas_preprocesamiento"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S0 <- 1981; S1 <- 2003     # años del satélite que se usan (periodo del proyecto)
# Años de entrenamiento de la partición del proyecto (val 1996-2000, prueba 2001-2003).
# Todo lo que se CALIBRA en la limpieza (umbrales, correlaciones, elección de vecinas) usa solo
# estos años, para que el preprocesamiento no use información de validación ni de prueba.
ENTRENA <- 1981:1995

C_MAIN <- "#2a6f97"; C_ALT <- "#e07a5f"; C_BAD <- "#b23a48"
theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()))
guardar_tabla <- function(df, nombre) write_xlsx(as.data.frame(df), file.path(OUT, paste0(nombre, ".xlsx")))
mostrar <- function(p) print(p)

# Trazabilidad: "foto" de los datos de estaciones después de cada paso (tabla antes/después)
traz <- list()
foto <- function(st, paso, cal = NULL) {
  # Faltantes: días sin dato dentro de la vida útil de cada estación (de su primer a su último día
  # con dato). Tras el criterio de cobertura se cuentan con el calendario de los años que quedan.
  vida <- st |> group_by(STATION_ID) |>
    summarise(dias_vida = as.numeric(max(FECHA) - min(FECHA)) + 1, n = n(), .groups = "drop")
  traz[[length(traz) + 1]] <<- tibble(
    paso = paso, estaciones = n_distinct(st$STATION_ID),
    pares_estacion_anio = nrow(distinct(st, STATION_ID, YEAR)),
    dias_con_dato = nrow(st),
    dias_faltantes = if (is.null(cal)) sum(vida$dias_vida - vida$n) else sum(!cal$presente),
    pct_dias_en_0 = round(100 * mean(st$VALUE == 0), 1), max_mm = max(st$VALUE),
    desde = min(st$YEAR), hasta = max(st$YEAR))
}

# ------------------------------------------------------------------ 1. Carga
# Los tres archivos de estaciones se unen en una sola tabla (1930-2003)
# (el histórico trae la columna UNITS; los CSV no, por eso se guarda de qué archivo viene cada fila)
hist <- read_excel(file.path(D, "P_1930_1979.xlsx")) |>
  select(STATION_ID, YEAR, MONTH, DAY, VALUE, UNITS)
st <- bind_rows(hist |> mutate(archivo = "1930-1979"),
                read.csv(file.path(D, "P_1980_1990.csv")) |> mutate(archivo = "1980-1990"),
                read.csv(file.path(D, "P_1991_2003.csv")) |> mutate(archivo = "1991-2003")) |>
  mutate(STATION_ID = as.numeric(STATION_ID), VALUE_txt = as.character(VALUE),
         VALUE = suppressWarnings(as.numeric(VALUE)),
         FECHA = suppressWarnings(make_date(YEAR, MONTH, DAY)))
rm(hist)

# El catálogo viene en latin-1 (tildes); se lee como bytes latin1 para no cortar filas
ll  <- read.csv(file.path(D, "CVC_Wheather_Stations_Long_Lat.csv"), encoding = "latin1") |>
  mutate(CODIGO = as.numeric(CODIGO))
stopifnot(nrow(ll) == 238)
# Catálogo de altitud: trae altitud, cuenca y municipio (algunas estaciones no están en el catálogo de coordenadas)
cat_alt <- read_excel(file.path(D, "CatalogoEstacionesElevacion.xlsx")) |>
  mutate(codigo = as.numeric(codigo)) |> distinct(codigo, .keep_all = TRUE)

# Máscara del Valle del Cauca: la grilla de altitud (misma malla del satélite) solo tiene valores
# dentro del departamento; fuera del Valle es NA.
dem   <- rast(file.path(D, "altitud_valle.tif"))                     # altitud de la grilla
valle <- ifel(is.na(dem), NA, 1)                                     # máscara: 1 = Valle, NA = fuera
borde <- as.data.frame(geom(as.lines(as.polygons(valle))))           # borde del Valle (para los mapas)

Y0 <- min(st$YEAR); Y1 <- max(st$YEAR)
estaciones <- sort(unique(st$STATION_ID))

# ------------------------------------------------- 2. Ubicación y altitud de TODAS las estaciones
# Solo con los datos disponibles (sin fuentes externas). Grupos según la información de cada estación:
#   A: coordenadas + altitud de catálogo     -> se dejan igual
#   B: coordenadas sin altitud               -> altitud de la grilla (interpolada y corregida)
#   C: altitud de catálogo sin coordenadas   -> ubicación por lluvia, ajustada con su altitud
#   D: sin coordenadas ni altitud            -> ubicación por lluvia, altitud de la grilla
# La ubicación por lluvia ("triangulación"): dos estaciones cercanas registran la lluvia de los
# mismos días. Se buscan las 3 estaciones con coordenadas cuya lluvia diaria más se correlaciona
# y se promedia su posición con peso = correlación^2. Validado escondiendo las coordenadas de las
# estaciones conocidas: error típico ~5 km (aprox. un píxel del satélite).

conocidas   <- sort(intersect(estaciones, ll$CODIGO))
desconocidas <- setdiff(estaciones, conocidas)

# Matriz diaria (filas = días, columnas = estaciones). Se excluyen valores >= 500 mm (errores)
lluvia_ancha <- st |> filter(VALUE < 500) |> select(FECHA, STATION_ID, VALUE) |>
  pivot_wider(names_from = STATION_ID, values_from = VALUE) |> arrange(FECHA)
M <- as.matrix(lluvia_ancha[, -1])
Mu <- M[, as.character(desconocidas), drop = FALSE]
Mk <- M[, as.character(conocidas), drop = FALSE]
comun <- crossprod(!is.na(Mu), !is.na(Mk))               # días en común de cada par
correl <- suppressWarnings(cor(Mu, Mk, use = "pairwise.complete.obs"))
correl[comun < 365] <- NA                                # se exige al menos un año en común

coord_k <- ll |> filter(CODIGO %in% conocidas) |> select(CODIGO, Lon, lat)
coord_k <- coord_k[match(conocidas, coord_k$CODIGO), ]
km_entre <- function(lo1, la1, lo2, la2) sqrt(((lo1 - lo2) * 111 * cos(3.8 * pi / 180))^2 + ((la1 - la2) * 111)^2)

# Respaldo por ID: centro de las estaciones conocidas con los mismos 5 primeros dígitos
grupo_id <- function(id) {
  g <- coord_k |> filter(substr(CODIGO, 1, 5) == substr(id, 1, 5))
  if (nrow(g) == 0) c(NA, NA) else c(mean(g$Lon), mean(g$lat))
}

# Celdas de la grilla (para el ajuste por altitud del grupo C)
celdas <- as.data.frame(dem, xy = TRUE) |> rename(alt_celda = 3)

estim <- lapply(desconocidas, function(id) {
  r <- correl[as.character(id), ]
  top <- sort(r[!is.na(r)], decreasing = TRUE)[1:3]
  top <- top[!is.na(top)]
  gid <- grupo_id(id)
  if (length(top) == 0) {                                   # sin datos en común: solo el ID
    return(tibble(STATION_ID = id, Lon_est = gid[1], lat_est = gid[2], corr_max = NA_real_,
                  fuente_coord = ifelse(is.na(gid[1]), "sin ubicación", "ID")))
  }
  w <- pmax(top, 0.01)^2
  k <- match(as.numeric(names(top)), coord_k$CODIGO)
  lo <- sum(coord_k$Lon[k] * w) / sum(w); la <- sum(coord_k$lat[k] * w) / sum(w)
  fuente <- "lluvia"
  a <- cat_alt$alt[match(id, cat_alt$codigo)]
  if (!is.na(a)) {                                          # grupo C: ajuste por altitud
    d <- km_entre(celdas$x, celdas$y, lo, la)
    cand <- celdas[d <= 10, ]; dc <- d[d <= 10]
    if (nrow(cand) > 0) {
      j <- which.min(abs(cand$alt_celda - a) / 200 + dc / 10)
      lo <- (lo + cand$x[j]) / 2; la <- (la + cand$y[j]) / 2
      fuente <- "lluvia + altitud"
    }
  }
  tibble(STATION_ID = id, Lon_est = lo, lat_est = la, corr_max = unname(top[1]), fuente_coord = fuente,
         dist_grupo_id_km = ifelse(is.na(gid[1]), NA, km_entre(lo, la, gid[1], gid[2])))
}) |> bind_rows()

# Confianza: alta si la correlación es >= 0,5 y la posición no cae a más de 20 km de su grupo de ID
estim <- estim |> mutate(confianza = ifelse(!is.na(corr_max) & corr_max >= 0.5 &
                                              (is.na(dist_grupo_id_km) | dist_grupo_id_km <= 20), "alta", "baja"))

ubic <- tibble(STATION_ID = estaciones) |>
  left_join(ll |> select(STATION_ID = CODIGO, Lon, lat), by = "STATION_ID") |>
  left_join(estim, by = "STATION_ID") |>
  mutate(fuente_coord = ifelse(!is.na(Lon), "catálogo", fuente_coord),
         confianza = ifelse(!is.na(Lon), "catálogo", confianza),
         Lon = coalesce(Lon, Lon_est), lat = coalesce(lat, lat_est)) |>
  select(STATION_ID, Lon, lat, fuente_coord, corr_max, confianza)

# Altitud: catálogo si existe; si no, grilla interpolada (bilineal) y corregida con una regresión
# ajustada en las estaciones que tienen catálogo y coordenadas (la grilla sobreestima en montaña)
xy_u <- as.matrix(ubic[, c("Lon", "lat")]); xy_u[is.na(xy_u)] <- 0
alt_bil <- terra::extract(dem, xy_u, method = "bilinear")[, 1]
alt_cel <- terra::extract(dem, xy_u)[, 1]
ubic <- ubic |> mutate(alt_grilla = ifelse(is.na(Lon), NA, coalesce(alt_bil, alt_cel)),
                       alt_catalogo = cat_alt$alt[match(STATION_ID, cat_alt$codigo)])
ajuste <- lm(alt_catalogo ~ alt_grilla, data = ubic |> filter(fuente_coord == "catálogo"))
ubic <- ubic |> mutate(alt_corregida = predict(ajuste, newdata = ubic),
                       altitud = coalesce(alt_catalogo, alt_corregida),
                       fuente_alt = case_when(!is.na(alt_catalogo) ~ "catálogo",
                                              !is.na(alt_corregida) ~ "grilla corregida",
                                              TRUE ~ "sin dato"),
                       en_valle = !is.na(terra::extract(valle, xy_u)[, 1]) & !is.na(Lon)) |>
  select(STATION_ID, Lon, lat, fuente_coord, corr_max, confianza, altitud, fuente_alt, en_valle)

print(count(ubic, fuente_coord, confianza))
print(count(ubic, fuente_alt))

# ------------------------------------------------- 3. Solo estaciones dentro del Valle del Cauca
# Se usan todas las estaciones con coordenadas como referencia para ubicar (paso 2); desde aquí
# solo se analizan las que quedan dentro del Valle.
ubic <- ubic |> left_join(ll |> select(STATION_ID = CODIGO, STATION_NAME), by = "STATION_ID")
base::cat("\nEstaciones con datos:", nrow(ubic),
          "| dentro del Valle:", sum(ubic$en_valle),
          "| fuera del Valle:", sum(!ubic$en_valle & !is.na(ubic$Lon)),
          "| sin ubicación:", sum(is.na(ubic$Lon)), "\n")
ubic <- ubic |> filter(en_valle)
# Tabla de verificación: todas las estaciones del Valle con su ubicación y de dónde salió
print(ubic |> select(STATION_ID, STATION_NAME, Lon, lat, fuente_coord, confianza, altitud, fuente_alt),
      n = Inf, width = Inf)
base::cat("Estaciones del Valle sin coordenadas:", sum(is.na(ubic$Lon)),
          "| sin altitud:", sum(is.na(ubic$altitud)), "\n")

foto(st, "0. Datos crudos (todas las estaciones)")
st <- st |> filter(STATION_ID %in% ubic$STATION_ID)
foto(st, "1. Solo estaciones dentro del Valle")
estaciones <- sort(unique(st$STATION_ID))

# ------------------------------------------------- 4. Verificación de formato: mm/día y mm/péntada
# Estaciones: lluvia en mm, un dato por día. Satélite: lluvia en mm, un archivo por péntada.
# Solo se cuenta cuántas estaciones / archivos NO cumplen cada condición.
fmt_est <- st |> group_by(STATION_ID) |>
  summarise(
    unidad_no_mm  = any(archivo == "1930-1979" & (is.na(UNITS) | UNITS != "mm")),
    tiene_hist    = any(archivo == "1930-1979"),
    fecha_mala    = any(is.na(FECHA)),
    dia_repetido  = any(duplicated(FECHA)),
    valor_malo    = any(is.na(VALUE) | VALUE < 0),
    .groups = "drop")
# Los CSV no traen unidad: se revisa que la escala no cambie entre archivos (si una estación
# pasara a décimas de mm o a cm, su total anual típico cambiaría ~10 veces). Se usan años con >= 330 días.
escala <- st |> group_by(STATION_ID, archivo, YEAR) |>
  summarise(total = sum(VALUE), n = n(), .groups = "drop") |> filter(n >= 330) |>
  group_by(STATION_ID, archivo) |> summarise(tipico = median(total), .groups = "drop") |>
  arrange(STATION_ID, archivo) |> group_by(STATION_ID) |>
  summarise(comparables = n() > 1,
            cambio_escala = any(tipico / lag(tipico) > 3 | tipico / lag(tipico) < 1/3, na.rm = TRUE))
# Dato diario: si una estación tiene 6 o más registros SEGUIDOS separados por 2 o más días
# (semanal, mensual, cada pocos días), en ese periodo puede que no haya reportado día a día.
# Aquí solo se cuenta; en la limpieza 1 se decide con las vecinas si son acumulados o días sueltos.
MIN_REG_RACHA <- 6
rachas_separadas <- function(st) {
  st |> arrange(STATION_ID, FECHA) |> group_by(STATION_ID) |>
    mutate(gap = as.numeric(FECHA - lag(FECHA)), sep = !is.na(gap) & gap >= 2,
           racha = cumsum(sep != lag(sep, default = FALSE))) |>
    group_by(STATION_ID, racha) |>
    mutate(n_racha = sum(sep), en_racha = sep & n_racha >= MIN_REG_RACHA) |> ungroup()
}
diaria <- rachas_separadas(st) |> group_by(STATION_ID) |> summarise(no_diaria = any(en_racha))

# ------------------------------------------------- 5. Completitud estación-año (1930-2003)
# Un día faltante = una fila que NO existe (no hay NA en los archivos).
# Por eso construimos el calendario completo y marcamos qué días tienen dato.

# Racha más larga de FALSE (días sin dato) en un vector lógico
hueco_max <- function(presente) {
  r <- rle(presente)
  if (any(!r$values)) max(r$lengths[!r$values]) else 0L
}

# Vida útil de cada estación: de su primer a su último día con dato.
# Antes y después de ese rango la estación no existía o ya se había retirado, así que esos días
# NO se cuentan como faltantes (serían nulos falsos).
# Calendario de cada estación solo dentro de su vida útil: una fila por estación y por día.
# Los días que no venían en los archivos quedan con VALUE = NA: estos son los faltantes reales.
# Cobertura estación-año (ea). Los años fuera de la vida útil quedan con 0 días (para el mapa de calor).
# Se escribe como función porque se calcula dos veces: 1930-2003 y, después del corte y la limpieza, 1981-2003.
completitud <- function(st, y0, y1) {
  dias_anio <- tibble(FECHA = seq(make_date(y0, 1, 1), make_date(y1, 12, 31), by = "day")) |>
    count(YEAR = year(FECHA), name = "dias")
  vida <- st |> group_by(STATION_ID) |> summarise(ini = min(FECHA), fin = max(FECHA))
  cal <- vida |> rowwise() |> reframe(STATION_ID, FECHA = seq(ini, fin, by = "day")) |>
    left_join(st |> select(STATION_ID, FECHA, VALUE), by = c("STATION_ID", "FECHA")) |>
    mutate(YEAR = year(FECHA), presente = !is.na(VALUE))
  ea <- cal |>
    group_by(STATION_ID, YEAR) |>
    summarise(n_dias = sum(presente), hueco = hueco_max(presente), .groups = "drop") |>
    complete(STATION_ID = sort(unique(st$STATION_ID)), YEAR = y0:y1, fill = list(n_dias = 0L)) |>
    left_join(dias_anio, by = "YEAR") |>
    mutate(hueco = coalesce(hueco, dias), cob = 100 * n_dias / dias)
  list(cal = cal, ea = ea)
}
ea <- completitud(st, Y0, Y1)$ea

# Estas dos gráficas se imprimen aquí porque son las únicas con todo 1930-2003.
# Orden de las estaciones: de más a menos años "buenos" (cobertura >= 90 %)
orden <- ea |> group_by(STATION_ID) |> summarise(a90 = sum(cob >= 90), cob = mean(cob)) |>
  arrange(desc(a90), desc(cob)) |> pull(STATION_ID)
eje_anios <- scale_x_continuous(breaks = seq(1930, Y1, 5))

# G1. Mapa de calor estación x año (1930-2003)
# Cada fila es una estación y cada columna un año. Arriba quedan las estaciones con más años
# completos y abajo las que menos tienen. Azul = año bien cubierto, rojo = año con muchos faltantes.
ea_plot <- ea |> mutate(est = factor(STATION_ID, levels = rev(orden)),
                        clase = cut(cob, c(-Inf, 0, 50, 80, 90, 99.999, Inf),
                                    labels = c("sin datos", "<50%", "50-80%", "80-90%", "90-<100%", "100%")))
mostrar(
  ggplot(ea_plot, aes(YEAR, est, fill = clase)) + geom_tile() +
    scale_fill_manual(values = c("#e8e8e8", "#b2182b", "#ef8a62", "#92c5de", "#4393c3", "#2166ac"),
                      name = "% de días\ncon dato", drop = FALSE) +
    geom_vline(xintercept = S0 - .5, linetype = "dashed") +
    annotate("text", x = S0 + 11, y = length(orden) + 8, label = "periodo con satélite 1981-2003", size = 3) +
    scale_x_continuous(breaks = seq(1930, Y1, 5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = "Cobertura diaria por estación y año",
         x = NULL, y = paste(length(orden), "estaciones")) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(), plot.margin = margin(15, 5, 5, 5))
)

# G2. Estaciones que reportan cada año
f2 <- ea |> group_by(YEAR) |>
  summarise(`con algún dato` = sum(cob > 0), `cobertura >= 80 %` = sum(cob >= 80),
            `cobertura >= 90 %` = sum(cob >= 90), `100 %` = sum(cob == 100)) |>
  pivot_longer(-YEAR, names_to = "serie", values_to = "n") |>
  mutate(serie = factor(serie, levels = c("con algún dato", "cobertura >= 80 %", "cobertura >= 90 %", "100 %")))
mostrar(
  ggplot(f2, aes(YEAR, n, colour = serie, linetype = serie)) + geom_line(linewidth = .8) +
    geom_vline(xintercept = S0 - .5, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("#000000", "#E69F00", "#009E73", "#CC79A7"), name = NULL) +
    scale_linetype_manual(values = c("solid", "solid", "longdash", "solid"), name = NULL) +
    eje_anios + expand_limits(y = 0) +
    labs(title = "¿Cuántas estaciones reportan cada año?", x = NULL, y = "n.º de estaciones") +
    theme(legend.position = "bottom")
)

# ------------------------------------------------- 6. Corte: desde aquí solo 1981-2003
# La partición del proyecto empieza en 1981 (train 1981-1995 / val 1996-2000 / test 2001-2003).
n_antes <- n_distinct(st$STATION_ID)
st <- st |> filter(YEAR >= S0, YEAR <= S1)
Y0 <- S0; Y1 <- S1
estaciones <- sort(unique(st$STATION_ID))
ubic <- ubic |> filter(STATION_ID %in% estaciones)
foto(st, "2. Corte 1981-2003")
base::cat("\nCorte 1981-2003: estaciones del Valle", n_antes, "->", length(estaciones),
          "(", n_antes - length(estaciones), "dejaron de medir antes de 1981 )\n")

# ------------------------------------------------- 7. Satélite (1981-2003), solo píxeles del Valle
fs  <- list.files(file.path(D, "ImagenesSatelitePentadal"), pattern = "^ppt\\d{6}\\.tif$", full.names = TRUE)
inf <- tibble(f = fs, anio = as.integer(substr(basename(fs), 4, 7)),
              pentada = as.integer(substr(basename(fs), 8, 9))) |>
  filter(anio >= S0, anio <= S1) |> arrange(anio, pentada)
sat_r <- rast(inf$f)                       # una capa por péntada; el NA (-9999) lo lee terra
cubo  <- as.array(sat_r)                   # filas x columnas x péntadas
fuera <- is.na(as.array(valle)[, , 1])     # píxeles fuera del Valle
cubo[rep(fuera, dim(cubo)[3])] <- NA       # se quitan del análisis

validos_px <- apply(!is.na(cubo), c(1, 2), mean) * 100
vmax <- max(cubo, na.rm = TRUE)

# Emparejamiento estación -> píxel (el píxel que contiene a cada estación)
celda <- cellFromXY(sat_r, as.matrix(ubic[, c("Lon", "lat")]))
mc <- ubic |> mutate(fila = rowFromCell(sat_r, celda), col = colFromCell(sat_r, celda))
mc$pixel_valido <- validos_px[cbind(mc$fila, mc$col)] == 100
mc <- mc |> group_by(fila, col) |> mutate(comparte_pixel = n() > 1) |> ungroup()
guardar_tabla(mc |> select(STATION_ID, STATION_NAME, Lon, lat, fuente_coord, fila, col, pixel_valido, comparte_pixel),
              "estacion_pixel")

# Verificación de formato del satélite y tabla de formato (estaciones de la sección 4 + satélite)
fmt_sat <- tibble(f = fs) |>
  mutate(anio = as.integer(substr(basename(f), 4, 7)), pentada = as.integer(substr(basename(f), 8, 9))) |>
  filter(anio >= S0, anio <= S1) |>
  rowwise() |>
  mutate(r = list(rast(f)), v = list(values(r)),
         grilla_distinta = !(nrow(r) == nrow(sat_r) && ncol(r) == ncol(sat_r) &&
                               isTRUE(all.equal(res(r), res(sat_r))) && ext(r) == ext(sat_r) && same.crs(r, sat_r)),
         varias_capas = nlyr(r) != 1,
         valor_malo   = any(v < 0, na.rm = TRUE) | any(v > 200, na.rm = TRUE)) |>
  ungroup()
pent_faltan <- sum(72 - table(factor(fmt_sat$anio, levels = S0:S1)))

formato <- tibble(
  fuente = c(rep("Estaciones", 6), rep("Satélite", 4)),
  verificacion = c(
    "Unidad = mm (columna UNITS del archivo 1930-1979)",
    "Misma escala entre archivos (los CSV no traen unidad)",
    "Fecha válida",
    "Un solo dato por día (sin días repetidos)",
    "Valor numérico y no negativo",
    "Dato diario, 1930-2003 (sin periodos de 6+ registros seguidos separados por 2+ días; en 1981-2003 se corrigen en la limpieza 1)",
    "Una péntada por archivo: 72 por año, 1981-2003 (péntadas que faltan)",
    "Misma grilla: 39 x 37 píxeles, 0,05°, WGS84",
    "Una sola capa por archivo",
    "Valor en mm por péntada: numérico, no negativo y <= 200"),
  revisados = c(sum(fmt_est$tiene_hist), sum(escala$comparables), rep(nrow(fmt_est), 3), nrow(diaria),
                S1 - S0 + 1, rep(nrow(fmt_sat), 3)),
  no_cumplen = c(sum(fmt_est$unidad_no_mm), sum(escala$cambio_escala), sum(fmt_est$fecha_mala),
                 sum(fmt_est$dia_repetido), sum(fmt_est$valor_malo), sum(diaria$no_diaria),
                 pent_faltan, sum(fmt_sat$grilla_distinta), sum(fmt_sat$varias_capas), sum(fmt_sat$valor_malo)))
print(formato, width = Inf)
guardar_tabla(formato, "verificacion_formato")
rm(fmt_sat)

# ------------------------------------------------- 8. Limpieza 1: registros acumulados (semanales, mensuales...)
# Algunas estaciones, en ciertos periodos, no reportaron día a día sino un valor cada semana, cada
# mes o cada pocos días. Si ese valor es el ACUMULADO desde el registro anterior, es información
# real: no se bota, se reparte en días (desagregación) con el patrón de la vecina más parecida.
#   - Periodo sospechoso: 6 o más registros seguidos separados por 2 o más días (ver sección 4).
#     Ventana de cada registro: los días desde el registro anterior hasta él.
#   - Donante: estación a <= 25 km, con >= 90 % de los días del periodo y la mayor correlación de
#     lluvia diaria con la estación (calculada con 1981-1995, sin los periodos sospechosos y con
#     >= 1 año de días en común; si ninguna cumple, la más cercana).
#   - ¿Acumulado o días sueltos? Se suma lo que reportó la estación (V) y se compara con lo que
#     llovió en la donante en las ventanas completas (Sw) y solo en los días con registro (S1).
#     Es acumulado si V se parece más a Sw que a S1 y además V >= 1,5 x S1. Si no, son días sueltos
#     con faltantes en medio y se dejan como están.
#     La decisión se toma POR ESTACIÓN (sumando todos sus periodos), no periodo por periodo:
#     en una prueba con 500 casos simulados (solo con años 1981-1995) de acumulados y de días sueltos
#     conocidos, acertó el 96 % por estación contra el 89 % por periodo (con V >= 3 x S1: 89 % y 83 %).
#   - Periodos mensuales que no pasan como acumulados (V muy por debajo de Sw, casi todo en 0) no
#     traen información diaria: se quitan.
#   - Día d de la ventana: valor x (lluvia de la donante el día d / lluvia de la donante en la ventana).
#     Se conserva el total y la forma real de la lluvia (días secos y aguaceros). No se promedian
#     varias vecinas porque eso suaviza: llovería poco todos los días.
#   - Si la donante no tuvo lluvia en la ventana y el valor es > 0, esa ventana queda faltante.
#   - Los días generados quedan marcados (origen = "desagregado") para no usarlos como verdad
#     al evaluar modelos.
RADIO_KM <- 25

matriz_diaria <- function(st) {
  st |> select(FECHA, STATION_ID, VALUE) |>
    pivot_wider(names_from = STATION_ID, values_from = VALUE) |> arrange(FECHA) |>
    complete(FECHA = seq(make_date(S0, 1, 1), make_date(S1, 12, 31), by = "day"))
}
distancias <- function(ids) {                                # km entre todas las estaciones
  xy <- ubic[match(ids, ubic$STATION_ID), c("Lon", "lat")]
  d <- outer(seq_along(ids), seq_along(ids),
             function(a, b) km_entre(xy$Lon[a], xy$lat[a], xy$Lon[b], xy$lat[b]))
  diag(d) <- Inf; d
}

st <- st |> mutate(origen = "medido", donante = NA_real_)
rs <- rachas_separadas(st)
sosp <- rs |> filter(en_racha) |>
  group_by(STATION_ID, racha) |>
  mutate(gap = ifelse(is.na(gap) | gap > 62, median(gap, na.rm = TRUE), gap)) |>   # 1er registro tras el corte
  ungroup() |> select(STATION_ID, racha, FECHA, VALUE, gap)
st_d <- st |> anti_join(sosp, by = c("STATION_ID", "FECHA"))            # sin los periodos sospechosos
ancha <- matriz_diaria(st_d)
M <- as.matrix(ancha[, -1]); ids_m <- as.numeric(colnames(M)); dist_m <- distancias(ids_m)
tr <- year(ancha$FECHA) %in% ENTRENA                         # filas de entrenamiento

# Paso 1: para cada periodo, su donante y las sumas V, Sw y S1
per <- list()
for (g in split(sosp, list(sosp$STATION_ID, sosp$racha), drop = TRUE)) {
  j <- match(g$STATION_ID[1], ids_m)
  fin <- match(g$FECHA, ancha$FECHA); ini <- pmax(fin - g$gap + 1, 1)
  vent <- lapply(seq_along(fin), function(k) ini[k]:fin[k])
  todas <- unique(unlist(vent))
  cand <- which(dist_m[j, ] <= RADIO_KM & colMeans(!is.na(M[todas, , drop = FALSE])) >= .9)
  r <- if (length(cand)) suppressWarnings(cor(M[tr, j], M[tr, cand, drop = FALSE], use = "pairwise.complete.obs")[1, ]) else numeric(0)
  comun <- if (length(cand)) colSums(!is.na(M[tr, cand, drop = FALSE]) & !is.na(M[tr, j])) else numeric(0)
  r[is.na(r) | comun < 365] <- -Inf                         # correlación solo con >= 1 año de días en común
  don <- if (!length(cand)) NA else if (any(is.finite(r))) cand[which.max(r)] else cand[which.min(dist_m[j, cand])]
  per[[length(per) + 1]] <- list(
    g = g, fin = fin, vent = vent, don = don,
    fila = tibble(
      STATION_ID = g$STATION_ID[1], desde = min(g$FECHA), hasta = max(g$FECHA), registros = nrow(g),
      dias_entre_registros = median(g$gap), V_estacion = sum(g$VALUE),
      donante = if (is.na(don)) NA else ids_m[don],
      corr_donante = if (is.na(don) || !is.finite(max(r))) NA else round(max(r), 2),
      km_donante = if (is.na(don)) NA else round(dist_m[j, don], 1),
      Sw_donante_ventanas = if (is.na(don)) NA else sum(sapply(vent, function(w) sum(M[w, don], na.rm = TRUE))),
      S1_donante_mismos_dias = if (is.na(don)) NA else sum(M[fin, don], na.rm = TRUE)))
}
traza <- bind_rows(lapply(per, `[[`, "fila"))

# Paso 2: decisión por estación (sumando todos sus periodos con donante)
es_acumulado <- function(V, Sw, S1) abs(log((V + 1) / (Sw + 1))) < abs(log((V + 1) / (S1 + 1))) & V >= 1.5 * S1
dec_est <- traza |> filter(!is.na(donante)) |> group_by(STATION_ID) |>
  summarise(acum_estacion = es_acumulado(sum(V_estacion), sum(Sw_donante_ventanas), sum(S1_donante_mismos_dias)))
traza <- traza |> left_join(dec_est, by = "STATION_ID") |>
  mutate(acum_estacion = coalesce(acum_estacion, FALSE),
         resultado = case_when(
           acum_estacion & !is.na(donante) ~ "acumulado: desagregado",
           !acum_estacion & !is.na(donante) & dias_entre_registros >= 28 &
             V_estacion < Sw_donante_ventanas / 3 ~ "mensual sin información diaria: se quita",
           TRUE ~ "días sueltos: se dejan"))

# Paso 3: desagregación de los periodos acumulados
desag <- list(); traza$registros_desagregados <- 0
for (i in seq_along(per)) {
  if (traza$resultado[i] != "acumulado: desagregado") next
  g <- per[[i]]$g; don <- per[[i]]$don; vent <- per[[i]]$vent
  for (k in seq_along(vent)) {
    v <- M[vent[[k]], don]; tv <- sum(v, na.rm = TRUE)
    if (tv > 0 || g$VALUE[k] == 0) {
      traza$registros_desagregados[i] <- traza$registros_desagregados[i] + 1
      desag[[length(desag) + 1]] <- tibble(STATION_ID = g$STATION_ID[1], FECHA = ancha$FECHA[vent[[k]]],
                                           VALUE = round(if (tv > 0) g$VALUE[k] * v / tv else 0 * v, 1)) |>
        filter(!is.na(VALUE)) |>
        mutate(YEAR = year(FECHA), MONTH = month(FECHA), DAY = day(FECHA),
               origen = "desagregado", donante = ids_m[don])
    }
  }
}
print(traza |> count(resultado, wt = registros, name = "registros"))
guardar_tabla(traza, "registros_acumulados")
# Se quitan los registros de periodos acumulados (se reemplazan por sus días) y los mensuales vacíos
fuera_p <- traza |> filter(resultado != "días sueltos: se dejan") |> select(STATION_ID, desde, hasta)
sacar <- sosp |> inner_join(fuera_p, by = "STATION_ID", relationship = "many-to-many") |>
  filter(FECHA >= desde, FECHA <= hasta) |> distinct(STATION_ID, FECHA)
st <- st |> anti_join(sacar, by = c("STATION_ID", "FECHA")) |>
  bind_rows(bind_rows(desag)) |> arrange(STATION_ID, FECHA)
stopifnot(!any(duplicated(st[, c("STATION_ID", "FECHA")])))
base::cat("Limpieza 1: periodos revisados:", nrow(traza), "| acumulados desagregados:",
          sum(traza$resultado == "acumulado: desagregado"), "en",
          n_distinct(traza$STATION_ID[traza$resultado == "acumulado: desagregado"]), "estaciones",
          "| registros reemplazados:", nrow(sacar), "| días generados:", nrow(bind_rows(desag)), "\n")
rm(st_d, rs)
foto(st, "3. Limpieza 1: acumulados repartidos en días")

# ------------------------------------------------- 9. Limpieza 2: lluvias de más de 300 mm y totales mensuales en un día
# Rango posible: 0 mm/día como mínimo. Para el máximo no se usa un tope fijo: todo día de más de
# 300 mm se contrasta con las estaciones vecinas, porque en la costa Pacífica un día así puede ser real.
#   - Vecinas: estaciones a <= 25 km. Se toma su lluvia máxima el mismo día, el anterior o el
#     siguiente (no todas leen el pluviómetro a la misma hora).
#   - Confirmado: alguna vecina registra al menos el 30 % de ese valor. Calibración (1981-1995): en
#     los días de 100-200 mm (casi todos reales) el 86 % pasa esta prueba.
#   - Satélite: se muestra la lluvia de la péntada en el píxel de la estación como apoyo. No decide,
#     porque un píxel de 5,5 km promedia y casi nunca alcanza el valor de un pluviómetro.
#   - No confirmado: el dato se quita (pasa a faltante, para imputarlo después).
UMBRAL_ALTO <- 300; PROP_VECINA <- 0.3

ancha <- matriz_diaria(st)
M <- as.matrix(ancha[, -1]); ids_m <- as.numeric(colnames(M))
# Máximo de cada estación en la ventana día anterior / mismo día / día siguiente
M3 <- pmax(M, rbind(NA, M[-nrow(M), ]), rbind(M[-1, ], NA), na.rm = TRUE)
dist_m <- distancias(ids_m)

pentada_de <- function(f) (month(f) - 1) * 6 + pmin(ceiling(day(f) / 5), 6)
altos <- st |> filter(VALUE > UMBRAL_ALTO) |>
  select(STATION_ID, FECHA, YEAR, MONTH, DAY, VALUE, origen, donante) |> arrange(desc(VALUE))
ia <- match(altos$FECHA, ancha$FECHA); ja <- match(altos$STATION_ID, ids_m)
# Un día desagregado no se puede confirmar con su propia donante (tiene su misma forma)
vecinas_de <- function(j, donante) which(dist_m[j, ] <= RADIO_KM & !(ids_m %in% donante))
altos$n_vecinas <- mapply(function(i, j, d) sum(!is.na(M3[i, vecinas_de(j, d)])), ia, ja, altos$donante)
altos$max_vecina_mm <- mapply(function(i, j, d) {
  v <- M3[i, vecinas_de(j, d)]; if (all(is.na(v))) NA else max(v, na.rm = TRUE) }, ia, ja, altos$donante)
altos$ultimo_dia_mes <- altos$DAY == days_in_month(altos$FECHA)
pix <- mc[match(altos$STATION_ID, mc$STATION_ID), c("fila", "col")]
capa <- match(paste(altos$YEAR, pentada_de(altos$FECHA)), paste(inf$anio, inf$pentada))
altos$satelite_pentada_mm <- cubo[cbind(pix$fila, pix$col, capa)]
altos <- altos |> mutate(
  decision = ifelse(!is.na(max_vecina_mm) & max_vecina_mm >= PROP_VECINA * VALUE,
                    "confirmado por vecinas", "no confirmado: se quita"))
print(altos |> count(decision, mas_de_500 = VALUE > 500, ultimo_dia_mes))
guardar_tabla(altos |> select(-FECHA, -donante), "dias_mayores_300mm")

# Se quitan los no confirmados (la fila desaparece = día faltante, igual que en los archivos)
n_filas <- nrow(st)
quitar <- altos |> filter(decision != "confirmado por vecinas") |> select(STATION_ID, FECHA)
# Mismo contraste para el "total mensual en un solo día" de menos de 300 mm: mes completo
# (>= 28 días con dato) en el que solo llovió un día, con >= 50 mm. Si ninguna vecina registra al
# menos el 30 % de ese valor (día anterior, mismo o siguiente), el valor se quita.
un_dia <- st |> group_by(STATION_ID, YEAR, MONTH) |>
  filter(n() >= 28, sum(VALUE > 0) == 1, max(VALUE) >= 50, max(VALUE) <= UMBRAL_ALTO) |>
  filter(VALUE > 0) |> ungroup() |> select(STATION_ID, FECHA, VALUE)
iu <- match(un_dia$FECHA, ancha$FECHA); ju <- match(un_dia$STATION_ID, ids_m)
un_dia$max_vecina_mm <- mapply(function(i, j) {
  v <- M3[i, vecinas_de(j, NA)]; if (all(is.na(v))) NA else max(v, na.rm = TRUE) }, iu, ju)
un_dia <- un_dia |> mutate(decision = ifelse(is.na(max_vecina_mm) | max_vecina_mm >= PROP_VECINA * VALUE,
                                             "se deja", "no confirmado: se quita"))
print(count(un_dia, decision))
quitar <- bind_rows(quitar, un_dia |> filter(decision != "se deja") |> select(STATION_ID, FECHA))
st <- st |> anti_join(quitar, by = c("STATION_ID", "FECHA"))
base::cat("Limpieza 2: se quitaron", n_filas - nrow(st), "días de", n_distinct(quitar$STATION_ID),
          "estaciones | máximo que queda:", max(st$VALUE), "mm\n")
foto(st, "4. Limpieza 2: > 300 mm y totales en un día no confirmados")

# ------------------------------------------------- 10. Limpieza 3: ceros falsos
# Un pluviómetro dañado o sin leer suele dejar ceros en vez de faltantes. Pero en el Valle hay
# sequías reales (dic-feb, jun-ago y años de El Niño: 1982, 1991-92, 1997, 2001-02), así que una
# racha larga de ceros no es falsa por sí sola: se contrasta con las vecinas.
#   - Racha: 15 o más días seguidos con 0 mm.
#   - Vecinas parecidas: a <= 25 km y con correlación de lluvia diaria >= 0,4 (calculada con
#     1981-1995, >= 1 año en común),
#     con datos en al menos el 90 % de los días de la racha.
#   - Cero falso: al menos la mitad de esas vecinas registró >= 50 mm por cada 30 días de racha.
#     Calibración (ventanas de 30 días, 1981-1995): si una vecina parecida recibe 50-100 mm, la
#     estación queda toda en 0 solo el 0,07 % de las veces; con 0 mm en la vecina, el 16 %.
#   - Sin vecinas parecidas se usa el satélite en el píxel de la estación, que es menos preciso:
#     falso solo si el píxel tiene >= 150 mm por cada 30 días (calibración 1981-1995: 0,07 % en 0).
#   - Sequía real: se deja. Cero falso: la racha pasa a faltante (para imputar después).
MIN_RACHA_CERO <- 15; CORR_VECINA <- 0.4

st_antes3 <- st                                               # para la gráfica antes/después
ancha <- matriz_diaria(st)
M <- as.matrix(ancha[, -1]); ids_m <- as.numeric(colnames(M)); dist_m <- distancias(ids_m)
tr <- year(ancha$FECHA) %in% ENTRENA                         # correlaciones solo con entrenamiento
P <- !is.na(M[tr, ]) * 1
comun_m <- crossprod(P)                                        # días en común entre estaciones
corr_m <- suppressWarnings(cor(M[tr, ], use = "pairwise.complete.obs"))
corr_m[comun_m < 365 | is.na(corr_m)] <- -1

# Satélite repartido por día en el píxel de cada estación (péntada / días de la péntada)
capa_dia <- match(paste(year(ancha$FECHA), pentada_de(ancha$FECHA)), paste(inf$anio, inf$pentada))
dias_pent <- tabulate(capa_dia, nbins = nrow(inf))
pix_m <- mc[match(ids_m, mc$STATION_ID), c("fila", "col")]

rachas0 <- list()
for (j in seq_along(ids_m)) {
  r <- rle(!is.na(M[, j]) & M[, j] == 0)
  fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  k <- which(r$values & r$lengths >= MIN_RACHA_CERO)
  if (!length(k)) next
  vec <- which(dist_m[j, ] <= RADIO_KM & corr_m[j, ] >= CORR_VECINA)
  sat_j <- if (is.na(pix_m$fila[j])) rep(NA, nrow(M)) else cubo[pix_m$fila[j], pix_m$col[j], capa_dia] / dias_pent[capa_dia]
  for (q in k) {
    w <- ini[q]:fin[q]; L <- length(w); umbral <- 50 * L / 30
    blk <- M[w, vec, drop = FALSE]; ok <- colMeans(is.na(blk)) <= .1
    tot <- colSums(blk[, ok, drop = FALSE], na.rm = TRUE)
    sat <- sum(sat_j[w])
    decision <- if (sum(ok) > 0) {
      ifelse(mean(tot >= umbral) >= .5, "cero falso: se quita", "sequía real: se deja")
    } else if (!is.na(sat)) {
      ifelse(sat >= 150 * L / 30, "cero falso: se quita", "sequía real: se deja")
    } else "sin información: se deja"
    rachas0[[length(rachas0) + 1]] <- tibble(
      STATION_ID = ids_m[j], desde = ancha$FECHA[ini[q]], hasta = ancha$FECHA[fin[q]], dias = L,
      vecinas = sum(ok), vecinas_con_lluvia = sum(tot >= umbral),
      mediana_vecinas_mm = if (sum(ok)) median(tot) else NA, satelite_mm = round(sat, 1),
      decidido_con = if (sum(ok)) "vecinas" else "satélite", decision = decision)
  }
}
rachas0 <- bind_rows(rachas0)
print(rachas0 |> count(decidido_con, decision, name = "rachas") |>
        left_join(rachas0 |> group_by(decidido_con, decision) |> summarise(dias = sum(dias), .groups = "drop"),
                  by = c("decidido_con", "decision")))
guardar_tabla(rachas0, "rachas_ceros_15dias")
falsos <- rachas0 |> filter(decision == "cero falso: se quita") |>
  rowwise() |> reframe(STATION_ID, FECHA = seq(desde, hasta, by = "day"))
st <- st |> anti_join(falsos, by = c("STATION_ID", "FECHA"))
base::cat("Limpieza 3: ceros falsos quitados:", nrow(falsos), "días en",
          n_distinct(falsos$STATION_ID), "estaciones\n")
rm(M, P, comun_m, corr_m)
foto(st, "5. Limpieza 3: ceros falsos")

# ------------------------------------------------- 11. Limpieza 4: meses inconsistentes con las vecinas
# Las limpiezas anteriores miran días sueltos. Un error sostenido (p. ej. 270 mm todos los días
# durante dos meses) no tiene ningún día imposible, pero su total mensual sí lo es.
#   - Se comparan los meses con >= 25 días con dato contra la vecina MÁS lluviosa (<= 25 km) con el
#     mes también completo: razón r = total de la estación / total de esa vecina.
#   - Cada estación tiene su razón típica R (mediana de r en 1981-1995): una estación de la costa
#     Pacífica puede ser normalmente más lluviosa que sus vecinas y eso no es un error.
#   - Mes inconsistente: total > 600 mm y r > 2 x max(R, 1). El mínimo de 600 mm evita marcar meses
#     secos (40 vs 15 mm también es "el doble").
#   - Estación con R > 3 (siempre más del triple que la vecina más lluviosa): la serie entera es
#     dudosa (otra unidad o mal ubicada) y se quita completa.
#   - Mes inconsistente -> todos sus días pasan a faltante.
MIN_MES <- 600; VECES_VECINA <- 2; MIN_DIAS_MES <- 25; R_MAX_ESTACION <- 3
mt <- st |> group_by(STATION_ID, YEAR, MONTH) |> summarise(total = sum(VALUE), n = n(), .groups = "drop") |>
  filter(n >= MIN_DIAS_MES)
piv <- mt |> select(-n) |> pivot_wider(names_from = STATION_ID, values_from = total) |> arrange(YEAR, MONTH)
Mm <- as.matrix(piv[, -(1:2)]); ids_mm <- as.numeric(colnames(Mm)); dist_mm <- distancias(ids_mm)
max_vec <- sapply(seq_along(ids_mm), function(j) {
  v <- Mm[, dist_mm[j, ] <= RADIO_KM, drop = FALSE]
  if (!ncol(v)) return(rep(NA_real_, nrow(Mm)))
  apply(v, 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE)) })
mt <- mt |> mutate(max_vecina_mm = max_vec[cbind(match(paste(YEAR, MONTH), paste(piv$YEAR, piv$MONTH)),
                                                 match(STATION_ID, ids_mm))],
                   r = total / pmax(max_vecina_mm, 1))
R_tipica <- mt |> filter(YEAR %in% ENTRENA, !is.na(r), max_vecina_mm >= 50) |>
  group_by(STATION_ID) |> summarise(R = median(r), meses_R = n())
mt <- mt |> left_join(R_tipica, by = "STATION_ID") |> mutate(R = coalesce(R, 1),
                                                             inconsistente = !is.na(max_vecina_mm) & total > MIN_MES & r > VECES_VECINA * pmax(R, 1))
base::cat("Limpieza 4 - calibración (1981-1995): meses comparables", sum(!is.na(mt$max_vecina_mm) & mt$YEAR %in% ENTRENA),
          "| inconsistentes", sum(mt$inconsistente & mt$YEAR %in% ENTRENA), "(",
          round(100 * mean(mt$inconsistente[!is.na(mt$max_vecina_mm) & mt$YEAR %in% ENTRENA]), 2), "% )\n")
est_mala <- R_tipica |> filter(R > R_MAX_ESTACION)
meses_malos <- mt |> filter(inconsistente | STATION_ID %in% est_mala$STATION_ID) |>
  mutate(decision = ifelse(STATION_ID %in% est_mala$STATION_ID, "estación completa: se quita",
                           "mes inconsistente: se quita"))
print(meses_malos |> count(decision))
print(est_mala)
guardar_tabla(meses_malos |> mutate(across(c(r, R), ~ round(.x, 2))), "meses_inconsistentes")
n_antes4 <- nrow(st)
st <- st |>
  anti_join(meses_malos |> filter(decision == "mes inconsistente: se quita"), by = c("STATION_ID", "YEAR", "MONTH")) |>
  filter(!STATION_ID %in% est_mala$STATION_ID)
base::cat("Limpieza 4: se quitaron", n_antes4 - nrow(st), "días |", sum(meses_malos$decision == "mes inconsistente: se quita"),
          "meses y", nrow(est_mala), "estaciones completas\n")
# Año sin vecinas para comparar: 5301100101 en 2001 (ninguna estación a <= 25 km tiene ese año).
# Se revisó con el satélite y con la propia estación, y falla en tres pruebas independientes:
#   - total de 14.240 mm, 48 % más que su año más lluvioso (9.636 mm en 1998);
#   - el satélite en su píxel marca un año normal (5.047 mm), y la razón estación/satélite pasa de
#     1,4 (su mediana de 1981-1995) a 2,8;
#   - solo el 7 % de días en 0 (en sus otros años, 14-46 %) y ningún mes por debajo de 940 mm.
# Se quita el año completo (no se imputa).
ANIO_DUDOSO <- tibble(STATION_ID = 5301100101, YEAR = 2001)
n_antes_a <- nrow(st)
st <- st |> anti_join(ANIO_DUDOSO, by = c("STATION_ID", "YEAR"))
base::cat("Año dudoso sin vecinas quitado:", n_antes_a - nrow(st), "días\n")
foto(st, "6. Limpieza 4: meses inconsistentes con las vecinas")

# ------------------------------------------------- 12. Completitud 1981-2003 (ya con las limpiezas 1 a 4)
comp <- completitud(st, Y0, Y1)
cal <- comp$cal; ea <- comp$ea; rm(comp)
base::cat("Filas del calendario (vida útil de cada estación):", nrow(cal), "\n",
          "Días con dato:", sum(cal$presente), "\n",
          "Días faltantes (NA):", sum(!cal$presente), "\n")

# Resumen por estación-mes con los datos limpios, antes de imputar (para las gráficas G4 y G7)
mes <- st |>
  group_by(STATION_ID, YEAR, MONTH) |>
  summarise(n = n(), total = sum(VALUE), .groups = "drop")
# Meses completos (>= 28 días con dato) con lluvia total = 0 que quedan tras la limpieza 3
ceros <- mes |> filter(n >= 28, total == 0)

# ------------------------------------------------- 13. Ranking por estación
# Cobertura sobre la vida útil de la estación (de su primer a su último año con dato) y
# cobertura en todo 1981-2003, para comparar estaciones de distinta antigüedad.
meta <- ea |> group_by(STATION_ID) |>
  summarise(primer_anio = min(YEAR[n_dias > 0]),
            ultimo_anio = max(YEAR[n_dias > 0]),
            anios_con_algun_dato = sum(n_dias > 0),
            anios_90 = sum(cob >= 90),
            cobertura_vida_util = 100 * sum(n_dias) / sum(dias[YEAR >= primer_anio & YEAR <= ultimo_anio]),
            cobertura_1981_2003 = 100 * sum(n_dias) / sum(dias),
            .groups = "drop") |>
  left_join(ll |> select(STATION_ID = CODIGO, CATEGORY, STATUS), by = "STATION_ID") |>
  left_join(ubic, by = "STATION_ID") |>
  arrange(desc(anios_90), desc(cobertura_vida_util))
guardar_tabla(meta, "ranking_estaciones")

# ------------------------------------------------- 14. Sensibilidad del criterio de cobertura (1981-2003)
sens <- expand_grid(cobertura_min = c(50, 60, 70, 75, 80, 85, 90, 95, 100),
                    hueco_max_dias = c(Inf, 31)) |>
  rowwise() |>
  mutate(ok = list(ea |> filter(cob >= cobertura_min, hueco <= hueco_max_dias)),
         pares_estacion_anio = nrow(ok),
         estaciones_con_1_anio = n_distinct(ok$STATION_ID),
         pares_entrenamiento_1981_1995 = sum(ok$YEAR %in% ENTRENA),
         estaciones_con_10_anios = sum(table(ok$STATION_ID) >= 10)) |>
  ungroup() |> select(-ok) |>
  mutate(hueco_max_dias = ifelse(is.infinite(hueco_max_dias), "sin limite", "31"))
guardar_tabla(sens, "sensibilidad_criterio_cobertura")

# ------------------------------------------------- 15. Criterio de cobertura por estación-año
# Un año de una estación entra al modelo si tiene >= 90 % de sus días con dato y ningún hueco de
# más de 31 días seguidos. La curva anual se construye con el año completo: imputar un año con
# medio año vacío sería inventar la curva, y un hueco de más de un mes puede borrar uno de los dos
# picos de lluvia. Los años que no cumplen se quitan; los huecos de los que sí cumplen se imputan
# después. (La tabla de sensibilidad muestra cuántos pares quedan con otros umbrales.)
COB_MIN <- 90; HUECO_MAX <- 31
pares_ok <- ea |> filter(cob >= COB_MIN, hueco <= HUECO_MAX) |> select(STATION_ID, YEAR)
st  <- st  |> semi_join(pares_ok, by = c("STATION_ID", "YEAR"))
# Calendario para la imputación: el AÑO COMPLETO de cada par aceptado (la curva necesita los 365/366
# días; si la estación empezó o terminó a mitad de un año aceptado, esos días también se imputan)
cal <- pares_ok |> rowwise() |>
  reframe(STATION_ID, FECHA = seq(make_date(YEAR, 1, 1), make_date(YEAR, 12, 31), by = "day")) |>
  left_join(st |> select(STATION_ID, FECHA, VALUE), by = c("STATION_ID", "FECHA")) |>
  mutate(YEAR = year(FECHA), presente = !is.na(VALUE))
base::cat("Criterio de cobertura: pares estación-año", nrow(distinct(ea |> filter(n_dias > 0), STATION_ID, YEAR)),
          "->", nrow(pares_ok), "| estaciones", n_distinct(ea$STATION_ID[ea$n_dias > 0]), "->",
          n_distinct(pares_ok$STATION_ID), "| días por imputar:", sum(!cal$presente), "\n")
foto(st, "7. Criterio de cobertura (>= 90 % y hueco <= 31 días)", cal)

# ------------------------------------------------- 16. Imputación de los días faltantes
# Solo se imputan los días faltantes de los años que pasaron el criterio de cobertura.
# Se comparan cuatro métodos que usan las estaciones vecinas (<= 50 km) que sí midieron ese día:
#   - IDW: promedio de las 8 más cercanas con peso 1/d^2 (solo distancia).
#   - Kriging ordinario: 8 más cercanas, pesos del sistema de kriging con un variograma
#     exponencial común ajustado con días de 1981-1995 (correlación modelada con la distancia).
#   - SBE (vecina donante + razón normal): la vecina con mayor correlación diaria (1981-1995).
#   - CCW: promedio de las 8 vecinas más correlacionadas con peso r^2.
# Todos trabajan sobre la lluvia RELATIVA (lluvia / promedio diario de la estación en ese mes,
# calculado con 1981-1995) para no mezclar sitios secos y húmedos (orografía), y luego se pasa a mm
# con el promedio de la estación a imputar. IDW, Kriging y CCW primero deciden si llovió (mayoría
# ponderada de vecinas con lluvia) para no llenar los días secos de llovizna.
# Si la estación no tiene años de entrenamiento, su promedio mensual sale del satélite en su píxel.
# Validación: se esconden 30.000 días conocidos de 1981-1995, se imputan y se comparan con lo medido.
# Se elige el método con mayor KGE (métrica del proyecto). Sin vecinas: satélite del píxel.
RADIO_IMP <- 50; K_VEC <- 8; N_VALIDA <- 30000
set.seed(2026)
ancha <- matriz_diaria(st)
M <- as.matrix(ancha[, -1]); ids_m <- as.numeric(colnames(M)); fechas <- ancha$FECHA
dist_m <- distancias(ids_m); mes_f <- month(fechas); tr <- year(fechas) %in% ENTRENA
capa_dia <- match(paste(year(fechas), pentada_de(fechas)), paste(inf$anio, inf$pentada))
dias_pent <- tabulate(capa_dia, nbins = nrow(inf))
pix_m <- mc[match(ids_m, mc$STATION_ID), c("fila", "col")]
sat_dia <- sapply(seq_along(ids_m), function(j) cubo[pix_m$fila[j], pix_m$col[j], capa_dia] / dias_pent[capa_dia])
clim_sat <- apply(sat_dia[tr, ], 2, function(x) tapply(x, mes_f[tr], mean))
clim <- sapply(seq_along(ids_m), function(j) tapply(M[tr, j], mes_f[tr], mean, na.rm = TRUE)[as.character(1:12)])
clim[is.na(clim)] <- clim_sat[is.na(clim)]
clim <- pmax(clim, 0.1)
idx_clim <- cbind(rep(mes_f, ncol(M)), rep(seq_along(ids_m), each = nrow(M)))
REL <- M / matrix(clim[idx_clim], nrow(M))                         # lluvia relativa
P <- !is.na(M[tr, ]) * 1; comun_m <- crossprod(P)
corr_m <- suppressWarnings(cor(M[tr, ], use = "pairwise.complete.obs"))
corr_m[comun_m < 365 | is.na(corr_m)] <- NA; diag(corr_m) <- NA

# Variograma experimental de la lluvia relativa (400 días de entrenamiento) y ajuste exponencial
dias_v <- sample(which(tr), 400)
vg <- do.call(rbind, lapply(dias_v, function(t) {
  o <- which(!is.na(REL[t, ])); if (length(o) < 5) return(NULL)
  d <- dist_m[o, o]; g <- outer(REL[t, o], REL[t, o], function(x, y) .5 * (x - y)^2)
  k <- upper.tri(d) & d <= 100; cbind(h = d[k], g = g[k])
}))
bins <- cut(vg[, "h"], seq(0, 100, 5))
variog <- tibble(h = as.vector(tapply(vg[, "h"], bins, mean)), gamma = as.vector(tapply(vg[, "g"], bins, mean)),
                 pares = as.vector(table(bins))) |> filter(!is.na(gamma))
ajuste_v <- optim(c(.5, .5, 20), function(p) sum(variog$pares * (variog$gamma - (p[1] + p[2] * (1 - exp(-variog$h / p[3]))))^2),
                  method = "L-BFGS-B", lower = c(0, .01, 1), upper = c(10, 10, 300))$par
gam <- function(h) ifelse(h == 0, 0, ajuste_v[1] + ajuste_v[2] * (1 - exp(-h / ajuste_v[3])))
base::cat("Variograma exponencial: pepita", round(ajuste_v[1], 2), "| meseta parcial", round(ajuste_v[2], 2),
          "| rango", round(ajuste_v[3], 1), "km\n")

imputa <- function(j, t) {
  o <- which(!is.na(M[t, ])); o <- o[o != j & dist_m[j, o] <= RADIO_IMP]
  out <- c(IDW = NA, Kriging = NA, SBE = NA, CCW = NA)
  if (!length(o)) return(out)
  ct <- clim[mes_f[t], j]; a <- REL[t, ]; llovio <- M[t, ] > 0
  n1 <- o[order(dist_m[j, o])][seq_len(min(K_VEC, length(o)))]
  w <- 1 / pmax(dist_m[j, n1], .5)^2; w <- w / sum(w)
  out["IDW"] <- if (sum(w * llovio[n1]) >= .5) sum(w * a[n1]) * ct else 0
  D1 <- dist_m[n1, n1, drop = FALSE]; diag(D1) <- 0; n <- length(n1)
  G <- matrix(1, n + 1, n + 1); G[1:n, 1:n] <- gam(D1); G[n + 1, n + 1] <- 0
  lam <- tryCatch(solve(G, c(gam(dist_m[j, n1]), 1))[1:n], error = function(e) w)
  out["Kriging"] <- if (sum(lam * llovio[n1]) >= .5) max(0, sum(lam * a[n1])) * ct else 0
  rc <- corr_m[j, o]; okc <- which(!is.na(rc) & rc > 0)
  if (length(okc)) {
    oc <- o[okc][order(-rc[okc])]
    out["SBE"] <- a[oc[1]] * ct
    n2 <- oc[seq_len(min(K_VEC, length(oc)))]; w2 <- corr_m[j, n2]^2; w2 <- w2 / sum(w2)
    out["CCW"] <- if (sum(w2 * llovio[n2]) >= .5) sum(w2 * a[n2]) * ct else 0
  } else out[c("SBE", "CCW")] <- out["IDW"]
  out
}

# Validación con días conocidos de entrenamiento
conocidos <- which(!is.na(M) & tr, arr.ind = TRUE)
muestra <- conocidos[sample(nrow(conocidos), N_VALIDA), ]
pred <- t(apply(muestra, 1, function(x) imputa(x[2], x[1])))
medido <- M[muestra]; okv <- complete.cases(pred)
kge <- function(s, o) { r <- cor(s, o); 1 - sqrt((r - 1)^2 + (sd(s) / sd(o) - 1)^2 + (mean(s) / mean(o) - 1)^2) }
boot <- replicate(300, { b <- sample(which(okv), sum(okv), TRUE); sapply(colnames(pred), function(m) kge(pred[b, m], medido[b])) })
comparacion <- bind_rows(lapply(colnames(pred), function(m) {
  s <- pred[okv, m]; o <- medido[okv]
  tibble(metodo = m, MAE_mm = mean(abs(s - o)), RMSE_mm = sqrt(mean((s - o)^2)),
         sesgo_pct = 100 * (sum(s) - sum(o)) / sum(o), correlacion = cor(s, o), KGE = kge(s, o),
         KGE_IC95_inf = quantile(boot[m, ], .025), KGE_IC95_sup = quantile(boot[m, ], .975),
         ratio_desv_estandar = sd(s) / sd(o), pct_dias_secos_medido = 100 * mean(o == 0),
         pct_dias_secos_imputado = 100 * mean(s == 0), acierto_lluvia_o_seco_pct = 100 * mean((s > 0) == (o > 0)),
         p95_medido = quantile(o, .95), p95_imputado = quantile(s, .95))
})) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
print(comparacion, width = Inf)
guardar_tabla(comparacion, "comparacion_metodos_imputacion")
MEJOR <- comparacion$metodo[which.max(comparacion$KGE)]
base::cat("Método elegido (mayor KGE):", MEJOR, "\n")

# Imputación de los faltantes con el método elegido
faltan <- cal |> filter(!presente) |> select(STATION_ID, FECHA)
fi <- cbind(match(faltan$FECHA, fechas), match(faltan$STATION_ID, ids_m))
val <- apply(fi, 1, function(x) imputa(x[2], x[1])[MEJOR])
sin_vec <- is.na(val)
val[sin_vec] <- (sat_dia[fi[sin_vec, , drop = FALSE]] / clim_sat[cbind(mes_f[fi[sin_vec, 1]], fi[sin_vec, 2])]) *
  clim[cbind(mes_f[fi[sin_vec, 1]], fi[sin_vec, 2])]
imputados <- faltan |> mutate(VALUE = round(val, 1), YEAR = year(FECHA), MONTH = month(FECHA), DAY = day(FECHA),
                              origen = ifelse(sin_vec, "imputado (satélite)", paste0("imputado (", MEJOR, ")")))
base::cat("Imputados:", nrow(imputados), "días |", "con vecinas:", sum(!sin_vec), "| con satélite:", sum(sin_vec), "\n")

st_antes_imp <- st
st <- bind_rows(st, imputados) |> arrange(STATION_ID, FECHA)
stopifnot(!any(duplicated(st[, c("STATION_ID", "FECHA")])))
cal <- cal |> left_join(imputados |> select(STATION_ID, FECHA, imputado = VALUE), by = c("STATION_ID", "FECHA")) |>
  mutate(VALUE = coalesce(VALUE, imputado), imputado = !is.na(imputado)) # calendario después (presente = antes)

# Antes / después: ¿la imputación cambia la distribución de la lluvia?
resumen_lluvia <- function(x, cuando) tibble(
  datos = cuando, dias = length(x), media_mm = mean(x), mediana_mm = median(x), desv_estandar = sd(x),
  pct_dias_secos = 100 * mean(x == 0), p90_mm = quantile(x, .9), p99_mm = quantile(x, .99), max_mm = max(x))
antes_despues <- bind_rows(resumen_lluvia(st_antes_imp$VALUE, "antes (solo medido)"),
                           resumen_lluvia(st$VALUE, "después (medido + imputado)"),
                           resumen_lluvia(imputados$VALUE, "solo los días imputados")) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
print(antes_despues, width = Inf)
guardar_tabla(antes_despues, "antes_despues_imputacion")
foto(st, "8. Imputación de faltantes", cal |> mutate(presente = TRUE))

# ------------------------------------------------- 17. Datos limpios: diario y por péntadas (CSV)
# Para construir los datos funcionales, la estación debe estar en la misma escala que el satélite:
# 72 péntadas por año (días 1-5, 6-10, 11-15, 16-20, 21-25 y 26-fin de cada mes; la sexta tiene
# de 3 a 6 días). La lluvia de la péntada es la SUMA de sus días (mm/péntada), igual que el satélite.
# Se agrega el valor del satélite en el píxel de la estación para la misma péntada.
diario <- st |> select(STATION_ID, FECHA, YEAR, MONTH, DAY, VALUE, origen) |> arrange(STATION_ID, FECHA)
pentadas <- diario |> mutate(PENTADA = pentada_de(FECHA)) |>
  group_by(STATION_ID, YEAR, PENTADA) |>
  summarise(lluvia_mm = sum(VALUE), dias = n(), dias_imputados = sum(grepl("imputado", origen)),
            dias_desagregados = sum(origen == "desagregado"), .groups = "drop")
pix_p <- mc[match(pentadas$STATION_ID, mc$STATION_ID), c("fila", "col")]
capa_p <- match(paste(pentadas$YEAR, pentadas$PENTADA), paste(inf$anio, inf$pentada))
pentadas$satelite_mm <- cubo[cbind(pix_p$fila, pix_p$col, capa_p)]
stopifnot(all(count(pentadas, STATION_ID, YEAR)$n == 72))           # cada año completo: 72 péntadas
base::cat("Datos limpios: diario", nrow(diario), "filas | péntadas", nrow(pentadas), "filas (",
          nrow(distinct(pentadas, STATION_ID, YEAR)), "curvas estación-año x 72 )\n")
write.csv(diario, file.path(OUT, "estaciones_diario_limpio.csv"), row.names = FALSE)
write.csv(pentadas, file.path(OUT, "estaciones_pentadas_limpio.csv"), row.names = FALSE)

trazabilidad <- bind_rows(traz)
print(trazabilidad, width = Inf)
guardar_tabla(trazabilidad, "trazabilidad")

# ------------------------------------------------- 18. Estadísticas descriptivas
# Estadísticas descriptivas de la lluvia: estaciones limpias e imputadas y satélite
moda <- function(x) { x <- x[!is.na(x)]; u <- unique(x); u[which.max(tabulate(match(x, u)))] }
descriptivas <- function(x, fuente) {
  tibble(fuente = fuente, media = mean(x, na.rm = TRUE), mediana = median(x, na.rm = TRUE),
         moda = moda(x), max = max(x, na.rm = TRUE), min = min(x, na.rm = TRUE))
}
estadisticas <- bind_rows(descriptivas(st$VALUE, "Estaciones (mm/día, 1981-2003, limpias e imputadas)"),
                          descriptivas(as.vector(cubo), "Satélite (mm/péntada, 1981-2003)"))
print(estadisticas)
guardar_tabla(estadisticas, "estadisticas_descriptivas")

# ================================================================== 19. GRÁFICAS (solo se imprimen)
eje_anios <- scale_x_continuous(breaks = seq(S0, S1, 2))

# G3. Número de estaciones según su cobertura, año por año (un panel por año)
# Para cada año se cuentan las estaciones que tienen algún dato ese año, según su cobertura:
# < 50 %, 50-80 %, 80-90 % y >= 90 %. La cobertura es el % de días del AÑO COMPLETO con dato
# (la misma del mapa de calor): para construir la curva anual de una estación se necesita el año
# entero, así que un año en que la estación solo funcionó unos meses no es un año completo.
ANIOS_PANEL <- Y0:Y1
cob_anio <- ea |> filter(cob > 0, YEAR %in% ANIOS_PANEL) |>
  mutate(clase = cut(cob, c(-Inf, 50, 80, 90, Inf), right = FALSE,
                     labels = c("< 50 %", "50-80 %", "80-90 %", ">= 90 %"))) |>
  count(YEAR, clase) |> complete(YEAR, clase, fill = list(n = 0))
print(cob_anio |> pivot_wider(names_from = clase, values_from = n) |>
        mutate(total = rowSums(across(-YEAR))), n = Inf)

mostrar(
  ggplot(cob_anio, aes(clase, n, fill = clase)) + geom_col(width = .8) +
    geom_text(aes(label = n), vjust = -.4, size = 2.3) +
    facet_wrap(~YEAR, ncol = 7) +
    scale_fill_manual(values = c("< 50 %" = "#c0392b", "50-80 %" = "#e8a33d",
                                 "80-90 %" = "#a6c8e6", ">= 90 %" = "#2166ac"), name = "Cobertura") +
    scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 50)) +
    labs(title = "Número de estaciones según su cobertura, año por año", x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(face = "bold"), plot.title = element_text(hjust = .5),
          legend.position = "right")
)

# G4. Serie de tiempo de la lluvia 1981-2003
# Total mensual de cada estación (solo meses con >= 28 días con dato, datos limpios sin imputar) y,
# para cada mes, la mediana entre estaciones (menos sensible a estaciones muy lluviosas que la media).
serie <- mes |> filter(n >= 28, YEAR >= S0, YEAR <= S1) |>
  group_by(YEAR, MONTH) |> summarise(lluvia = median(total), estaciones = n(), .groups = "drop") |>
  mutate(fecha = make_date(YEAR, MONTH, 1))
mostrar(
  ggplot(serie, aes(fecha, lluvia)) + geom_line(colour = "#2166ac", linewidth = .4) +
    scale_x_date(breaks = make_date(seq(S0, S1, 2), 1, 1), date_labels = "%Y") +
    labs(title = "Serie de tiempo de la lluvia mensual 1981-2003", x = NULL,
         y = "lluvia mensual (mm, mediana entre estaciones)")
)

# G5. Hueco más largo en los pares estación-año con algún dato
f5 <- ea |> filter(n_dias > 0) |>
  mutate(rango = cut(hueco, c(-Inf, 0, 1, 7, 31, 91, Inf),
                     labels = c("0 (año completo)", "1 día", "2-7 días", "8-31 días", "32-91 días", "> 3 meses"))) |>
  count(rango)
mostrar(
  ggplot(f5, aes(rango, n, fill = as.integer(rango) > 3)) + geom_col() + geom_text(aes(label = n), vjust = -.3) +
    scale_fill_manual(values = c("#2166ac", "#b2182b"), guide = "none") +
    labs(title = "¿Los faltantes son días sueltos o bloques largos?",
         x = "hueco más largo dentro del año (días consecutivos sin dato)", y = "pares estación-año")
)

# G6. Días con más de 300 mm por año: confirmados por las vecinas (azul) y quitados (rojo)
mostrar(
  altos |> count(YEAR, decision) |> complete(YEAR = Y0:Y1, decision, fill = list(n = 0)) |>
    ggplot(aes(YEAR, n, fill = decision)) + geom_col() + eje_anios +
    scale_fill_manual(values = c("confirmado por vecinas" = "#2166ac", "no confirmado: se quita" = "#b2182b"),
                      name = NULL) +
    labs(title = "Días con más de 300 mm, por año", x = NULL, y = "cantidad") +
    theme(legend.position = "bottom")
)

# G6b. Días secos por mes, antes y después de quitar los ceros falsos (limpieza 3)
# % de días con 0 mm sobre los días con dato. Si la limpieza funciona, bajan los días secos de más
# en los meses lluviosos y las temporadas secas reales (dic-feb, jun-ago) se mantienen.
secos <- bind_rows(st_antes3 |> mutate(cuando = "antes"), st |> mutate(cuando = "después")) |>
  group_by(cuando, MONTH) |> summarise(pct = 100 * mean(VALUE == 0), .groups = "drop") |>
  mutate(cuando = factor(cuando, levels = c("antes", "después")))
mostrar(
  ggplot(secos, aes(factor(MONTH), pct, fill = cuando)) + geom_col(position = position_dodge(.8), width = .75) +
    geom_text(aes(label = round(pct, 1)), position = position_dodge(.8), vjust = -.3, size = 2.5) +
    scale_x_discrete(labels = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")) +
    scale_fill_manual(values = c("antes" = "#bdbdbd", "después" = "#2166ac"), name = NULL) +
    labs(title = "Días secos por mes, antes y después de quitar los ceros falsos", x = NULL,
         y = "% de días con 0 mm") +
    theme(legend.position = "bottom")
)

# G6c. Lluvia media por mes del año, antes y después de imputar
# Si la imputación es buena, el régimen bimodal (picos en abr-may y oct-nov) no debe cambiar.
clim_ad <- bind_rows(st_antes_imp |> mutate(cuando = "antes (solo medido)"),
                     st |> mutate(cuando = "después (medido + imputado)")) |>
  group_by(cuando, MONTH) |> summarise(mm = mean(VALUE), .groups = "drop")
mostrar(
  ggplot(clim_ad, aes(MONTH, mm, colour = cuando, linetype = cuando)) + geom_line(linewidth = .9) + geom_point() +
    scale_x_continuous(breaks = 1:12, labels = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")) +
    scale_colour_manual(values = c("#bdbdbd", "#2166ac"), name = NULL) +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    labs(title = "Lluvia media diaria por mes, antes y después de imputar", x = NULL, y = "mm/día") +
    theme(legend.position = "bottom")
)

# G7. Satélite en los meses donde la estación marca 0 (después de la limpieza 3)
# Lluvia del satélite (mismo mes y píxel) en los meses completos que la estación marca en 0.
# La péntada k pertenece al mes ceiling(k/6).
mes_idx <- (inf$anio - S0) * 12 + ceiling(inf$pentada / 6)
cz <- ceros |> filter(YEAR >= S0, YEAR <= S1) |>
  inner_join(mc |> filter(pixel_valido) |> select(STATION_ID, fila, col), by = "STATION_ID")
cz$sat_mm_mes <- mapply(function(y, m, f, c) sum(cubo[f, c, mes_idx == (y - S0) * 12 + m]),
                        cz$YEAR, cz$MONTH, cz$fila, cz$col)
mostrar(
  ggplot(cz, aes(sat_mm_mes)) + geom_histogram(breaks = seq(0, 400, 20), fill = "#2166ac", colour = "white") +
    geom_vline(xintercept = c(50, 100), linetype = "dashed", colour = "#b2182b") +
    labs(title = "Satélite en los meses donde la estación marca 0 (1981-2003)",
         x = "lluvia satelital en ese mes y píxel (mm)", y = "n.º de meses")
)

# Capa común para los mapas: borde del Valle del Cauca
capa_borde <- geom_path(data = borde, aes(x, y, group = interaction(geom, part)),
                        inherit.aes = FALSE, colour = "black", linewidth = .5)

# G8. Estaciones con coordenadas de catálogo y estaciones ubicadas con la propuesta
valle_df <- as.data.frame(valle, xy = TRUE)                          # solo los píxeles del Valle
mostrar(
  ggplot() + geom_raster(data = valle_df, aes(x, y), fill = "grey92") + capa_borde +
    geom_point(data = ubic |> mutate(origen = ifelse(fuente_coord == "catálogo",
                                                     "coordenadas de catálogo", "ubicación estimada")),
               aes(Lon, lat, colour = origen), size = 1.3) +
    scale_colour_manual(values = c("coordenadas de catálogo" = "#2166ac", "ubicación estimada" = "#e66101"),
                        labels = function(x) paste0(x, " (", table(ifelse(ubic$fuente_coord == "catálogo",
                                                                          "coordenadas de catálogo", "ubicación estimada"))[x], ")"),
                        name = NULL) +
    coord_equal() + labs(title = "Estaciones del Valle: con coordenadas y ubicadas", x = "Lon", y = "Lat") +
    theme(legend.position = "bottom")
)

# G9-G10. Mapas del satélite con el borde del Valle y las estaciones encima
# Cada cuadro es un píxel de 0,05° (~5,5 km). Solo píxeles del Valle. Puntos negros = estaciones.
xy <- as.data.frame(sat_r[[1]], xy = TRUE)[, 1:2]                  # centros de píxel en el orden de terra
mapa <- xy |> mutate(media_anual = as.vector(t(apply(cubo, c(1, 2), mean, na.rm = TRUE) * 72)),
                     pct_tope    = as.vector(t(apply(cubo == vmax, c(1, 2), mean, na.rm = TRUE) * 100)))
ok <- mc |> filter(pixel_valido)
mapa <- mapa |> filter(!is.na(media_anual))
mostrar(
  ggplot(mapa, aes(x, y, fill = media_anual)) + geom_raster() +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "mm/año") + capa_borde +
    geom_point(data = ok, aes(Lon, lat), inherit.aes = FALSE, size = .6) + coord_equal() +
    labs(title = "Satélite: lluvia media anual 1981-2003", x = "Lon", y = "Lat")
)
mostrar(
  ggplot(mapa, aes(x, y, fill = pct_tope)) + geom_raster() +
    scale_fill_distiller(palette = "Reds", direction = 1, name = "% de\npéntadas") + capa_borde +
    geom_point(data = ok, aes(Lon, lat), inherit.aes = FALSE, size = .6) + coord_equal() +
    labs(title = sprintf("Satélite: %% de péntadas que llegan al tope de %d mm", vmax), x = "Lon", y = "Lat")
)

base::cat("Tablas guardadas en", OUT, "\n")