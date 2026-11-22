# =============================================================================
# Preprocesamiento
# =============================================================================
paquetes <- c("dplyr", "tidyr", "lubridate", "ggplot2", "readxl", "writexl", "terra")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan) > 0) install.packages(faltan, repos = "https://cloud.r-project.org")
invisible(lapply(paquetes, library, character.only = TRUE))

# Generador aleatorio y semilla fijos (el muestreo de la calibración y el bootstrap de la Tabla 7)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
SEMILLA <- 2026
set.seed(SEMILLA)
options(stringsAsFactors = FALSE, warn = 1)

# Versiones usadas (para reportar en el informe)
cat("R:", R.version.string, "\n")
for (pq in paquetes) cat(sprintf("  %-10s %s\n", pq, as.character(packageVersion(pq))))

args <- commandArgs(trailingOnly = TRUE)
D   <- if (length(args) >= 1) args[1] else "./Datos"
OUT <- if (length(args) >= 2) args[2] else "./salidas_preprocesamiento"
if (!dir.exists(D)) stop("No se encuentra la carpeta de datos '", D, "' desde ", getwd(),
                         ". Ponga el directorio de trabajo en la carpeta del proyecto o pase la ruta como argumento.")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S0 <- 1981; S1 <- 2003     # años del satélite que se usan (periodo del proyecto)
# Partición del proyecto (la misma de la fase 4, PARTICION = "ferro"):
# entrenamiento 1981-1993, validación 1994-1996, prueba 1997-2003.
# Todo lo que se CALIBRA en la limpieza (umbrales, correlaciones, elección de vecinas) usa solo
# los años de entrenamiento, para que el preprocesamiento no use información de validación ni de prueba.
ENTRENA <- 1981:1993; VALIDA <- 1994:1996; PRUEBA <- 1997:2003

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

# Coordenadas reales de las estaciones (solo fuentes de la CVC; no se estima ninguna ubicación).
# Orden de prioridad si una estación aparece en varias:
#   1) CVC_Wheather_Stations_Long_Lat.csv (archivo original del proyecto, ya en lon/lat)
#   2) Catálogo_de_estaciones_hidroclimatológicas_de_la_CVC_20261120.csv (datos.gov.co)
#   3) Registros_de_precipitación,_caudal_y_nivel_de_estaciones_de_monitoreo_automático_de_la_CVC_20261120.csv
# Los archivos 2 y 3 traen coordenadas planas MAGNA-SIRGAS con comas de miles: ESTE ~ 1.000.000
# es origen Oeste (EPSG:3115) y ESTE ~ 4.600.000 es Origen Nacional CTM12 (EPSG:9377).
numero <- function(x) as.numeric(gsub(",", "", x))                  # "2,612,000,109" -> 2612000109
a_lonlat <- function(este, norte) {
  sis <- ifelse(este > 3e6, "EPSG:9377", "EPSG:3115"); out <- matrix(NA_real_, length(este), 2)
  for (s in unique(sis)) {
    i <- sis == s
    out[i, ] <- crds(project(vect(cbind(este[i], norte[i]), crs = s), "EPSG:4326"))
  }
  out
}
plano <- function(df, fuente) {                                     # una fila por estación
  xy <- a_lonlat(df$este, df$norte)
  df |> mutate(Lon = xy[, 1], lat = xy[, 2]) |> group_by(CODIGO) |>
    summarise(STATION_NAME = first(STATION_NAME), Lon = mean(Lon), lat = mean(lat), .groups = "drop") |>
    mutate(fuente_coord = fuente)
}

# El archivo original viene en latin-1 (tildes); se lee como bytes latin1 para no cortar filas
ll  <- read.csv(file.path(D, "CVC_Wheather_Stations_Long_Lat.csv"), encoding = "latin1") |>
  mutate(CODIGO = as.numeric(CODIGO))
stopifnot(nrow(ll) == 238)
# Los dos archivos de la CVC se buscan por el inicio de su nombre largo (el de la descarga de
# datos.gov.co), sin escribir las tildes, para que no falle por la codificación en Windows.
buscar_archivo <- function(patron) {
  f <- list.files(D, pattern = patron, full.names = TRUE)
  if (length(f) != 1) stop("No se encontró un único archivo en ", D, " que coincida con: ", patron)
  f
}
cvc_cat <- read.csv(buscar_archivo("^Cat.*logo_de_estaciones_hidroclim.*\\.csv$"), encoding = "UTF-8",
                    colClasses = "character", check.names = FALSE)
names(cvc_cat) <- trimws(gsub("\\s+", " ", names(cvc_cat)))
cvc_cat <- cvc_cat |> transmute(CODIGO = numero(`CODIGO ESTACION`), STATION_NAME = `NOMBRE ESTACION`,
                                este = numero(`ESTE X`), norte = numero(`NORTE Y`)) |>
  distinct() |> plano("catálogo CVC")
cvc_aut <- read.csv(buscar_archivo("^Registros_de_precipitaci.*autom.*tico.*\\.csv$"), encoding = "UTF-8",
                    colClasses = "character", check.names = FALSE) |>
  transmute(CODIGO = numero(`CODIGO ESTACION`), STATION_NAME = `NOMBRE ESTACION`,
            este = numero(`COORDENADA X`), norte = numero(`COORDENADA Y`)) |>
  distinct(CODIGO, .keep_all = TRUE) |> plano("estaciones automáticas CVC")

coords <- ll |> transmute(CODIGO, STATION_NAME, Lon, lat, fuente_coord = "archivo original")
coords <- bind_rows(coords, cvc_cat |> filter(!CODIGO %in% coords$CODIGO))
coords <- bind_rows(coords, cvc_aut |> filter(!CODIGO %in% coords$CODIGO))
# 4) coordenadas_encontradas_pendientes.csv: estaciones halladas en documentos PDF de la CVC
#    (búsqueda "código" filetype:pdf), ya en lon/lat; trae la altitud del documento (elev_tesis)
pdf_c <- read.csv(file.path(D, "coordenadas_encontradas_pendientes.csv")) |>
  transmute(CODIGO = as.numeric(STATION_ID), STATION_NAME = NA_character_, Lon = lon, lat = lat,
            fuente_coord = "documentos PDF", alt_pdf = elev_tesis)
coords <- bind_rows(coords, pdf_c |> select(-alt_pdf) |> filter(!CODIGO %in% coords$CODIGO))
print(count(coords, fuente_coord))
# Catálogo de altitud: trae altitud, cuenca y municipio (algunas estaciones no están en el catálogo de coordenadas)
cat_alt <- read_excel(file.path(D, "CatalogoEstacionesElevacion.xlsx")) |>
  mutate(codigo = as.numeric(codigo)) |> distinct(codigo, .keep_all = TRUE)
# Las estaciones halladas en documentos PDF traen su altitud: se agrega como altitud de catálogo
cat_alt <- bind_rows(cat_alt, pdf_c |> filter(!CODIGO %in% cat_alt$codigo) |> transmute(codigo = CODIGO, alt = alt_pdf))

# Máscara del Valle del Cauca: la grilla de altitud (misma malla del satélite) solo tiene valores
# dentro del departamento; fuera del Valle es NA.
dem   <- rast(file.path(D, "altitud_valle.tif"))                     # altitud de la grilla
valle <- ifel(is.na(dem), NA, 1)                                     # máscara: 1 = Valle, NA = fuera
borde <- as.data.frame(geom(as.lines(as.polygons(valle))))           # borde del Valle (para los mapas)

Y0 <- min(st$YEAR); Y1 <- max(st$YEAR)
estaciones <- sort(unique(st$STATION_ID))

# ------------------------------------------------- 2. Ubicación y altitud de TODAS las estaciones
# Solo se usan estaciones con coordenadas reales (fuentes de la CVC, paso 1). Las que no aparecen en
# ninguna fuente se quitan: su ubicación no se estima.
# Altitud: la de catálogo si existe; si no, la de la grilla en la posición real (corregida abajo).
# Estaciones del Valle cuya posición real cae en un píxel de borde (costa o cordillera) que la grilla
# de altitud (0,05°) no cubre: su municipio es del Valle, así que se incluyen igual y su píxel del
# satélite se añade al análisis (sección 7).
EN_VALLE_BORDE <- c(5300000202,   # Colpuertos, Buenaventura (HIMAT)
                    5301000104,   # San Francisco, Buenaventura
                    2612800103,   # Bella Vista, Palmira
                    2612720106,   # Bolivia, Pradera
                    2612710104,   # Los Andes, Florida
                    2615400110,   # El Alambrado, Sevilla
                    2615400103,   # Santa Ana, Cartago
                    2620000111)   # Tablanca, Jamundí
km_entre <- function(lo1, la1, lo2, la2) sqrt(((lo1 - lo2) * 111 * cos(3.8 * pi / 180))^2 + ((la1 - la2) * 111)^2)
sin_coord <- setdiff(estaciones, coords$CODIGO)
base::cat("Estaciones con datos:", length(estaciones), "| con coordenadas reales:",
          length(estaciones) - length(sin_coord), "| sin coordenadas (se quitan):", length(sin_coord), "\n")
ubic <- tibble(STATION_ID = setdiff(estaciones, sin_coord)) |>
  left_join(coords |> select(STATION_ID = CODIGO, STATION_NAME, Lon, lat, fuente_coord), by = "STATION_ID")

# Altitud: catálogo si existe; si no, grilla interpolada (bilineal) y corregida con una regresión
# ajustada en las estaciones que tienen catálogo y coordenadas (la grilla sobreestima en montaña)
xy_u <- as.matrix(ubic[, c("Lon", "lat")]); xy_u[is.na(xy_u)] <- 0
alt_bil <- terra::extract(dem, xy_u, method = "bilinear")[, 1]
alt_cel <- terra::extract(dem, xy_u)[, 1]
ubic <- ubic |> mutate(alt_grilla = ifelse(is.na(Lon), NA, coalesce(alt_bil, alt_cel)),
                       alt_catalogo = cat_alt$alt[match(STATION_ID, cat_alt$codigo)])
ajuste <- lm(alt_catalogo ~ alt_grilla, data = ubic)
ubic <- ubic |> mutate(alt_corregida = predict(ajuste, newdata = ubic),
                       altitud = coalesce(alt_catalogo, alt_corregida),
                       fuente_alt = case_when(!is.na(alt_catalogo) ~ "catálogo",
                                              !is.na(alt_corregida) ~ "grilla corregida",
                                              TRUE ~ "sin dato"),
                       en_valle = (!is.na(terra::extract(valle, xy_u)[, 1]) | STATION_ID %in% EN_VALLE_BORDE) & !is.na(Lon)) |>
  select(STATION_ID, STATION_NAME, Lon, lat, fuente_coord, altitud, fuente_alt, en_valle)

print(count(ubic, fuente_coord))
print(count(ubic, fuente_alt))

# ------------------------------------------------- 3. Solo estaciones dentro del Valle del Cauca
# Se usan todas las estaciones con coordenadas como referencia para ubicar (paso 2); desde aquí
# solo se analizan las que quedan dentro del Valle.
base::cat("\nEstaciones con datos:", nrow(ubic),
          "| dentro del Valle:", sum(ubic$en_valle),
          "| fuera del Valle:", sum(!ubic$en_valle), "\n")
ubic <- ubic |> filter(en_valle)
# Tabla de verificación: todas las estaciones del Valle con su ubicación y de dónde salió
print(ubic |> select(STATION_ID, STATION_NAME, Lon, lat, fuente_coord, altitud, fuente_alt),
      n = Inf, width = Inf)
base::cat("Estaciones del Valle sin altitud:", sum(is.na(ubic$altitud)), "\n")

foto(st, "0. Datos crudos (todas las estaciones)")
st_crudo <- st |> select(STATION_ID, YEAR, DAY, FECHA, VALUE)   # para la Tabla 1
st <- st |> filter(STATION_ID %in% ubic$STATION_ID)
foto(st, "1. Solo estaciones dentro del Valle")
estaciones <- sort(unique(st$STATION_ID))

# ------------------------------------------------- 4. Registros no diarios
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

# Figura 1. Mapa de calor estación x año (1930-2003)
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

# Figura 2. Estaciones que reportan cada año
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
# La partición del proyecto empieza en 1981 (train 1981-1993 / val 1994-1996 / test 1997-2003).
n_antes <- n_distinct(st$STATION_ID)
st <- st |> filter(YEAR >= S0, YEAR <= S1)
Y0 <- S0; Y1 <- S1
estaciones <- sort(unique(st$STATION_ID))
ubic <- ubic |> filter(STATION_ID %in% estaciones)
foto(st, "2. Corte 1981-2003")
base::cat("\nCorte 1981-2003: estaciones del Valle", n_antes, "->", length(estaciones),
          "(", n_antes - length(estaciones), "dejaron de medir antes de 1981 )\n")

# Tabla 1. Qué contienen las estaciones al inicio (antes de limpiar): todas vs. Valle 1981-2003
fmt <- function(x, d = 0) format(round(x, d), big.mark = ".", decimal.mark = ",", nsmall = d,
                                 scientific = FALSE, trim = TRUE)
inventario <- function(x) {
  ids  <- unique(x$STATION_ID)
  vida <- x |> group_by(STATION_ID) |>
    summarise(dias_vida = as.numeric(max(FECHA) - min(FECHA)) + 1, n = n(), .groups = "drop")
  no_diaria <- rachas_separadas(x) |> filter(en_racha) |> distinct(STATION_ID)
  con_coord <- ids %in% coords$CODIGO; con_alt <- ids %in% cat_alt$codigo
  c(fmt(length(ids)), paste0(min(x$YEAR), "-", max(x$YEAR)), fmt(nrow(x)),
    fmt(sum(vida$dias_vida - vida$n)), paste0(fmt(100 * sum(vida$n) / sum(vida$dias_vida), 1), "%"),
    paste0(fmt(sum(con_coord)), " / ", fmt(sum(!con_coord))),
    paste0(fmt(sum(con_alt)), " / ", fmt(sum(!con_alt))), fmt(sum(!con_coord & !con_alt)),
    fmt(nrow(no_diaria)), paste0(fmt(100 * mean(x$VALUE == 0), 1), "%"),
    paste(fmt(min(x$VALUE)), "/", fmt(max(x$VALUE))))
}
tabla1 <- tibble(
  Aspecto = c("Estaciones con datos", "Periodo", "Días con dato", "Días faltantes", "Cobertura de días",
              "Con / sin coordenadas", "Con / sin altitud", "Sin coordenadas, ni altitud",
              "Con periodos no diarios (semanal/mensual)", "Días en 0 mm", "Mínimo / máximo (mm/día)"),
  Datos = inventario(st_crudo), `1981-2003` = inventario(st))
print(tabla1, width = Inf)
guardar_tabla(tabla1, "tabla1_estadisticas_basicas")
rm(st_crudo)

# ------------------------------------------------- 7. Satélite (1981-2003), solo píxeles del Valle
fs  <- list.files(file.path(D, "ImagenesSatelitePentadal"), pattern = "^ppt\\d{6}\\.tif$", full.names = TRUE)
inf <- tibble(f = fs, anio = as.integer(substr(basename(fs), 4, 7)),
              pentada = as.integer(substr(basename(fs), 8, 9))) |>
  filter(anio >= S0, anio <= S1) |> arrange(anio, pentada)
sat_r <- rast(inf$f)                       # una capa por péntada; el NA (-9999) lo lee terra
cubo  <- as.array(sat_r)                   # filas x columnas x péntadas
fuera <- is.na(as.array(valle)[, , 1])     # píxeles fuera del Valle
costa <- ubic |> filter(STATION_ID %in% EN_VALLE_BORDE)  # píxeles de borde de estaciones del Valle
if (nrow(costa)) {
  cc <- cellFromXY(sat_r, as.matrix(costa[, c("Lon", "lat")]))
  fuera[cbind(rowFromCell(sat_r, cc), colFromCell(sat_r, cc))] <- FALSE
}
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

# Tabla 2. Verificación del satélite: formato y tope de 200 mm/péntada en los píxeles del Valle
pct_px <- apply(cubo == vmax, c(1, 2), mean, na.rm = TRUE) * 100
pct_px <- pct_px[!is.na(validos_px) & validos_px > 0]
cumple <- function(n) ifelse(n == 0, "Cumple", paste(n, "no cumplen"))
tabla2 <- tibble(
  Verificacion = c("72 péntadas por año, 1981-2003", "Misma grilla (39 × 37, 0,05°, WGS84)",
                   "Una sola capa por archivo", "mm/péntada, no negativo y ≤ 200", "Píxeles del Valle / péntadas",
                   paste("Péntadas en el tope de", vmax, "mm"), "Píxeles que alguna vez tocan el tope",
                   "Píxel más saturado"),
  Resultado = c(cumple(pent_faltan), cumple(sum(fmt_sat$grilla_distinta)), cumple(sum(fmt_sat$varias_capas)),
                cumple(sum(fmt_sat$valor_malo)), paste(length(pct_px), "/", fmt(dim(cubo)[3])),
                paste0(fmt(sum(cubo == vmax, na.rm = TRUE)), " (", fmt(100 * mean(cubo == vmax, na.rm = TRUE), 2), " %)"),
                paste(sum(pct_px > 0), "de", length(pct_px)), paste0(fmt(max(pct_px), 1), " % de sus péntadas")))
print(tabla2, width = Inf)
guardar_tabla(tabla2, "tabla2_verificacion_satelite")
rm(fmt_sat)

# ------------------------------------------------- 8. Limpieza 1: registros acumulados (semanales, mensuales...)
# Algunas estaciones, en ciertos periodos, no reportaron día a día sino un valor cada semana, cada
# mes o cada pocos días. Si ese valor es el ACUMULADO desde el registro anterior, es información
# real: no se bota, se reparte en días (desagregación) con el patrón de la vecina más parecida.
#   - Periodo sospechoso: 6 o más registros seguidos separados por 2 o más días (ver sección 4).
#     Ventana de cada registro: los días desde el registro anterior hasta él.
#   - Donante: estación a <= 25 km, con >= 90 % de los días del periodo y la mayor correlación de
#     lluvia diaria con la estación (calculada con 1981-1993, sin los periodos sospechosos y con
#     >= 1 año de días en común; si ninguna cumple, la más cercana).
#   - ¿Acumulado o días sueltos? Se suma lo que reportó la estación (V) y se compara con lo que
#     llovió en la donante en las ventanas completas (Sw) y solo en los días con registro (S1).
#     Es acumulado si V se parece más a Sw que a S1 y además V >= 1,5 x S1. Si no, son días sueltos
#     con faltantes en medio y se dejan como están.
#     La decisión se toma POR ESTACIÓN (sumando todos sus periodos), no periodo por periodo:
#     en una prueba con 500 casos simulados (solo con años 1981-1993) de acumulados y de días sueltos
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
# Estaciones con reporte mensual (un registro cada ~31 días): su dato original, para la Figura 8
mensuales <- traza |> filter(dias_entre_registros >= 28)
crudo_mens <- st |> filter(STATION_ID %in% c(mensuales$STATION_ID, traza$STATION_ID[traza$resultado == "acumulado: desagregado"])) |>
  select(STATION_ID, FECHA, YEAR, VALUE)                      # dato original, para la Figura 8
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
#   - Confirmado: alguna vecina registra al menos el 30 % de ese valor. Calibración (1981-1993): en
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
#     1981-1993, >= 1 año en común),
#     con datos en al menos el 90 % de los días de la racha.
#   - Cero falso: al menos la mitad de esas vecinas registró >= 50 mm por cada 30 días de racha.
#     Calibración (ventanas de 30 días, 1981-1993): si una vecina parecida recibe 50-100 mm, la
#     estación queda toda en 0 solo el 0,07 % de las veces; con 0 mm en la vecina, el 16 %.
#   - Sin vecinas parecidas se usa el satélite en el píxel de la estación, que es menos preciso:
#     falso solo si el píxel tiene >= 150 mm por cada 30 días (calibración 1981-1993: 0,07 % en 0).
#   - Sequía real: se deja. Cero falso: la racha pasa a faltante (para imputar después).
MIN_RACHA_CERO <- 15; CORR_VECINA <- 0.4

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
#   - Cada estación tiene su razón típica R (mediana de r en 1981-1993): una estación de la costa
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
base::cat("Limpieza 4 - calibración (1981-1993): meses comparables", sum(!is.na(mt$max_vecina_mm) & mt$YEAR %in% ENTRENA),
          "| inconsistentes", sum(mt$inconsistente & mt$YEAR %in% ENTRENA), "(",
          round(100 * mean(mt$inconsistente[!is.na(mt$max_vecina_mm) & mt$YEAR %in% ENTRENA]), 2), "% )\n")
est_mala <- R_tipica |> filter(R > R_MAX_ESTACION)
meses_malos <- mt |> filter(inconsistente | STATION_ID %in% est_mala$STATION_ID) |>
  mutate(decision = ifelse(STATION_ID %in% est_mala$STATION_ID, "estación completa: se quita",
                           "mes inconsistente: se quita"))
print(meses_malos |> count(decision))
print(est_mala)
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
#     1,4 (su mediana de 1981-1993) a 2,8;
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

# Ranking por estación: cobertura sobre su vida útil (de su primer a su último año con dato) y
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

# ------------------------------------------------- 13. Criterio de cobertura por estación-año
# Un año de una estación entra al modelo si tiene >= 90 % de sus días con dato y ningún hueco de
# más de 31 días seguidos. La curva anual se construye con el año completo: imputar un año con
# medio año vacío sería inventar la curva, y un hueco de más de un mes puede borrar uno de los dos
# picos de lluvia. Los años que no cumplen se quitan; los huecos de los que sí cumplen se imputan
# después.
COB_MIN <- 90; HUECO_MAX <- 31
pares_ok <- ea |> filter(cob >= COB_MIN, hueco <= HUECO_MAX) |> select(STATION_ID, YEAR)

# Tabla 10. Sensibilidad del criterio: cuántas curvas quedan (por bloque de la partición) y cuántos
# días habría que imputar con otros umbrales de cobertura y de hueco máximo. Más laxo = más curvas pero
# más días inventados y huecos que pueden borrar un pico; más estricto = menos curvas y estaciones.
# Solo cuenta pares; no ajusta nada con validación ni prueba (los umbrales se fijan antes de mirar el modelo).
bloque_de <- function(y) ifelse(y %in% ENTRENA, "entrenamiento", ifelse(y %in% VALIDA, "validacion", "prueba"))
sens <- expand_grid(cobertura_min_pct = c(80, 85, 90, 95, 100), hueco_max_dias = c(7, 15, 31, 61, 366)) |>
  rowwise() |>
  mutate(ok = list(ea |> filter(n_dias > 0, cob >= cobertura_min_pct, hueco <= hueco_max_dias)),
         curvas = nrow(ok), estaciones = n_distinct(ok$STATION_ID),
         curvas_entrenamiento = sum(bloque_de(ok$YEAR) == "entrenamiento"),
         curvas_validacion = sum(bloque_de(ok$YEAR) == "validacion"),
         curvas_prueba = sum(bloque_de(ok$YEAR) == "prueba"),
         dias_por_imputar = sum(ok$dias - ok$n_dias),
         pct_dias_imputados = round(100 * dias_por_imputar / sum(ok$dias), 2),
         max_dias_imputados_en_un_anio = if (curvas > 0) max(ok$dias - ok$n_dias) else NA_integer_,
         criterio_usado = cobertura_min_pct == COB_MIN & hueco_max_dias == HUECO_MAX) |>
  ungroup() |> select(-ok) |>
  mutate(hueco_max_dias = ifelse(hueco_max_dias == 366, "sin límite", as.character(hueco_max_dias)))
print(sens, n = Inf, width = Inf)
guardar_tabla(sens, "tabla10_sensibilidad_cobertura")

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

# ------------------------------------------------- 14. Imputación de los días faltantes
# Solo se imputan los días faltantes de los años que pasaron el criterio de cobertura.
# Se comparan dos métodos que usan las estaciones vecinas (<= 50 km) que sí midieron ese día:
#   - CCW: promedio de las 8 vecinas más correlacionadas con peso r^2.
#   - IDW_altitud: promedio de las 8 vecinas más cercanas con peso 1/d^2 (IDW), pero (1) solo vecinas de la misma vertiente (código 26 = interior, cuenca del río Cauca,
#     53/54 = Pacífico; ambas dentro del Valle; si no hay ninguna, se usan todas) y (2) con una distancia "efectiva" que
#     suma la diferencia de altitud: d_ef = sqrt(d^2 + (LAMBDA * dif_altitud / 100)^2), es decir,
#     cada 100 m de diferencia de altitud cuentan como LAMBDA km. LAMBDA se calibra con días de
#     1981-1993 distintos de los de la validación.
# Todos trabajan sobre la lluvia RELATIVA (lluvia / promedio diario de la estación en ese mes,
# calculado con 1981-1993) para no mezclar sitios secos y húmedos (orografía), y luego se pasa a mm
# con el promedio de la estación a imputar. Los dos métodos primero deciden si llovió (mayoría
# ponderada de vecinas con lluvia) para no llenar los días secos de llovizna.
# Si la estación no tiene años de entrenamiento, su promedio mensual sale del satélite en su píxel.
# Temporadas: la lluvia del Valle es bimodal y las vecinas se parecen distinto en época lluviosa y
# seca. Un mes es "Lluviosa" si su lluvia media diaria en 1981-1993 (promedio de todas las estaciones
# del Valle) supera el promedio de los 12 meses; si no, "Seca". Las correlaciones (CCW) y la
# validación se hacen por separado en cada temporada.
# Validación: se esconden 30.000 días conocidos de 1981-1993, se imputan y se comparan con lo medido.
# Se usa IDW_altitud en las dos temporadas (ver la justificación junto a la Tabla 7). Sin vecinas: satélite del píxel.
RADIO_IMP <- 50; K_VEC <- 8; N_VALIDA <- 30000
set.seed(SEMILLA)
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

# Temporada de cada mes (solo con 1981-1993) y de cada día
lluvia_mes <- tapply(rowMeans(M[tr, ], na.rm = TRUE), mes_f[tr], mean, na.rm = TRUE)[as.character(1:12)]
TEMP_MES <- ifelse(lluvia_mes > mean(lluvia_mes), "Lluviosa", "Seca")
TEMPS <- c("Lluviosa", "Seca")
tem_f <- TEMP_MES[mes_f]
print(tibble(mes = 1:12, lluvia_media_mm_dia = round(lluvia_mes, 2), temporada = TEMP_MES))

# Correlaciones entre estaciones, por temporada (solo días de entrenamiento)
corr_t <- list()
for (tp in TEMPS) {
  dt <- tr & tem_f == tp
  P <- !is.na(M[dt, ]) * 1; comun_m <- crossprod(P)
  cm <- suppressWarnings(cor(M[dt, ], use = "pairwise.complete.obs"))
  cm[comun_m < 180 | is.na(cm)] <- NA; diag(cm) <- NA            # al menos 180 días en común
  corr_t[[tp]] <- cm
}

kge <- function(s, o) { r <- cor(s, o); 1 - sqrt((r - 1)^2 + (sd(s) / sd(o) - 1)^2 + (mean(s) / mean(o) - 1)^2) }
# Altitud y vertiente de cada estación de la matriz
alt_m  <- ubic$altitud[match(ids_m, ubic$STATION_ID)]
vert_m <- ifelse(substr(ids_m, 1, 2) %in% c("53", "54"), "Pacífico", "Interior (río Cauca)")
dz_m   <- abs(outer(alt_m, alt_m, "-"))

idw_alt <- function(j, t, o, lambda) {
  ov <- o[vert_m[o] == vert_m[j]]; if (length(ov)) o <- ov
  de <- sqrt(dist_m[j, o]^2 + (lambda * dz_m[j, o] / 100)^2)
  k <- order(de)[seq_len(min(K_VEC, length(o)))]; n3 <- o[k]
  w <- 1 / pmax(de[k], .5)^2; w <- w / sum(w)
  if (sum(w * (M[t, n3] > 0)) >= .5) sum(w * REL[t, n3]) * clim[mes_f[t], j] else 0
}

# Calibración de LAMBDA (km por cada 100 m): 5.000 días conocidos de 1981-1993 (no se usan en la validación)
conocidos <- which(!is.na(M) & tr, arr.ind = TRUE)
id_cal <- sample(nrow(conocidos), 5000); cal_d <- conocidos[id_cal, ]
LAMBDAS <- c(0, 0.5, 1, 2, 5, 10)
kge_lambda <- sapply(LAMBDAS, function(l) {
  s <- apply(cal_d, 1, function(x) {
    o <- which(!is.na(M[x[1], ])); o <- o[o != x[2] & dist_m[x[2], o] <= RADIO_IMP]
    if (length(o)) idw_alt(x[2], x[1], o, l) else NA })
  ok <- !is.na(s); kge(s[ok], M[cal_d][ok])
})
LAMBDA <- LAMBDAS[which.max(kge_lambda)]
print(tibble(lambda_km_por_100m = LAMBDAS, KGE = round(kge_lambda, 3)))
base::cat("LAMBDA elegido:", LAMBDA, "km por cada 100 m de diferencia de altitud\n")

imputa <- function(j, t) {
  o <- which(!is.na(M[t, ])); o <- o[o != j & dist_m[j, o] <= RADIO_IMP]
  out <- c(CCW = NA, IDW_altitud = NA)
  if (!length(o)) return(out)
  ct <- clim[mes_f[t], j]; a <- REL[t, ]; llovio <- M[t, ] > 0
  corr_m <- corr_t[[tem_f[t]]]
  n1 <- o[order(dist_m[j, o])][seq_len(min(K_VEC, length(o)))]
  w <- 1 / pmax(dist_m[j, n1], .5)^2; w <- w / sum(w)
  idw <- if (sum(w * llovio[n1]) >= .5) sum(w * a[n1]) * ct else 0   # respaldo de CCW
  rc <- corr_m[j, o]; okc <- which(!is.na(rc) & rc > 0)
  if (length(okc)) {
    oc <- o[okc][order(-rc[okc])]
    n2 <- oc[seq_len(min(K_VEC, length(oc)))]; w2 <- corr_m[j, n2]^2; w2 <- w2 / sum(w2)
    out["CCW"] <- if (sum(w2 * llovio[n2]) >= .5) sum(w2 * a[n2]) * ct else 0
  } else out["CCW"] <- idw
  out["IDW_altitud"] <- idw_alt(j, t, o, LAMBDA)
  out
}

# Validación con días conocidos de entrenamiento
muestra <- conocidos[sample(setdiff(seq_len(nrow(conocidos)), id_cal), N_VALIDA), ]
pred <- t(apply(muestra, 1, function(x) imputa(x[2], x[1])))
medido <- M[muestra]; okv <- complete.cases(pred); tem_v <- tem_f[muestra[, 1]]
comparacion <- bind_rows(lapply(TEMPS, function(tp) {
  ok_t <- okv & tem_v == tp
  boot <- replicate(300, { b <- sample(which(ok_t), sum(ok_t), TRUE); sapply(colnames(pred), function(m) kge(pred[b, m], medido[b])) })
  bind_rows(lapply(colnames(pred), function(m) {
    s <- pred[ok_t, m]; o <- medido[ok_t]
    tibble(temporada = tp, dias_validados = sum(ok_t), metodo = m, MAE_mm = mean(abs(s - o)), RMSE_mm = sqrt(mean((s - o)^2)),
           sesgo_pct = 100 * (sum(s) - sum(o)) / sum(o), correlacion = cor(s, o), KGE = kge(s, o),
           KGE_IC95_inf = quantile(boot[m, ], .025), KGE_IC95_sup = quantile(boot[m, ], .975),
           ratio_desv_estandar = sd(s) / sd(o), pct_dias_secos_medido = 100 * mean(o == 0),
           pct_dias_secos_imputado = 100 * mean(s == 0), acierto_lluvia_o_seco_pct = 100 * mean((s > 0) == (o > 0)),
           p95_medido = quantile(o, .95), p95_imputado = quantile(s, .95))
  }))
})) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
print(comparacion, width = Inf)
guardar_tabla(comparacion, "tabla7_resultados_imputacion")             # Tabla 7
# Método elegido: IDW_altitud en las dos temporadas. En validación su KGE empata con CCW
# (intervalos de confianza solapados) y su MAE es apenas mayor que el de CCW. Se prefiere porque
# tiene en cuenta explícitamente la vertiente (solo vecinas del mismo lado de la cordillera) y la
# diferencia de altitud, y porque es un solo método para todas las estaciones: no depende de que
# existan correlaciones (CCW necesita IDW de respaldo cuando no las hay).
MEJOR <- c(Lluviosa = "IDW_altitud", Seca = "IDW_altitud")
mayor_kge <- sapply(TEMPS, function(tp) { x <- comparacion |> filter(temporada == tp); x$metodo[which.max(x$KGE)] })
base::cat("Mayor KGE | Lluviosa:", mayor_kge["Lluviosa"], "| Seca:", mayor_kge["Seca"],
          "|| Método usado: IDW_altitud en ambas temporadas\n")
# KGE de todo el año con el método elegido en cada temporada (comparable con un solo método para todo el año)
s_comb <- pred[cbind(which(okv), match(MEJOR[tem_v[okv]], colnames(pred)))]
base::cat("KGE de todo el año | método usado:", round(kge(s_comb, medido[okv]), 3),
          "|", paste(colnames(pred), round(apply(pred[okv, ], 2, kge, o = medido[okv]), 3), collapse = " | "), "\n")

# Imputación de los faltantes con el método elegido
faltan <- cal |> filter(!presente) |> select(STATION_ID, FECHA)
fi <- cbind(match(faltan$FECHA, fechas), match(faltan$STATION_ID, ids_m))
met <- MEJOR[tem_f[fi[, 1]]]                                       # método de la temporada de cada día
val <- sapply(seq_len(nrow(fi)), function(i) imputa(fi[i, 2], fi[i, 1])[met[i]])
sin_vec <- is.na(val)
val[sin_vec] <- (sat_dia[fi[sin_vec, , drop = FALSE]] / clim_sat[cbind(mes_f[fi[sin_vec, 1]], fi[sin_vec, 2])]) *
  clim[cbind(mes_f[fi[sin_vec, 1]], fi[sin_vec, 2])]
imputados <- faltan |> mutate(VALUE = round(val, 1), YEAR = year(FECHA), MONTH = month(FECHA), DAY = day(FECHA),
                              origen = ifelse(sin_vec, "imputado (satélite)", paste0("imputado (", met, ")")))
base::cat("Imputados:", nrow(imputados), "días |", "con vecinas:", sum(!sin_vec), "| con satélite:", sum(sin_vec), "\n")

st_antes_imp <- st
st <- bind_rows(st, imputados) |> arrange(STATION_ID, FECHA)
stopifnot(!any(duplicated(st[, c("STATION_ID", "FECHA")])))
cal <- cal |> left_join(imputados |> select(STATION_ID, FECHA, imputado = VALUE), by = c("STATION_ID", "FECHA")) |>
  mutate(VALUE = coalesce(VALUE, imputado), imputado = !is.na(imputado)) # calendario después (presente = antes)

foto(st, "8. Imputación de faltantes", cal |> mutate(presente = TRUE))

# Tabla 9. Estadísticas básicas de la lluvia diaria antes y después de imputar (mismas 3.428 curvas
# estación-año). "Solo imputados" muestra los días que llenó el método, para ver si se parecen a los medidos.
resumen_lluvia <- function(x, etiqueta) {
  anual <- x |> group_by(STATION_ID, YEAR) |> summarise(t = sum(VALUE), .groups = "drop")
  tibble(Datos = etiqueta, `Días` = fmt(nrow(x)),
         `Días faltantes` = fmt(if (etiqueta == "Antes de imputar") nrow(imputados) else 0),
         `Media (mm/día)` = fmt(mean(x$VALUE), 2), `Desv. estándar (mm/día)` = fmt(sd(x$VALUE), 2),
         `Mediana (mm/día)` = fmt(median(x$VALUE), 1), `Percentil 95 (mm/día)` = fmt(quantile(x$VALUE, .95), 1),
         `Máximo (mm/día)` = fmt(max(x$VALUE), 1), `Días secos (%)` = fmt(100 * mean(x$VALUE == 0), 1),
         `Lluvia anual media (mm)` = if (etiqueta == "Solo días imputados") "-" else fmt(mean(anual$t)))
}
tabla9 <- bind_rows(resumen_lluvia(st_antes_imp, "Antes de imputar"),
                    resumen_lluvia(st, "Después de imputar"),
                    resumen_lluvia(imputados, "Solo días imputados"))
print(tabla9, width = Inf)
guardar_tabla(tabla9, "tabla9_antes_despues_imputacion")

# ------------------------------------------------- 15. Datos limpios: diario y por péntadas (CSV)
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

# ------------------------------------------------- 16. Tablas 3, 4, 5 y 8 de la presentación
# Tabla 3. Pasos del preprocesamiento y su resultado
n_crudo <- traz[[1]]$estaciones; n_corte <- traz[[3]]$estaciones
tabla3 <- tibble(
  Paso = c("1. Ubicación", "2. Valle y 1981-2003", "3. Acumulados", "4. Días > 300 mm", "5. Ceros falsos",
           "6. Meses imposibles", "7. Cobertura"),
  `Qué se hizo` = c("Coordenadas reales de la CVC; sin coordenada se quita",
                    "Solo el Valle y el periodo del satélite",
                    "Totales semanales o mensuales repartidos en días",
                    "Contraste con las estaciones vecinas",
                    "Rachas de 0 mm desmentidas por las vecinas",
                    "Totales > 600 mm y más del doble de la vecina",
                    "Año con ≥ 90 % de días y sin huecos > 31 días"),
  Resultado = c(paste(length(sin_coord), "estaciones sin coordenadas quitadas"),
                paste(n_crudo, "→", n_corte, "estaciones"),
                paste(sum(traza$resultado == "acumulado: desagregado"), "periodos corregidos"),
                paste(sum(altos$decision != "confirmado por vecinas"), "días quitados"),
                paste(fmt(nrow(falsos)), "días quitados"),
                paste(sum(meses_malos$decision == "mes inconsistente: se quita"), "meses y", nrow(est_mala),
                      "estaciones quitados"),
                paste0(fmt(nrow(pares_ok)), " estación-año, ", n_distinct(pares_ok$STATION_ID), " estaciones")))
print(tabla3, width = Inf)
guardar_tabla(tabla3, "tabla3_preprocesamiento")

# Tabla 4. Ceros reales vs. falsos por temporada (cada racha va al mes de su día central)
# Temporadas: las mismas de la imputación (TEMP_MES, calculadas con 1981-1993).
MESES <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")
temporada <- function(m) factor(TEMP_MES[m], levels = c("Seca", "Lluviosa"))   # misma definición de la imputación
r0 <- rachas0 |> mutate(MONTH = month(desde + floor((hasta - desde) / 2)), falso = decision == "cero falso: se quita")
resumen_ceros <- function(x) summarise(x, `Rachas revisadas` = n(), `Sequía real` = sum(!falso),
                                       `Cero falso` = sum(falso), `% falsos` = round(100 * mean(falso), 1),
                                       `Días quitados` = sum(dias[falso]))
tabla4 <- bind_rows(r0 |> mutate(Temporada = as.character(temporada(MONTH))) |> group_by(Temporada) |>
                      resumen_ceros() |> arrange(desc(Temporada)),
                    r0 |> resumen_ceros() |> mutate(Temporada = "Total"))
print(tabla4)
guardar_tabla(tabla4, "tabla4_ceros_por_temporada")

# Tabla 5. Estaciones con reporte mensual: su valor se parece al total del MES de la vecina
# (Sw), no a la lluvia de la vecina ese mismo día (S1); por eso es un total mensual.
tabla5 <- mensuales |> transmute(Estación = STATION_ID, Desde = desde, Hasta = hasta, Registros = registros,
                                 `Total estación (mm)` = V_estacion, `Vecina, todo el mes (mm)` = Sw_donante_ventanas,
                                 `Vecina, mismo día (mm)` = S1_donante_mismos_dias,
                                 Resultado = ifelse(resultado == "acumulado: desagregado", "Repartido en días",
                                                    "Quitado (sin información)"))
print(tabla5, width = Inf)
guardar_tabla(tabla5, "tabla5_estaciones_mensuales")

# Tabla 8. Resultado final
tabla8 <- tibble(
  `Resultado final` = c("Estaciones", "Curvas estación-año completas", "Días con dato", "Días faltantes",
                        "Péntadas (escala del satélite)", "Entrenamiento 1981-1993", "Validación 1994-1996",
                        "Prueba 1997-2003"),
  Valor = c(paste(n_distinct(diario$STATION_ID), "de", n_crudo),
            fmt(nrow(distinct(diario, STATION_ID, YEAR))), fmt(nrow(diario)), fmt(sum(is.na(diario$VALUE))),
            paste(fmt(nrow(pentadas)), "(72 por año), cada una con su píxel del satélite"),
            sapply(list(ENTRENA, VALIDA, PRUEBA), function(a) {
              x <- distinct(diario, STATION_ID, YEAR) |> filter(YEAR %in% a)
              paste(fmt(nrow(x)), "curvas,", n_distinct(x$STATION_ID), "estaciones") })))
print(tabla8, width = Inf)
guardar_tabla(tabla8, "tabla8_resultado_final")

# ================================================================== 17. FIGURAS (solo se imprimen)
eje_anios <- scale_x_continuous(breaks = seq(S0, S1, 2))

# Figura 3. Número de estaciones según su cobertura, año por año (un panel por año)
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

# Figura 4. Hueco más largo en los pares estación-año con algún dato
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

# Capa común para los mapas: borde del Valle del Cauca
capa_borde <- geom_path(data = borde, aes(x, y, group = interaction(geom, part)),
                        inherit.aes = FALSE, colour = "black", linewidth = .5)

# Figura 5. Mapa del satélite (lluvia media anual) con el borde del Valle y las estaciones encima
# Cada cuadro es un píxel de 0,05° (~5,5 km). Solo píxeles del Valle. Puntos negros = estaciones.
xy <- as.data.frame(sat_r[[1]], xy = TRUE)[, 1:2]                  # centros de píxel en el orden de terra
mapa <- xy |> mutate(media_anual = as.vector(t(apply(cubo, c(1, 2), mean, na.rm = TRUE) * 72)))
ok <- mc |> filter(pixel_valido)
mapa <- mapa |> filter(!is.na(media_anual))
mostrar(
  ggplot(mapa, aes(x, y, fill = media_anual)) + geom_raster() +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "mm/año") + capa_borde +
    geom_point(data = ok, aes(Lon, lat), inherit.aes = FALSE, size = .6) + coord_equal() +
    labs(title = "Satélite: lluvia media anual 1981-2003", x = "Lon", y = "Lat")
)
# Figura 6. Estaciones del Valle y la fuente de sus coordenadas (todas reales, de la CVC)
valle_df <- as.data.frame(valle, xy = TRUE)                          # solo los píxeles del Valle
col_fuente <- c("archivo original" = "#2166ac", "catálogo CVC" = "#67a9cf", "estaciones automáticas CVC" = "#053061",
                "documentos PDF" = "#92c5de")
mostrar(
  ggplot() + geom_raster(data = valle_df, aes(x, y), fill = "grey92") + capa_borde +
    geom_point(data = ubic, aes(Lon, lat, colour = fuente_coord), size = 1.3) +
    scale_colour_manual(values = col_fuente,
                        labels = function(x) paste0(x, " (", table(ubic$fuente_coord)[x], ")"), name = NULL) +
    coord_equal() + labs(title = "Estaciones del Valle con coordenadas reales", x = "Lon", y = "Lat") +
    theme(legend.position = "bottom")
)

# Figura 7. Temporadas secas y lluviosas: lluvia media por mes y % de rachas de 15+ días en 0
# que resultaron ser ceros falsos (limpieza 3)
g_ceros <- bind_rows(
  diario |> group_by(MONTH) |> summarise(valor = mean(VALUE)) |>
    mutate(panel = "Lluvia media (mm/día)", tipo = as.character(temporada(MONTH))),
  r0 |> group_by(MONTH) |> summarise(valor = 100 * mean(falso)) |>
    mutate(panel = "% de rachas de 15+ días en 0 que eran ceros falsos", tipo = "Cero falso")) |>
  mutate(panel = factor(panel, levels = c("Lluvia media (mm/día)", "% de rachas de 15+ días en 0 que eran ceros falsos")),
         tipo = factor(tipo, levels = c("Seca", "Lluviosa", "Cero falso")))
mostrar(
  ggplot(g_ceros, aes(factor(MONTH), valor, fill = tipo)) + geom_col(width = .75) +
    geom_text(aes(label = round(valor, 1)), vjust = -.3, size = 3) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_x_discrete(labels = MESES) + scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    scale_fill_manual(values = c("Seca" = "#9ecae1", "Lluviosa" = "#2166ac", "Cero falso" = "#b2182b"), name = NULL) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 11),
          panel.grid.major.x = element_blank())
)

# Figura 8. Estación con registros acumulados: dato original y repartido en días.
# Ejemplo: el periodo desagregado con mayor total reportado cuyo año pasó el criterio de cobertura
# (para que el panel "Después" tenga datos).
ej <- traza |> filter(resultado == "acumulado: desagregado") |> mutate(YEAR = year(desde)) |>
  semi_join(diario |> distinct(STATION_ID, YEAR), by = c("STATION_ID", "YEAR")) |>
  slice_max(V_estacion, n = 1, with_ties = FALSE)
EJ <- ej$STATION_ID; ANIO_EJ <- year(ej$desde)
g_mens <- bind_rows(
  crudo_mens |> filter(STATION_ID == EJ, YEAR == ANIO_EJ) |> mutate(panel = "Antes: registros acumulados"),
  diario |> filter(STATION_ID == EJ, YEAR == ANIO_EJ) |> mutate(panel = "Después: repartido en días")) |>
  mutate(panel = factor(panel, levels = c("Antes: registros acumulados", "Después: repartido en días")))
mostrar(
  ggplot(g_mens, aes(FECHA, VALUE, fill = panel)) + geom_col(width = 1) +
    facet_wrap(~ panel, ncol = 1) +
    scale_x_date(date_breaks = "1 month", labels = function(x) MESES[month(x)]) +
    scale_fill_manual(values = c("#b2182b", "#2166ac"), guide = "none") +
    labs(x = NULL, y = "mm/día") +
    theme(strip.text = element_text(face = "bold", size = 11))
)

# Figura 9. Serie de tiempo de la lluvia mensual, antes y después de imputar
# Lluvia media diaria de cada mes sobre todas las estaciones (antes: solo días medidos; después:
# medidos + imputados). Si la imputación es buena, las dos series casi son iguales.
serie_ad <- bind_rows(st_antes_imp |> mutate(cuando = "Antes de imputar"),
                      st |> mutate(cuando = "Después de imputar")) |>
  group_by(cuando, YEAR, MONTH) |> summarise(mm = mean(VALUE), .groups = "drop") |>
  mutate(fecha = make_date(YEAR, MONTH, 1))
mostrar(
  ggplot(serie_ad, aes(fecha, mm, colour = cuando)) + geom_line(linewidth = .5) +
    facet_wrap(~ factor(cuando, levels = c("Antes de imputar", "Después de imputar")), ncol = 1) +
    scale_x_date(breaks = make_date(seq(S0, S1, 2), 1, 1), date_labels = "%Y") +
    scale_colour_manual(values = c("#737373", "#2166ac"), guide = "none") +
    labs(x = NULL, y = "lluvia media (mm/día)") +
    theme(strip.text = element_text(face = "bold", size = 11))
)

base::cat("Tablas guardadas en", OUT, "\n")