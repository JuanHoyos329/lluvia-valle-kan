# =============================================================================
# Proyecto final de Analítica de Datos 2026-2S
# =============================================================================

# =============================================================================
# PARTE 0. Configuración común
# =============================================================================

# ---- Paquetes: se instalan solo si faltan --------------------------------------
# fda y fda.usc se cargan primero para que dplyr quede encima (select, filter) como en el
# preprocesamiento original.
paquetes <- c("fda", "fda.usc", "dplyr", "tidyr", "lubridate", "ggplot2", "patchwork", "readxl",
              "writexl", "terra")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan) > 0) install.packages(faltan, repos = "https://cloud.r-project.org", dependencies = TRUE)
invisible(lapply(paquetes, library, character.only = TRUE))

# ---- Reproducibilidad: generador y semilla fijos -------------------------------
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
SEMILLA <- 2026
set.seed(SEMILLA)
options(stringsAsFactors = FALSE, warn = 1)

# Versiones usadas (para reportar en el informe)
cat("R:", R.version.string, "\n")
for (pq in paquetes) cat(sprintf("  %-10s %s\n", pq, as.character(packageVersion(pq))))

# ---- Partición del proyecto (profesor Ferro) -----------------------------------
# entrenamiento 1981-1993, validación 1994-1996, prueba 1997-2003. Sin estaciones reservadas.
# Todo lo que se calibra o ajusta usa solo los años de entrenamiento.
ENTRENA <- 1981:1993; VALIDA <- 1994:1996; PRUEBA <- 1997:2003

# ---- Rutas -----------------------------------------------------------------------
DATOS <- "./Datos"
if (!dir.exists(DATOS)) stop("No se encuentra la carpeta de datos '", DATOS, "' desde ", getwd(),
                             ". Ponga el directorio de trabajo en la carpeta del proyecto.")

# ---- Parte en Python -----------------------------------------------------------
# TRUE: en la PARTE 4 se corre proyecto_final.py (entrena las redes; es la parte más lenta).
# FALSE: se usan las salidas que ya están en salidas_fase4/ (de una corrida anterior del .py).
CORRER_PYTHON <- TRUE
# Intérprete de Python: el del entorno virtual del proyecto si existe; si no, el del sistema.
PYTHON <- c(".venv/Scripts/python.exe", ".venv/bin/python")
PYTHON <- if (any(file.exists(PYTHON))) PYTHON[file.exists(PYTHON)][1] else Sys.which(c("python3", "python"))[Sys.which(c("python3", "python")) != ""][1]

# ---- Limpieza entre partes -------------------------------------------------------
# Cada parte empieza con el área de trabajo limpia (como si fuera un script aparte), pero se
# conserva esta configuración. Así ninguna parte depende de objetos de otra: se comunican
# solo por los archivos que guardan.
CONSERVAR <- c("paquetes", "SEMILLA", "ENTRENA", "VALIDA", "PRUEBA", "DATOS", "CORRER_PYTHON",
               "PYTHON", "CONSERVAR", "limpiar")
limpiar <- function() {
  rm(list = setdiff(ls(envir = .GlobalEnv, all.names = TRUE), c(CONSERVAR, ".Random.seed")),
     envir = .GlobalEnv)
  par(mfrow = c(1, 1))
  invisible(NULL)
}


# ===================================================================================
# PARTE 1. Preprocesamiento de las estaciones  (ítem 4.1)
# ===================================================================================

D   <- DATOS
OUT <- "./salidas_preprocesamiento"
if (!dir.exists(D)) stop("No se encuentra la carpeta de datos '", D, "' desde ", getwd(),
                         ". Ponga el directorio de trabajo en la carpeta del proyecto o pase la ruta como argumento.")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

S0 <- 1981; S1 <- 2003     # años del satélite que se usan (periodo del proyecto)
# Partición del proyecto (la misma de la fase 4, PARTICION = "ferro"):
# entrenamiento 1981-1993, validación 1994-1996, prueba 1997-2003.
# Todo lo que se CALIBRA en la limpieza (umbrales, correlaciones, elección de vecinas) usa solo
# los años de entrenamiento, para que el preprocesamiento no use información de validación ni de prueba.
# (ENTRENA, VALIDA y PRUEBA se fijan en la PARTE 0)

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
# Se imprimen juntas, lado a lado, en una sola figura.
# Orden de las estaciones: de más a menos años "buenos" (cobertura >= 90 %)
orden <- ea |> group_by(STATION_ID) |> summarise(a90 = sum(cob >= 90), cob = mean(cob)) |>
  arrange(desc(a90), desc(cob)) |> pull(STATION_ID)

# Figura 1. Mapa de calor estación x año (1930-2003)
# Cada fila es una estación y cada columna un año. Arriba quedan las estaciones con más años
# completos y abajo las que menos tienen. Azul = año bien cubierto, rojo = año con muchos faltantes.
ea_plot <- ea |> mutate(est = factor(STATION_ID, levels = rev(orden)),
                        clase = cut(cob, c(-Inf, 0, 50, 80, 90, 99.999, Inf),
                                    labels = c("sin datos", "<50%", "50-80%", "80-90%", "90-<100%", "100%")))
fig1_calor <- (
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
fig2_reportan <- (
  ggplot(f2, aes(YEAR, n, colour = serie, linetype = serie)) + geom_line(linewidth = .8) +
    geom_vline(xintercept = S0 - .5, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("#000000", "#E69F00", "#009E73", "#CC79A7"), name = NULL) +
    scale_linetype_manual(values = c("solid", "solid", "longdash", "solid"), name = NULL) +
    scale_x_continuous(breaks = seq(1930, Y1, 5), expand = c(0, 0)) + expand_limits(y = 0) +
    labs(title = "¿Cuántas estaciones reportan cada año?", x = NULL, y = "n.º de estaciones") +
    theme(legend.position = "bottom")
)

# Figuras 1 y 2 juntas en una sola figura: a la izquierda la cobertura diaria, a la derecha las
# estaciones que reportan cada año
mostrar(fig1_calor + fig2_reportan)

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

fig3_cobertura <- (
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
fig4_huecos <- (
  ggplot(f5, aes(rango, n, fill = as.integer(rango) > 3)) + geom_col() + geom_text(aes(label = n), vjust = -.3) +
    scale_fill_manual(values = c("#2166ac", "#b2182b"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +              # espacio para las etiquetas
    labs(title = "¿Los faltantes son días sueltos o bloques largos?",
         x = "hueco más largo dentro del año (días consecutivos sin dato)", y = "pares estación-año")
)

# Figuras 3 y 4 juntas en una sola figura: a la izquierda la cobertura año por año, a la derecha
# el hueco más largo (la leyenda de la cobertura va abajo para no quedar entre las dos gráficas)
mostrar(
  (fig3_cobertura + theme(legend.position = "bottom")) +
    (fig4_huecos + scale_x_discrete(labels = label_wrap_gen(10))) +
    plot_layout(widths = c(1.8, 1))
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


# ===================================================================================
# PARTE 2. Construcción de datos funcionales (fase 1)  (ítem 1)
# ===================================================================================

# =============================================================================
# Fase 1 - Construcción de datos funcionales (ítem 1 de la rúbrica)
# Proyecto final de Analítica de Datos: lluvia del Valle del Cauca, 1981-2003
#
# Qué hace este script:
#   1. Lee los datos limpios del preprocesamiento (estaciones diarias y satélite por péntadas).
#   2. Pone las dos fuentes sobre el MISMO dominio anual: t en días, de 0 a 365.
#   3. Construye la base B-spline A MANO (recursión de Cox-de Boor) y la matriz Phi.
#   4. Calcula los coeficientes por mínimos cuadrados: c = (Phi'Phi)^-1 Phi' y.
#   5. Verifica la matriz suavizadora S = Phi (Phi'Phi)^-1 Phi': simétrica, idempotente, tr(S) = K.
#   6. Elige K con ECM, VC y VCG, para cada fuente, sin usar los años de prueba.
#   7. Calcula sigma^2, la covarianza residual y las bandas de confianza PUNTO A PUNTO.
#   8. Guarda los coeficientes de todas las curvas para las fases 2, 3 y 4.
#
# Las gráficas se imprimen en pantalla (panel Plots de RStudio) y no se guardan como archivos.
#
# La librería fda se usa SOLO para verificar que lo programado a mano da lo mismo.
#
# Decisiones:
#   - Escala: raíz cuadrada de la lluvia en mm/día (sqrt(mm/día)), como pide el modelo de la fase 4.
#   - Años bisiestos: la lluvia del 29 de febrero se suma al 28 de febrero. Todos los años
#     quedan con 365 días y se conserva el total anual.
#   - Estaciones: cada día se ubica en su centro (día 1 -> t = 0,5; día 365 -> t = 364,5).
#   - Satélite: cada péntada se ubica en el centro de sus días (calendario no bisiesto).
#   - Base: B-splines cúbicos (orden 4) con nudos equiespaciados en [0, 365].
#   - Partición del proyecto (PARTE 0): entrenamiento 1981-1993, validación 1994-1996 y prueba
#     1997-2003. K se elige solo con los años de entrenamiento (sin validación ni prueba).
#
# Paquetes: fda (solo para verificar)
# =============================================================================

limpiar()

ENT <- "./salidas_preprocesamiento"    # datos limpios de la fase 0 (Hoyos)
OUT <- "./salidas_fase1"               # salidas de esta fase
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
guardar <- function(df, nombre) write_xlsx(as.data.frame(df), file.path(OUT, paste0(nombre, ".xlsx")))

# ENTRENA (1981-1993) viene de la PARTE 0: K se elige solo con esos años
ORDEN   <- 4                           # B-splines cúbicos
RANGO   <- c(0, 365)                   # dominio anual común



#####################################################################################
# 1. Carga de datos
#####################################################################################

# ---- Estaciones: lluvia diaria limpia (incluye días imputados y desagregados)
d <- read.csv(file.path(ENT, "estaciones_diario_limpio.csv"))
d <- d[order(d$STATION_ID, d$YEAR, d$MONTH, d$DAY), ]

# 29 de febrero -> se suma al 28 de febrero de la misma estación y año
es29    <- d$MONTH == 2 & d$DAY == 29
clave   <- paste(d$STATION_ID, d$YEAR, d$MONTH, d$DAY)
pos28   <- match(paste(d$STATION_ID[es29], d$YEAR[es29], 2, 28), clave)
stopifnot(!anyNA(pos28))
d$VALUE[pos28] <- d$VALUE[pos28] + d$VALUE[es29]
d <- d[!es29, ]

# Claves de las curvas: una curva = una estación en un año
claves <- unique(d[, c("STATION_ID", "YEAR")])
rownames(claves) <- NULL
stopifnot(nrow(d) == 365 * nrow(claves))         # todos los años con 365 días

# Matriz Y: filas = días (365), columnas = curvas estación-año. Escala sqrt(mm/día)
Y <- matrix(sqrt(d$VALUE), nrow = 365)
N <- ncol(Y)
cat("Curvas estación-año:", N, "| estaciones:", length(unique(claves$STATION_ID)), "\n")

# ---- Satélite: lluvia por péntada en el píxel de cada estación
p <- read.csv(file.path(ENT, "estaciones_pentadas_limpio.csv"))
p <- p[order(p$STATION_ID, p$YEAR, p$PENTADA), ]
claves_sat <- unique(p[, c("STATION_ID", "YEAR")])
rownames(claves_sat) <- NULL
stopifnot(identical(claves$STATION_ID, claves_sat$STATION_ID),
          identical(claves$YEAR, claves_sat$YEAR))  # misma curva en la misma columna
stopifnot(nrow(p) == 72 * N, !anyNA(p$satelite_mm), all(p$dias %in% 3:6))

#Matriz X: filas = péntadas (72), columnas = las mismas curvas de Y
#El satélite trae el total de la péntada; se divide por sus días (3 a 6) para que quede en
#mm/día, la misma escala de las estaciones. Si no, la sexta péntada de cada mes (hasta 6 días)
#sale más alta solo por tener más días, y estación - satélite no se podría comparar
X <- matrix(sqrt(p$satelite_mm / p$dias), nrow = 72)

# Años fuera de la prueba
tr <- claves$YEAR %in% ENTRENA
cat("Curvas para elegir K (", paste(range(ENTRENA), collapse = "-"), "):", sum(tr), "\n")


#####################################################################################
# 2. Dominio común: t en días [0, 365]
#####################################################################################

# Estaciones: centro de cada día
t_est <- (1:365) - 0.5

# Satélite: centro de cada péntada (días 1-5, 6-10, 11-15, 16-20, 21-25 y 26-fin de mes)
dias_mes <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
inicio   <- cumsum(c(0, dias_mes[-12]))
t_sat <- unlist(lapply(1:12, function(m) sapply(1:6, function(k) {
  a <- 5 * (k - 1) + 1
  b <- if (k < 6) 5 * k else dias_mes[m]
  inicio[m] + (a + b) / 2 - 0.5
})))
stopifnot(length(t_sat) == 72)

n_est <- length(t_est)   # 365 puntos por curva
n_sat <- length(t_sat)   # 72 puntos por curva

# Malla fina para dibujar las curvas
t_fino <- seq(RANGO[1], RANGO[2], length.out = 731)


#####################################################################################
# 3. Base B-spline programada a mano (recursión de Cox-de Boor)
#####################################################################################
# Nudos: los extremos se repiten ORDEN veces y los nudos interiores se reparten igual.
# Orden 1: "escalones" que valen 1 en su intervalo y 0 fuera.
# Orden m: cada función se arma combinando dos funciones de orden m-1:
#   B_{j,m}(t) = (t - u_j)/(u_{j+m-1} - u_j) B_{j,m-1}(t) + (u_{j+m} - t)/(u_{j+m} - u_{j+1}) B_{j+1,m-1}(t)

nudos_bspline <- function(K, orden = ORDEN, rango = RANGO) {
  internos <- seq(rango[1], rango[2], length.out = K - orden + 2)
  internos <- internos[-c(1, length(internos))]
  c(rep(rango[1], orden), internos, rep(rango[2], orden))
}

base_bspline <- function(t, K, orden = ORDEN, rango = RANGO) {
  u <- nudos_bspline(K, orden, rango)
  nt <- length(t)
  # Orden 1
  B <- matrix(0, nt, length(u) - 1)
  for (j in 1:(length(u) - 1)) B[, j] <- as.numeric(t >= u[j] & t < u[j + 1])
  B[t == rango[2], max(which(u < rango[2]))] <- 1      # el extremo derecho entra en el último tramo
  # Órdenes 2, ..., orden
  if (orden > 1) for (m in 2:orden) {
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
  B                                                    # nt x K
}

# Ejemplo: base con K = 8 en la malla fina
Phi_ej <- base_bspline(t_fino, 8)
stopifnot(max(abs(rowSums(Phi_ej) - 1)) < 1e-12)       # las piezas suman 1 en todo t

matplot(t_fino, Phi_ej, type = "l", lty = 1, lwd = 2, col = rainbow(8),
        xlab = "Día del año", ylab = "Valor de la función base",
        main = "")                                     # el título va en el pie de figura del informe
abline(v = nudos_bspline(8), lty = 3, col = "grey50")
par(mfrow = c(1, 1))


#####################################################################################
# 4. Mínimos cuadrados paso a paso y matriz suavizadora S
#####################################################################################
# c_hat = (Phi'Phi)^-1 Phi' y      (una columna de coeficientes por curva)
# S     = Phi (Phi'Phi)^-1 Phi'    (y_hat = S y)

coeficientes <- function(Phi, Ymat) solve(crossprod(Phi), crossprod(Phi, Ymat))
matriz_S     <- function(Phi) Phi %*% solve(crossprod(Phi), t(Phi))

propiedades_S <- function(S, K) data.frame(
  K = K,
  simetria_max_dif     = max(abs(S - t(S))),          # S = S'
  idempotencia_max_dif = max(abs(S %*% S - S)),       # S S = S
  traza                = sum(diag(S)),                # tr(S) = K
  traza_menos_K        = sum(diag(S)) - K)


#####################################################################################
# 5. Selección de K con ECM, VC y VCG (solo entrenamiento, sin validación ni prueba)
#####################################################################################
# ECM(K) = promedio de (y - y_hat)^2                    -> siempre baja con K; no sirve solo
# VC(K)  = promedio de ((y - y_hat) / (1 - s_ii))^2     -> atajo de dejar-uno-fuera con la diagonal de S
# VCG(K) = promedio sobre curvas de n * SSE / (n - K)^2 -> validación cruzada generalizada

criterios_K <- function(t, Ymat, Ks) {
  n <- length(t)
  res <- NULL
  for (K in Ks) {
    Phi  <- base_bspline(t, K)
    S    <- matriz_S(Phi)
    E    <- Ymat - S %*% Ymat
    s_ii <- diag(S)
    res <- rbind(res, data.frame(
      K   = K,
      ECM = mean(E^2),
      VC  = mean((E / (1 - s_ii))^2),
      VCG = mean(n * colSums(E^2) / (n - K)^2)))
  }
  res
}

# Estaciones: curvas de entrenamiento
Ks_est    <- seq(6, 60, by = 2)
tabla_est <- criterios_K(t_est, Y[, tr], Ks_est)

# Satélite: varias estaciones comparten píxel; para elegir K se usa cada píxel-año una sola vez
X_tr_unicas <- t(unique(t(X[, tr])))
cat("Curvas de satélite (píxel-año) de entrenamiento:", ncol(X_tr_unicas), "\n")
Ks_sat    <- seq(6, 40, by = 2)
tabla_sat <- criterios_K(t_sat, X_tr_unicas, Ks_sat)

K_est <- tabla_est$K[which.min(tabla_est$VCG)]
K_sat <- tabla_sat$K[which.min(tabla_sat$VCG)]
cat("\nK elegido por VCG -> estaciones:", K_est, "| satélite:", K_sat, "\n")
cat("K elegido por VC  -> estaciones:", tabla_est$K[which.min(tabla_est$VC)],
    "| satélite:", tabla_sat$K[which.min(tabla_sat$VC)], "\n")
print(round(tabla_est, 4)); print(round(tabla_sat, 4))

grafica_K <- function(tabla, K_elegido, fuente) {
  par(mfrow = c(1, 3))
  plot(tabla$K, tabla$ECM, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "ECM",
       main = "(a) ECM")
  grid(5, 5, lwd = 2)
  plot(tabla$K, tabla$VC, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VC",
       main = "(b) VC")
  abline(v = tabla$K[which.min(tabla$VC)], lty = 3, col = "red")
  grid(5, 5, lwd = 2)
  plot(tabla$K, tabla$VCG, type = "b", pch = 19, xlab = "Número de bases (K)", ylab = "VCG",
       main = "(c) VCG")
  abline(v = K_elegido, lty = 3, col = "red")
  abline(h = min(tabla$VCG), lty = 2, col = "red")
  grid(5, 5, lwd = 2)
}
grafica_K(tabla_est, K_est, "estaciones")
par(mfrow = c(1, 1))
grafica_K(tabla_sat, K_sat, "satélite")
par(mfrow = c(1, 1))
# Nota: el ECM baja siempre al subir K (sobreajuste). La VC y la VCG penalizan la complejidad
# y tienen un mínimo interior. Con datos diarios la lluvia es muy intermitente: la VCG prefiere
# curvas suaves (K pequeño) y la curva conserva los dos picos del régimen bimodal (ver sección 8).


#####################################################################################
# 6. Ajuste final con el K elegido y verificación contra la librería fda
#####################################################################################

Phi_est <- base_bspline(t_est, K_est)              # 365 x K_est
Phi_sat <- base_bspline(t_sat, K_sat)              # 72  x K_sat
C_est   <- coeficientes(Phi_est, Y)                # K_est x N
C_sat   <- coeficientes(Phi_sat, X)                # K_sat x N
S_est   <- matriz_S(Phi_est)
S_sat   <- matriz_S(Phi_sat)

# Propiedades de S (deben dar ~0, ~0 y K)
verif_S <- rbind(cbind(fuente = "estaciones", propiedades_S(S_est, K_est)),
                 cbind(fuente = "satélite",   propiedades_S(S_sat, K_sat)))
print(verif_S)

# Verificación con fda: misma base y mismos coeficientes (diferencias del orden de 1e-12)
base_fda_est <- create.bspline.basis(RANGO, nbasis = K_est, norder = ORDEN)
base_fda_sat <- create.bspline.basis(RANGO, nbasis = K_sat, norder = ORDEN)
muestra <- 1:50
verif_fda <- data.frame(
  fuente = c("estaciones", "satélite"),
  max_dif_Phi = c(max(abs(Phi_est - eval.basis(t_est, base_fda_est))),
                  max(abs(Phi_sat - eval.basis(t_sat, base_fda_sat)))),
  max_dif_coeficientes = c(
    max(abs(C_est[, muestra] - smooth.basis(t_est, Y[, muestra], base_fda_est)$fd$coefs)),
    max(abs(C_sat[, muestra] - smooth.basis(t_sat, X[, muestra], base_fda_sat)$fd$coefs))))
print(verif_fda)

# Las dos fuentes evaluadas en la misma malla fina (mismo dominio anual)
Phi_fino_est <- base_bspline(t_fino, K_est)
Phi_fino_sat <- base_bspline(t_fino, K_sat)
curvas_est <- Phi_fino_est %*% C_est               # 731 x N
curvas_sat <- Phi_fino_sat %*% C_sat

par(mfrow = c(1, 2))
set.seed(2026)
m <- sample(which(tr), 200)
# Estilo de Ospina-Galindez et al. (2019): curvas en gris y la media en negro
matplot(t_fino, curvas_est[, m], type = "l", lty = 1, col = adjustcolor("grey55", 0.35),
        xlab = "Día del año", ylab = "Precipitación, √(mm/día)", main = "(a) Estaciones")
lines(t_fino, rowMeans(curvas_est[, tr]), lwd = 2.5, col = "black")
legend("topleft", c("Curvas individuales", paste("Media", paste(range(ENTRENA), collapse = "-"))), col = c("grey55", "black"),
       lwd = c(1, 2.5), bty = "n", cex = 0.8)
matplot(t_fino, curvas_sat[, m], type = "l", lty = 1, col = adjustcolor("grey55", 0.35),
        xlab = "Día del año", ylab = "Precipitación, √(mm/día)", main = "(b) Satélite")
lines(t_fino, rowMeans(curvas_sat[, tr]), lwd = 2.5, col = "black")
par(mfrow = c(1, 1))


#####################################################################################
# 7. Varianza residual, covarianza residual y bandas de confianza punto a punto
#####################################################################################
# Bajo Var(y) = sigma^2 I:
#   sigma^2 (insesgado)     = SSE / (n - K)            (un valor por curva)
#   Var(c_hat)              = sigma^2 (Phi'Phi)^-1
#   Var(x_hat(t))           = sigma^2 phi(t)' (Phi'Phi)^-1 phi(t)
#   Banda punto a punto     = x_hat(t) +/- 2 EE(t)
# Es PUNTO A PUNTO: en cada t, por separado, cubre el valor de la curva con ~95 %.
# No es una banda simultánea para toda la curva a la vez.

sigma2 <- function(Ymat, S, K) colSums((Ymat - S %*% Ymat)^2) / (nrow(Ymat) - K)
sig2_est <- sigma2(Y, S_est, K_est)
sig2_sat <- sigma2(X, S_sat, K_sat)
cat("\nsigma^2 mediana -> estaciones:", round(median(sig2_est), 3),
    "| satélite:", round(median(sig2_sat), 3), "\n")

# h(t) = phi(t)' (Phi'Phi)^-1 phi(t) en la malla fina
h_fino <- function(Phi_fino, Phi) rowSums((Phi_fino %*% solve(crossprod(Phi))) * Phi_fino)
h_est  <- h_fino(Phi_fino_est, Phi_est)
h_sat  <- h_fino(Phi_fino_sat, Phi_sat)

# Covarianza residual: ¿se cumple Var(y) = sigma^2 I?  Si los residuos de días vecinos
# están poco correlacionados, el supuesto es razonable y las bandas son válidas.
E_est   <- (Y - S_est %*% Y)[, tr]
Cor_res <- cor(t(E_est))                                   # 365 x 365
cor_lag1 <- mean(diag(Cor_res[-1, -365]))
cat("Correlación media entre residuos de días consecutivos (estaciones):", round(cor_lag1, 3), "\n")

# Correlación media según la distancia en días (rezago): resume el mapa en una sola curva
rezagos <- 1:30
cor_rezago <- sapply(rezagos, function(L) mean(diag(Cor_res[-(1:L), -((366 - L):365)])))
par(mfrow = c(1, 2))
image(t_est, t_est, Cor_res, col = hcl.colors(50, "Blue-Red 3"), zlim = c(-1, 1),
      xlab = "Día s", ylab = "Día t", main = "(a)")
plot(rezagos, cor_rezago, type = "h", lwd = 3, col = "grey30", ylim = c(-0.05, 0.3),
     xlab = "Distancia entre días (rezago)", ylab = "Correlación media", main = "(b)")
abline(h = 0, col = "grey40")
abline(h = c(-0.1, 0.1), lty = 3, col = "red")
par(mfrow = c(1, 1))

# Ejemplo de bandas: una estación-año de entrenamiento, ambas fuentes
j <- which(tr)[1]
grafica_banda <- function(t_obs, y_obs, curva, sig2, h, titulo, unidad) {
  ee <- sqrt(sig2 * h)
  plot(t_obs, y_obs, pch = 19, cex = 0.5, col = "grey60", xlab = "Día del año", ylab = unidad,
       ylim = range(y_obs, curva + 2 * ee, curva - 2 * ee), main = titulo)
  lines(t_fino, curva, lwd = 2, col = "black")
  lines(t_fino, curva + 2 * ee, lty = 2, col = "black")
  lines(t_fino, curva - 2 * ee, lty = 2, col = "black")
  legend("topleft", c("Observado", "Curva estimada", "Banda +/- 2 EE (punto a punto)"),
         pch = c(19, NA, NA), lty = c(NA, 1, 2), col = c("grey60", "black", "black"), bty = "n", cex = 0.8)
}
par(mfrow = c(1, 2))
grafica_banda(t_est, Y[, j], curvas_est[, j], sig2_est[j], h_est,
              paste0("(a) Estación ", claves$STATION_ID[j], ", ", claves$YEAR[j]), "Precipitación, √(mm/día)")
grafica_banda(t_sat, X[, j], curvas_sat[, j], sig2_sat[j], h_sat,
              paste0("(b) Satélite en su píxel, ", claves$YEAR[j]), "Precipitación, √(mm/día)")
par(mfrow = c(1, 1))


#####################################################################################
# 8. Control: ¿la curva media conserva el régimen bimodal?
#####################################################################################
media_est <- rowMeans(curvas_est[, tr])
picos <- which(diff(sign(diff(media_est))) == -2) + 1
cat("Picos de la curva media de estaciones (día del año):", round(t_fino[picos]), "\n")


#####################################################################################
# 9. Salidas para las fases 2, 3 y 4
#####################################################################################
coef_est <- data.frame(claves, t(C_est));  names(coef_est)[-(1:2)] <- paste0("c", 1:K_est)
coef_sat <- data.frame(claves, t(C_sat));  names(coef_sat)[-(1:2)] <- paste0("c", 1:K_sat)
# Los coeficientes y base_info quedan en CSV: son datos de entrada de la parte en Python (fase 4)
coef_est$particion <- coef_sat$particion <- ifelse(claves$YEAR %in% ENTRENA, "entrenamiento",
                                             ifelse(claves$YEAR %in% VALIDA, "validacion", "prueba"))
write.csv(coef_est, file.path(OUT, "coeficientes_estaciones.csv"), row.names = FALSE)
write.csv(coef_sat, file.path(OUT, "coeficientes_satelite.csv"),   row.names = FALSE)

# Tablas (Excel)
guardar(tabla_est, "seleccion_K_estaciones")
guardar(tabla_sat, "seleccion_K_satelite")
guardar(verif_S,   "verificacion_matriz_S")
guardar(verif_fda, "verificacion_vs_fda")

# Cómo reconstruir la base en las otras fases (mismos nudos, mismo orden)
base_info <- data.frame(fuente = c("estaciones", "satélite"), K = c(K_est, K_sat), orden = ORDEN,
                        rango = "0-365", escala = "sqrt(mm/día)",
                        nudos = c(paste(round(nudos_bspline(K_est), 4), collapse = " "),
                                  paste(round(nudos_bspline(K_sat), 4), collapse = " ")))
write.csv(base_info, file.path(OUT, "base_info.csv"), row.names = FALSE)

cat("\nListo. Salidas en", OUT, "\n")


#####################################################################################
# PARTE 2b. Complemento de la fase 1: ¿por qué B-spline, mínimos cuadrados ordinarios y K = 8?
#
# Tres revisiones que respaldan las decisiones de la construcción funcional (PARTE 2):
#   A. Base de Fourier frente a B-spline (VCG y condición de Phi'Phi)
#   B. Varianza de los residuos por día del año: ¿hace falta mínimos cuadrados ponderados?
#   C. Suavizado penalizado (segunda derivada) con lambda elegido por VCG
#
# Usa solo las curvas de entrenamiento (las mismas con que se eligió K). No modifica nada de la
# fase 1 ni guarda archivos: imprime los números en la consola y muestra las figuras.
#####################################################################################

limpiar()
ENT     <- "./salidas_preprocesamiento"
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
Ytr <- Y[, tr]                       # estaciones, años de entrenamiento
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


# ===================================================================================
# PARTE 3. Análisis exploratorio de datos funcionales (fase 2)  (ítem 2)
# ===================================================================================

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
#     regla de la fase 4 (Python, PARTE 4); si ya existe salidas_fase4/fase_enso_por_anio.csv,
#     se comprueba que las dos clasificaciones coinciden.

#Limpiar área de trabajo (se conserva la configuración de la PARTE 0)
limpiar()

ENT0 <- "./salidas_preprocesamiento"  
ENT1 <- "./salidas_fase1"              
ENT4 <- "./salidas_fase4"              
# DATOS (carpeta de datos, con el ONI de la NOAA CPC) viene de la PARTE 0

#Carpeta para guardar las tablas (Excel). Las figuras solo se imprimen.
OUT  <- "./salidas_fase2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

set.seed(SEMILLA)

# Imprime la figura en pantalla (panel Plots de RStudio); no se guarda como archivo.
# El nombre sirve para ubicar la figura del informe en el código.
figura <- function(nombre, codigo, ancho = 10, alto = 4.5) {
  codigo <- substitute(codigo); entorno <- parent.frame()
  eval(codigo, entorno); par(mfrow = c(1, 1))
  invisible(nombre)
}
guardar <- function(df, nombre) write_xlsx(as.data.frame(df), file.path(OUT, paste0(nombre, ".xlsx")))

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

cat("\nListo. Tablas en", OUT, "\n")


# ===================================================================================
# PARTE 4. Red funcional, ablación e interpretación de la KAN (fase 4, en Python)  (ítems 4.2 y 4.3)
# ===================================================================================
# El código está en proyecto_final.py. Lee salidas_fase1/ y salidas_preprocesamiento/ y guarda
# en salidas_fase4/ los coeficientes predichos, las tablas de la ablación y los modelos.

limpiar()
if (CORRER_PYTHON) {
  if (is.na(PYTHON) || !nzchar(PYTHON)) stop("No se encontró Python. Instálelo (o cree .venv) o ponga CORRER_PYTHON <- FALSE.")
  cat("Corriendo proyecto_final.py con", PYTHON, "(puede tardar)...\n")
  estado <- system2(PYTHON, c("proyecto_final.py", "--sin-graficas"))
  if (estado != 0) stop("proyecto_final.py terminó con error (código ", estado, ").")
}
necesarios <- file.path("./salidas_fase4", c("pred_coef_lineal_semilla0.csv", "fase_enso_por_anio.csv",
                                             "ablacion_resultados.csv", "reconstruccion_mm.csv"))
if (!all(file.exists(necesarios)))
  stop("Faltan salidas de la fase 4: corra proyecto_final.py (o ponga CORRER_PYTHON <- TRUE).")


# ===================================================================================
# PARTE 5. Comparación funcional por fuente y regresión funcional (fase 3)  (ítem 3)
# ===================================================================================

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
# Necesita las salidas de la fase 4 (PARTE 4, Python): de ahí toma la partición y las predicciones
# de las redes para comparar con la regresión funcional (sección 6).
# Paquetes: fda, readxl, terra (terra solo para verificar la extracción del píxel)
# =============================================================================

limpiar()

ENT0 <- "./salidas_preprocesamiento"     # datos limpios (Hoyos)
ENT1 <- "./salidas_fase1"                # coeficientes de la fase 1 (Daniela): orden de las curvas
DSAT <- file.path(DATOS, "ImagenesSatelitePentadal")
OUT  <- "./salidas_fase3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

set.seed(SEMILLA)
# ENTRENA (años de entrenamiento) se toma de la partición de la fase 1, más abajo (sección 0.1)

# Nombres de las zonas (acuerdo del grupo, 21 nov): "Cauca" y "Pacífico" solos se confunden con otros
# departamentos. Internamente se usan las claves cortas; en tablas y figuras se muestran los nombres completos.
ZONA <- c("Cauca" = "Cuenca del río Cauca", "Pacífico" = "Vertiente del Pacífico")
nombre_zona <- function(x) {
  x <- gsub("Pacífico", ZONA[["Pacífico"]], x, fixed = TRUE)
  gsub("(?<!del )(?<!río )\\bCauca\\b", ZONA[["Cauca"]], x, perl = TRUE)
}
guardar <- function(df, nombre) {
  for (v in names(df)) if (is.character(df[[v]]) || is.factor(df[[v]])) df[[v]] <- nombre_zona(as.character(df[[v]]))
  write_xlsx(as.data.frame(df), file.path(OUT, paste0(nombre, ".xlsx")))
}

# Imprime la figura en pantalla (panel Plots de RStudio); no se guarda como archivo.
# El nombre sirve para ubicar la figura del informe en el código.
figura <- function(nombre, codigo, ancho = 10, alto = 4.5, res = 150) {
  codigo <- substitute(codigo); entorno <- parent.frame()
  eval(codigo, entorno); par(mfrow = c(1, 1))
  invisible(nombre)
}


#####################################################################################
# 0. Curvas comparables
#####################################################################################

# ---- 0.1 Curvas estación-año, en el mismo orden que la fase 1 (y que la fase 4)
ce     <- read.csv(file.path(ENT1, "coeficientes_estaciones.csv"))
claves <- ce[, c("STATION_ID", "YEAR", "particion")]
# Partición oficial del grupo (Fabián, 21 nov): la columna "particion" de salidas_fase4/pred_coef_*.csv
# (entrenamiento 1981-1993, validación 1994-1996, prueba 1997-2003). Se toma de ahí para que todos
# evaluemos con el mismo examen; el orden de las curvas es el mismo de la fase 1.
part4 <- read.csv(file.path(".", "salidas_fase4", "pred_coef_lineal_semilla0.csv"))
stopifnot(nrow(part4) == nrow(claves), all(part4$STATION_ID == claves$STATION_ID), all(part4$YEAR == claves$YEAR))
stopifnot("La partición de la fase 4 no es la de la PARTE 0" = identical(part4$particion, claves$particion))
claves$particion <- part4$particion
ENTRENA <- sort(unique(claves$YEAR[claves$particion == "entrenamiento"]))   # misma partición que las demás fases
ANIOS   <- function(part) paste(range(claves$YEAR[claves$particion == part]), collapse = "-")
cat("Partición: entrenamiento", ANIOS("entrenamiento"), "| validación", ANIOS("validacion"), "| prueba", ANIOS("prueba"), "\n")
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
          ylab = "VCG / mínimo", main = paste0("Selección de una base común (péntadas, ", ANIOS("entrenamiento"), ")"))
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
# (la columna "confianza" ya no existe en el preprocesamiento del 21 nov: todas las coordenadas son reales)
emp$curvas     <- tabulate(match(claves$STATION_ID, emp$STATION_ID), nbins = nrow(emp))

# ---- 1.1 Regla de asignación: el píxel que CONTIENE la estación
# Malla de 0,05° (37 x 39), esquina superior izquierda en (-77,55; 5,00); fila 1 = norte
col_r  <- floor((emp$Lon + 77.55) / 0.05) + 1
fila_r <- floor((5.00 - emp$lat) / 0.05) + 1
# Una estación que cae EXACTAMENTE sobre el borde entre píxeles (coordenadas redondeadas al minuto, como
# 5300000202 Colpuertos: 3°48' N, 77°06' W = esquina de 4 píxeles) puede asignarse a cualquiera de ellos según
# la convención. Se acepta la asignación de Hoyos (con ella se extrajo el satélite) y se deja anotada.
en_borde <- function(v) abs(v - round(v)) < 1e-6
emp$en_borde_pixel <- en_borde((emp$Lon + 77.55) / 0.05) | en_borde((5.00 - emp$lat) / 0.05)
dif_regla <- col_r != emp$col | fila_r != emp$fila
if (any(dif_regla)) {
  cat("Estaciones sobre el borde entre píxeles (se conserva el píxel del preprocesamiento):\n")
  print(data.frame(emp[dif_regla, c("STATION_ID", "Lon", "lat", "fila", "col")], fila_regla = fila_r[dif_regla], col_regla = col_r[dif_regla]))
}
stopifnot(all(!dif_regla | emp$en_borde_pixel))   # cualquier otra diferencia sí es un error
emp$lon_pixel      <- -77.55 + (emp$col - 0.5) * 0.05
emp$lat_pixel      <-   5.00 - (emp$fila - 0.5) * 0.05
emp$dist_centro_km <- round(111.32 * sqrt(((emp$Lon - emp$lon_pixel) * cos(emp$lat * pi / 180))^2 +
                                          (emp$lat - emp$lat_pixel)^2), 2)

# ---- 1.2 Estaciones que comparten píxel (reciben la MISMA curva satelital)
emp$pixel        <- paste0("F", emp$fila, "-C", emp$col)
emp$est_en_pixel <- ave(rep(1, nrow(emp)), emp$pixel, FUN = sum)
# (En la corrida del 17 nov había 11 estaciones con coordenadas trianguladas "no confiables". Desde el
# preprocesamiento del 21 nov todas las coordenadas son reales y las estaciones sin coordenadas se
# excluyeron, así que esa lista se eliminó del script.)
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
               paste("Coordenadas:", names(table(emp$fuente_coord))),
               "Altitud de la grilla (no de catálogo)"),
  valor = c(nrow(emp), length(unique(emp$pixel)), sum(emp$est_en_pixel > 1),
            max(emp$est_en_pixel), median(emp$dist_centro_km),
            as.vector(table(emp$fuente_coord)),
            sum(emp$fuente_alt != "catálogo")))
print(resumen_emp)
guardar(resumen_emp, "tabla_f3_05_resumen_emparejamiento")

figura("fig_f3_03_mapa_emparejamiento", {
  fuentes <- sort(unique(emp$fuente_coord))
  cols <- setNames(c("black", "steelblue", "orange", "red", "purple")[seq_along(fuentes)], fuentes)
  par(mar = c(4.2, 4.2, 2.5, 1))
  plot(emp$Lon, emp$lat, type = "n", xlim = c(-77.55, -75.70), ylim = c(3.05, 5.00),
       asp = 1 / cos(4 * pi / 180), xlab = "Longitud", ylab = "Latitud",
       main = "Estaciones y malla del satélite (píxeles de 0,05°)")
  abline(v = -77.55 + 0.05 * (0:37), h = 5 - 0.05 * (0:39), col = "grey88")
  o <- order(emp$fuente_coord != "archivo original")                    # las coordenadas rescatadas, encima
  points(emp$Lon[o], emp$lat[o], pch = ifelse(emp$est_en_pixel[o] > 1, 17, 19),
         cex = ifelse(emp$fuente_coord[o] == "archivo original", 0.9, 1.1), col = cols[emp$fuente_coord[o]])
  legend("topleft", c(names(cols), "comparte píxel (triángulo)"), col = c(cols, "grey40"),
         pch = c(rep(19, length(cols)), 17), bg = "white", cex = 0.8)
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
#     MEDIANA por estación y el % de estaciones con d > 0, y se repite solo con las coordenadas del
#     archivo original de la CVC (sin las coordenadas rescatadas el 21 nov), como prueba de sensibilidad.
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
                       coord  = ifelse(emp$fuente_coord[ie] == "archivo original", "original", "rescatada"))

grupos <- list(
  "Todas"               = rep(TRUE, n_est),
  "Solo originales"     = info_est$coord == "original",
  "Cauca"               = info_est$region == "Cauca",
  "Cauca - originales"  = info_est$region == "Cauca" & info_est$coord == "original",
  "Pacífico"            = info_est$region == "Pacífico",
  "Pacífico - originales" = info_est$region == "Pacífico" & info_est$coord == "original")
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
                    fuente_coord = emp$fuente_coord[ie],
                    altitud = emp$altitud[ie], Lon = round(emp$Lon[ie], 3), lat = round(emp$lat[ie], 3),
                    pixel = emp$pixel[ie], curvas = emp$curvas[ie],
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
etiq <- function(g) nombre_zona(paste0(g, " (n = ", sum(grupos[[g]]), ")"))

figura("fig_f3_04_diferencia_mm", {
  par(mfrow = c(1, 2))
  col_reg <- ifelse(info_est$region == "Cauca", "grey75", adjustcolor("orange", 0.7))
  matplot(t_eval, G_dmm, type = "l", lty = 1, col = col_reg, xlab = "Día del año", ylab = "mm/día",
          main = "d(t) por estación (promedio de sus años)", cex.main = 0.95)
  lines(t_eval, bandas$Todas$media, lwd = 3, col = "darkred")
  lines(t_eval, bandas$Todas$mediana, lwd = 2, lty = 2)
  abline(h = 0, lty = 3)
  legend("topright", c("cuenca del río Cauca", "vertiente del Pacífico", "media", "mediana"),
         col = c("grey60", "orange", "darkred", "black"), lty = c(1, 1, 1, 2), lwd = c(1, 1, 3, 2), bty = "n", cex = 0.75)
  dibujar_banda(bandas$Todas, "mm/día", paste("Todas las estaciones (n =", n_est, ")"), espacio = 0.6)
  leyenda_banda()
})

figura("fig_f3_05_diferencia_por_region", {
  par(mfrow = c(2, 2))
  dibujar_banda(bandas$Cauca, "mm/día", etiq("Cauca"), espacio = 0.6)
  leyenda_banda()
  dibujar_banda(bandas$`Pacífico`, "mm/día", etiq("Pacífico"))
  comparar_bandas(bandas$Cauca, bandas$`Cauca - originales`, etiq("Cauca"), etiq("Cauca - originales"),
                  "Río Cauca: todas vs. coord. originales")
  comparar_bandas(bandas$`Pacífico`, bandas$`Pacífico - originales`, etiq("Pacífico"), etiq("Pacífico - originales"),
                  "V. del Pacífico: todas vs. coord. originales")
}, ancho = 12, alto = 8.5)

figura("fig_f3_06_diferencia_sqrt_y_pct", {
  par(mfrow = c(1, 2))
  dibujar_banda(b_sq, "sqrt(mm/día)", "Todas, en sqrt(mm/día)")
  plot(t_eval, bandas$Cauca$pct_subestima, type = "l", lwd = 3, col = "grey40", ylim = c(0, 100),
       xlab = "Día del año", ylab = "% de estaciones", main = "Estaciones donde el satélite subestima", cex.main = 0.95)
  lines(t_eval, bandas$`Pacífico`$pct_subestima, lwd = 3, col = "orange")
  lines(t_eval, bandas$`Solo originales`$pct_subestima, lwd = 2, lty = 2, col = "darkgreen")
  abline(h = 50, lty = 3)
  legend("topright", c(etiq("Cauca"), etiq("Pacífico"), etiq("Solo originales")),
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
          main = paste0("r(s,t) río Cauca, ", sub("Cauca - ", "", a, fixed = TRUE), " (", n_cc[a], " curvas)"), cex.main = 0.95)
    contour(t_eval, t_eval, analisis_cc[[a]], levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
    abline(0, 1, lty = 3, col = "black")
  }
  cols_g <- c("Todas" = "black", "Cauca" = "grey50", "Pacífico" = "orange")
  matplot(t_eval, diagonales, type = "l", ylim = c(-0.2, 1), xlab = "Día del año", ylab = "r(t, t)",
          lty = rep(c(1, 2), 3), lwd = 2, col = rep(cols_g, each = 2),
          main = "Diagonal de la correlación cruzada", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -2, c(BORDE, 365), 2, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("bottomright", nombre_zona(colnames(diagonales)), col = rep(cols_g, each = 2), lty = rep(c(1, 2), 3),
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
  for (q in list(list(CY, "Covarianza estación - río Cauca, total (mm²/día²)"),
                 list(CX, "Covarianza satélite - río Cauca, total (misma escala)"))) {
    image(t_eval, t_eval, q[[1]], zlim = range(niveles), col = hcl.colors(40, "YlOrRd", rev = TRUE),
          xlab = "s (día)", ylab = "t (día)", main = q[[2]], cex.main = 0.95)
    contour(t_eval, t_eval, q[[1]], levels = niveles, labcex = 0.6, add = TRUE)
  }
  de <- sapply(c("Cauca - total", "Cauca - anomalías"), function(a)
    cbind(sqrt(diag(cov_grupo[[a]]$Y)), sqrt(diag(cov_grupo[[a]]$X))))
  de <- matrix(de, ncol = 4)
  matplot(t_eval, de, type = "l", lty = c(1, 2, 1, 2), lwd = c(3, 3, 1.5, 1.5),
          col = c("slategray", "darkturquoise", "slategray", "darkturquoise"), ylim = c(0, max(de)),
          xlab = "Día del año", ylab = "mm/día", main = "Desviación estándar funcional - río Cauca", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = adjustcolor("grey90", 0.6), border = NA)
  legend("topright", c("estación total", "satélite total", "estación anomalías", "satélite anomalías"),
         col = c("slategray", "darkturquoise", "slategray", "darkturquoise"), lty = c(1, 2, 1, 2),
         lwd = c(3, 3, 1.5, 1.5), bty = "n", cex = 0.75)
  pc <- rbind(pct_componentes(cov_grupo[["Cauca - anomalías"]]$Y)[1:6],
              pct_componentes(cov_grupo[["Cauca - anomalías"]]$X)[1:6])
  barplot(pc, beside = TRUE, names.arg = 1:6, col = c("slategray", "darkturquoise"),
          xlab = "Componente", ylab = "% de varianza", main = "Componentes de la covarianza - río Cauca, anomalías",
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
cor_alt <- cor.test(info_est$altitud[cauca], info_est$sesgo_pct[cauca], method = "spearman", exact = FALSE)
# Compresión espacial: ¿el satélite reproduce las diferencias de lluvia entre sitios?
# Se ajusta log(satélite) ~ log(estación) con la lluvia media de cada estación (todas).
# Pendiente = 1: el satélite reproduce las diferencias; pendiente < 1: las comprime (sobreestima los
# sitios secos y subestima los lluviosos). No se correlaciona la razón est/sat con la lluvia de la
# estación porque esa correlación es en parte espuria (la estación está en los dos lados).
fit_comp <- lm(log(x_anual) ~ log(y_anual))
ic_comp  <- confint(fit_comp)[2, ]
tabla_cor_alt <- data.frame(
  relacion   = c("Cauca: altitud vs. sesgo anual (%) - rho de Spearman",
                 "Todas: pendiente de log(satélite) ~ log(estación)"),
  estaciones = c(sum(cauca), length(y_anual)),
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
          ylim = yl, main = "Cuenca del río Cauca: d(t) por franja de altitud", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("topright", paste0(franjas_cauca, " (n = ", sapply(franjas_cauca, function(fr) sum(cauca & info_est$franja == fr)), ")"),
         col = cols_fr, lwd = 2.5, bty = "n", cex = 0.75, ncol = 2)
  plot(t_int, an_alt$F_t, type = "l", lwd = 2, xlab = "Día del año", ylab = "F(t)",
       main = paste0("ANOVA funcional por franja (p global ", texto_p(an_alt$p_global), ")"), cex.main = 0.95,
       ylim = c(0, max(an_alt$F_t, an_alt$c_max) * 1.45))
  abline(h = an_alt$c_max, lty = 2, col = "red")
  legend("topright", c("F(t) observado", "valor crítico máximo (5 %)"), col = c("black", "red"), lty = 1:2, bty = "n", cex = 0.75)
  plot(info_est$altitud, info_est$sesgo_pct, pch = 19, cex = 0.8,
       col = ifelse(cauca, "grey40", "orange"), xlab = "Altitud (m)", ylab = "Sesgo anual (%)",
       main = paste0("Sesgo anual por estación (río Cauca: rho = ", round(cor_alt$estimate, 2), ")"), cex.main = 0.95)
  abline(h = 0, lty = 3)
  legend("bottomright", c("cuenca del río Cauca", "vertiente del Pacífico"),
         col = c("grey40", "orange"), pch = c(19, 19), bg = "white", cex = 0.75)
  plot(y_anual, x_anual, log = "xy", pch = 19, cex = 0.8,
       col = ifelse(cauca, "grey40", "orange"), xlab = "Lluvia media de la estación (mm/día, log)",
       ylab = "Lluvia media del satélite en su píxel (mm/día, log)",
       main = paste0("Compresión espacial: pendiente = ", round(coef(fit_comp)[2], 2), " (1 = sin compresión)"),
       cex.main = 0.95)
  abline(0, 1, lty = 3)
  xs <- range(y_anual); lines(xs, exp(coef(fit_comp)[1]) * xs^coef(fit_comp)[2], lwd = 2, col = "darkred")
  legend("topleft", c("1:1 (sin sesgo)", "ajuste log-log"), lty = c(3, 1), lwd = c(1, 2),
         col = c("black", "darkred"), bty = "n", cex = 0.75)
}, ancho = 12, alto = 9)

cols_fase <- c("Niño" = "firebrick", "Neutro" = "grey40", "Niña" = "steelblue")
figura("fig_f3_10_sesgo_enso", {
  par(mfrow = c(1, 3))
  med_fase <- sapply(fases, function(f) rowMeans(G_fase_d[[f]][, es_cauca(G_fase_d[[f]]), drop = FALSE]))
  matplot(t_eval, med_fase, type = "l", lty = 1, lwd = 2.5, col = cols_fase, xlab = "Día del año", ylab = "mm/día",
          main = "Cuenca del río Cauca: d(t) por fase ENSO", cex.main = 0.95)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = adjustcolor("grey90", 0.6), border = NA)
  abline(h = 0, lty = 3)
  legend("topright", fases, col = cols_fase, lwd = 2.5, bty = "n", cex = 0.8)
  dibujar_banda(b_nino, "mm/día", paste0("Río Cauca: Niño - Neutro, pareado (n = ", ncol(par_nino$D), ")"), col = "firebrick")
  dibujar_banda(b_nina, "mm/día", paste0("Río Cauca: Niña - Neutro, pareado (n = ", ncol(par_nina$D), ")"), col = "steelblue")
}, ancho = 15, alto = 5)


#####################################################################################
# 6. Regresión funcional función-función (modelo de referencia de la fase 4)
#####################################################################################
# Y_i(t) = beta0(t) + ∫ beta1(s, t) X_i(s) ds + e_i(t)       (script de clase: linmod)
#
# Para que sea comparable con las redes de la fase 4 se usan EXACTAMENTE sus entradas y su métrica:
#   - X = curva satelital de la fase 1; Y = curva de la estación de la fase 1 (sqrt, diaria). Los K se
#     leen de base_info.csv. Mismo orden de curvas y misma partición.
#   - Error L² relativo = 100 * sqrt( suma ||y_hat - y||² / suma ||y||² ), con la matriz de Gram de la
#     base de las estaciones (igual que error_L2 y L2_relativo en fase4_red_funcional.py).
# beta1 en base B-spline (K_s1 en s x K_e1 en t), penalización de la segunda derivada en s y en t con un
# lambda común elegido con VALIDACIÓN. La PRUEBA solo se usa para reportar (años: ver ANIOS()).
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
# Control de consistencia con la fase 4: el satélite sin corregir debe dar lo mismo que en
# salidas_fase4/ablacion_resultados.csv (si no, las curvas o la partición no son las mismas)
abl4 <- read.csv(file.path(".", "salidas_fase4", "ablacion_resultados.csv"), fileEncoding = "UTF-8")
cat("Satélite sin corregir: validación", round(L2rel(err_sat(i_va), i_va), 2),
    "% | prueba", round(L2rel(err_sat(i_te), i_te), 2), "%  (fase 4:",
    round(abl4$val_media[1], 2), "% |", round(abl4$prueba_media[1], 2), "%)\n")

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

# Error en prueba por bloque ENSO (misma división que la fase 4: Niño 1997, Niña 1998-2000, Neutro 2001-2003)
bloques <- list("Niño 1997" = 1997, "Niña 1998-2000" = 1998:2000, "Neutro 2001-2003" = 2001:2003)
tabla_bloq <- do.call(rbind, lapply(names(bloques), function(b) {
  idx <- i_te[claves$YEAR[i_te] %in% bloques[[b]]]
  data.frame(bloque = b, curvas = length(idx), estaciones = length(unique(claves$STATION_ID[idx])),
             t(sapply(err_todos, function(e) round(L2rel(e[idx], idx), 2))), check.names = FALSE)
}))
print(tabla_bloq, row.names = FALSE)
guardar(tabla_bloq, "tabla_f3_19b_regresion_por_bloque_enso")

# IC bootstrap por estación (prueba), con el mismo esquema de la fase 4
est_te <- claves$STATION_ID[i_te]
sum_est <- function(e) tapply(e[i_te], est_te, sum)
E_est  <- lapply(err_todos, sum_est)
nor_te <- tapply(norma_y(i_te), est_te, sum)
cont   <- t(replicate(2000, tabulate(sample(length(nor_te), replace = TRUE), nbins = length(nor_te))))
l2b    <- function(E) 100 * sqrt(as.vector(cont %*% E) / as.vector(cont %*% nor_te))
# Además del agregado (bootstrap) se reporta lo que pasa en la estación TÍPICA (mediana de la diferencia
# de error relativo por estación y % de estaciones donde gana el modelo a), y se repite por subconjunto:
# así se ve si una ganancia viene de toda la zona o de unas pocas estaciones.
ids_te <- as.numeric(names(nor_te))
subconj <- list("Todas"            = rep(TRUE, length(nor_te)),
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
             pct_estaciones_gana_a = round(100 * mean(ra < rb), 1), row.names = NULL)
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
write.csv(pred_ff, file.path(OUT, "pred_coef_regresion_funcional.csv"), row.names = FALSE)   # datos (coeficientes predichos), no tabla del informe

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
          xlab = "Día del año", ylab = expression(R^2*(t)), main = paste0("R²(t) en prueba (", ANIOS("prueba"), ")"), cex.main = 0.95)
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


#####################################################################################
# Figuras para el informe (leyenda dentro del panel del río Cauca, en la esquina superior derecha)
#####################################################################################
# fig_informe_6_1: bandas de d(t) de las dos zonas (los dos paneles de arriba de fig_f3_05).
# Se guardan a 300 dpi para pegarlas en el informe.

# 'espacio' agranda el eje y por arriba (fracción del rango) para que quepa la leyenda sin tapar las curvas
panel_banda <- function(b, main, col = "darkred", espacio = 0) {
  yl <- range(b$inf, b$sup, b$mediana, 0)
  yl[1] <- yl[1] - 0.08 * diff(yl)                           # espacio abajo para las barras de tramos
  yl[2] <- yl[2] + espacio * diff(yl)                        # espacio arriba para la leyenda
  plot(b$t, b$media, type = "n", ylim = yl, xaxs = "i", xlab = "Día del año", ylab = "mm/día",
       main = main, cex.main = 1)
  rect(c(0, 365 - BORDE), -1e6, c(BORDE, 365), 1e6, col = "grey90", border = NA)
  polygon(c(b$t, rev(b$t)), c(b$inf, rev(b$sup)), col = adjustcolor(col, 0.25), border = NA)
  lines(b$t, b$media, lwd = 2.5, col = col)
  lines(b$t, b$mediana, lwd = 2, lty = 2)
  abline(h = 0, col = "grey60")
  y0 <- par("usr")[3] + 0.04 * diff(par("usr")[3:4])
  lines(b$t, ifelse(b$inf > 0, y0, NA), lwd = 5, col = "steelblue", lend = 1)    # subestima
  lines(b$t, ifelse(b$sup < 0, y0, NA), lwd = 5, col = "orange", lend = 1)       # sobreestima
  box()
}
leyenda_informe <- function() {
  legend("topright", c("mediana", "subestima", "sobreestima"),
         col = c("black", "steelblue", "orange"), lty = c(2, 1, 1), lwd = c(2, 5, 5), seg.len = 2.2,
         bty = "n", cex = 0.8, inset = c(0.01, 0.01))
}
panel_diagonal <- function() {
  sel <- c("Cauca - total", "Cauca - anomalías", "Pacífico - total", "Pacífico - anomalías")
  cols <- c("grey20", "grey20", "#1b9e77", "#1b9e77")         # gris: río Cauca; verde: Pacífico
  matplot(t_eval, diagonales[, sel], type = "n", ylim = c(-0.2, 1.75), xaxs = "i", yaxt = "n",
          xlab = "Día del año", ylab = "r(t, t)", main = "a) Diagonal de la correlación cruzada", cex.main = 1)
  axis(2, at = seq(-0.2, 1, 0.2))
  rect(c(0, 365 - BORDE), -2, c(BORDE, 365), 3, col = "grey90", border = NA)
  abline(h = 0, col = "grey60")
  matlines(t_eval, diagonales[, sel], lty = c(1, 2, 1, 2), lwd = 2.5, col = cols)
  legend("topright", c("río Cauca, total", "río Cauca, anomalías", "Pacífico, total", "Pacífico, anomalías"),
         col = cols, lty = c(1, 2, 1, 2), lwd = 2.5, seg.len = 2.2, bty = "n", cex = 0.8,
         y.intersp = 0.9, inset = c(0.01, 0.01))
  box()
}

figura("fig_informe_6_1", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  panel_banda(bandas$Cauca, etiq("Cauca"), espacio = 0.45); leyenda_informe()
  panel_banda(bandas$`Pacífico`, etiq("Pacífico"))
}, ancho = 12, alto = 4.2, res = 300)

# fig_informe_6_2: (a) diagonal de la correlación cruzada r(t, t) por zona; (b) desviación estándar
# funcional de cada fuente en la cuenca del río Cauca. Juntas responden si el satélite sigue a la
# estación (correlación) y con qué amplitud (desviación estándar).
panel_de <- function() {
  de <- sapply(c("Cauca - total", "Cauca - anomalías"), function(a)
    cbind(sqrt(diag(cov_grupo[[a]]$Y)), sqrt(diag(cov_grupo[[a]]$X))))
  de <- matrix(de, ncol = 4)                    # estación total, satélite total, estación anom., satélite anom.
  cols <- c("grey20", "#7570b3", "grey20", "#7570b3")          # gris: estación; morado: satélite
  matplot(t_eval, de, type = "n", ylim = c(0, 1.6 * max(de[t_eval >= BORDE & t_eval <= 365 - BORDE, ])), xaxs = "i",
          xlab = "Día del año", ylab = "mm/día", main = "b) Desviación estándar - cuenca del río Cauca", cex.main = 1)
  rect(c(0, 365 - BORDE), -1e3, c(BORDE, 365), 1e3, col = "grey90", border = NA)
  matlines(t_eval, de, lty = c(1, 1, 2, 2), lwd = 2.5, col = cols)
  legend("topright", c("estación, total", "satélite, total", "estación, anomalías", "satélite, anomalías"),
         col = cols, lty = c(1, 1, 2, 2), lwd = 2.5, seg.len = 2.2, bty = "n", cex = 0.8,
         y.intersp = 0.9, inset = c(0.01, 0.01))
  box()
}

figura("fig_informe_6_2", {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  panel_diagonal()
  panel_de()
}, ancho = 12, alto = 4.2, res = 300)


# ===================================================================================
# PARTE 6. Evaluación  (ítem 4.4)
# ===================================================================================

limpiar()
set.seed(SEMILLA)


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


# =============================================================================
# Fin
# =============================================================================
cat("\nProyecto completo. Tablas en salidas_preprocesamiento, salidas_fase1, salidas_fase2,",
    "salidas_fase3, salidas_fase4 y salidas_evaluacion.\n")
print(sessionInfo())
