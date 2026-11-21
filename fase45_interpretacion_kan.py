# Fase 4.5 - Interpretación climatológica de la KAN
# Proyecto final de Analítica de Datos: lluvia del Valle del Cauca, 1981-2003
#
# Lee los modelos guardados por la fase 4 (salidas_fase4/modelos/) y produce:
#   fig07  w_k(t) de la KAN y de la cabeza lineal (ventana del año que usa cada neurona)
#   fig08  activaciones de la 1.a capa para ONI y altitud, con el ajuste simbólico
#   fig09  efectos parciales: ONI por época y región, altitud y escala del satélite (KAN vs lineal)
#   fig10  peores casos de prueba con sus curvas
#   tablas interp_*.csv en salidas_fase4/
#
# Los .pt se leen con numpy (no hace falta torch para esta parte). La red se reproduce con las
# mismas operaciones de fase4_red_funcional.py; la función verificar_reproduccion() compara contra
# pred_coef_*.csv y debe dar una diferencia del orden de 1e-14.
#
# Uso: python fase45_interpretacion_kan.py   (desde la carpeta del repositorio)
# o, dentro del script de la fase 4: interpretacion_kan()

import collections
import pickle
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DIR_F1 = Path("./salidas_fase1")
DIR_F4 = Path("./salidas_fase4")
DIR_FIG = DIR_F4 / "figuras"
SEMILLAS_45 = (0, 1, 2)
NOMBRES_Z = [f"z{k+1}" for k in range(8)] + ["altitud", "ONI"]
NO_CONFIABLES_45 = [5300800101, 5300900102, 5300600101, 5300000203, 2621500105, 5300000202,
                    5301300101, 5301000103, 5400000201, 2621500104, 5300900101]
MESES = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
NOM_MES = ["ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"]
EPOCAS = {"Dic-Feb": np.r_[0:59, 334:365], "Mar-May": np.r_[59:151],
          "Jun-Ago": np.r_[151:243], "Sep-Nov": np.r_[243:334]}


# 1. Lectura de los modelos y reproducción de la red en numpy

def cargar_pt(ruta):
    """state_dict de PyTorch (.pt) como diccionario de arreglos numpy."""
    z = zipfile.ZipFile(ruta)
    raiz = z.namelist()[0].split("/")[0]
    tipos = {"DoubleStorage": np.float64, "FloatStorage": np.float32, "LongStorage": np.int64}

    class Lector(pickle.Unpickler):
        def find_class(self, modulo, nombre):
            if nombre == "_rebuild_tensor_v2":
                def tensor(storage, offset, size, stride, *otros):
                    if len(size) == 0:
                        return storage[offset].copy()
                    pasos = [s * storage.itemsize for s in stride]
                    return np.lib.stride_tricks.as_strided(storage[offset:], tuple(size), pasos).copy()
                return tensor
            if nombre == "OrderedDict":
                return collections.OrderedDict
            if modulo.startswith("torch") and nombre.endswith("Storage"):
                return nombre
            return super().find_class(modulo, nombre)

        def persistent_load(self, pid):
            _, tipo, clave, _, _ = pid
            return np.frombuffer(z.read(f"{raiz}/data/{clave}"), dtype=tipos[tipo])

    return Lector(z.open(f"{raiz}/data.pkl")).load()


def base_bspline_45(t, K, orden=4, rango=(0.0, 365.0)):
    """Base B-spline con nudos equiespaciados (igual que la fase 1), recursión de Cox-de Boor."""
    interiores = np.linspace(rango[0], rango[1], K - orden + 2)[1:-1]
    u = np.r_[[rango[0]] * orden, interiores, [rango[1]] * orden]
    B = np.zeros((len(t), len(u) - 1))
    for j in range(len(u) - 1):
        B[:, j] = (t >= u[j]) & (t < u[j + 1])
    B[t == rango[1], K - 1] = 1.0
    for p in range(1, orden):
        nueva = np.zeros((len(t), len(u) - 1 - p))
        for j in range(len(u) - 1 - p):
            a = (t - u[j]) / (u[j + p] - u[j]) * B[:, j] if u[j + p] > u[j] else 0
            b = (u[j + p + 1] - t) / (u[j + p + 1] - u[j + 1]) * B[:, j + 1] if u[j + p + 1] > u[j + 1] else 0
            nueva[:, j] = a + b
        B = nueva
    return B


def silu(x):
    return x / (1 + np.exp(-x))


def bases_kan(x, nudos, k=3):
    u = nudos
    x = x[..., None]
    B = ((x >= u[:-1]) & (x < u[1:])).astype(float)
    for p in range(1, k + 1):
        B = ((x - u[:-(p + 1)]) / (u[p:-1] - u[:-(p + 1)]) * B[..., :-1]
             + (u[p + 1:] - x) / (u[p + 1:] - u[1:-p]) * B[..., 1:])
    return B


def capa_kan(sd, pre, x, layernorm=False):
    if layernorm:
        m, v = x.mean(1, keepdims=True), x.var(1, keepdims=True)
        x = (x - m) / np.sqrt(v + 1e-5) * sd[pre + "norm.weight"] + sd[pre + "norm.bias"]
    base = silu(x) @ sd[pre + "w_b"].T
    spline = np.einsum("big,oig->bo", bases_kan(x, sd[pre + "nudos"]), sd[pre + "c"] * sd[pre + "w_s"][..., None])
    return base + spline


def arista(sd, o, i, x):
    """phi_{o,i}(x) de la 1.a capa de la KAN (entrada i hacia la neurona oculta o)."""
    pre = "cabeza.capa1."
    return sd[pre + "w_b"][o, i] * silu(x) + bases_kan(x, sd[pre + "nudos"]) @ (sd[pre + "c"][o, i] * sd[pre + "w_s"][o, i])


def codificar(sd, C_sat, extra_std):
    z = np.column_stack([C_sat @ (sd["cod.A"] @ sd["cod.W"]).T, extra_std])
    return (z - sd["cod.norm.running_mean"]) / np.sqrt(sd["cod.norm.running_var"] + 1e-5)


def predecir(tipo, sd, C_sat, extra_std):
    """Coeficientes de la curva de la estación (base K = 8) que predice la red."""
    z = codificar(sd, C_sat, extra_std)
    if tipo == "kan":
        s = capa_kan(sd, "cabeza.capa2.", capa_kan(sd, "cabeza.capa1.", z), layernorm=True)
    else:
        s = z @ sd["cabeza.weight"].T + sd["cabeza.bias"]
    return sd["dec.mu"] + s @ sd["dec.Psi"].T


def cargar_datos_45():
    ce = pd.read_csv(DIR_F1 / "coeficientes_estaciones.csv")
    cs = pd.read_csv(DIR_F1 / "coeficientes_satelite.csv")
    ranking = pd.read_excel("./salidas_preprocesamiento/ranking_estaciones.xlsx")
    enso = pd.read_csv(DIR_F4 / "fase_enso_por_anio.csv").set_index("YR")
    altitud = ce["STATION_ID"].map(ranking.set_index("STATION_ID")["altitud"]).to_numpy(float)
    oni = ce["YEAR"].map(enso["ONI_medio"]).to_numpy(float)
    extra = np.column_stack([altitud, oni])
    tr = (ce["particion"] == "entrenamiento").to_numpy()
    d = dict(ce=ce, enso=enso, extra=extra, tr=tr, te=(ce["particion"] == "prueba").to_numpy(),
             media_extra=extra[tr].mean(0), sd_extra=extra[tr].std(0),
             C_est=ce[[c for c in ce.columns if c.startswith("c")]].to_numpy(float),
             C_sat=cs[[c for c in cs.columns if c.startswith("c")]].to_numpy(float),
             region=np.where(ce["STATION_ID"].astype(str).str[:2] == "26", "Cauca", "Pacífico"))
    d["extra_std"] = (extra - d["media_extra"]) / d["sd_extra"]
    d["t"] = np.arange(365) + 0.5
    d["P_est"] = base_bspline_45(d["t"], d["C_est"].shape[1])
    d["P_sat"] = base_bspline_45(d["t"], d["C_sat"].shape[1])
    d["modelos"] = {(tp, s): cargar_pt(DIR_F4 / "modelos" / f"red_{tp}_semilla{s}.pt")
                    for tp in ("kan", "lineal") for s in SEMILLAS_45}
    return d


def verificar_reproduccion(d):
    filas = []
    for tp in ("kan", "lineal"):
        for s in SEMILLAS_45:
            guardado = pd.read_csv(DIR_F4 / f"pred_coef_{tp}_semilla{s}.csv")
            guardado = guardado[[c for c in guardado.columns if c.startswith("c")]].to_numpy(float)
            dif = np.abs(predecir(tp, d["modelos"][(tp, s)], d["C_sat"], d["extra_std"]) - guardado).max()
            filas.append({"cabeza": tp, "semilla": s, "max_dif_vs_pred_coef": dif})
    tabla = pd.DataFrame(filas)
    print("Reproducción de la red en numpy vs pred_coef_*.csv:\n", tabla.to_string(index=False))
    assert tabla["max_dif_vs_pred_coef"].max() < 1e-8
    return tabla


# 2. Ventanas w_k(t)

def fecha(dia):
    m = np.searchsorted(MESES, dia, side="right") - 1
    return f"{int(dia - MESES[m]) + 1} {NOM_MES[m]}"


def tramos(mascara, t):
    """Intervalos contiguos donde mascara es True, como texto '1 ene-20 feb; ...'."""
    idx = np.where(mascara)[0]
    if len(idx) == 0:
        return ""
    cortes = np.where(np.diff(idx) > 1)[0]
    ini = np.r_[idx[0], idx[cortes + 1]]
    fin = np.r_[idx[cortes], idx[-1]]
    return "; ".join(f"{fecha(t[a])}-{fecha(t[b])}" for a, b in zip(ini, fin))


def ventanas_wk(d):
    filas = []
    for tp in ("kan", "lineal"):
        for s in SEMILLAS_45:
            W = d["P_sat"] @ d["modelos"][(tp, s)]["cod.A"].T          # 365 x 8
            for k in range(W.shape[1]):
                w = W[:, k]
                i = np.argmax(np.abs(w))
                signo = np.sign(w[i])
                filas.append({"cabeza": tp, "semilla": s, "neurona": f"w{k+1}",
                              "pico": fecha(d["t"][i]), "signo_pico": "+" if signo > 0 else "-",
                              "ventana_50pct": tramos(signo * w >= 0.5 * np.abs(w[i]), d["t"]),
                              "dias_ventana": int((signo * w >= 0.5 * np.abs(w[i])).sum())})
    return pd.DataFrame(filas)


def figura_wk(d, semilla=0):
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.2), sharey=True)
    for a, tp in zip(ax, ("kan", "lineal")):
        W = d["P_sat"] @ d["modelos"][(tp, semilla)]["cod.A"].T
        for (ini, fin) in [(90, 151), (273, 334)]:
            a.axvspan(ini, fin, color="#9ecae1", alpha=0.35, lw=0)
        for k in range(W.shape[1]):
            a.plot(d["t"], W[:, k] / np.abs(W[:, k]).max(), lw=1.4, label=f"w{k+1}")
        a.axhline(0, color="k", lw=0.5)
        a.set_xticks(MESES[:-1] + 15, NOM_MES)
        a.set_title(f"w_k(t), cabeza {tp} (semilla {semilla})")
    ax[0].set_ylabel("w_k(t) escalada a máx. 1")
    ax[1].legend(fontsize=7, ncol=2)
    fig.text(0.5, -0.02, "Sombra: picos de lluvia abr-may y oct-nov", ha="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(DIR_FIG / "fase4_fig07_wk_kan_vs_lineal.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 3. Activaciones de la 1.a capa: importancia y ajuste simbólico

CANDIDATAS = {
    "lineal": lambda x, a, b: a * x + b,
    "cuadrática": lambda x, a, b, c: a * x ** 2 + b * x + c,
    "tanh": lambda x, a, b, c, e: a * np.tanh(b * (x - c)) + e,
    "exponencial": lambda x, a, b, c: a * np.exp(b * x) + c,
    "seno": lambda x, a, b, c, e: a * np.sin(b * x + c) + e,
}


def gauss_newton(f, x, y, p0, iteraciones=200):
    """Mínimos cuadrados no lineales por Gauss-Newton con amortiguamiento (Levenberg-Marquardt).
    El jacobiano se aproxima por diferencias finitas."""
    p = np.asarray(p0, float)
    lam = 1e-3
    sce = np.sum((y - f(x, *p)) ** 2)
    for _ in range(iteraciones):
        r = y - f(x, *p)
        J = np.empty((len(x), len(p)))
        for j in range(len(p)):
            h = 1e-6 * max(1.0, abs(p[j]))
            dp = p.copy()
            dp[j] += h
            J[:, j] = (f(x, *dp) - f(x, *p)) / h
        A = J.T @ J
        paso = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-12), J.T @ r)
        nuevo = p + paso
        sce_nuevo = np.sum((y - f(x, *nuevo)) ** 2)
        if np.isfinite(sce_nuevo) and sce_nuevo < sce:
            p, lam = nuevo, lam / 10
            if sce - sce_nuevo < 1e-12 * max(sce, 1e-12):
                break
            sce = sce_nuevo
        else:
            lam *= 10
            if lam > 1e10:
                break
    return p


def inicios(nombre, x, y):
    """Valores iniciales tomados de los datos (varios, para no quedar en un mínimo local)."""
    amp, centro, nivel = (y.max() - y.min()) / 2, np.median(x), y.mean()
    escala = 1 / (x.std() + 1e-12)
    if nombre == "tanh":
        return [[s * amp, b * escala, centro, nivel] for s in (1, -1) for b in (0.5, 1, 2, 4)]
    if nombre == "exponencial":
        return [[s * amp, b * escala, nivel] for s in (1, -1) for b in (-1, -0.3, 0.3, 1)]
    if nombre == "seno":
        return [[amp, w * escala, fase, nivel] for w in (0.5, 1, 2) for fase in (0, np.pi / 2, np.pi, -np.pi / 2)]
    return []


def ajuste_simbolico(x, y):
    """Ajusta cada expresión candidata por mínimos cuadrados y devuelve R² y parámetros.
    Lineal y cuadrática tienen solución cerrada; las demás se ajustan con gauss_newton()."""
    sct = np.sum((y - y.mean()) ** 2)
    filas = []
    for nombre, f in CANDIDATAS.items():
        if nombre in ("lineal", "cuadrática"):
            grado = 1 if nombre == "lineal" else 2
            X = np.vander(x, grado + 1)
            p = np.linalg.lstsq(X, y, rcond=None)[0]
            candidatos = [p]
        else:
            candidatos = [gauss_newton(f, x, y, p0) for p0 in inicios(nombre, x, y)]
        r2 = [1 - np.sum((y - f(x, *p)) ** 2) / sct for p in candidatos]
        k = int(np.nanargmax(r2))
        filas.append({"expresion": nombre, "R2": r2[k], "parametros": np.round(candidatos[k], 4).tolist()})
    return sorted(filas, key=lambda r: -r["R2"])


def activaciones(d):
    """Importancia de cada arista (sd de phi sobre las z de entrenamiento) y ajuste simbólico de
    las aristas más fuertes de ONI, altitud y de la neurona funcional principal."""
    imp, simb = [], []
    for s in SEMILLAS_45:
        sd = d["modelos"][("kan", s)]
        Z = codificar(sd, d["C_sat"], d["extra_std"])[d["tr"]]
        M = np.array([[np.std(arista(sd, o, i, Z[:, i])) for i in range(10)] for o in range(4)])
        for i in range(10):
            imp.append({"semilla": s, "entrada": NOMBRES_Z[i], "importancia_total": M[:, i].sum(),
                        "oculta_mas_fuerte": int(np.argmax(M[:, i])) + 1, "importancia_max": M[:, i].max()})
        o_f, i_f = np.unravel_index(np.argmax(M[:, :8]), M[:, :8].shape)
        for i, o in [(9, int(np.argmax(M[:, 9]))), (8, int(np.argmax(M[:, 8]))), (int(i_f), int(o_f))]:
            x = Z[:, i]
            for r in ajuste_simbolico(x, arista(sd, o, i, x)):
                simb.append({"semilla": s, "arista": f"{NOMBRES_Z[i]} -> oculta {o+1}", **r})
    imp = pd.DataFrame(imp)
    imp["rango_en_semilla"] = imp.groupby("semilla")["importancia_total"].rank(ascending=False).astype(int)
    return imp, pd.DataFrame(simb)


def oni_en_grados(d, sd, z):
    """Pasa un valor de la entrada ONI ya normalizada (z) a grados de anomalía ONI."""
    m, v = sd["cod.norm.running_mean"][9], sd["cod.norm.running_var"][9]
    return (z * np.sqrt(v + 1e-5) + m) * d["sd_extra"][1] + d["media_extra"][1]


def figura_activaciones(d, simb):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    f_tanh = CANDIDATAS["tanh"]
    colores = ["#1b9e77", "#d95f02", "#7570b3"]
    for s, c in zip(SEMILLAS_45, colores):
        sd = d["modelos"][("kan", s)]
        Z = codificar(sd, d["C_sat"], d["extra_std"])[d["tr"]]
        for a, i in zip(ax[:2], (9, 8)):
            lo, hi = np.percentile(Z[:, i], [1, 99])
            xm = np.linspace(lo, hi, 200)
            fila = simb[(simb.semilla == s) & simb.arista.str.startswith(NOMBRES_Z[i] + " ")]
            o = int(fila.arista.iloc[0].split()[-1]) - 1
            eje_x = oni_en_grados(d, sd, xm) if i == 9 else xm
            a.plot(eje_x, arista(sd, o, i, xm), color=c, lw=2, label=f"semilla {s}, oculta {o+1}")
            if i == 9:
                ft = fila[fila.expresion == "tanh"].iloc[0]
                a.plot(eje_x, f_tanh(xm, *ft.parametros), "--", color="k", lw=0.9)
                a.text(0.02, 0.05 + 0.08 * s, f"s{s}: tanh R² = {ft.R2:.3f}", transform=a.transAxes, fontsize=8, color=c)
    ax[0].set(xlabel="ONI medio del año (°C)", ylabel="φ(x)", title="Arista más fuerte del ONI (- - ajuste tanh)")
    ax[1].set(xlabel="altitud normalizada", title="Arista más fuerte de la altitud")
    ax[0].legend(fontsize=7, loc="upper right")
    ax[1].legend(fontsize=7)
    imp = simb.drop_duplicates(["semilla", "arista"])
    ax[2].axis("off")
    texto = "Mejor expresión por arista (R²)\n\n" + "\n".join(
        f"s{r.semilla}  {r.arista:22s} {r.expresion:11s} {r.R2:.3f}" for r in imp.itertuples())
    ax[2].text(0, 1, texto, va="top", family="monospace", fontsize=8)
    plt.tight_layout()
    plt.savefig(DIR_FIG / "fase4_fig08_activaciones_ajuste_simbolico.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 4. Efectos parciales (lo que la red le hace a la lluvia) y contraste con la cabeza lineal
# El signo de una arista de la 1.a capa no es el efecto sobre la lluvia, porque pasa por la 2.a capa
# y por el decodificador. Por eso se mueve una entrada y se mira el cambio en la curva predicha.

def curva_media(d, tipo, C_sat, extra_std):
    return np.mean([predecir(tipo, d["modelos"][(tipo, s)], C_sat, extra_std) @ d["P_est"].T
                    for s in SEMILLAS_45], axis=0)


def efectos_parciales(d):
    tr, filas = d["tr"], []
    curvas_oni = {}
    for tp in ("kan", "lineal"):
        for rg in ("Cauca", "Pacífico"):
            m = tr & (d["region"] == rg)
            alto, bajo = d["extra_std"][m].copy(), d["extra_std"][m].copy()
            alto[:, 1], bajo[:, 1] = 1.0, -1.0
            y_alto = curva_media(d, tp, d["C_sat"][m], alto)
            y_bajo = curva_media(d, tp, d["C_sat"][m], bajo)
            efecto = ((y_alto - y_bajo) / 2).mean(0)                      # por +1 sd de ONI
            curvas_oni[(tp, rg)] = efecto
            base = curva_media(d, tp, d["C_sat"][m], d["extra_std"][m]).mean(0)
            for ep, dias in {"Año": np.arange(365), **EPOCAS}.items():
                y0 = base[dias]
                filas.append({"efecto": "ONI +1 sd", "cabeza": tp, "region": rg, "epoca": ep,
                              "cambio_sqrt": efecto[dias].mean(),
                              "cambio_pct_aprox": 100 * (np.mean((y0 + efecto[dias]) ** 2) / np.mean(y0 ** 2) - 1)})
    # altitud: respuesta media anual (sqrt) en una malla de altitudes
    alturas = np.array([100, 250, 500, 750, 1000, 1500, 2000, 2500])
    for tp in ("kan", "lineal"):
        for h in alturas:
            e = d["extra_std"][tr].copy()
            e[:, 0] = (h - d["media_extra"][0]) / d["sd_extra"][0]
            filas.append({"efecto": f"altitud {h} m", "cabeza": tp, "region": "todas", "epoca": "Año",
                          "cambio_sqrt": curva_media(d, tp, d["C_sat"][tr], e).mean(), "cambio_pct_aprox": np.nan})
    # escala del satélite: se multiplica la curva satelital por f
    factores = np.array([0.6, 0.8, 1.0, 1.2, 1.4])
    for tp in ("kan", "lineal"):
        for rg in ("Cauca", "Pacífico"):
            m = tr & (d["region"] == rg)
            for f in factores:
                filas.append({"efecto": f"satélite x{f}", "cabeza": tp, "region": rg, "epoca": "Año",
                              "cambio_sqrt": curva_media(d, tp, d["C_sat"][m] * f, d["extra_std"][m]).mean(),
                              "cambio_pct_aprox": np.nan})
    return pd.DataFrame(filas), curvas_oni, alturas, factores


def figura_efectos(d, tabla, curvas_oni, alturas, factores):
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.3))
    estilos = {("kan", "Cauca"): ("#d95f02", "-"), ("kan", "Pacífico"): ("#1b9e77", "-"),
               ("lineal", "Cauca"): ("k", "--")}
    for clave, (c, ls) in estilos.items():
        etiqueta = "lineal (igual en las dos regiones)" if clave[0] == "lineal" else f"KAN {clave[1]}"
        ax[0].plot(d["t"], curvas_oni[clave], color=c, ls=ls, lw=2, label=etiqueta)
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].set_xticks(MESES[:-1] + 15, NOM_MES)
    ax[0].set(ylabel="cambio en sqrt(mm/día)", title="Efecto de +1 sd de ONI (≈ +0,6 °C)")
    ax[0].legend(fontsize=8)
    for tp, c in (("kan", "#7570b3"), ("lineal", "k")):
        v = tabla[(tabla.cabeza == tp) & tabla.efecto.str.startswith("altitud")]["cambio_sqrt"]
        ax[1].plot(alturas, v, "o-", color=c, label=tp)
    alt_tr = d["extra"][d["tr"], 0]
    ax[1].plot(alt_tr, np.full(len(alt_tr), ax[1].get_ylim()[0]), "|", color="gray", alpha=0.3)
    ax[1].set(xlabel="altitud (m)", ylabel="media anual predicha, sqrt(mm/día)", title="Respuesta a la altitud")
    ax[1].legend(fontsize=8)
    for (tp, rg), (c, ls) in {("kan", "Cauca"): ("#d95f02", "-"), ("kan", "Pacífico"): ("#1b9e77", "-"),
                              ("lineal", "Cauca"): ("#d95f02", "--"), ("lineal", "Pacífico"): ("#1b9e77", "--")}.items():
        v = tabla[(tabla.cabeza == tp) & (tabla.region == rg) & tabla.efecto.str.startswith("satélite")]["cambio_sqrt"]
        ax[2].plot(factores, v, "o", color=c, ls=ls, label=f"{tp} {rg}")
    ax[2].set(xlabel="factor sobre la curva satelital", ylabel="media anual predicha, sqrt(mm/día)",
              title="Respuesta a la intensidad del satélite")
    ax[2].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(DIR_FIG / "fase4_fig09_efectos_parciales_kan_vs_lineal.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 5. Peores casos de prueba

def peores_casos(d, n=8):
    G = d["P_est"].T @ d["P_est"]
    def l2rel(C_hat):
        e = C_hat - d["C_est"]
        return 100 * np.sqrt(np.einsum("bi,ij,bj->b", e, G, e) / np.einsum("bi,ij,bj->b", d["C_est"], G, d["C_est"]))
    err = {tp: np.mean([l2rel(predecir(tp, d["modelos"][(tp, s)], d["C_sat"], d["extra_std"])) for s in SEMILLAS_45], 0)
           for tp in ("kan", "lineal")}
    obs = (d["C_est"] @ d["P_est"].T).mean(1)
    sat = (d["C_sat"] @ d["P_sat"].T).mean(1)
    kan = curva_media(d, "kan", d["C_sat"], d["extra_std"]).mean(1)
    ce = d["ce"]
    tabla = pd.DataFrame({"STATION_ID": ce.STATION_ID, "YEAR": ce.YEAR, "region": d["region"],
                          "altitud": d["extra"][:, 0], "fase_enso": ce.YEAR.map(d["enso"]["fase_enso"]),
                          "no_confiable": ce.STATION_ID.isin(NO_CONFIABLES_45),
                          "L2rel_kan": err["kan"], "L2rel_lineal": err["lineal"],
                          "media_obs_sqrt": obs, "media_kan_sqrt": kan, "media_sat_sqrt": sat})[d["te"]]
    tabla = tabla.sort_values("L2rel_kan", ascending=False)
    print(f"Mediana del L2 relativo de la KAN en prueba: {tabla.L2rel_kan.median():.1f} %")
    print(f"Lluvia media (sqrt) de las {n} peores: {tabla.media_obs_sqrt.head(n).mean():.2f}"
          f" | resto de prueba: {tabla.media_obs_sqrt.iloc[n:].mean():.2f}")
    return tabla


def figura_peores(d, tabla, n=6):
    fig, ax = plt.subplots(2, 3, figsize=(15, 6.5), sharex=True)
    for a, fila in zip(ax.flat, tabla.head(n).itertuples()):
        i = np.where((d["ce"].STATION_ID == fila.STATION_ID) & (d["ce"].YEAR == fila.YEAR))[0]
        a.plot(d["t"], (d["C_est"][i] @ d["P_est"].T)[0], "k", lw=2, label="estación")
        a.plot(d["t"], curva_media(d, "kan", d["C_sat"][i], d["extra_std"][i])[0], color="#7570b3", label="KAN")
        a.plot(d["t"], curva_media(d, "lineal", d["C_sat"][i], d["extra_std"][i])[0], color="#d95f02", ls="--", label="lineal")
        a.plot(d["t"], (d["C_sat"][i] @ d["P_sat"].T)[0], color="gray", ls=":", label="satélite")
        a.set_title(f"{fila.STATION_ID} {fila.YEAR} ({fila.region}, {fila.altitud:.0f} m, {fila.fase_enso})\n"
                    f"L2 KAN {fila.L2rel_kan:.0f} %, lineal {fila.L2rel_lineal:.0f} %", fontsize=8)
        a.set_xticks(MESES[:-1:2] + 15, NOM_MES[::2])
    ax[0, 0].set_ylabel("sqrt(mm/día)")
    ax[1, 0].set_ylabel("sqrt(mm/día)")
    ax[0, 0].legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(DIR_FIG / "fase4_fig10_peores_casos.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 6. Todo junto

def interpretacion_kan():
    DIR_FIG.mkdir(parents=True, exist_ok=True)
    d = cargar_datos_45()
    verificar_reproduccion(d)

    tabla_w = ventanas_wk(d)
    tabla_w.to_csv(DIR_F4 / "interp_ventanas_wk.csv", index=False)
    figura_wk(d)

    imp, simb = activaciones(d)
    imp.to_csv(DIR_F4 / "interp_importancia_aristas.csv", index=False)
    simb.to_csv(DIR_F4 / "interp_ajuste_simbolico.csv", index=False)
    figura_activaciones(d, simb)
    print("\nAjuste simbólico, mejor expresión por arista:\n",
          simb.drop_duplicates(["semilla", "arista"])[["semilla", "arista", "expresion", "R2"]].round(3).to_string(index=False))
    print("\nAjuste tanh de la arista del ONI:\n",
          simb[simb.arista.str.startswith("ONI") & (simb.expresion == "tanh")][["semilla", "arista", "R2", "parametros"]].to_string(index=False))

    efectos, curvas_oni, alturas, factores = efectos_parciales(d)
    efectos.to_csv(DIR_F4 / "interp_efectos_parciales.csv", index=False)
    figura_efectos(d, efectos, curvas_oni, alturas, factores)
    print("\nEfecto de +1 sd de ONI:\n",
          efectos[efectos.efecto == "ONI +1 sd"].round(3).to_string(index=False))

    peores = peores_casos(d)
    peores.to_csv(DIR_F4 / "interp_peores_casos.csv", index=False)
    figura_peores(d, peores)
    print(peores.head(8).round(2).to_string(index=False))
    print("\nFase 4.5 lista. Figuras fase4_fig07 a fig10 y tablas interp_*.csv en", DIR_F4)


if __name__ == "__main__":
    interpretacion_kan()
