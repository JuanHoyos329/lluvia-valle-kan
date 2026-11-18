# =============================================================================
# Fase 4 - Red neuronal funcional y ablación de la cabeza (ítems 4.2 y 4.3)
# Proyecto final de Analítica de Datos: lluvia del Valle del Cauca, 1981-2003
#
# Qué hace este script (se irá completando por bloques):
#   Bloque 1. Bases B-spline, matrices de Gram (W) y de rugosidad (R) por cuadratura,
#             y FPCA de las curvas de estaciones con la métrica de W (solo entrenamiento).
#   Bloque 2. Red funcional en PyTorch: codificador z = A W c, cabeza y decodificador FPCA.
#   Bloque 3. Cabeza lineal: demostración y verificación numérica de la equivalencia con
#             la regresión funcional función a función.
#   Bloque 4. KAN (B-splines en las aristas) verificada con una función sintética.
#   Bloque 5. Ablación: lineal vs KAN vs MLP, parámetros igualados, 3 semillas.
#
# Entradas (salidas de la fase 1, Daniela):
#   salidas_fase1/coeficientes_estaciones.csv  -> c de las estaciones (K = 8), escala sqrt(mm/día)
#   salidas_fase1/coeficientes_satelite.csv    -> c del satélite (K = 12), misma escala
#   salidas_fase1/base_info.csv                -> K, orden y nudos de cada base
# La fila i de los dos archivos de coeficientes es la misma estación en el mismo año.
#
# Partición: entrenamiento 1981-1995 / validación 1996-2000 / prueba 2001-2003.
# Todo lo que se AJUSTA (FPCA, estandarizaciones, sigma^2 del paso a mm) usa solo entrenamiento.
#
# Uso (desde la carpeta del repositorio, con el .venv activo):
#   python fase4_red_funcional.py
# Las gráficas se muestran en pantalla y no se guardan (misma convención del grupo).
# =============================================================================

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SEED = 2026
np.random.seed(SEED)

# True: muestra cada gráfica (hay que cerrar la ventana para que el script siga).
# False: no las muestra; útil mientras se itera porque el script corre de corrido.
MOSTRAR_GRAFICAS = False

# True: entrena todo desde cero (lo que se entrega). False: si ya existen los modelos guardados en
# salidas_fase4/modelos/, los carga en vez de reentrenar (para iterar rápido sobre los bloques 6+).
REENTRENAR = True

def mostrar(nombre):
    """Guarda la figura actual en salidas_fase4/figuras/<nombre>.png (para el informe) y la muestra
    en pantalla si MOSTRAR_GRAFICAS es True."""
    carpeta = Path("./salidas_fase4/figuras")
    carpeta.mkdir(parents=True, exist_ok=True)
    plt.savefig(carpeta / f"{nombre}.png", dpi=150, bbox_inches="tight")
    if MOSTRAR_GRAFICAS:
        plt.show()
    else:
        plt.close("all")

ENT1 = Path("./salidas_fase1")          # coeficientes de la fase 1
OUT = Path("./salidas_fase4")           # salidas de esta fase
OUT.mkdir(exist_ok=True)

ORDEN = 4                               # B-splines cúbicos (igual que la fase 1)
RANGO = (0.0, 365.0)                    # dominio anual común


###############################################################################
# BLOQUE 1. Bases, matrices de Gram y de rugosidad, y FPCA con la métrica W
###############################################################################

# -----------------------------------------------------------------------------
# 1.1 Carga de coeficientes
# -----------------------------------------------------------------------------
coef_est = pd.read_csv(ENT1 / "coeficientes_estaciones.csv")
coef_sat = pd.read_csv(ENT1 / "coeficientes_satelite.csv")
base_info = pd.read_csv(ENT1 / "base_info.csv")

# Misma curva en la misma fila (lo garantiza la fase 1; se vuelve a comprobar)
assert (coef_est["STATION_ID"].values == coef_sat["STATION_ID"].values).all()
assert (coef_est["YEAR"].values == coef_sat["YEAR"].values).all()

cols_est = [c for c in coef_est.columns if c.startswith("c")]
cols_sat = [c for c in coef_sat.columns if c.startswith("c")]
C_est = coef_est[cols_est].to_numpy(np.float64)          # N x 8   (salida que se quiere reconstruir)
C_sat = coef_sat[cols_sat].to_numpy(np.float64)          # N x 12  (entrada de la red)
K_est, K_sat = C_est.shape[1], C_sat.shape[1]
particion = coef_est["particion"].to_numpy()
tr = particion == "entrenamiento"
va = particion == "validacion"
te = particion == "prueba"
print(f"Curvas: {len(particion)} | entrenamiento {tr.sum()} | validación {va.sum()} | prueba {te.sum()}")
print(f"K estaciones = {K_est} | K satélite = {K_sat}")


# -----------------------------------------------------------------------------
# 1.2 Base B-spline y sus derivadas (Cox-de Boor, la misma construcción de la fase 1)
# -----------------------------------------------------------------------------
# Nudos: extremos repetidos ORDEN veces y nudos interiores equiespaciados en [0, 365].
def nudos_bspline(K, orden=ORDEN, rango=RANGO):
    internos = np.linspace(rango[0], rango[1], K - orden + 2)[1:-1]
    return np.r_[[rango[0]] * orden, internos, [rango[1]] * orden]


def _bspline_orden(t, u, m):
    """Todas las B-splines de orden m sobre el vector de nudos u, evaluadas en t (n x (len(u)-m))."""
    t = np.asarray(t, dtype=np.float64)
    B = np.zeros((len(t), len(u) - 1))
    for j in range(len(u) - 1):
        B[:, j] = (t >= u[j]) & (t < u[j + 1])
    B[t == u[-1], np.max(np.where(u < u[-1]))] = 1.0      # el extremo derecho entra en el último tramo
    for k in range(2, m + 1):
        Bn = np.zeros((len(t), len(u) - k))
        for j in range(len(u) - k):
            d1, d2 = u[j + k - 1] - u[j], u[j + k] - u[j + 1]
            izq = (t - u[j]) / d1 * B[:, j] if d1 > 0 else 0.0
            der = (u[j + k] - t) / d2 * B[:, j + 1] if d2 > 0 else 0.0
            Bn[:, j] = izq + der
        B = Bn
    return B


def base_bspline(t, K, deriv=0, orden=ORDEN):
    """Matriz n x K de la base (deriv = 0) o de su derivada de orden 'deriv'.

    Derivada de una B-spline de orden m (de Boor):
        B'_{j,m}(t) = (m-1) [ B_{j,m-1}(t)/(u_{j+m-1}-u_j) - B_{j+1,m-1}(t)/(u_{j+m}-u_{j+1}) ]
    Se aplica 'deriv' veces bajando el orden en cada paso.
    """
    u = nudos_bspline(K, orden)
    if deriv == 0:
        return _bspline_orden(t, u, orden)
    # Derivada recursiva: D^d B_{j,m} se arma con D^{d-1} B_{.,m-1}
    def D(d, m):
        if d == 0:
            return _bspline_orden(t, u, m)
        Bm1 = D(d - 1, m - 1)                              # n x (len(u) - m + 1)
        out = np.zeros((len(t), len(u) - m))
        for j in range(len(u) - m):
            d1, d2 = u[j + m - 1] - u[j], u[j + m] - u[j + 1]
            a = Bm1[:, j] / d1 if d1 > 0 else 0.0
            b = Bm1[:, j + 1] / d2 if d2 > 0 else 0.0
            out[:, j] = (m - 1) * (a - b)
        return out
    return D(deriv, orden)


# Los nudos tienen que ser exactamente los de la fase 1 (base_info.csv)
for fuente, K in [("estaciones", K_est), ("satélite", K_sat)]:
    nudos_f1 = np.array(base_info.loc[base_info["fuente"] == fuente, "nudos"].iloc[0].split(), float)
    assert np.allclose(nudos_bspline(K), nudos_f1, atol=1e-3), f"nudos distintos a la fase 1 ({fuente})"


# -----------------------------------------------------------------------------
# 1.3 Matrices de Gram W = ∫ φ φᵀ dt y de rugosidad R = ∫ φ'' φ''ᵀ dt por cuadratura
# -----------------------------------------------------------------------------
# Entre dos nudos consecutivos cada B-spline cúbica es un polinomio de grado 3, así que:
#   φ_i φ_j      es un polinomio de grado 6 en cada tramo
#   φ_i'' φ_j''  es un polinomio de grado 2 en cada tramo
# Gauss-Legendre con q puntos integra EXACTO polinomios de grado <= 2q - 1. Con q = 4 por tramo
# (grado <= 7) las dos integrales son exactas, salvo el redondeo del computador.
def cuadratura_gauss(K, q=4, orden=ORDEN):
    """Nodos y pesos de Gauss-Legendre en cada tramo entre nudos distintos."""
    u = np.unique(nudos_bspline(K, orden))
    x, w = np.polynomial.legendre.leggauss(q)              # nodos y pesos en [-1, 1]
    nodos, pesos = [], []
    for a, b in zip(u[:-1], u[1:]):
        nodos.append((b - a) / 2 * x + (a + b) / 2)        # cambio de variable [-1, 1] -> [a, b]
        pesos.append((b - a) / 2 * w)
    return np.concatenate(nodos), np.concatenate(pesos)


def matriz_gram(K, deriv=0):
    t_q, w_q = cuadratura_gauss(K)
    P = base_bspline(t_q, K, deriv)                        # n_q x K
    return P.T @ (w_q[:, None] * P)                        # Σ_q w_q φ(t_q) φ(t_q)ᵀ


W_est, W_sat = matriz_gram(K_est), matriz_gram(K_sat)
R_est, R_sat = matriz_gram(K_est, deriv=2), matriz_gram(K_sat, deriv=2)

# Verificaciones de W y R (todas deben dar ~0)
def coef_greville(K, orden=ORDEN):
    """Coeficientes que representan exactamente f(t) = t en la base (abscisas de Greville)."""
    u = nudos_bspline(K, orden)
    return np.array([u[j + 1:j + orden].mean() for j in range(K)])

t_fino = np.linspace(*RANGO, 20001)
verif_W = []
for fuente, K, W, R in [("estaciones", K_est, W_est, R_est), ("satélite", K_sat, W_sat, R_sat)]:
    uno = np.ones(K)                                       # f(t) = 1  (las B-splines suman 1)
    g = coef_greville(K)                                   # f(t) = t
    P_f = base_bspline(t_fino, K)
    W_trap = np.trapezoid(P_f[:, :, None] * P_f[:, None, :], t_fino, axis=0)
    verif_W.append({
        "fuente": fuente, "K": K,
        "W_simetrica": np.abs(W - W.T).max(),
        "W_min_autovalor": np.linalg.eigvalsh(W).min(),            # > 0: definida positiva
        "∫1dt - 365": uno @ W @ uno - 365,
        "∫t²dt - 365³/3": g @ W @ g - 365 ** 3 / 3,
        "W - trapecio(malla fina)": np.abs(W - W_trap).max(),
        "rugosidad de f(t)=1": uno @ R @ uno,                       # φ'' de una constante = 0
        "rugosidad de f(t)=t": g @ R @ g,                           # φ'' de una recta = 0
    })
verif_W = pd.DataFrame(verif_W)
with pd.option_context("display.float_format", "{:.2e}".format, "display.width", 200):
    print("\nVerificación de W y R:\n", verif_W.T)


# -----------------------------------------------------------------------------
# 1.4 FPCA de las curvas de estaciones con la métrica de W (SOLO entrenamiento)
# -----------------------------------------------------------------------------
# Cada curva es y_i(t) = φ(t)ᵀ c_i. Con c centrados y Σ_c = cov(c):
#   operador de covarianza:  v(s,t) = φ(s)ᵀ Σ_c φ(t)
#   autofunción ψ(t) = φ(t)ᵀ b  cumple  ∫ v(s,t) ψ(t) dt = λ ψ(s)  <=>  Σ_c W b = λ b
# Con u = W^{1/2} b el problema queda simétrico:  (W^{1/2} Σ_c W^{1/2}) u = λ u,
# y la normalización ∫ψ_m ψ_l = b_mᵀ W b_l = δ_ml sale automática.
def raiz_W(W):
    val, vec = np.linalg.eigh(W)
    return vec @ np.diag(np.sqrt(val)) @ vec.T, vec @ np.diag(1 / np.sqrt(val)) @ vec.T


def fpca_coeficientes(C, W):
    mu = C.mean(0)
    Cc = C - mu
    Sigma = Cc.T @ Cc / (len(C) - 1)
    W12, Wm12 = raiz_W(W)
    lam, U = np.linalg.eigh(W12 @ Sigma @ W12)
    orden_desc = np.argsort(lam)[::-1]
    lam, U = lam[orden_desc], U[:, orden_desc]
    Bpsi = Wm12 @ U                                        # columnas = coeficientes de ψ_m en la base
    # Signo fijo (igual que fit_fpca del curso): ∫ψ_m > 0
    signo = np.sign(np.ones(len(mu)) @ W @ Bpsi)
    signo[signo == 0] = 1
    Bpsi = Bpsi * signo
    return {"mu": mu, "psi": Bpsi, "lam": lam, "explicada": lam / lam.sum()}


fpca = fpca_coeficientes(C_est[tr], W_est)
var_acum = np.cumsum(fpca["explicada"])
tabla_fpca = pd.DataFrame({"componente": np.arange(1, K_est + 1), "autovalor": fpca["lam"],
                           "var_explicada_pct": 100 * fpca["explicada"],
                           "var_acumulada_pct": 100 * var_acum})
print("\nFPCA de estaciones (entrenamiento):\n", tabla_fpca.round(2).to_string(index=False))

# Verificación 1: ortonormalidad en L²  ->  Bψᵀ W Bψ = I
orto = np.abs(fpca["psi"].T @ W_est @ fpca["psi"] - np.eye(K_est)).max()

# Verificación 2: misma FPCA por el camino discreto del curso (SVD sobre una malla fina).
# En una malla con paso dt, ∫ f g ≈ Σ f g dt, así que los autovalores del curso (sing²/(N-1))
# multiplicados por dt deben coincidir con λ.
t_m = np.linspace(*RANGO, 3651)
dt = t_m[1] - t_m[0]
Ym = C_est[tr] @ base_bspline(t_m, K_est).T                        # N_tr x 3651
_, sing, Vt = np.linalg.svd(Ym - Ym.mean(0), full_matrices=False)
lam_malla = sing[:K_est] ** 2 / (tr.sum() - 1) * dt
dif_lam_rel = np.abs(lam_malla - fpca["lam"]).max() / fpca["lam"][0]
# Las autofunciones deben coincidir salvo signo: ψ_m(t) vs Vt_m / sqrt(dt)
Psi_m = base_bspline(t_m, K_est) @ fpca["psi"]                     # 3651 x K
dif_psi = max(min(np.abs(Psi_m[:, m] - Vt[m] / np.sqrt(dt)).max(),
                  np.abs(Psi_m[:, m] + Vt[m] / np.sqrt(dt)).max()) for m in range(4))
print(f"\nVerificación FPCA | máx |BᵀWB - I| = {orto:.1e} | dif. relativa de λ vs malla = {dif_lam_rel:.1e}"
      f" | máx dif ψ_1..4 vs malla = {dif_psi:.1e}")

# Error de reconstrucción L² con M componentes (entrenamiento y validación).
# El FPCA se ajustó con entrenamiento; la validación solo se proyecta (no se reajusta nada).
def error_reconstruccion(C, M):
    S = (C - fpca["mu"]) @ W_est @ fpca["psi"][:, :M]              # puntajes s_m = ∫(y-μ)ψ_m
    E = C - (fpca["mu"] + S @ fpca["psi"][:, :M].T)
    return np.einsum("ij,jk,ik->i", E, W_est, E)                   # ∫ (y - ŷ)² dt por curva

def norma2(C):
    return np.einsum("ij,jk,ik->i", C, W_est, C)                   # ∫ y² dt por curva

tabla_M = pd.DataFrame([{
    "M": M, "var_acumulada_pct": 100 * var_acum[M - 1],
    "L2_rel_entrenamiento_pct": 100 * np.sqrt(error_reconstruccion(C_est[tr], M).sum() / norma2(C_est[tr]).sum()),
    "L2_rel_validacion_pct": 100 * np.sqrt(error_reconstruccion(C_est[va], M).sum() / norma2(C_est[va]).sum()),
} for M in range(1, K_est + 1)])
print("\nElección de M (decodificador):\n", tabla_M.round(2).to_string(index=False))

# Gráficas: varianza explicada y autofunciones (efecto de cada componente sobre la media)
t_g = np.linspace(*RANGO, 731)
P_g = base_bspline(t_g, K_est)
mu_g = P_g @ fpca["mu"]
sd_s = np.sqrt(fpca["lam"])
fig, ax = plt.subplots(1, 3, figsize=(16, 4))
ax[0].bar(range(1, K_est + 1), 100 * fpca["explicada"], color="#2a6f97")
ax[0].plot(range(1, K_est + 1), 100 * var_acum, "o-", color="#b23a48")
ax[0].set(xlabel="Componente m", ylabel="%", title="FPCA estaciones: varianza explicada (entrenamiento)")
for m in range(4):
    ax[1].plot(t_g, P_g @ fpca["psi"][:, m], label=f"ψ{m+1} ({100*fpca['explicada'][m]:.1f} %)")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set(xlabel="Día del año", title="Autofunciones ψ_m(t)")
ax[1].legend(fontsize=8)
m_show = 1
ax[2].plot(t_g, mu_g, "k", lw=2, label="μ(t)")
ax[2].plot(t_g, mu_g + 2 * sd_s[m_show] * (P_g @ fpca["psi"][:, m_show]), "g--", label=f"μ + 2·sd·ψ{m_show+1}")
ax[2].plot(t_g, mu_g - 2 * sd_s[m_show] * (P_g @ fpca["psi"][:, m_show]), "m--", label=f"μ − 2·sd·ψ{m_show+1}")
ax[2].set(xlabel="Día del año", ylabel="sqrt(mm/día)", title=f"Efecto de la componente {m_show+1}")
ax[2].legend(fontsize=8)
plt.tight_layout()
mostrar("fase4_fig01_fpca_estaciones")

# Salidas del bloque 1
verif_W.to_csv(OUT / "verificacion_W_R.csv", index=False)
tabla_fpca.to_csv(OUT / "fpca_varianza_explicada.csv", index=False)
tabla_M.to_csv(OUT / "fpca_eleccion_M.csv", index=False)
print("\nBloque 1 listo. Salidas en", OUT)


###############################################################################
# BLOQUE 2. Red neuronal funcional en PyTorch
###############################################################################
# Arquitectura (igual para las tres cabezas; solo cambia la cabeza):
#   Codificador funcional: z_k = (1/365) ∫ w_k(t) X(t) dt = a_kᵀ (W_sat/365) c_sat,  k = 1..N_NEUR
#                          + altitud y ONI del año (estandarizados con entrenamiento)
#                          + BatchNorm sin parámetros (affine=False): en evaluación es una
#                            transformación afín FIJA, así que no rompe la equivalencia lineal.
#   Cabeza:                z (N_NEUR + 2) -> s (M puntajes FPCA)
#   Decodificador FPCA:    ĉ = μ_c + Ψ s   (μ_c y Ψ del bloque 1, congelados)
#   Pérdida:               (1/365)(ĉ - c)ᵀ W_est (ĉ - c)  +  λ_rug · rugosidad relativa de las w_k
import copy
import math
import torch
import torch.nn as nn

torch.set_default_dtype(torch.float64)        # float64: la verificación de la equivalencia pide <= 1e-4
torch.set_num_threads(max(1, torch.get_num_threads()))

N_NEUR = 8                 # neuronas funcionales (el enunciado pide entre 6 y 8)
M = 6                      # componentes FPCA (bloque 1: primer M con >= 95 % de varianza)
LAMBDA_RUG = 1e-2          # peso de la penalización de rugosidad (igual para todas las cabezas)


def fijar_semilla(s):
    np.random.seed(s)
    torch.manual_seed(s)


# -----------------------------------------------------------------------------
# 2.1 Covariables escalares: altitud de la estación y ONI del año
# -----------------------------------------------------------------------------
# Altitud: tabla del preprocesamiento (catálogo si existe; si no, grilla de altitud corregida).
ranking = pd.read_excel("./salidas_preprocesamiento/ranking_estaciones.xlsx")
altitud = coef_est["STATION_ID"].map(ranking.set_index("STATION_ID")["altitud"]).to_numpy(float)
assert not np.isnan(altitud).any(), "hay estaciones sin altitud"

# ONI (NOAA CPC, oni.ascii.txt, ERSST.v5): media 3-mensual de la anomalía de TSM en Niño 3.4.
# Covariable de la red: promedio de las 12 temporadas del año.
# Fase ENSO del año (para reportar resultados): regla de episodios de la CPC (umbral ±0,5 °C
# durante >= 5 temporadas solapadas seguidas); el año toma la fase si el episodio cubre >= 6 de
# sus 12 temporadas; si no, es Neutro.
oni = pd.read_csv("./Datos/oni_cpc.csv")

def episodios(mask, minimo=5):
    out = np.zeros(len(mask), bool)
    i = 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j < len(mask) and mask[j]:
                j += 1
            if j - i >= minimo:
                out[i:j] = True
            i = j
        else:
            i += 1
    return out

oni["nino"] = episodios((oni["ANOM"] >= 0.5).to_numpy())
oni["nina"] = episodios((oni["ANOM"] <= -0.5).to_numpy())
tabla_enso = oni.groupby("YR").agg(ONI_medio=("ANOM", "mean"), temporadas_nino=("nino", "sum"),
                                   temporadas_nina=("nina", "sum"))
tabla_enso["fase_enso"] = np.where(tabla_enso["temporadas_nino"] >= 6, "Niño",
                                   np.where(tabla_enso["temporadas_nina"] >= 6, "Niña", "Neutro"))
tabla_enso = tabla_enso.loc[1981:2003]
print("\nFase ENSO por año (ONI, NOAA CPC):\n", tabla_enso.round(2).to_string())
tabla_enso.to_csv(OUT / "fase_enso_por_anio.csv")

oni_anio = coef_est["YEAR"].map(tabla_enso["ONI_medio"]).to_numpy(float)
fase_enso = coef_est["YEAR"].map(tabla_enso["fase_enso"]).to_numpy()

# Estandarización con media y desviación SOLO de entrenamiento
extra = np.column_stack([altitud, oni_anio])
media_extra, sd_extra = extra[tr].mean(0), extra[tr].std(0)
extra_std = (extra - media_extra) / sd_extra


# -----------------------------------------------------------------------------
# 2.2 Módulos de la red
# -----------------------------------------------------------------------------
def a_tensor(x):
    return torch.as_tensor(np.asarray(x), dtype=torch.float64)


def init_ventanas(n_neur, K, semilla):
    """a_k iniciales: 'campanas' sobre los índices de la base, repartidas en el año.
    Así cada neurona empieza mirando una ventana distinta del año (ayuda a interpretar las w_k)."""
    g = np.random.default_rng(semilla)
    centros = np.linspace(0, K - 1, n_neur)
    j = np.arange(K)
    A = np.exp(-0.5 * ((j[None, :] - centros[:, None]) / 1.0) ** 2)
    return A + 0.05 * g.standard_normal(A.shape)


class CodificadorFuncional(nn.Module):
    """z_k = a_kᵀ (W/365) c  (integral de w_k(t) X(t) con la matriz de Gram) + covariables + BatchNorm."""

    def __init__(self, W_sat, R_sat, n_neur, n_extra, semilla):
        super().__init__()
        K = W_sat.shape[0]
        self.register_buffer("W", a_tensor(W_sat / 365.0))
        self.register_buffer("W_full", a_tensor(W_sat))
        self.register_buffer("R", a_tensor(R_sat))
        self.A = nn.Parameter(a_tensor(init_ventanas(n_neur, K, semilla)))     # n_neur x K
        self.norm = nn.BatchNorm1d(n_neur + n_extra, affine=False)
        # Máxima rugosidad relativa posible en esta base (autovalor máx. de W^{-1} R): sirve para que
        # la penalización quede entre 0 y 1 por neurona y no dependa de las unidades de t.
        self.register_buffer("rug_max", a_tensor(np.linalg.eigvals(np.linalg.solve(W_sat, R_sat)).real.max()))

    def integrales(self, c_sat):
        return c_sat @ (self.A @ self.W).T                                      # (B, n_neur)

    def forward(self, c_sat, extra):
        return self.norm(torch.cat([self.integrales(c_sat), extra], dim=1))

    def rugosidad(self):
        """Rugosidad relativa de cada w_k: ∫(w_k'')² / ∫ w_k², dividida por la máxima posible.
        Es invariante a la escala de a_k (la BatchNorm hace que la escala de a_k no importe,
        así que una penalización absoluta aᵀRa se podría 'hacer trampa' encogiendo a_k)."""
        num = torch.einsum("ki,ij,kj->k", self.A, self.R, self.A)
        den = torch.einsum("ki,ij,kj->k", self.A, self.W_full, self.A)
        return (num / den / self.rug_max).sum()


class DecodificadorFPCA(nn.Module):
    """ĉ = μ_c + Ψ s : de M puntajes a los coeficientes de la curva de la estación. Congelado."""

    def __init__(self, fpca, M):
        super().__init__()
        self.register_buffer("mu", a_tensor(fpca["mu"]))                        # K_est
        self.register_buffer("Psi", a_tensor(fpca["psi"][:, :M]))               # K_est x M

    def forward(self, s):
        return self.mu + s @ self.Psi.T


def crear_cabeza(tipo, d_in, M, **kw):
    if tipo == "lineal":
        return nn.Linear(d_in, M)
    if tipo == "kan":
        return CabezaKAN(d_in, M, **kw)          # bloque 4
    if tipo == "mlp":
        return CabezaMLP(d_in, M, **kw)          # bloque 5
    raise ValueError(f"cabeza desconocida: {tipo}")


class RedFuncional(nn.Module):
    def __init__(self, tipo_cabeza, semilla, **kw_cabeza):
        super().__init__()
        fijar_semilla(semilla)
        self.tipo = tipo_cabeza
        if tipo_cabeza in ("kan", "mlp"):
            kw_cabeza = {**kw_cabeza, "semilla": semilla}
        self.cod = CodificadorFuncional(W_sat, R_sat, N_NEUR, extra.shape[1], semilla)
        self.cabeza = crear_cabeza(tipo_cabeza, N_NEUR + extra.shape[1], M, **kw_cabeza)
        self.dec = DecodificadorFPCA(fpca, M)

    def forward(self, c_sat, extra):
        return self.dec(self.cabeza(self.cod(c_sat, extra)))


W_est_t = a_tensor(W_est)

def error_L2(c_hat, c):
    """∫ (ŷ - y)² dt por curva, con la matriz de Gram de las estaciones."""
    e = c_hat - c
    return torch.einsum("bi,ij,bj->b", e, W_est_t, e)


# -----------------------------------------------------------------------------
# 2.3 Entrenamiento (mismo presupuesto para todas las cabezas)
# -----------------------------------------------------------------------------
# Basado en train_model del curso (S14): AdamW, recorte de gradiente y parada temprana con
# validación. La pérdida no es el MSE de PyTorch sino la norma L² funcional + rugosidad.
PRESUPUESTO = dict(epocas=400, paciencia=30, lote=128, lr=3e-3, weight_decay=1e-4, clip=1.0)

Xs = a_tensor(C_sat); Xe = a_tensor(extra_std); Yc = a_tensor(C_est)
idx_tr, idx_va, idx_te = np.where(tr)[0], np.where(va)[0], np.where(te)[0]


def recalcular_bn(modelo):
    """Estadísticas de la BatchNorm EXACTAS sobre todo el entrenamiento con las a_k actuales.
    Las estadísticas 'móviles' por defecto se quedan atrás cuando las a_k cambian y hacen que la
    pérdida de validación salte de una época a otra. Solo usa entrenamiento (no hay fuga)."""
    bn = modelo.cod.norm
    bn.reset_running_stats()
    bn.momentum = None                               # promedio acumulado: con un solo paso = exacto
    modelo.train()
    with torch.no_grad():
        modelo.cod(Xs[idx_tr], Xe[idx_tr])
    modelo.eval()


def entrenar(modelo, semilla, epocas, paciencia, lote, lr, weight_decay, clip, lam_rug=LAMBDA_RUG,
             verbose=False):
    fijar_semilla(semilla)
    opt = torch.optim.AdamW([p for p in modelo.parameters() if p.requires_grad], lr=lr,
                            weight_decay=weight_decay)
    gen = torch.Generator().manual_seed(semilla)
    hist = {"train": [], "val": []}
    mejor, mejor_estado, malas = np.inf, None, 0
    for ep in range(1, epocas + 1):
        modelo.train()
        perm = idx_tr[torch.randperm(len(idx_tr), generator=gen).numpy()]
        acum = 0.0
        for i in range(0, len(perm), lote):
            b = perm[i:i + lote]
            if len(b) < 2:                           # BatchNorm necesita al menos 2 curvas
                continue
            c_hat = modelo(Xs[b], Xe[b])
            perdida = error_L2(c_hat, Yc[b]).mean() / 365 + lam_rug * modelo.cod.rugosidad()
            opt.zero_grad(set_to_none=True)
            perdida.backward()
            nn.utils.clip_grad_norm_(modelo.parameters(), clip)
            opt.step()
            acum += perdida.item() * len(b)
        recalcular_bn(modelo)
        with torch.no_grad():
            val = (error_L2(modelo(Xs[idx_va], Xe[idx_va]), Yc[idx_va]).mean() / 365).item()
        hist["train"].append(acum / len(idx_tr)); hist["val"].append(val)
        if val < mejor - 1e-7:
            mejor, mejor_estado, malas = val, copy.deepcopy(modelo.state_dict()), 0
        else:
            malas += 1
        if verbose and (ep == 1 or ep % 25 == 0):
            print(f"  época {ep:3d} | train {hist['train'][-1]:.4f} | val {val:.4f}")
        if malas >= paciencia:
            if verbose:
                print(f"  parada temprana en la época {ep} (mejor val = {mejor:.4f})")
            break
    modelo.load_state_dict(mejor_estado)
    modelo.eval()
    return hist


# -----------------------------------------------------------------------------
# 2.4 Primera prueba: cabeza lineal, semilla 0
# -----------------------------------------------------------------------------
# Referencia obligatoria: la curva satelital SIN corregir comparada con la de la estación.
# Las dos bases son distintas (K = 12 vs K = 8), así que el error usa la Gram cruzada
# W_cruz = ∫ φ_sat φ_estᵀ dt, integrada exacto con Gauss-Legendre en la unión de los nudos.
def gram_cruzada(K1, K2, q=4):
    u = np.union1d(np.unique(nudos_bspline(K1)), np.unique(nudos_bspline(K2)))
    x, w = np.polynomial.legendre.leggauss(q)
    t_q = np.concatenate([(b - a) / 2 * x + (a + b) / 2 for a, b in zip(u[:-1], u[1:])])
    w_q = np.concatenate([(b - a) / 2 * w for a, b in zip(u[:-1], u[1:])])
    return base_bspline(t_q, K1).T @ (w_q[:, None] * base_bspline(t_q, K2))

W_cruz = gram_cruzada(K_sat, K_est)                                   # K_sat x K_est
W_sat_exacta = gram_cruzada(K_sat, K_sat)
assert np.abs(W_sat_exacta - W_sat).max() < 1e-9                      # control de la Gram cruzada

def L2_satelite(idx):
    cs, ce = C_sat[idx], C_est[idx]
    return (np.einsum("bi,ij,bj->b", cs, W_sat, cs) - 2 * np.einsum("bi,ij,bj->b", cs, W_cruz, ce)
            + np.einsum("bi,ij,bj->b", ce, W_est, ce))

def L2_relativo(err, idx):
    return 100 * np.sqrt(err.sum() / norma2(C_est[idx]).sum())

print("\nEntrenando la red funcional con cabeza lineal (semilla 0)...")
red_lin = RedFuncional("lineal", semilla=0)
hist_lin = entrenar(red_lin, semilla=0, verbose=True, **PRESUPUESTO)
with torch.no_grad():
    err_va_lin = error_L2(red_lin(Xs[idx_va], Xe[idx_va]), Yc[idx_va]).numpy()

print("\nError L² relativo en VALIDACIÓN (1996-2000), escala sqrt(mm/día):")
print(f"  Satélite sin corregir         : {L2_relativo(L2_satelite(idx_va), idx_va):6.2f} %")
print(f"  Red funcional, cabeza lineal  : {L2_relativo(err_va_lin, idx_va):6.2f} %")
print(f"  Techo del decodificador (M={M}) : {tabla_M.loc[M - 1, 'L2_rel_validacion_pct']:6.2f} %")

# Control: el mejor modelo LINEAL posible con las mismas entradas, por mínimos cuadrados cerrados
# (puntajes FPCA ~ 1 + c_sat + altitud + ONI). Si la red lineal queda cerca, el entrenamiento
# funcionó y la diferencia con el techo es falta de información en el satélite, no de optimización.
X_ols = np.column_stack([np.ones(len(C_sat)), C_sat, extra_std])
S_obj = (C_est - fpca["mu"]) @ W_est @ fpca["psi"][:, :M]
B_ols = np.linalg.lstsq(X_ols[idx_tr], S_obj[idx_tr], rcond=None)[0]
C_ols = fpca["mu"] + X_ols @ B_ols @ fpca["psi"][:, :M].T
E_ols = C_ols - C_est
err_ols = np.einsum("bi,ij,bj->b", E_ols, W_est, E_ols)
print(f"  Mínimos cuadrados (control)   : {L2_relativo(err_ols[idx_va], idx_va):6.2f} %")

# Curvas de aprendizaje y funciones de pesos w_k(t)
P_sat_g = base_bspline(t_g, K_sat)
fig, ax = plt.subplots(1, 2, figsize=(14, 4))
ax[0].plot(hist_lin["train"], label="entrenamiento (con rugosidad)")
ax[0].plot(hist_lin["val"], label="validación")
ax[0].set(xlabel="Época", ylabel="L² / 365", title="Curvas de aprendizaje - cabeza lineal")
ax[0].legend()
Wk = P_sat_g @ red_lin.cod.A.detach().numpy().T                        # 731 x N_NEUR
for k in range(N_NEUR):
    ax[1].plot(t_g, Wk[:, k] / np.abs(Wk[:, k]).max(), label=f"w{k+1}")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set(xlabel="Día del año", ylabel="w_k(t) (escalada a máx. 1)",
          title="Funciones de pesos del codificador (cabeza lineal, semilla 0)")
ax[1].legend(fontsize=7, ncol=2)
plt.tight_layout()
mostrar("fase4_fig02_aprendizaje_y_wk_lineal")

print("\nBloque 2 listo.")


###############################################################################
# BLOQUE 3. Cabeza lineal = regresión funcional función a función (demostración y verificación)
###############################################################################
# DEMOSTRACIÓN (va escrita en el informe):
# En evaluación, para una curva satelital X(u) = φ_sat(u)ᵀ c y covariables e = (altitud, ONI):
#   (1) Codificador:  z̃_k = (1/365) ∫ w_k(u) X(u) du,  w_k(u) = φ_sat(u)ᵀ a_k      (k = 1..8)
#                     z̃_{8+j} = e_j                                                  (j = 1, 2)
#   (2) BatchNorm en evaluación (fija):  z_k = d_k (z̃_k − m_k),  d_k = 1/sqrt(v_k + ε)
#   (3) Cabeza lineal:  s_m = b_m + Σ_k B_mk z_k                                     (m = 1..M)
#   (4) Decodificador:  ŷ(t) = μ(t) + Σ_m s_m ψ_m(t)
# Sustituyendo (1)-(3) en (4) e intercambiando la suma finita con la integral:
#   ŷ(t) = β0(t) + ∫ β(u,t) X(u) du + Σ_j γ_j(t) e_j
# con
#   β(u,t)  = (1/365) Σ_{k<=8} Σ_m  w_k(u) d_k B_mk ψ_m(t)
#   γ_j(t)  = Σ_m d_{8+j} B_{m,8+j} ψ_m(t)
#   β0(t)   = μ(t) + Σ_m [ b_m − Σ_k B_mk d_k m_k ] ψ_m(t)
# Es el modelo de regresión funcional función a función Y(t) = β0(t) + ∫ β(u,t) X(u) du (+ efectos
# de covariables escalares), con β(u,t) de rango <= min(8, M) = 6: "regresión funcional de rango bajo".
# En forma matricial:  β(u,t) = φ_sat(u)ᵀ Bβ φ_est(t),  Bβ = (1/365) Aᵀ diag(d) Bᵀ Ψᵀ   (12 x 8).
# Con la BatchNorm el factor d_k se absorbe en la escala de w_k; el enunciado lo escribe sin él.

def superficie_beta(modelo):
    """Coeficientes de β(u,t), β0(t) y γ_j(t) a partir de los parámetros de la red entrenada."""
    with torch.no_grad():
        A = modelo.cod.A.numpy()                                       # 8 x K_sat
        bn = modelo.cod.norm
        d = 1 / np.sqrt(bn.running_var.numpy() + bn.eps)
        m_bn = bn.running_mean.numpy()
        Bh = modelo.cabeza.weight.numpy()                               # M x 10
        bh = modelo.cabeza.bias.numpy()                                 # M
        Psi = modelo.dec.Psi.numpy()                                    # K_est x M
        mu = modelo.dec.mu.numpy()
    n = A.shape[0]
    B_beta = (A.T * (d[:n] / 365)) @ Bh[:, :n].T @ Psi.T                # K_sat x K_est
    c_gamma = (Psi @ (Bh[:, n:] * d[n:])).T                             # 2 x K_est
    c_beta0 = mu + Psi @ (bh - Bh @ (d * m_bn))                         # K_est
    return B_beta, c_beta0, c_gamma

B_beta, c_beta0, c_gamma = superficie_beta(red_lin)

# VERIFICACIÓN NUMÉRICA: la salida de la red vs la regresión funcional evaluada por separado.
# La integral ∫ β(u,t) X(u) du se calcula con cuadratura de Gauss en los nudos del satélite
# (β(·,t) y X son cúbicos por tramos en esos nudos: el producto es de grado 6 y la cuadratura
# con 4 puntos es exacta). Se evalúa en una malla de 731 días para TODAS las curvas.
t_q, w_q = cuadratura_gauss(K_sat)
P_q = base_bspline(t_q, K_sat)                                          # n_q x K_sat
beta_ut = P_q @ B_beta @ P_g.T                                          # β(u_q, t_g): n_q x 731
X_q = C_sat @ P_q.T                                                     # X_i(u_q): N x n_q
y_regfun = (c_beta0 @ P_g.T)[None, :] + (X_q * w_q) @ beta_ut + extra_std @ (c_gamma @ P_g.T)
with torch.no_grad():
    y_red = red_lin(Xs, Xe).numpy() @ P_g.T                             # salida de la red en la malla
dif_equiv = np.abs(y_red - y_regfun).max()
print(f"\nEquivalencia cabeza lineal = regresión funcional | máx |ŷ_red − ŷ_regfun| = {dif_equiv:.2e}"
      f"  (criterio <= 1e-4: {'CUMPLE' if dif_equiv <= 1e-4 else 'NO CUMPLE'})")
pd.DataFrame({"verificacion": ["max_dif_red_vs_regresion_funcional"], "valor": [dif_equiv],
              "curvas": [len(C_sat)], "puntos_t": [len(t_g)]}).to_csv(OUT / "verificacion_equivalencia_lineal.csv", index=False)

# Superficie β(u,t): qué días del satélite (u) influyen en cada día de la estación (t)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
beta_gg = P_sat_g @ B_beta @ P_g.T
lim = np.abs(beta_gg).max()
im = ax[0].contourf(t_g, t_g, beta_gg.T, levels=21, cmap="RdBu_r", vmin=-lim, vmax=lim)
ax[0].plot([0, 365], [0, 365], "k--", lw=0.6)
ax[0].set(xlabel="u: día del satélite", ylabel="t: día de la estación", title="β(u,t) de la cabeza lineal")
plt.colorbar(im, ax=ax[0])
ax[1].plot(t_g, c_beta0 @ P_g.T, label="β0(t)")
ax[1].plot(t_g, c_gamma[0] @ P_g.T, label="γ_altitud(t) (por 1 sd)")
ax[1].plot(t_g, c_gamma[1] @ P_g.T, label="γ_ONI(t) (por 1 sd)")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set(xlabel="Día del año", ylabel="sqrt(mm/día)", title="Intercepto funcional y efectos escalares")
ax[1].legend(fontsize=8)
plt.tight_layout()
mostrar("fase4_fig03_beta_cabeza_lineal")

print("\nBloque 3 listo.")


###############################################################################
# BLOQUE 4. KAN con B-splines en las aristas, verificada con una función sintética
###############################################################################
# Basada en la KANLayer del curso (S14), con un cambio: el curso usa bases gaussianas (FastKAN) y
# el enunciado pide B-splines:  φ(x) = w_b·SiLU(x) + w_s·Σ_i c_i B_i(x).
# Las B_i son B-splines cúbicas sobre una rejilla UNIFORME en [x_min, x_max] con G intervalos,
# extendida k nudos a cada lado (como pykan) para que en todo el rango sumen 1. Hay G + k bases.
# Fuera del rango las B-splines valen 0 y solo queda la parte SiLU (la función no se "corta").

class CapaKAN(nn.Module):
    def __init__(self, n_in, n_out, G=5, k=3, rango=(-3.0, 3.0), layernorm=False, semilla=0):
        super().__init__()
        gen = torch.Generator().manual_seed(semilla)
        self.n_in, self.n_out, self.G, self.k = n_in, n_out, G, k
        h = (rango[1] - rango[0]) / G
        self.register_buffer("nudos", rango[0] + h * torch.arange(-k, G + k + 1, dtype=torch.float64))
        self.norm = nn.LayerNorm(n_in) if layernorm else nn.Identity()
        self.w_b = nn.Parameter(torch.empty(n_out, n_in))                         # peso de SiLU
        nn.init.kaiming_uniform_(self.w_b, a=math.sqrt(5), generator=gen)
        self.w_s = nn.Parameter(torch.ones(n_out, n_in))                          # escala de la spline
        self.c = nn.Parameter(0.1 * torch.randn(n_out, n_in, G + k, generator=gen))  # coeficientes c_i

    def bases(self, x):
        """B-splines de grado k en x: (B, n_in) -> (B, n_in, G + k). Recursión de Cox-de Boor."""
        u = self.nudos
        x = x.unsqueeze(-1)
        B = ((x >= u[:-1]) & (x < u[1:])).to(x.dtype)
        for p in range(1, self.k + 1):
            B = ((x - u[:-(p + 1)]) / (u[p:-1] - u[:-(p + 1)]) * B[..., :-1]
                 + (u[p + 1:] - x) / (u[p + 1:] - u[1:-p]) * B[..., 1:])
        return B

    def forward(self, x):
        x = self.norm(x)
        base = nn.functional.silu(x) @ self.w_b.T                                   # (B, out)
        spline = torch.einsum("big,oig->bo", self.bases(x), self.c * self.w_s[..., None])
        return base + spline

    @torch.no_grad()
    def funcion_arista(self, o, i, x):
        """φ_{o,i}(x) sobre una malla 1D (sin la LayerNorm), para graficar e interpretar."""
        Bx = self.bases(x[:, None])[:, 0, :]
        return self.w_b[o, i] * nn.functional.silu(x) + Bx @ (self.c[o, i] * self.w_s[o, i])


class CabezaKAN(nn.Module):
    """Cabeza KAN [d_in -> H -> M]. La 1.ª capa actúa DIRECTO sobre z (ya estandarizado por la
    BatchNorm del codificador), así sus activaciones se leen sobre cada neurona funcional / altitud /
    ONI. La 2.ª capa lleva LayerNorm a la entrada, como la KANLayer del curso."""

    def __init__(self, d_in, M, H=4, G=5, k=3, semilla=0):
        super().__init__()
        self.capa1 = CapaKAN(d_in, H, G, k, rango=(-3.0, 3.0), layernorm=False, semilla=semilla)
        self.capa2 = CapaKAN(H, M, G, k, rango=(-3.0, 3.0), layernorm=True, semilla=semilla + 1)

    def forward(self, z):
        return self.capa2(self.capa1(z))


# --- Verificación con una función sintética (antes de usar datos reales) -----------------------
# y = sin(π x1) + x2² + 0·x3, sin ruido, x ~ U(-1, 1). Es aditiva, así que una sola capa KAN
# (y = Σ_i φ_i(x_i)) la puede representar: cada arista debe aprender sin(πx), x² y una constante.
# Criterio de la rúbrica: RMSE <= 1e-2 en puntos de prueba que no se usaron para entrenar.
def prueba_sintetica_kan(semilla=0, pasos=3000):
    fijar_semilla(semilla)
    g = torch.Generator().manual_seed(semilla)
    x = torch.rand(3000, 3, generator=g) * 2 - 1
    f = lambda x: torch.sin(math.pi * x[:, 0]) + x[:, 1] ** 2
    y = f(x)
    x_tr, y_tr, x_te, y_te = x[:2000], y[:2000], x[2000:], y[2000:]
    kan = CapaKAN(3, 1, G=10, k=3, rango=(-1.0, 1.0), semilla=semilla)
    opt = torch.optim.Adam(kan.parameters(), lr=1e-2)
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=1000, gamma=0.3)
    for _ in range(pasos):
        perdida = ((kan(x_tr).squeeze(-1) - y_tr) ** 2).mean()
        opt.zero_grad(); perdida.backward(); opt.step(); sched.step()
    with torch.no_grad():
        rmse = torch.sqrt(((kan(x_te).squeeze(-1) - y_te) ** 2).mean()).item()
    return kan, rmse

kan_sint, rmse_sint = prueba_sintetica_kan()
print(f"\nKAN sintética: RMSE en prueba = {rmse_sint:.2e}  (criterio <= 1e-2: "
      f"{'CUMPLE' if rmse_sint <= 1e-2 else 'NO CUMPLE'})")
pd.DataFrame({"funcion": ["sin(pi x1) + x2^2 + 0 x3"], "G": [10], "k": [3],
              "rmse_prueba": [rmse_sint]}).to_csv(OUT / "verificacion_kan_sintetica.csv", index=False)

malla = torch.linspace(-1, 1, 200)
verdad = [torch.sin(math.pi * malla), malla ** 2, torch.zeros_like(malla)]
fig, ax = plt.subplots(1, 3, figsize=(13, 3.2))
for j in range(3):
    centrar = lambda v: v - v.mean()          # las funciones aditivas se definen salvo una constante
    ax[j].plot(malla, centrar(verdad[j]), "k--", lw=2, label="verdadera")
    ax[j].plot(malla, centrar(kan_sint.funcion_arista(0, j, malla)), label="KAN φ_j")
    ax[j].set_title(f"x{j+1}"); ax[j].legend(fontsize=8)
plt.suptitle(f"Prueba sintética de la KAN con B-splines (RMSE = {rmse_sint:.1e})")
plt.tight_layout()
mostrar("fase4_fig04_kan_sintetica")
print("\nBloque 4 listo.")


###############################################################################
# BLOQUE 5. Ablación de la cabeza: lineal vs KAN vs MLP
###############################################################################
# Diseño (rúbrica 4.3):
#   - Mismo codificador, mismo decodificador FPCA, misma pérdida: solo cambia la cabeza.
#   - KAN y MLP con el mismo número de parámetros (±10 %), contados con sum(p.numel()).
#   - Mismo presupuesto: épocas máximas, parada temprana con validación (1996-2000) y el MISMO
#     número de pruebas de hiperparámetros (4 por cabeza: lr x λ_rug), elegidas con la semilla 0.
#   - Con la configuración elegida, 3 semillas (0, 1, 2): media y desviación estándar.
#   - La prueba (2001-2003) NO se usa para elegir nada; solo se reporta al final.
#   - Una cabeza supera a otra solo si el IC bootstrap del 95 % de la diferencia excluye el 0
#     (remuestreo de ESTACIONES, 2000 réplicas).
from scipy.stats import wilcoxon


class CabezaMLP(nn.Module):
    """Red densa [d_in -> h -> M] con GELU (como make_head del curso), sin dropout (igual que la KAN)."""

    def __init__(self, d_in, M, h=32, semilla=0):
        super().__init__()
        self.red = nn.Sequential(nn.Linear(d_in, h), nn.GELU(), nn.Linear(h, M))

    def forward(self, z):
        return self.red(z)


def n_parametros(modulo):
    return sum(p.numel() for p in modulo.parameters())


# --- Igualar parámetros: se fija la KAN y se elige el ancho h del MLP más cercano ---------------
D_IN = N_NEUR + extra.shape[1]
KW_KAN = dict(H=4, G=5, k=3)
p_kan = n_parametros(CabezaKAN(D_IN, M, **KW_KAN))
h_mlp = min(range(4, 200), key=lambda h: abs(n_parametros(CabezaMLP(D_IN, M, h)) - p_kan))
KW_MLP = dict(h=h_mlp)
p_mlp = n_parametros(CabezaMLP(D_IN, M, **KW_MLP))
p_lin = n_parametros(nn.Linear(D_IN, M))
p_cod = n_parametros(CodificadorFuncional(W_sat, R_sat, N_NEUR, extra.shape[1], 0))
tabla_param = pd.DataFrame({
    "cabeza": ["lineal", "kan", "mlp"],
    "arquitectura": [f"Linear({D_IN}->{M})", f"KAN[{D_IN}->{KW_KAN['H']}->{M}], G={KW_KAN['G']}, k={KW_KAN['k']}",
                     f"MLP[{D_IN}->{h_mlp}->{M}], GELU"],
    "parametros_cabeza": [p_lin, p_kan, p_mlp],
    "parametros_codificador": [p_cod] * 3,
    "parametros_decodificador_entrenables": [0] * 3})
tabla_param["dif_vs_kan_pct"] = 100 * (tabla_param["parametros_cabeza"] - p_kan) / p_kan
print("\nParámetros por cabeza (sum(p.numel())):\n", tabla_param.round(1).to_string(index=False))
assert abs(p_mlp - p_kan) / p_kan <= 0.10, "KAN y MLP no quedaron dentro del ±10 %"
tabla_param.to_csv(OUT / "ablacion_parametros.csv", index=False)

KW_CABEZA = {"lineal": {}, "kan": KW_KAN, "mlp": KW_MLP}


def construir(tipo, semilla):
    return RedFuncional(tipo, semilla, **KW_CABEZA[tipo])


# --- Búsqueda de hiperparámetros: mismo número de pruebas por cabeza (semilla 0, validación) ----
REJILLA = [dict(lr=lr, lam_rug=lam) for lr in (1e-3, 3e-3) for lam in (1e-2, 1e-1)]
RUTA_BUSQ = OUT / "ablacion_busqueda_hiperparametros.csv"
busqueda = []
for tipo in ("lineal", "kan", "mlp"):
    if not REENTRENAR and RUTA_BUSQ.exists():
        break
    for cfg in REJILLA:
        modelo = construir(tipo, 0)
        pres = dict(PRESUPUESTO); pres.update(lr=cfg["lr"])
        h = entrenar(modelo, semilla=0, lam_rug=cfg["lam_rug"], **pres)
        busqueda.append({"cabeza": tipo, **cfg, "val_L2": min(h["val"]), "epocas": len(h["val"])})
        print(f"  {tipo:6s} lr={cfg['lr']:.0e} λ_rug={cfg['lam_rug']:.0e} -> val {min(h['val']):.4f}"
              f" ({len(h['val'])} épocas)")
if busqueda:
    busqueda = pd.DataFrame(busqueda)
    busqueda.to_csv(RUTA_BUSQ, index=False)
else:
    busqueda = pd.read_csv(RUTA_BUSQ)
    print("Búsqueda de hiperparámetros cargada de", RUTA_BUSQ)
mejor_cfg = {t: busqueda[busqueda.cabeza == t].sort_values("val_L2").iloc[0][["lr", "lam_rug"]].to_dict()
             for t in ("lineal", "kan", "mlp")}
print("\nConfiguración elegida (solo con validación):", mejor_cfg)


# --- 3 semillas por cabeza con la configuración elegida ------------------------------------------
DIR_MOD = OUT / "modelos"; DIR_MOD.mkdir(exist_ok=True)
SEMILLAS = (0, 1, 2)
modelos, err_curva = {}, {}
for tipo in ("lineal", "kan", "mlp"):
    for sem in SEMILLAS:
        modelo = construir(tipo, sem)
        ruta_mod = DIR_MOD / f"red_{tipo}_semilla{sem}.pt"
        if not REENTRENAR and ruta_mod.exists():
            modelo.load_state_dict(torch.load(ruta_mod))
            modelo.eval()
        else:
            pres = dict(PRESUPUESTO); pres.update(lr=mejor_cfg[tipo]["lr"])
            entrenar(modelo, semilla=sem, lam_rug=mejor_cfg[tipo]["lam_rug"], **pres)
            torch.save(modelo.state_dict(), ruta_mod)
        modelos[(tipo, sem)] = modelo
        with torch.no_grad():
            C_hat = modelo(Xs, Xe).numpy()
        err_curva[(tipo, sem)] = np.einsum("bi,ij,bj->b", C_hat - C_est, W_est, C_hat - C_est)
        # Predicciones (coeficientes de la base de estaciones, escala sqrt) para la evaluación 4.4
        pred = coef_est[["STATION_ID", "YEAR", "particion"]].copy()
        pred[[f"c{j+1}" for j in range(K_est)]] = C_hat
        pred.to_csv(OUT / f"pred_coef_{tipo}_semilla{sem}.csv", index=False)
        print(f"  {tipo:6s} semilla {sem}: val {L2_relativo(err_curva[(tipo, sem)][idx_va], idx_va):.2f} %"
              f" | prueba {L2_relativo(err_curva[(tipo, sem)][idx_te], idx_te):.2f} %")

err_sat = L2_satelite(np.arange(len(C_sat)))
filas = [{"modelo": "satélite sin corregir", "val_media": L2_relativo(err_sat[idx_va], idx_va), "val_sd": 0.0,
          "prueba_media": L2_relativo(err_sat[idx_te], idx_te), "prueba_sd": 0.0}]
for tipo in ("lineal", "kan", "mlp"):
    v = [L2_relativo(err_curva[(tipo, s)][idx_va], idx_va) for s in SEMILLAS]
    p = [L2_relativo(err_curva[(tipo, s)][idx_te], idx_te) for s in SEMILLAS]
    filas.append({"modelo": f"red funcional - cabeza {tipo}", "val_media": np.mean(v), "val_sd": np.std(v, ddof=1),
                  "prueba_media": np.mean(p), "prueba_sd": np.std(p, ddof=1)})
tabla_abl = pd.DataFrame(filas)
print("\nAblación - error L² relativo (%), media y desviación de 3 semillas:\n",
      tabla_abl.round(2).to_string(index=False))
tabla_abl.to_csv(OUT / "ablacion_resultados.csv", index=False)


# --- Diferencias entre cabezas en PRUEBA: IC bootstrap por estación, Wilcoxon y lectura por estación
# Error de cada cabeza por curva = promedio de las 3 semillas. Se remuestrean estaciones (no curvas)
# porque las curvas de una misma estación no son independientes.
# Además del agregado (lo que pide la regla del IC) se reporta la estación típica: mediana de la
# diferencia por estación y % de estaciones donde gana la primera cabeza. Y se separa por región
# (zona hidrográfica del código IDEAM: 26 = Cauca, 53/54 = Pacífico, igual que la fase 3) y sin las
# estaciones con coordenadas no confiables (su píxel, y por tanto su curva satelital, puede ser otro).
NO_CONFIABLES = [5300800101, 5300900102, 5300600101, 5300000203, 2621500105, 5300000202,
                 5301300101, 5301000103, 5400000201, 2621500104, 5300900101]
ids_all = coef_est["STATION_ID"].to_numpy()
region_curva = np.where(pd.Series(ids_all).astype(str).str[:2] == "26", "Cauca", "Pacífico")
SUBCONJUNTOS = {"Todas": np.ones(len(ids_all), bool),
                "Sin no confiables": ~np.isin(ids_all, NO_CONFIABLES),
                "Cauca": region_curva == "Cauca",
                "Pacífico": region_curva == "Pacífico",
                "Pacífico sin no confiables": (region_curva == "Pacífico") & ~np.isin(ids_all, NO_CONFIABLES)}

err_prom = {t: np.mean([err_curva[(t, s)] for s in SEMILLAS], axis=0) for t in ("lineal", "kan", "mlp")}
err_prom["satelite"] = L2_satelite(np.arange(len(C_sat)))
norma_all = norma2(C_est)
rng = np.random.default_rng(SEED)
n_boot = 2000
comparaciones = []
for nombre_sub, mascara in SUBCONJUNTOS.items():
    idx = idx_te[mascara[idx_te]]
    est, pos = np.unique(ids_all[idx], return_inverse=True)
    nor = np.bincount(pos, weights=norma_all[idx])
    E = {t: np.bincount(pos, weights=e[idx]) for t, e in err_prom.items()}
    conteos = np.stack([np.bincount(rng.integers(0, len(est), len(est)), minlength=len(est))
                        for _ in range(n_boot)])                             # n_boot x n_estaciones
    l2b = lambda Et: 100 * np.sqrt((conteos @ Et) / (conteos @ nor))
    for a, b in [("kan", "lineal"), ("mlp", "lineal"), ("kan", "mlp"), ("lineal", "satelite")]:
        dif_obs = 100 * (np.sqrt(E[a].sum() / nor.sum()) - np.sqrt(E[b].sum() / nor.sum()))
        li, ls = np.percentile(l2b(E[a]) - l2b(E[b]), [2.5, 97.5])
        dif_est = 100 * (np.sqrt(E[a] / nor) - np.sqrt(E[b] / nor))         # L² relativo por estación
        comparaciones.append({
            "subconjunto": nombre_sub, "comparacion": f"{a} - {b}", "n_estaciones": len(est),
            "dif_L2rel_pp": dif_obs, "IC95_inf": li, "IC95_sup": ls, "IC_excluye_0": (li > 0) or (ls < 0),
            "mediana_dif_por_estacion_pp": np.median(dif_est),
            "pct_estaciones_gana_primera": 100 * np.mean(dif_est < 0),
            "wilcoxon_p": wilcoxon(dif_est).pvalue if len(est) > 5 else np.nan})
tabla_comp = pd.DataFrame(comparaciones)
with pd.option_context("display.width", 220):
    print("\nDiferencias en PRUEBA (puntos de L² relativo; negativo = la primera cabeza es mejor):\n",
          tabla_comp.round(3).to_string(index=False))
tabla_comp.to_csv(OUT / "ablacion_comparaciones_bootstrap.csv", index=False)


# --- Primer vistazo a las activaciones de la 1.ª capa de la KAN (semilla 0) ----------------------
kan0 = modelos[("kan", 0)]
with torch.no_grad():
    Z = kan0.cod(Xs[idx_tr], Xe[idx_tr]).numpy()                     # z de entrenamiento (rango real)
nombres_z = [f"z{k+1} (neurona {k+1})" for k in range(N_NEUR)] + ["altitud", "ONI"]
fig, ax = plt.subplots(2, 5, figsize=(18, 6), sharey=False)
for i in range(D_IN):
    lo, hi = np.percentile(Z[:, i], [1, 99])
    xm = torch.linspace(lo, hi, 200)
    a_ = ax.flat[i]
    for o in range(KW_KAN["H"]):
        a_.plot(xm, kan0.cabeza.capa1.funcion_arista(o, i, xm), label=f"→ oculta {o+1}")
    a_.set_title(nombres_z[i], fontsize=9)
ax.flat[0].legend(fontsize=7)
plt.suptitle("Activaciones de la 1.ª capa de la KAN en el rango de los datos (percentiles 1-99)")
plt.tight_layout()
mostrar("fase4_fig05_activaciones_kan")
print("\nBloque 5 listo.")


###############################################################################
# BLOQUE 6. Reconstrucción en mm: no negativa y sin sesgo de retransformación
###############################################################################
# La red reconstruye ŷ(t) en sqrt(mm/día). Para volver a mm:
#   1) Recorte:  ŷ⁺(t) = max(ŷ(t), 0)   -> se reporta el % de valores recortados.
#   2) Retransformación con corrección de varianza (smearing de Duan, 1983, por niveles):
#        E[P | x] = E[(ŷ⁺ + r)²]  con r = sqrt(P_obs) − ŷ⁺ el residuo diario.
#      Elevar solo ŷ⁺² pierde la varianza de r: con lluvia diaria intermitente eso subestima el
#      acumulado anual en ~65 %. Se estima, SOLO con entrenamiento y para cada modelo,
#      la media r̄_b y la media de r² (m2_b) dentro de 10 niveles (deciles) de ŷ⁺, y
#        P̂(t) = ŷ⁺² + 2 ŷ⁺ r̄_b + m2_b  = promedio de (ŷ⁺ + r)² en el nivel b   (>= 0 siempre).
# Fuente de la observación diaria: estaciones_diario_limpio.csv (el 29 de febrero se suma al 28,
# como en la fase 1).

diario = pd.read_csv("./salidas_preprocesamiento/estaciones_diario_limpio.csv",
                     usecols=["STATION_ID", "YEAR", "MONTH", "DAY", "VALUE"])
cum_mes = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
diario["doy"] = np.where((diario["MONTH"] == 2) & (diario["DAY"] == 29), 58,
                         cum_mes[diario["MONTH"] - 1] + diario["DAY"] - 1)
fila = pd.Series(np.arange(len(coef_est)),
                 index=pd.MultiIndex.from_frame(coef_est[["STATION_ID", "YEAR"]]))
diario["fila"] = fila.reindex(pd.MultiIndex.from_frame(diario[["STATION_ID", "YEAR"]])).to_numpy()
assert not diario["fila"].isna().any()
Y_mm = np.zeros((len(coef_est), 365))
np.add.at(Y_mm, (diario["fila"].to_numpy(int), diario["doy"].to_numpy(int)), diario["VALUE"].to_numpy(float))
del diario

t_dia = np.arange(365) + 0.5                       # centro de cada día (igual que la fase 1)
P_dia = base_bspline(t_dia, K_est)                 # 365 x K_est
obs_anual = Y_mm.sum(1)


def calibrar_smearing(Y_raiz_tr, Yobs_tr, n_niveles=10):
    """Niveles (deciles de ŷ⁺) y, en cada uno, media de r y de r² con residuos de ENTRENAMIENTO."""
    yp = np.clip(Y_raiz_tr, 0, None)
    r = np.sqrt(Yobs_tr) - yp
    cortes = np.quantile(yp, np.linspace(0, 1, n_niveles + 1)[1:-1])
    nivel = np.digitize(yp, cortes)
    rbar = np.array([r[nivel == b].mean() for b in range(n_niveles)])
    m2 = np.array([(r[nivel == b] ** 2).mean() for b in range(n_niveles)])
    return {"cortes": cortes, "rbar": rbar, "m2": m2}


def a_mm(Y_raiz, cal):
    yp = np.clip(Y_raiz, 0, None)
    b = np.digitize(yp, cal["cortes"])
    return yp ** 2 + 2 * yp * cal["rbar"][b] + cal["m2"][b]


def mm_de_modelo(tipo):
    """Curvas diarias en mm (promedio de las 3 semillas en escala sqrt) + calibración del smearing."""
    with torch.no_grad():
        C_hat = np.mean([modelos[(tipo, s)](Xs, Xe).numpy() for s in SEMILLAS], axis=0)
    Y_raiz = C_hat @ P_dia.T                                            # N x 365, sqrt(mm/día)
    cal = calibrar_smearing(Y_raiz[idx_tr], Y_mm[idx_tr])
    return Y_raiz, a_mm(Y_raiz, cal), cal


# Satélite sin corregir en mm: el total de sus péntadas (sin pasar por curvas ni transformaciones)
pent = pd.read_csv("./salidas_preprocesamiento/estaciones_pentadas_limpio.csv",
                   usecols=["STATION_ID", "YEAR", "satelite_mm"])
sat_anual = pent.groupby(["STATION_ID", "YEAR"])["satelite_mm"].sum().reindex(
    pd.MultiIndex.from_frame(coef_est[["STATION_ID", "YEAR"]])).to_numpy()

def sesgo_anual(total):
    return 100 * (total - obs_anual) / obs_anual   # total del mismo largo que obs_anual (todas las curvas)

filas_mm = [{"modelo": "satélite sin corregir (suma de péntadas)", "particion": n,
             "pct_recortados": np.nan, "sesgo_anual_mediana_sin_correccion": np.nan,
             "sesgo_anual_mediana": np.median(sesgo_anual(sat_anual)[i]),
             "sesgo_anual_agregado": 100 * (sat_anual[i].sum() - obs_anual[i].sum()) / obs_anual[i].sum()}
            for n, i in [("entrenamiento", idx_tr), ("validacion", idx_va), ("prueba", idx_te)]]
Y_mm_pred = {}
for tipo in ("lineal", "kan", "mlp"):
    Y_raiz, Y_hat_mm, cal = mm_de_modelo(tipo)
    Y_mm_pred[tipo] = Y_hat_mm
    assert (Y_hat_mm >= 0).all()
    for n, i in [("entrenamiento", idx_tr), ("validacion", idx_va), ("prueba", idx_te)]:
        filas_mm.append({
            "modelo": f"red - cabeza {tipo}", "particion": n,
            "pct_recortados": 100 * (Y_raiz[i] < 0).mean(),
            "sesgo_anual_mediana_sin_correccion":
                np.median(100 * ((np.clip(Y_raiz[i], 0, None) ** 2).sum(1) - obs_anual[i]) / obs_anual[i]),
            "sesgo_anual_mediana": np.median(100 * (Y_hat_mm[i].sum(1) - obs_anual[i]) / obs_anual[i]),
            "sesgo_anual_agregado": 100 * (Y_hat_mm[i].sum() - obs_anual[i].sum()) / obs_anual[i].sum()})
tabla_mm = pd.DataFrame(filas_mm)
print("\nReconstrucción en mm: % de valores recortados a 0 y sesgo del acumulado anual (%):\n",
      tabla_mm.round(2).to_string(index=False))
tabla_mm.to_csv(OUT / "reconstruccion_mm.csv", index=False)


def exportar_csv_mm(tipo, idx, ruta):
    """Formato exigido por el docente: station_id, year, day, value_mm (día 1..365)."""
    Yh = Y_mm_pred[tipo][idx]
    out = pd.DataFrame({
        "station_id": np.repeat(coef_est["STATION_ID"].to_numpy()[idx], 365),
        "year": np.repeat(coef_est["YEAR"].to_numpy()[idx], 365),
        "day": np.tile(np.arange(1, 366), len(idx)),
        "value_mm": np.round(Yh.ravel(), 3)})
    out.to_csv(ruta, index=False)
    return out

# Curvas diarias en mm de TODAS las curvas, una fila por estación-año (d1..d365), para la evaluación
# 4.4 (KGE, fechas de picos, bootstrap, fases ENSO). Incluye la observación y el satélite crudo.
base_ancha = coef_est[["STATION_ID", "YEAR", "particion"]].assign(fase_enso=fase_enso)
dias_cols = [f"d{d}" for d in range(1, 366)]
pd.concat([base_ancha, pd.DataFrame(Y_mm, columns=dias_cols)], axis=1).to_csv(
    OUT / "mm_diario_observado.csv", index=False)
for tipo in ("lineal", "kan", "mlp"):
    pd.concat([base_ancha, pd.DataFrame(np.round(Y_mm_pred[tipo], 3), columns=dias_cols)], axis=1).to_csv(
        OUT / f"mm_diario_pred_{tipo}.csv", index=False)

# Mientras se conocen las estaciones reservadas, se deja el CSV de la prueba (2001-2003) con la KAN
# (el objetivo general del proyecto). El mismo exportar_csv_mm sirve para las reservadas.
exportar_csv_mm("kan", idx_te, OUT / "predicciones_prueba_kan_mm.csv")

# Ejemplo en mm, en una escala comparable: todo se agrega a PÉNTADAS (mm/día promedio de cada
# péntada), que es la resolución del satélite. Comparar barras diarias contra curvas suaves no es justo:
# la red predice el valor esperado, no el aguacero de un día. Se muestran la estación de prueba con
# error mediano de la KAN en el Cauca y en el Pacífico (región por código IDEAM: 26 = Cauca).
dias_mes = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
pent_de_dia = np.concatenate([np.minimum(np.arange(d) // 5, 5) + 6 * m for m, d in enumerate(dias_mes)])  # 365 -> 0..71
def a_pentadas(Y):                                                    # N x 365 -> N x 72 (mm/día)
    return np.stack([Y[:, pent_de_dia == k].mean(1) for k in range(72)], axis=1)

pent_sat = pd.read_csv("./salidas_preprocesamiento/estaciones_pentadas_limpio.csv",
                       usecols=["STATION_ID", "YEAR", "PENTADA", "dias", "satelite_mm"])
pent_sat = pent_sat.sort_values(["STATION_ID", "YEAR", "PENTADA"])
Sat_p = (pent_sat["satelite_mm"] / pent_sat["dias"]).to_numpy().reshape(-1, 72)   # mismo orden que coef_est
assert (pent_sat["STATION_ID"].to_numpy()[::72] == coef_est["STATION_ID"].to_numpy()).all()

region_c = np.where(coef_est["STATION_ID"].astype(str).str[:2] == "26", "Cauca", "Pacífico")
err_kan_rel = np.sqrt(np.mean([err_curva[("kan", s)] for s in SEMILLAS], axis=0) / norma2(C_est))
x_p = np.arange(1, 73)
fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
for a_, reg in zip(ax, ["Cauca", "Pacífico"]):
    cand = idx_te[region_c[idx_te] == reg]
    j = cand[np.argsort(err_kan_rel[cand])[len(cand) // 2]]           # caso típico (error mediano)
    a_.bar(x_p, a_pentadas(Y_mm[j:j + 1])[0], color="lightgray", width=0.9, label="estación observada")
    a_.plot(x_p, Sat_p[j], "k--", lw=1.5, label="satélite sin corregir")
    for tipo, col in [("lineal", "#2a6f97"), ("kan", "#b23a48"), ("mlp", "#e07a5f")]:
        a_.plot(x_p, a_pentadas(Y_mm_pred[tipo][j:j + 1])[0], color=col, lw=1.8, label=f"red - {tipo}")
    a_.set(ylabel="mm/día (promedio de la péntada)",
           title=f"{reg}: estación {coef_est['STATION_ID'].iloc[j]}, {coef_est['YEAR'].iloc[j]} (prueba, error mediano de la KAN)")
    a_.legend(fontsize=7, ncol=3)
ax[1].set_xlabel("Péntada del año (1-72)")
plt.tight_layout()
mostrar("fase4_fig06_ejemplo_pentadas")
print("\nBloque 6 listo.")
