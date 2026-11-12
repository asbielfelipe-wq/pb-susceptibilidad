#!/usr/bin/env python3
"""
Calibracion de la cola alta del FST nulo bajo error de asignacion.

POR QUE
-------
La simulacion del script 27 compara el FST medio estimado con el verdadero y
muestra que el error de asignacion lo atenua. Eso no basta para el analisis de
este trabajo, que no usa medias sino percentiles: identifica variantes
atipicas por su posicion en la cola alta de la distribucion de FST. Que la
media no se desplace no demuestra que la cola no se ensanche.

Este script lo comprueba directamente. Simula un conjunto grande de variantes
sin diferenciacion verdadera, con y sin error de asignacion, y compara los
percentiles altos de la distribucion de FST estimado. Si el error ensanchara
la cola, el percentil 99 del escenario con error quedaria por encima del
percentil 99 sin error, y el umbral del analisis seleccionaria variantes por
ruido de asignacion.

QUE NO CUBRE
------------
Los cromosomas se simulan de forma independiente para cada variante, de modo
que el error de asignacion no esta correlacionado entre variantes proximas ni
entre regiones. El error tampoco es asimetrico entre componentes: ambos
reciben el mismo esquema. Un error sistematico que afectara de distinta manera
a un componente, o que estuviera correlacionado a lo largo de un haplotipo, no
queda evaluado aqui y sigue siendo una limitacion declarada.

Correr: python3 scripts/31_calibracion_cola.py
Salida: results/suplementario/Tabla_S8_calibracion_cola.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
SAL = RAIZ / "results" / "suplementario"

SEMILLA = 20261111
CONCENTRACION = 20.0
PROP_IMX = 0.61
N_VARIANTES = 20000
TROZO = 500                      # variantes por trozo
TAMANOS = [2000, 20000]          # cromosomas; enmarcan el estrato de wgs
EXACTITUDES = [1.00, 0.95, 0.85]
MAF_MIN, MAF_MAX = 0.01, 0.50    # mismo rango que el barrido


def fst_hudson(p1, n1, p2, n2):
    with np.errstate(invalid="ignore", divide="ignore"):
        num = ((p1 - p2) ** 2
               - p1 * (1 - p1) / np.maximum(n1 - 1, 1e-12)
               - p2 * (1 - p2) / np.maximum(n2 - 1, 1e-12))
        den = p1 * (1 - p2) + p2 * (1 - p1)
        return np.where(den > 0, num / den, np.nan)


def simular(rng, p, n_crom, exactitud):
    """FST estimado de cada variante, sin diferenciacion verdadera.

    p es un vector de frecuencias, una por variante; las dos componentes
    comparten la misma, de modo que el FST verdadero es cero en todas.
    """
    m = len(p)
    if exactitud < 1.0:
        a = exactitud * CONCENTRACION
        b = (1 - exactitud) * CONCENTRACION
    ac_i = np.zeros(m); an_i = np.zeros(m)
    ac_e = np.zeros(m); an_e = np.zeros(m)
    paso = max(1, int(2e7 // max(m, 1)))
    hechos = 0
    while hechos < n_crom:
        k = min(paso, n_crom - hechos)
        hechos += k
        es_imx = rng.random((m, k)) < PROP_IMX
        g = (rng.random((m, k)) < p[:, None]).astype(np.float64)
        if exactitud >= 1.0:
            q = es_imx.astype(np.float64)
        else:
            q = np.where(es_imx, rng.beta(a, b, (m, k)), rng.beta(b, a, (m, k)))
        ac_i += (q * g).sum(1); an_i += q.sum(1)
        ac_e += ((1 - q) * g).sum(1); an_e += (1 - q).sum(1)
    f_i = np.divide(ac_i, an_i, out=np.full(m, np.nan), where=an_i > 0)
    f_e = np.divide(ac_e, an_e, out=np.full(m, np.nan), where=an_e > 0)
    return fst_hudson(f_i, an_i, f_e, an_e)


def main():
    rng = np.random.default_rng(SEMILLA)
    print("=" * 78)
    print("COLA ALTA DEL FST NULO CON Y SIN ERROR DE ASIGNACION")
    print("=" * 78)
    print(f"\n{N_VARIANTES:,} variantes por condicion, sin diferenciacion "
          f"verdadera, semilla {SEMILLA}")

    filas = []
    for n_crom in TAMANOS:
        print(f"\n{n_crom:,} cromosomas")
        print(f"  {'exactitud':>10}{'mediana':>10}{'p95':>10}{'p99':>10}"
              f"{'p99.9':>10}{'maxima':>10}{'sobre p99 sin error':>22}")
        ref99 = None
        for exactitud in EXACTITUDES:
            est = np.empty(N_VARIANTES)
            for i in range(0, N_VARIANTES, TROZO):
                p = rng.uniform(MAF_MIN, MAF_MAX, min(TROZO, N_VARIANTES - i))
                est[i:i + len(p)] = simular(rng, p, n_crom, exactitud)
            v = est[~np.isnan(est)]
            q50, q95, q99, q999 = np.percentile(v, [50, 95, 99, 99.9])
            if ref99 is None:
                ref99 = q99
                exceso = ""
            else:
                exceso = f"{100.0 * np.mean(v > ref99):.2f}%"
            filas.append({
                "cromosomas": n_crom, "exactitud_asignacion": exactitud,
                "variantes": len(v),
                "FST_mediana": round(float(q50), 6),
                "FST_p95": round(float(q95), 6),
                "FST_p99": round(float(q99), 6),
                "FST_p999": round(float(q999), 6),
                "FST_maxima": round(float(v.max()), 6),
                "proporcion_sobre_p99_sin_error": exceso,
            })
            print(f"  {exactitud:>10.2f}{q50:>10.5f}{q95:>10.5f}{q99:>10.5f}"
                  f"{q999:>10.5f}{v.max():>10.5f}{exceso:>22}")

    SAL.mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame(filas)
    t.to_csv(SAL / "Tabla_S8_calibracion_cola.tsv", sep="\t", index=False)

    print("\n" + "=" * 78)
    print("LECTURA")
    print("=" * 78)
    peor = 0.0
    for n_crom in TAMANOS:
        s = t[t["cromosomas"] == n_crom]
        base = float(s[s["exactitud_asignacion"] == 1.0]["FST_p99"].iloc[0])
        for _, r in s[s["exactitud_asignacion"] < 1.0].iterrows():
            peor = max(peor, float(r["FST_p99"]) / base if base > 0 else np.inf)
    if peor <= 1.0:
        print("\nCon error de asignacion, el percentil 99 de la distribucion "
              "nula de FST no\nsupera al de la asignacion perfecta en ningun "
              "tamano evaluado. El error no\nensancha la cola alta: la "
              "comprime, igual que comprime la media. Un umbral\nde percentil "
              "aplicado a esta distribucion no selecciona variantes por ruido "
              "de\nasignacion.")
    else:
        print(f"\nCon error de asignacion, el percentil 99 llega a ser "
              f"{peor:.2f} veces el de la\nasignacion perfecta. El error SI "
              f"ensancha la cola alta y el umbral del analisis\npuede estar "
              f"seleccionando variantes por ruido. Hay que decirlo en "
              f"Limitaciones.")
    print(f"\nEscrito en {SAL / 'Tabla_S8_calibracion_cola.tsv'}")


if __name__ == "__main__":
    main()
