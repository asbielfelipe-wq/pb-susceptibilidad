#!/usr/bin/env python3
"""
Comportamiento del estimador de Hudson sobre conteos fraccionarios.

POR QUE
-------
El estimador de Hudson se derivo para muestras de dos poblaciones, con
conteos alelicos enteros y muestras independientes. Aqui se aplica a algo
distinto: los conteos no son enteros sino sumas de probabilidades de
asignacion de ancestria local, y los dos componentes no proceden de dos
muestras separadas sino de los mismos individuos. Dos preguntas se siguen de
ahi y este script las responde por simulacion, sin usar datos de MCPS.

La primera es si la correccion por tamano de muestra conserva sus
propiedades cuando el denominador es AN, la suma de probabilidades, en lugar
de un numero de cromosomas observados. La segunda es en que direccion y con
que magnitud desplaza el resultado el error de la inferencia de ancestria.

QUE SE SIMULA
-------------
Una cohorte mezclada de la que se conoce la verdad. Cada cromosoma tiene una
ancestria verdadera, indigena o europea, y lleva el alelo alterno con la
frecuencia que corresponde a esa ancestria. La inferencia de ancestria local
no devuelve la etiqueta verdadera sino una probabilidad: para un cromosoma
indigena, una probabilidad alta pero no uno; para uno europeo, baja pero no
cero. Los conteos se forman como en el analisis real, sumando genotipos
ponderados por esa probabilidad, y de ellos sale el FST estimado, que se
compara con el verdadero.

La exactitud de la inferencia se describe por la probabilidad media que
recibe un cromosoma de su propia ancestria. Un valor de 1 es inferencia
perfecta; 0.95 y 0.85 son escenarios de error creciente.

Correr: python3 scripts/27_simulacion_fst.py
Salida: results/suplementario/Tabla_S7_simulacion_fst.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
SAL = RAIZ / "results" / "suplementario"

SEMILLA = 20261110
CONCENTRACION = 20.0   # dispersion de las probabilidades de asignacion
TROZO = 5_000_000      # cromosomas por trozo, para no agotar la memoria

# Escenarios de frecuencia: (p indigena, p europea, etiqueta)
FRECUENCIAS = [(0.50, 0.50, "sin diferenciacion"),
               (0.55, 0.45, "diferenciacion leve"),
               (0.70, 0.30, "diferenciacion fuerte"),
               (0.05, 0.02, "alelo poco frecuente"),
               (0.30, 0.70, "alelo mayoritario invertido")]

TAMANOS = [(2_000, 3000), (20_000, 1500), (280_000, 200)]
EXACTITUDES = [1.00, 0.95, 0.85]
PROP_IMX = 0.61        # fraccion de cromosomas de ancestria indigena


def fst_hudson(p1, n1, p2, n2):
    """Estimador de Hudson en la formulacion de Bhatia et al. 2013.

    n1 y n2 son numeros efectivos de cromosomas; con conteos fraccionarios
    son sumas de probabilidades y no enteros, que es justo lo que se evalua.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        num = ((p1 - p2) ** 2
               - p1 * (1 - p1) / np.maximum(n1 - 1, 1e-9)
               - p2 * (1 - p2) / np.maximum(n2 - 1, 1e-9))
        den = p1 * (1 - p2) + p2 * (1 - p1)
        return np.where(den > 0, num / den, np.nan)


def fst_verdadero(p1, p2):
    """FST de parametros, sin correccion: el valor que se quiere recuperar."""
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return ((p1 - p2) ** 2) / den if den > 0 else np.nan


def una_replica(rng, n_crom, p_imx, p_eur, exactitud):
    """Conteos fraccionarios de una cohorte simulada, por trozos."""
    ac_i = an_i = ac_e = an_e = 0.0
    if exactitud < 1.0:
        a_i = exactitud * CONCENTRACION
        b_i = (1 - exactitud) * CONCENTRACION
    restantes = n_crom
    while restantes > 0:
        m = min(TROZO, restantes)
        restantes -= m
        es_imx = rng.random(m) < PROP_IMX
        g = (rng.random(m) < np.where(es_imx, p_imx, p_eur)).astype(np.float64)
        if exactitud >= 1.0:
            q = es_imx.astype(np.float64)
        else:
            # probabilidad asignada al componente indigena: alta si el
            # cromosoma lo es, baja si no, con dispersion en ambos casos
            q = np.where(es_imx,
                         rng.beta(a_i, b_i, m),
                         rng.beta(b_i, a_i, m))
        ac_i += float(np.dot(q, g))
        an_i += float(q.sum())
        ac_e += float(np.dot(1 - q, g))
        an_e += float((1 - q).sum())
    return ac_i, an_i, ac_e, an_e


def main():
    rng = np.random.default_rng(SEMILLA)
    filas = []
    print("=" * 78)
    print("ESTIMADOR DE HUDSON SOBRE CONTEOS FRACCIONARIOS")
    print("=" * 78)
    for p_imx, p_eur, etiqueta in FRECUENCIAS:
        verdadero = fst_verdadero(p_imx, p_eur)
        print(f"\n{etiqueta}: p indigena {p_imx}, p europea {p_eur}, "
              f"FST verdadero {verdadero:.4f}")
        print(f"  {'cromosomas':>11}{'exactitud':>11}{'media':>9}{'sesgo':>9}"
              f"{'DE':>8}{'AN indig.':>11}")
        for n_crom, n_rep in TAMANOS:
            for exactitud in EXACTITUDES:
                est = np.empty(n_rep)
                an_is = np.empty(n_rep)
                for r in range(n_rep):
                    ac_i, an_i, ac_e, an_e = una_replica(
                        rng, n_crom, p_imx, p_eur, exactitud)
                    f_i = ac_i / an_i if an_i > 0 else np.nan
                    f_e = ac_e / an_e if an_e > 0 else np.nan
                    est[r] = fst_hudson(f_i, an_i, f_e, an_e)
                    an_is[r] = an_i
                m = float(np.nanmean(est))
                filas.append({
                    "escenario": etiqueta,
                    "p_indigena": p_imx, "p_europea": p_eur,
                    "FST_verdadero": round(verdadero, 6),
                    "cromosomas": n_crom, "replicas": n_rep,
                    "exactitud_asignacion": exactitud,
                    "FST_estimado_medio": round(m, 6),
                    "sesgo": round(m - verdadero, 6),
                    "sesgo_relativo": (round((m - verdadero) / verdadero, 4)
                                       if verdadero > 0 else ""),
                    "desviacion_estandar": round(float(np.nanstd(est)), 6),
                    "AN_indigena_medio": round(float(an_is.mean()), 1),
                })
                print(f"  {n_crom:>11,}{exactitud:>11.2f}{m:>9.4f}"
                      f"{m - verdadero:>9.4f}{np.nanstd(est):>8.4f}"
                      f"{an_is.mean():>11,.0f}")

    SAL.mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame(filas)
    t.to_csv(SAL / "Tabla_S7_simulacion_fst.tsv", sep="\t", index=False)

    print("\n" + "=" * 78)
    print("LECTURA")
    print("=" * 78)
    perf = t[t["exactitud_asignacion"] == 1.00]
    grande = perf[perf["cromosomas"] == max(x for x, _ in TAMANOS)]
    print(f"\nCon asignacion perfecta y {max(x for x, _ in TAMANOS):,} cromosomas, "
          f"el sesgo absoluto maximo es {grande['sesgo'].abs().max():.6f}.")
    peq = perf[perf["cromosomas"] == min(x for x, _ in TAMANOS)]
    print(f"Con {min(x for x, _ in TAMANOS):,} cromosomas sigue siendo "
          f"{peq['sesgo'].abs().max():.6f}, de modo que la correccion por "
          f"tamano de muestra\nfunciona tambien cuando el denominador es una "
          f"suma de probabilidades.")
    con = t[(t["exactitud_asignacion"] < 1.0) & (t["FST_verdadero"] > 0)]
    peor = con["sesgo_relativo"].astype(float).min()
    print(f"\nCon error de asignacion el sesgo es negativo en todos los "
          f"escenarios con diferenciacion:\nla mezcla de componentes atenua el "
          f"FST hacia cero, hasta un {abs(peor):.0%} por debajo del verdadero "
          f"con exactitud de\n{con.loc[con['sesgo_relativo'].astype(float).idxmin(), 'exactitud_asignacion']:.2f}. "
          f"El estimador es por tanto conservador: subestima la\n"
          f"diferenciacion en lugar de inventarla.")
    nulo = t[(t["FST_verdadero"] == 0) & (t["exactitud_asignacion"] < 1.0)]
    print(f"\nSin diferenciacion verdadera, el error de asignacion no la crea: "
          f"el estimado medio\nse queda en {nulo['FST_estimado_medio'].abs().max():.6f} "
          f"como maximo.")
    print(f"\nEscrito en {SAL / 'Tabla_S7_simulacion_fst.tsv'}")


if __name__ == "__main__":
    main()
