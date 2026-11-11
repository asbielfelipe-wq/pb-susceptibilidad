#!/usr/bin/env python3
"""
Barrido de todas las variantes comunes de los genes candidatos: busca valores
atipicos de diferenciacion IMX vs EUR sin partir de una lista curada.

POR QUE UN BARRIDO
------------------
El analisis dirigido (script 11) depende de una lista de 18 variantes tomada de
la literatura. Un barrido no selecciona: entran todas las variantes comunes de
las regiones, de modo que el resultado no depende de que variantes se decidio
mirar.

METODO
------
1. Estratos separados (alta = Array+WES, wgs), por la diferencia de un orden de
   magnitud en tamano de muestra efectivo.

2. Para cada variante candidata se calcula su percentil de FST contra las
   N_VECINOS variantes de genes CONTROL con MAF media mas cercana, dentro del
   mismo estrato. El fondo control es independiente de los genes evaluados.

3. AGRUPAMIENTO POR BLOQUE. Variantes fisicamente proximas comparten haplotipo
   y no son pruebas independientes: reportar veinte atipicos que son la misma
   senal seria un error. Se agrupan las variantes atipicas contiguas separadas
   por menos de VENTANA_LD pb y se conserva solo el maximo de cada bloque.
   Esto es una aproximacion fisica al desequilibrio de ligamiento; sin
   genotipos individuales no se puede calcular r2 directamente, y conviene
   declararlo como limitacion.

4. CONTROL DE FALSOS POSITIVOS. Se estima la tasa empirica de falsos positivos
   aplicando el mismo procedimiento a las variantes CONTROL (dejando fuera, en
   cada caso, las de su propio gen). Si el barrido en genes candidatos produce
   la misma proporcion de atipicos que en genes control, no hay senal.

Correr: python3 scripts/12_barrido.py
Salidas:
  results/barrido_todas.tsv      todas las variantes con su percentil
  results/barrido_atipicos.tsv   bloques atipicos, uno por bloque
  results/barrido_resumen.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
ENTRADA = RAIZ / "data" / "processed" / "mcps_frecuencias.tsv"
RES = RAIZ / "results"

VERSION = "PLANA"
A, B = "IMX", "EUR"
N_VECINOS = 300
MAF_MIN = 0.01
PCT_UMBRAL = 99.0        # percentil a partir del cual se llama atipico
VENTANA_LD = 50000       # pb: variantes atipicas mas proximas se agrupan

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return np.where(den > 0, num / den, np.nan)


def percentiles_contra_fondo(objetivo, fondo):
    """Percentil de FST de cada fila de `objetivo` contra sus vecinos en MAF
    dentro de `fondo`. Vectorizado por busqueda ordenada."""
    f = fondo.sort_values("maf_media")
    maf_f = f["maf_media"].values
    fst_f = f["FST"].values
    n = len(f)
    k = min(N_VECINOS, n)
    out = np.empty(len(objetivo))
    for i, (m, s) in enumerate(zip(objetivo["maf_media"].values,
                                   objetivo["FST"].values)):
        j = np.searchsorted(maf_f, m)
        lo, hi = max(0, j - k // 2), min(n, max(0, j - k // 2) + k)
        lo = max(0, hi - k)
        out[i] = (fst_f[lo:hi] < s).mean() * 100
    return out


def agrupar_bloques(df_at):
    """Un bloque por racimo de atipicos separados por < VENTANA_LD pb en el
    mismo cromosoma. Devuelve el maximo de cada bloque."""
    if not len(df_at):
        return df_at.assign(bloque=[])
    d = df_at.sort_values(["chrom", "pos"]).copy()
    bloque, actual = [], 0
    prev_c, prev_p = None, None
    for c, p in zip(d["chrom"], d["pos"]):
        if prev_c != c or (p - prev_p) > VENTANA_LD:
            actual += 1
        bloque.append(actual)
        prev_c, prev_p = c, p
    d["bloque"] = bloque
    return d.loc[d.groupby("bloque")["pct"].idxmax()]


def main():
    if not ENTRADA.exists():
        sys.exit(f"Falta {ENTRADA}. Correr antes el script 10.")

    df = pd.read_csv(ENTRADA, sep="\t", low_memory=False)
    af_a, an_a = f"AF_{A}_{VERSION}", f"AN_{A}_{VERSION}"
    af_b, an_b = f"AF_{B}_{VERSION}", f"AN_{B}_{VERSION}"

    n0 = len(df)
    df = df[(df["reparto_default"] == 0) & (df["es_indel"] == 0)]
    df = df.dropna(subset=[af_a, af_b, an_a, an_b])
    df = df[(df[an_a] > 0) & (df[an_b] > 0)]
    df["maf_media"] = ((np.minimum(df[af_a], 1 - df[af_a]) +
                        np.minimum(df[af_b], 1 - df[af_b])) / 2)
    df = df[df["maf_media"] >= MAF_MIN]
    df["estrato"] = np.where(df["source"] == "WGS", "wgs", "alta")
    df["FST"] = fst_hudson(df[af_a].values, df[an_a].values,
                           df[af_b].values, df[an_b].values)
    df["dif_abs"] = (df[af_a] - df[af_b]).abs()
    df = df.dropna(subset=["FST"]).reset_index(drop=True)

    di("=" * 76)
    di("BARRIDO: TODAS LAS VARIANTES COMUNES DE LOS GENES CANDIDATOS")
    di("=" * 76)
    di(f"\n{n0:,} variantes -> {len(df):,} tras filtros")

    partes = []
    for est in ["alta", "wgs"]:
        sub = df[df["estrato"] == est]
        ctrl = sub[sub["clase"] == "control"]
        cand = sub[sub["clase"] == "candidato"].copy()
        if len(ctrl) < N_VECINOS or not len(cand):
            di(f"\nEstrato {est}: fondo insuficiente ({len(ctrl)} control)")
            continue

        cand["pct"] = percentiles_contra_fondo(cand, ctrl)
        cand["grupo"] = "candidato"

        # Control evaluado contra control de OTROS genes: tasa de fondo
        ctrl_ev = ctrl.copy()
        pcts = []
        for gen, g in ctrl.groupby("gen"):
            otros = ctrl[ctrl["gen"] != gen]
            pcts.append(pd.Series(percentiles_contra_fondo(g, otros), index=g.index))
        ctrl_ev["pct"] = pd.concat(pcts).reindex(ctrl_ev.index)
        ctrl_ev["grupo"] = "control"

        partes.append(pd.concat([cand, ctrl_ev]))

        di(f"\nEstrato {est}: {len(cand):,} candidatas, {len(ctrl):,} control")

    if not partes:
        sys.exit("Ningun estrato evaluable.")

    todo = pd.concat(partes, ignore_index=True)
    RES.mkdir(parents=True, exist_ok=True)
    cols = ["variant_id", "chrom", "pos", "gen", "clase", "grupo", "estrato",
            "source", af_a, af_b, "maf_media", "dif_abs", "FST", "pct"]
    todo[cols].sort_values("pct", ascending=False).to_csv(
        RES / "barrido_todas.tsv", sep="\t", index=False)

    # --- tasa de atipicos: candidatos vs control -------------------------
    di("\n" + "=" * 76)
    di(f"PROPORCION DE VARIANTES EN EL PERCENTIL >= {PCT_UMBRAL:.0f}")
    di("=" * 76)
    di("\nSi los genes candidatos no difieren de los control, las dos")
    di("proporciones deben ser parecidas y cercanas al valor esperado por azar.")
    di(f"\n{'estrato':<8}{'grupo':<12}{'n':>8}{'atipicos':>10}{'%':>8}"
       f"{'esperado %':>12}")
    esperado = 100 - PCT_UMBRAL
    for est in ["alta", "wgs"]:
        for grp in ["candidato", "control"]:
            s = todo[(todo["estrato"] == est) & (todo["grupo"] == grp)]
            if not len(s):
                continue
            n_at = (s["pct"] >= PCT_UMBRAL).sum()
            di(f"{est:<8}{grp:<12}{len(s):>8,}{n_at:>10}"
               f"{100*n_at/len(s):>8.2f}{esperado:>12.2f}")

    # --- bloques atipicos en genes candidatos ----------------------------
    at = todo[(todo["grupo"] == "candidato") & (todo["pct"] >= PCT_UMBRAL)].copy()
    bloques = agrupar_bloques(at)
    bloques = bloques.sort_values("pct", ascending=False)
    bloques[cols + ["bloque"]].to_csv(RES / "barrido_atipicos.tsv",
                                      sep="\t", index=False)

    di("\n" + "=" * 76)
    di("BLOQUES ATIPICOS EN GENES CANDIDATOS")
    di("=" * 76)
    di(f"\n{len(at)} variantes atipicas agrupadas en {len(bloques)} bloques "
       f"(separacion < {VENTANA_LD/1000:.0f} kb)")
    if len(bloques):
        di(f"\n{'gen':<10}{'variante':<20}{'estrato':<7}{'MAF':>7}"
           f"{'AF_IMX':>8}{'AF_EUR':>8}{'FST':>8}{'pct':>7}{'n blq':>7}")
        # tamano de cada bloque, para reportar cuantas variantes agrupa
        d2 = at.sort_values(["chrom", "pos"]).copy()
        b, actual, pc, pp = [], 0, None, None
        for c, p in zip(d2["chrom"], d2["pos"]):
            if pc != c or (p - pp) > VENTANA_LD:
                actual += 1
            b.append(actual); pc, pp = c, p
        d2["bloque"] = b
        tam = d2.groupby("bloque").size()
        for _, r in bloques.iterrows():
            di(f"{r['gen']:<10}{r['variant_id']:<20}{r['estrato']:<7}"
               f"{r['maf_media']:>7.3f}{r[af_a]:>8.3f}{r[af_b]:>8.3f}"
               f"{r['FST']:>8.4f}{r['pct']:>7.1f}{tam.get(r['bloque'],1):>7}")

    di("\n" + "=" * 76)
    di("RESUMEN")
    di("=" * 76)
    for est in ["alta", "wgs"]:
        c = todo[(todo["estrato"] == est) & (todo["grupo"] == "candidato")]
        k = todo[(todo["estrato"] == est) & (todo["grupo"] == "control")]
        if not len(c):
            continue
        pc = 100 * (c["pct"] >= PCT_UMBRAL).mean()
        pk = 100 * (k["pct"] >= PCT_UMBRAL).mean() if len(k) else float("nan")
        di(f"\nEstrato {est}: {pc:.2f}% de candidatas atipicas frente a "
           f"{pk:.2f}% de control")
        di(f"  mediana de percentil: candidatas {c['pct'].median():.1f}, "
           f"control {k['pct'].median():.1f}")

    di("\nLectura: una proporcion de atipicos en genes candidatos similar a la")
    di("de genes control indica que no hay senal de diferenciacion particular")
    di("en estas regiones. Un exceso claro, concentrado en pocos bloques,")
    di("senalaria variantes que merecen seguimiento.")

    (RES / "barrido_resumen.txt").write_text("\n".join(salida))
    print(f"\n\nTodas     -> {RES / 'barrido_todas.tsv'}")
    print(f"Atipicos  -> {RES / 'barrido_atipicos.tsv'}")
    print(f"Resumen   -> {RES / 'barrido_resumen.txt'}")


if __name__ == "__main__":
    main()
