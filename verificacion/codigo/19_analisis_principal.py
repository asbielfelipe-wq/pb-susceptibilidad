#!/usr/bin/env python3
"""
Analisis principal definitivo: percentiles de diferenciacion IMX vs EUR con la
frecuencia de apareamiento definida sobre la frecuencia conjunta.

CAMBIO RESPECTO A VERSIONES ANTERIORES
--------------------------------------
La frecuencia usada para aparear cada locus con su fondo pasa a ser

    maf = min(pbar, 1 - pbar),    pbar = (p_IMX + p_EUR) / 2

en lugar del promedio de las frecuencias menores de cada componente,

    maf_alt = [ min(p_IMX, 1-p_IMX) + min(p_EUR, 1-p_EUR) ] / 2

El proposito del apareamiento es comparar cada candidata con variantes que
tengan una capacidad similar de alcanzar un FST dado, y esa capacidad depende
de la frecuencia conjunta. Un locus a 0.90 en un componente y 0.10 en el otro
tiene un margen de diferenciacion amplio; la segunda definicion le asigna 0.10
y lo aparea contra variantes de frecuencia baja en ambos componentes, que no
pueden alcanzar ese FST. La primera le asigna 0.50, que refleja su margen real.
Las dos divergen solo cuando el alelo mayoritario cambia entre componentes, lo
que ocurre en siete de los 24 loci.

La definicion alternativa se conserva como analisis de sensibilidad y ambas se
reportan.

El resto del procedimiento no cambia: estratificacion por fuente, apareamiento
acotado por distancia maxima de frecuencia, remuestreo por bloques fisicos.

Correr: python3 scripts/19_analisis_principal.py
Salidas:
  results/principal_percentiles.tsv
  results/principal_resumen.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"

VERSION = "PLANA"
MAF_MIN, D_MAX = 0.01, 0.02
N_OBJETIVO, N_MINIMO = 300, 50
BLOQUE_KB, N_BOOT = 50_000, 2000

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
    with np.errstate(invalid="ignore", divide="ignore"):
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        return np.where(den > 0, num / den, np.nan)


def asignar_bloques(df):
    d = df.sort_values(["chrom", "pos"]).copy()
    b, actual, pc, pp = [], 0, None, None
    for c, p in zip(d["chrom"].astype(str), d["pos"]):
        if pc != c or (p - pp) > BLOQUE_KB:
            actual += 1
        b.append(actual)
        pc, pp = c, p
    d["bloque"] = b
    return d


def main():
    f_mcps = PROC / "mcps_frecuencias.tsv"
    f_cand = RES / "percentiles_v2.tsv"
    if not f_mcps.exists() or not f_cand.exists():
        sys.exit("Faltan mcps_frecuencias.tsv o percentiles_v2.tsv")

    mc = pd.read_csv(f_mcps, sep="\t", low_memory=False)
    cand = pd.read_csv(f_cand, sep="\t")
    af_a, an_a = f"AF_IMX_{VERSION}", f"AN_IMX_{VERSION}"
    af_b, an_b = f"AF_EUR_{VERSION}", f"AN_EUR_{VERSION}"

    d = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)].copy()
    d = d.dropna(subset=[af_a, af_b, an_a, an_b])
    d = d[(d[an_a] > 0) & (d[an_b] > 0)]
    p1, p2 = d[af_a].values, d[af_b].values
    pbar = (p1 + p2) / 2
    d["maf"] = np.minimum(pbar, 1 - pbar)
    d["maf_alt"] = (np.minimum(p1, 1 - p1) + np.minimum(p2, 1 - p2)) / 2
    d["FST"] = fst_hudson(p1, d[an_a].values, p2, d[an_b].values)
    d["estrato"] = np.where(d["source"] == "WGS", "wgs", "alta")
    d = d[d["maf"] >= MAF_MIN].dropna(subset=["FST"])
    fondo = asignar_bloques(d[d["clase"] == "control"])

    di("=" * 78)
    di("ANALISIS PRINCIPAL: PERCENTILES CON FRECUENCIA CONJUNTA")
    di("=" * 78)
    di(f"\nFrecuencia de apareamiento: min(pbar, 1-pbar), pbar = (p_IMX+p_EUR)/2")
    di(f"Distancia maxima {D_MAX}, hasta {N_OBJETIVO} vecinos, minimo {N_MINIMO}")
    di(f"Remuestreo: {N_BOOT} replicas sobre bloques de {BLOQUE_KB//1000} kb")
    for est, g in fondo.groupby("estrato"):
        di(f"  fondo {est}: {len(g):,} variantes en {g['bloque'].nunique()} bloques")

    rng = np.random.default_rng(20261006)
    filas = []
    for _, c in cand.iterrows():
        f = d[d["variant_id"] == c["variant_id"]]
        if not len(f):
            continue
        r = f.iloc[0]
        g = fondo[fondo["estrato"] == r["estrato"]]
        dd = (g["maf"] - r["maf"]).abs()
        vec = g[dd <= D_MAX]
        if len(vec) > N_OBJETIVO:
            vec = vec.loc[(vec["maf"] - r["maf"]).abs().nsmallest(N_OBJETIVO).index]
        base = {"gen": c["gen"], "rsid": c["rsid"], "variant_id": c["variant_id"],
                "estrato": r["estrato"], "AF_IMX": round(float(r[af_a]), 4),
                "AF_EUR": round(float(r[af_b]), 4),
                "maf_conjunta": round(float(r["maf"]), 4),
                "maf_alternativa": round(float(r["maf_alt"]), 4),
                "invierte_mayoritario": "si" if (r[af_a] > 0.5) != (r[af_b] > 0.5) else "no",
                "FST": round(float(r["FST"]), 4), "n_vecinos": len(vec)}
        if len(vec) < N_MINIMO:
            base.update({"n_bloques": 0, "d_max": None, "pct": None,
                         "IC95_inf": None, "IC95_sup": None, "evaluable": "no"})
            filas.append(base)
            continue
        fst_v = vec["FST"].values
        pct = (fst_v < r["FST"]).mean() * 100
        grupos = [x["FST"].values for _, x in vec.groupby("bloque")]
        nb = len(grupos)
        pcts = np.empty(N_BOOT)
        for i in range(N_BOOT):
            sel = rng.integers(0, nb, nb)
            pcts[i] = (np.concatenate([grupos[j] for j in sel]) < r["FST"]).mean() * 100
        lo, hi = np.percentile(pcts, [2.5, 97.5])
        base.update({"n_bloques": nb,
                     "d_max": round(float((vec["maf"] - r["maf"]).abs().max()), 5),
                     "pct": round(pct, 1), "IC95_inf": round(lo, 1),
                     "IC95_sup": round(hi, 1), "evaluable": "si"})
        filas.append(base)

    t = pd.DataFrame(filas).sort_values("pct", ascending=False, na_position="last")
    RES.mkdir(parents=True, exist_ok=True)
    t.to_csv(RES / "principal_percentiles.tsv", sep="\t", index=False)

    di("\n" + "=" * 78)
    di("RESULTADOS")
    di("=" * 78)
    di(f"\n{'gen':<8}{'rsid':<13}{'est':<6}{'AF_IMX':>8}{'AF_EUR':>8}{'maf':>8}"
       f"{'inv':>5}{'vec':>5}{'blq':>5}{'FST':>8}{'pct':>7}{'IC95':>13}")
    for _, r in t.iterrows():
        if r["evaluable"] == "no":
            di(f"{r['gen']:<8}{r['rsid']:<13}{r['estrato']:<6}{r['AF_IMX']:>8.4f}"
               f"{r['AF_EUR']:>8.4f}{r['maf_conjunta']:>8.4f}"
               f"{r['invierte_mayoritario']:>5}{r['n_vecinos']:>5}"
               f"{'-':>5}{r['FST']:>8.4f}{'no eval':>7}")
            continue
        ic = f"[{r['IC95_inf']:.0f}-{r['IC95_sup']:.0f}]"
        di(f"{r['gen']:<8}{r['rsid']:<13}{r['estrato']:<6}{r['AF_IMX']:>8.4f}"
           f"{r['AF_EUR']:>8.4f}{r['maf_conjunta']:>8.4f}"
           f"{r['invierte_mayoritario']:>5}{r['n_vecinos']:>5}{r['n_bloques']:>5}"
           f"{r['FST']:>8.4f}{r['pct']:>7.1f}{ic:>13}")

    ev = t[t["evaluable"] == "si"]
    di("\n" + "=" * 78)
    di("RESUMEN")
    di("=" * 78)
    di(f"\nLoci evaluables: {len(ev)} de {len(t)}")
    di(f"Percentil >= 95: {int((ev['pct'] >= 95).sum())}")
    di(f"Con IC95 completo por encima de 95: {int((ev['IC95_inf'] >= 95).sum())}")
    di(f"Rango {ev['pct'].min():.1f} a {ev['pct'].max():.1f}, "
       f"mediana {ev['pct'].median():.1f}")
    di(f"Amplitud media del IC95: "
       f"{(ev['IC95_sup']-ev['IC95_inf']).mean():.1f} puntos")

    altos = ev[ev["pct"] >= 95]
    if len(altos):
        di(f"\nLoci en el percentil 95 o superior:")
        for _, r in altos.iterrows():
            di(f"  {r['gen']} {r['rsid']}: percentil {r['pct']:.1f} "
               f"(IC95 {r['IC95_inf']:.0f}-{r['IC95_sup']:.0f}), "
               f"FST {r['FST']:.4f}")
        di(f"\nCon {len(ev)} comparaciones se esperaria aproximadamente "
           f"{0.05*len(ev):.1f} locus por azar a ese umbral.")
        di("Ningun intervalo queda enteramente por encima de 95."
           if not (altos["IC95_inf"] >= 95).any() else
           "Algun intervalo queda enteramente por encima de 95.")

    inv = ev[ev["invierte_mayoritario"] == "si"]
    di(f"\n{len(inv)} loci tienen el alelo mayoritario invertido entre")
    di("componentes; son los unicos en los que esta definicion difiere de la")
    di("alternativa. Sus percentiles aqui:")
    for _, r in inv.sort_values("pct", ascending=False).iterrows():
        di(f"  {r['gen']:<8}{r['rsid']:<13}maf conjunta {r['maf_conjunta']:.4f}"
           f"  alternativa {r['maf_alternativa']:.4f}  percentil {r['pct']:.1f}")

    (RES / "principal_resumen.txt").write_text("\n".join(salida))
    print(f"\n\nTabla  -> {RES / 'principal_percentiles.tsv'}")
    print(f"Resumen-> {RES / 'principal_resumen.txt'}")


if __name__ == "__main__":
    main()
