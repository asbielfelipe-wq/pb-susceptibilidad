#!/usr/bin/env python3
"""
Percentiles con apareamiento acotado por distancia de frecuencia y bootstrap
por bloques, en respuesta a la revision.

QUE CAMBIA RESPECTO AL SCRIPT 11
--------------------------------
1. APAREAMIENTO ACOTADO. Antes se tomaban siempre 300 vecinos, lo que en el
   estrato de alta precision equivale al 26.7% del fondo y produce vecinos a
   distancias de frecuencia de hasta 0.19. Ahora se exige ademas que la
   distancia no supere D_MAX, y se reporta cuantos vecinos quedaron. Un locus
   con menos de N_MINIMO vecinos dentro de esa ventana se marca como no
   evaluable con precision suficiente en lugar de recibir un percentil
   aparentemente valido.

2. BOOTSTRAP POR BLOQUES. Antes se remuestreaban variantes control
   individuales. Variantes proximas comparten haplotipo, de modo que 300
   variantes no son 300 observaciones independientes y el intervalo resultaba
   demasiado estrecho. Ahora el remuestreo es por bloque fisico (ventanas de
   BLOQUE_KB), que es la unidad aproximada de independencia.

3. DEFINICION EXPLICITA DE LA FRECUENCIA DE APAREAMIENTO. Se calcula
       maf = ( min(p_IMX, 1-p_IMX) + min(p_EUR, 1-p_EUR) ) / 2
   y se reporta tambien min(pbar, 1-pbar) para los loci donde el alelo
   mayoritario cambia entre ancestrias, porque ambas definiciones divergen ahi.

4. SENSIBILIDAD AL NUMERO DE VECINOS. El percentil se recalcula con varios
   tamanos de vecindario para mostrar que no depende de esa eleccion.

Correr: python3 scripts/16_percentiles_v2.py
Salidas:
  results/percentiles_v2.tsv
  results/percentiles_v2_sensibilidad.tsv
  results/percentiles_v2.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
PROC = RAIZ / "data" / "processed"

VERSION = "PLANA"
MAF_MIN = 0.01
D_MAX = 0.02          # distancia maxima de frecuencia admitida al fondo
N_OBJETIVO = 300      # vecinos buscados si caben dentro de D_MAX
N_MINIMO = 50         # por debajo de esto el percentil no se reporta
BLOQUE_KB = 50_000    # unidad de remuestreo
N_BOOT = 2000
VECINDARIOS = [100, 200, 300, 500]

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
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
    f_pct = RES / "ding_percentiles.tsv"
    f_mcps = PROC / "mcps_frecuencias.tsv"
    if not f_pct.exists() or not f_mcps.exists():
        sys.exit("Faltan ding_percentiles.tsv o mcps_frecuencias.tsv.")

    cand = pd.read_csv(f_pct, sep="\t")
    mc = pd.read_csv(f_mcps, sep="\t", low_memory=False)
    af_a, an_a = f"AF_IMX_{VERSION}", f"AN_IMX_{VERSION}"
    af_b, an_b = f"AF_EUR_{VERSION}", f"AN_EUR_{VERSION}"

    mc = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)]
    mc = mc.dropna(subset=[af_a, af_b, an_a, an_b])
    mc = mc[(mc[an_a] > 0) & (mc[an_b] > 0)].copy()
    p1, p2 = mc[af_a].values, mc[af_b].values
    mc["maf"] = (np.minimum(p1, 1 - p1) + np.minimum(p2, 1 - p2)) / 2
    mc["FST"] = fst_hudson(p1, mc[an_a].values, p2, mc[an_b].values)
    mc["estrato"] = np.where(mc["source"] == "WGS", "wgs", "alta")
    mc = mc[(mc["maf"] >= MAF_MIN)].dropna(subset=["FST"])
    fondo = asignar_bloques(mc[mc["clase"] == "control"])

    di("=" * 78)
    di("PERCENTILES CON APAREAMIENTO ACOTADO Y BOOTSTRAP POR BLOQUES")
    di("=" * 78)
    di(f"\nDistancia maxima de frecuencia admitida: {D_MAX}")
    di(f"Vecinos buscados: {N_OBJETIVO}   minimo exigido: {N_MINIMO}")
    di(f"Unidad de remuestreo: bloques de {BLOQUE_KB//1000} kb")
    for est, g in fondo.groupby("estrato"):
        di(f"  fondo {est}: {len(g):,} variantes en {g['bloque'].nunique():,} bloques")

    rng = np.random.default_rng(20261005)
    filas, sens = [], []

    for _, r in cand.iterrows():
        f = fondo[fondo["estrato"] == r["estrato"]]
        d = (f["maf"] - r["maf_media"]).abs()
        dentro = f[d <= D_MAX]
        if len(dentro) > N_OBJETIVO:
            idx = (dentro["maf"] - r["maf_media"]).abs().nsmallest(N_OBJETIVO).index
            dentro = dentro.loc[idx]

        base = {
            "gen": r["gen"], "rsid": r["rsid"], "variant_id": r["variant_id"],
            "estrato": r["estrato"], "maf": r["maf_media"],
            "AF_IMX": r["AF_IMX"], "AF_EUR": r["AF_EUR"], "FST": r["FST"],
            "pct_v1_300vecinos": r["pct_FST"],
            "IC95_v1_inf": r["IC95_inf"], "IC95_v1_sup": r["IC95_sup"],
            "n_vecinos_v2": len(dentro),
            "n_bloques_v2": int(dentro["bloque"].nunique()) if len(dentro) else 0,
            "d_max_v2": round(float((dentro["maf"] - r["maf_media"]).abs().max()), 5)
                        if len(dentro) else None,
        }

        if len(dentro) < N_MINIMO:
            base.update({"pct_v2": None, "IC95_v2_inf": None,
                         "IC95_v2_sup": None, "evaluable": "no"})
            filas.append(base)
            continue

        fst_v = dentro["FST"].values
        pct = (fst_v < r["FST"]).mean() * 100

        # bootstrap por bloque: se remuestrean bloques, no variantes
        grupos = [g["FST"].values for _, g in dentro.groupby("bloque")]
        nb = len(grupos)
        pcts = np.empty(N_BOOT)
        for i in range(N_BOOT):
            sel = rng.integers(0, nb, nb)
            v = np.concatenate([grupos[j] for j in sel])
            pcts[i] = (v < r["FST"]).mean() * 100
        lo, hi = np.percentile(pcts, [2.5, 97.5])

        base.update({"pct_v2": round(pct, 1), "IC95_v2_inf": round(lo, 1),
                     "IC95_v2_sup": round(hi, 1), "evaluable": "si"})
        filas.append(base)

        # sensibilidad al tamano del vecindario, siempre dentro de D_MAX
        for n in VECINDARIOS:
            sub = f[(f["maf"] - r["maf_media"]).abs() <= D_MAX]
            if len(sub) > n:
                i2 = (sub["maf"] - r["maf_media"]).abs().nsmallest(n).index
                sub = sub.loc[i2]
            if len(sub) >= N_MINIMO:
                sens.append({"gen": r["gen"], "rsid": r["rsid"],
                             "estrato": r["estrato"], "n_vecinos": n,
                             "n_usados": len(sub),
                             "pct": round((sub["FST"].values < r["FST"]).mean() * 100, 1)})

    t = pd.DataFrame(filas)
    RES.mkdir(parents=True, exist_ok=True)
    t.to_csv(RES / "percentiles_v2.tsv", sep="\t", index=False)
    pd.DataFrame(sens).to_csv(RES / "percentiles_v2_sensibilidad.tsv",
                              sep="\t", index=False)

    di("\n" + "=" * 78)
    di("COMPARACION CON EL ANALISIS ORIGINAL")
    di("=" * 78)
    di(f"\n{'gen':<9}{'rsid':<13}{'est':<5}{'vec':>5}{'blq':>5}{'dmax':>8}"
       f"{'pct v1':>8}{'pct v2':>8}{'IC95 v2':>14}{'camb':>7}")
    for _, r in t.sort_values("pct_v2", ascending=False, na_position="last").iterrows():
        if r["evaluable"] == "no":
            di(f"{r['gen']:<9}{r['rsid']:<13}{r['estrato']:<5}"
               f"{r['n_vecinos_v2']:>5}{'-':>5}{'-':>8}"
               f"{r['pct_v1_300vecinos']:>8.1f}{'no eval':>8}")
            continue
        ic = f"[{r['IC95_v2_inf']:.0f}-{r['IC95_v2_sup']:.0f}]"
        camb = r["pct_v2"] - r["pct_v1_300vecinos"]
        di(f"{r['gen']:<9}{r['rsid']:<13}{r['estrato']:<5}"
           f"{r['n_vecinos_v2']:>5}{r['n_bloques_v2']:>5}{r['d_max_v2']:>8.4f}"
           f"{r['pct_v1_300vecinos']:>8.1f}{r['pct_v2']:>8.1f}{ic:>14}{camb:>+7.1f}")

    ev = t[t["evaluable"] == "si"]
    di("\n" + "=" * 78)
    di("RESUMEN")
    di("=" * 78)
    di(f"\nLoci evaluables con apareamiento acotado: {len(ev)} de {len(t)}")
    noev = t[t["evaluable"] == "no"]
    if len(noev):
        di(f"No evaluables (menos de {N_MINIMO} vecinos dentro de "
           f"{D_MAX} de distancia):")
        for _, r in noev.iterrows():
            di(f"  {r['gen']:<9}{r['rsid']:<13}maf {r['maf']:.4f}  "
               f"vecinos {r['n_vecinos_v2']}")
    if len(ev):
        di(f"\nPercentil >= 95: {(ev['pct_v2'] >= 95).sum()}")
        di(f"Con IC95 completo por encima de 95: {(ev['IC95_v2_inf'] >= 95).sum()}")
        amp_v1 = (ev["IC95_v1_sup"] - ev["IC95_v1_inf"]).mean()
        amp_v2 = (ev["IC95_v2_sup"] - ev["IC95_v2_inf"]).mean()
        di(f"\nAmplitud media del IC95: {amp_v1:.1f} puntos con bootstrap por "
           f"variante, {amp_v2:.1f} con bootstrap por bloque")
        di("El bootstrap por bloque reconoce que variantes proximas no son")
        di("observaciones independientes y produce intervalos mas anchos.")
        cam = (ev["pct_v2"] - ev["pct_v1_300vecinos"]).abs()
        di(f"\nCambio absoluto medio de percentil: {cam.mean():.1f} puntos "
           f"(maximo {cam.max():.1f})")

    (RES / "percentiles_v2.txt").write_text("\n".join(salida))
    print(f"\n\nTabla       -> {RES / 'percentiles_v2.tsv'}")
    print(f"Sensibilidad-> {RES / 'percentiles_v2_sensibilidad.tsv'}")
    print(f"Resumen     -> {RES / 'percentiles_v2.txt'}")


if __name__ == "__main__":
    main()
