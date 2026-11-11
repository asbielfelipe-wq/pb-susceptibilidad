#!/usr/bin/env python3
"""
Recalcula los percentiles de las candidatas contra el fondo de genes control, y
compara las dos distribuciones de fondo entre si.

Dos preguntas:

A. Las regiones de los genes candidatos tienen una distribucion de FST distinta
   de la de genes sin relacion con metales? Si la respuesta es no, el fondo
   original era valido y el resultado nulo del script 05 se sostiene.

B. Como quedan los percentiles de las candidatas usando el fondo control, que
   es independiente por construccion?

Entradas:
  results/tabla3_fondo.tsv            fondo de genes candidatos (script 05)
  results/tabla4_percentiles.tsv      candidatas (script 05)
  data/raw/gnomad_cache_control/      fondo control (script 06)

Salidas:
  results/tabla5_fondo_control.tsv
  results/tabla6_percentiles_control.tsv
  results/comparacion_fondos.txt

Correr: python3 scripts/07_comparar_fondos.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ks_2samp

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
CACHE_CTRL = RAIZ / "data" / "raw" / "gnomad_cache_control"

POB_A, POB_B = "amr", "nfe"
AN_MINIMO_A, AN_MINIMO_B, AC_MINIMO = 10000, 20000, 10
N_VECINOS = 500
N_BOOT = 2000

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return np.where(den > 0, num / den, np.nan)


def conteos(v, pob):
    tot_ac = tot_an = 0
    for bloque in ("genome", "exome"):
        b = v.get(bloque)
        if not b:
            continue
        for p in b.get("populations") or []:
            if p["id"] == pob:
                tot_ac += p["ac"] or 0
                tot_an += p["an"] or 0
    return tot_ac, tot_an


def cargar_control():
    if not CACHE_CTRL.exists():
        sys.exit(f"Falta {CACHE_CTRL}. Correr antes el script 06.")
    vistas, filas = set(), []
    for f in sorted(CACHE_CTRL.glob("*.json")):
        d = json.loads(f.read_text())
        for v in d.get("variants") or []:
            vid = v.get("variant_id")
            if not vid or vid in vistas:
                continue
            vistas.add(vid)
            ac_a, an_a = conteos(v, POB_A)
            ac_b, an_b = conteos(v, POB_B)
            filas.append({"variant_id": vid, "origen": d.get("simbolo", f.stem),
                          "consecuencia": v.get("consequence") or "",
                          "AC_a": ac_a, "AN_a": an_a, "AC_b": ac_b, "AN_b": an_b})
    if not filas:
        sys.exit("El cache control esta vacio.")
    return pd.DataFrame(filas)


def preparar(df):
    df = df[(df["AN_a"] >= AN_MINIMO_A) & (df["AN_b"] >= AN_MINIMO_B)]
    df = df[(df["AC_a"] + df["AC_b"]) >= AC_MINIMO].copy()
    df["AF_a"] = df["AC_a"] / df["AN_a"]
    df["AF_b"] = df["AC_b"] / df["AN_b"]
    df = df[(df["AF_a"] > 0) | (df["AF_b"] > 0)].copy()
    df["dif_abs"] = (df["AF_a"] - df["AF_b"]).abs()
    df["FST"] = fst_hudson(df["AF_a"].values, df["AN_a"].values,
                           df["AF_b"].values, df["AN_b"].values)
    df["MAF_media"] = (np.minimum(df["AF_a"], 1 - df["AF_a"]) +
                       np.minimum(df["AF_b"], 1 - df["AF_b"])) / 2
    return df.dropna(subset=["FST"])


def main():
    ctrl = preparar(cargar_control())
    di("=" * 72)
    di("FONDO CONTROL")
    di("=" * 72)
    di(f"\n{len(ctrl)} variantes de {ctrl['origen'].nunique()} genes control "
       f"tras filtros")

    t3 = RES / "tabla3_fondo.tsv"
    if not t3.exists():
        sys.exit(f"Falta {t3}. Correr antes el script 05.")
    cand_fondo = pd.read_csv(t3, sep="\t").dropna(subset=["FST"])
    di(f"{len(cand_fondo)} variantes de las regiones de genes candidatos")

    # ---- A. comparacion de las dos distribuciones ------------------------
    di("\n" + "=" * 72)
    di("A. LAS DOS DISTRIBUCIONES DE FONDO SON COMPARABLES?")
    di("=" * 72)

    bins = [0, 0.01, 0.05, 0.10, 0.20, 0.35, 0.50]
    cand_fondo["bin"] = pd.cut(cand_fondo["MAF_media"], bins=bins, include_lowest=True)
    ctrl["bin"] = pd.cut(ctrl["MAF_media"], bins=bins, include_lowest=True)

    di(f"\n{'intervalo':<16}{'n cand':>8}{'n ctrl':>8}"
       f"{'FST med cand':>14}{'FST med ctrl':>14}{'p (MW)':>10}")
    for b in cand_fondo["bin"].cat.categories:
        a = cand_fondo[cand_fondo["bin"] == b]["FST"].dropna()
        c = ctrl[ctrl["bin"] == b]["FST"].dropna()
        if len(a) < 30 or len(c) < 30:
            continue
        _, p = mannwhitneyu(a, c)
        di(f"{str(b):<16}{len(a):>8}{len(c):>8}"
           f"{a.median():>14.5f}{c.median():>14.5f}{p:>10.3g}")

    ks = ks_2samp(cand_fondo["FST"].dropna(), ctrl["FST"].dropna())
    di(f"\nKolmogorov-Smirnov global: D = {ks.statistic:.4f}, p = {ks.pvalue:.3g}")
    di("Con decenas de miles de variantes, cualquier diferencia minima alcanza")
    di("significancia. Lo informativo es la comparacion de medianas por")
    di("intervalo: si son del mismo orden, los fondos son intercambiables.")

    # ---- B. percentiles contra el fondo control --------------------------
    di("\n" + "=" * 72)
    di("B. PERCENTILES DE LAS CANDIDATAS CONTRA EL FONDO CONTROL")
    di("=" * 72)

    t4 = RES / "tabla4_percentiles.tsv"
    cand = pd.read_csv(t4, sep="\t")
    rng = np.random.default_rng(20261004)
    filas = []

    for _, r in cand.iterrows():
        maf_c = (min(r["AF_amr"], 1 - r["AF_amr"]) +
                 min(r["AF_nfe"], 1 - r["AF_nfe"])) / 2
        d = (ctrl["MAF_media"] - maf_c).abs()
        vec = ctrl.loc[d.nsmallest(min(N_VECINOS, len(ctrl))).index]
        fst_v = vec["FST"].values
        pct = (fst_v < r["FST"]).mean() * 100
        idx = rng.integers(0, len(vec), size=(N_BOOT, len(vec)))
        pcts = (fst_v[idx] < r["FST"]).mean(axis=1) * 100
        lo, hi = np.percentile(pcts, [2.5, 97.5])
        filas.append({
            "gen": r["gen"], "rsid": r["rsid"], "MAF_media": round(maf_c, 4),
            "AF_amr": r["AF_amr"], "AF_nfe": r["AF_nfe"], "FST": r["FST"],
            "pct_fondo_candidato": r["pct_FST"],
            "pct_fondo_control": round(pct, 1),
            "IC95_inf": round(lo, 1), "IC95_sup": round(hi, 1),
            "n_vecinos": len(vec),
        })

    t6 = pd.DataFrame(filas).sort_values("pct_fondo_control", ascending=False)

    di(f"\n{'gen':<10}{'rsid':<12}{'MAF':>7}{'FST':>8}"
       f"{'f.candidato':>12}{'f.control':>11}{'IC95':>14}")
    for _, r in t6.iterrows():
        ic = f"[{r['IC95_inf']:.0f}-{r['IC95_sup']:.0f}]"
        di(f"{r['gen']:<10}{r['rsid']:<12}{r['MAF_media']:>7.3f}{r['FST']:>8.4f}"
           f"{r['pct_fondo_candidato']:>12.1f}{r['pct_fondo_control']:>11.1f}{ic:>14}")

    RES.mkdir(parents=True, exist_ok=True)
    ctrl.to_csv(RES / "tabla5_fondo_control.tsv", sep="\t", index=False)
    t6.to_csv(RES / "tabla6_percentiles_control.tsv", sep="\t", index=False)

    di("\n" + "=" * 72)
    di("RESUMEN")
    di("=" * 72)
    n95 = (t6["pct_fondo_control"] >= 95).sum()
    n95ic = (t6["IC95_inf"] >= 95).sum()
    di(f"\nCandidatas con percentil >=95 contra el fondo control: {n95} de {len(t6)}")
    di(f"Con IC95 completo por encima de 95:                     {n95ic} de {len(t6)}")
    dif = (t6["pct_fondo_control"] - t6["pct_fondo_candidato"]).abs()
    di(f"\nCambio medio de percentil al usar el fondo control: {dif.mean():.1f} puntos")
    di(f"Cambio maximo: {dif.max():.1f} puntos "
       f"({t6.loc[dif.idxmax(), 'gen']} {t6.loc[dif.idxmax(), 'rsid']})")

    (RES / "comparacion_fondos.txt").write_text("\n".join(salida))
    print(f"\n\nGuardado en {RES / 'comparacion_fondos.txt'}")


if __name__ == "__main__":
    main()
