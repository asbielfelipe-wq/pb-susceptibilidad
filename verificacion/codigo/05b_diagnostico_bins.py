#!/usr/bin/env python3
"""
Diagnostico de la estratificacion por frecuencia usada en el script 05.

Tres preguntas:

1. El intervalo mas bajo concentra 15 mil variantes con FST casi nulo. Que hay
   ahi y como se comporta segun el umbral de AC.

2. Dentro de un intervalo, la frecuencia no es homogenea. Una variante situada
   cerca del borde superior se compara contra variantes mucho mas raras, que
   por construccion no pueden alcanzar FST alto. Eso infla su percentil.

3. Que tan estable es el percentil de cada candidata si en vez de intervalos
   fijos se usa apareamiento por vecinos mas cercanos en frecuencia.

Salida: results/diagnostico_bins.txt y results/tabla4b_percentiles_apareados.tsv

Correr: python3 scripts/05b_diagnostico_bins.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RES = RAIZ / "results"
FONDO = RES / "tabla3_fondo.tsv"
CAND = RES / "tabla4_percentiles.tsv"

N_VECINOS = 500      # tamano del fondo apareado por candidata
N_BOOT = 2000        # remuestreos para el intervalo del percentil

salida = []


def di(txt=""):
    print(txt)
    salida.append(txt)


def main():
    if not FONDO.exists() or not CAND.exists():
        sys.exit("Faltan tabla3_fondo.tsv o tabla4_percentiles.tsv. Correr el 05 antes.")

    fondo = pd.read_csv(FONDO, sep="\t")
    cand = pd.read_csv(CAND, sep="\t")

    di("=" * 72)
    di("1. COMPOSICION DEL FONDO POR FRECUENCIA")
    di("=" * 72)

    fondo["AC_total"] = (fondo["AF_a"] * fondo["AN_a"] +
                         fondo["AF_b"] * fondo["AN_b"]).round()

    bins_finos = [0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50]
    fondo["bin_fino"] = pd.cut(fondo["MAF_media"], bins=bins_finos,
                               include_lowest=True)

    di(f"\n{'intervalo':<18}{'n':>8}{'% del total':>12}{'FST med':>10}"
       f"{'FST p95':>10}{'FST max':>10}")
    for b, g in fondo.groupby("bin_fino", observed=True):
        di(f"{str(b):<18}{len(g):>8}{100*len(g)/len(fondo):>11.1f}%"
           f"{g['FST'].median():>10.4f}{g['FST'].quantile(.95):>10.4f}"
           f"{g['FST'].max():>10.4f}")

    muy_raras = fondo[fondo["MAF_media"] < 0.001]
    di(f"\nVariantes con MAF media < 0.001: {len(muy_raras)} "
       f"({100*len(muy_raras)/len(fondo):.1f}% del fondo)")
    di(f"  FST mediano: {muy_raras['FST'].median():.5f}")
    di(f"  Son el {100*len(muy_raras)/len(fondo):.0f}% del fondo pero no pueden "
       f"alcanzar FST alto: con AF < 0.001 en ambas poblaciones")
    di(f"  la diferencia maxima posible es minuscula. Inflan el percentil de")
    di(f"  cualquier candidata con la que compartan intervalo.")

    di("\n" + "=" * 72)
    di("2. HETEROGENEIDAD DENTRO DEL INTERVALO")
    di("=" * 72)
    di("\nSi una candidata esta cerca del borde superior de su intervalo, se")
    di("compara contra variantes mas raras que no pueden igualar su FST.")

    bins_05 = [0, 0.01, 0.05, 0.10, 0.20, 0.35, 0.50]
    fondo["bin_05"] = pd.cut(fondo["MAF_media"], bins=bins_05, include_lowest=True)

    di(f"\n{'gen':<10}{'rsid':<12}{'MAF':>7}{'intervalo':<16}"
       f"{'MAF med fondo':>14}{'posicion':>10}")
    for _, r in cand.iterrows():
        maf_c = min(r["AF_amr"], 1 - r["AF_amr"],
                    r["AF_nfe"], 1 - r["AF_nfe"])
        maf_c = (min(r["AF_amr"], 1 - r["AF_amr"]) +
                 min(r["AF_nfe"], 1 - r["AF_nfe"])) / 2
        b = pd.cut([maf_c], bins=bins_05, include_lowest=True)[0]
        g = fondo[fondo["bin_05"] == b]
        if not len(g):
            continue
        pos = (g["MAF_media"] < maf_c).mean() * 100
        di(f"{r['gen']:<10}{r['rsid']:<12}{maf_c:>7.3f}{str(b):<16}"
           f"{g['MAF_media'].median():>14.3f}{pos:>9.0f}%")

    di("\nLa columna 'posicion' es el percentil de frecuencia de la candidata")
    di("dentro de su propio intervalo. Valores altos indican que se esta")
    di("comparando contra un fondo sistematicamente mas raro.")

    di("\n" + "=" * 72)
    di("3. PERCENTILES CON FONDO APAREADO POR FRECUENCIA")
    di("=" * 72)
    di(f"\nPara cada candidata se toman las {N_VECINOS} variantes de fondo con")
    di("MAF media mas cercana, en lugar de un intervalo fijo. Elimina el")
    di("problema de los bordes.")

    rng = np.random.default_rng(20261004)
    filas = []

    for _, r in cand.iterrows():
        maf_c = (min(r["AF_amr"], 1 - r["AF_amr"]) +
                 min(r["AF_nfe"], 1 - r["AF_nfe"])) / 2
        d = (fondo["MAF_media"] - maf_c).abs()
        vecinos = fondo.loc[d.nsmallest(N_VECINOS).index]

        pct = (vecinos["FST"] < r["FST"]).mean() * 100

        # Intervalo del percentil por bootstrap sobre el fondo apareado
        idx = rng.integers(0, len(vecinos), size=(N_BOOT, len(vecinos)))
        fst_v = vecinos["FST"].values
        pcts = (fst_v[idx] < r["FST"]).mean(axis=1) * 100
        lo, hi = np.percentile(pcts, [2.5, 97.5])

        filas.append({
            "gen": r["gen"], "rsid": r["rsid"],
            "MAF_media": round(maf_c, 4),
            "AF_amr": r["AF_amr"], "AF_nfe": r["AF_nfe"],
            "FST": r["FST"],
            "pct_bin_fijo": r["pct_FST"],
            "pct_apareado": round(pct, 1),
            "IC95_inf": round(lo, 1),
            "IC95_sup": round(hi, 1),
            "MAF_fondo_min": round(vecinos["MAF_media"].min(), 4),
            "MAF_fondo_max": round(vecinos["MAF_media"].max(), 4),
            "cambio": round(pct - r["pct_FST"], 1),
        })

    t4b = pd.DataFrame(filas).sort_values("pct_apareado", ascending=False)
    t4b.to_csv(RES / "tabla4b_percentiles_apareados.tsv", sep="\t", index=False)

    di(f"\n{'gen':<10}{'rsid':<12}{'MAF':>7}{'FST':>8}"
       f"{'bin fijo':>10}{'apareado':>10}{'IC95':>16}{'cambio':>8}")
    for _, r in t4b.iterrows():
        ic = f"[{r['IC95_inf']:.0f}-{r['IC95_sup']:.0f}]"
        di(f"{r['gen']:<10}{r['rsid']:<12}{r['MAF_media']:>7.3f}{r['FST']:>8.4f}"
           f"{r['pct_bin_fijo']:>10.1f}{r['pct_apareado']:>10.1f}{ic:>16}"
           f"{r['cambio']:>+8.1f}")

    di("\n" + "=" * 72)
    di("RESUMEN")
    di("=" * 72)
    altos_fijo = (cand["pct_FST"] >= 95).sum()
    altos_apar = (t4b["pct_apareado"] >= 95).sum()
    robustas = ((t4b["IC95_inf"] >= 95)).sum()
    di(f"\nCandidatas en percentil >=95:")
    di(f"  con intervalos fijos:      {altos_fijo} de {len(cand)}")
    di(f"  con fondo apareado:        {altos_apar} de {len(t4b)}")
    di(f"  con IC95 completo >=95:    {robustas} de {len(t4b)}")

    cambio_max = t4b.reindex(t4b["cambio"].abs().sort_values(ascending=False).index)
    di(f"\nMayores cambios al aparear:")
    for _, r in cambio_max.head(5).iterrows():
        di(f"  {r['gen']:<9} {r['rsid']:<12} "
           f"{r['pct_bin_fijo']:.1f} -> {r['pct_apareado']:.1f} "
           f"({r['cambio']:+.1f})")

    (RES / "diagnostico_bins.txt").write_text("\n".join(salida))
    print(f"\n\nGuardado en {RES / 'diagnostico_bins.txt'}")
    print(f"Tabla   -> {RES / 'tabla4b_percentiles_apareados.tsv'}")


if __name__ == "__main__":
    main()
