#!/usr/bin/env python3
"""
Auditoria del filtro de frecuencia y del denominador del analisis dirigido.

MOTIVO
------
El texto del manuscrito afirmaba que cinco loci quedaron fuera del analisis de
percentiles por tener frecuencia menor alelica media inferior a 0.01. Esa
afirmacion debe comprobarse contra la salida real, porque la definicion de
"frecuencia menor alelica media" admite dos lecturas que no son equivalentes:

  definicion A (la implementada)
      ( min(p_IMX, 1-p_IMX) + min(p_EUR, 1-p_EUR) ) / 2

  definicion B
      min(pbar, 1-pbar),  con  pbar = (p_IMX + p_EUR) / 2

Cuando el alelo mayoritario cambia entre ancestrias, las dos difieren mucho.
Este script calcula ambas para los 24 loci, indica cuales quedarian excluidos
bajo cada una, y contrasta ese resultado con el numero de loci efectivamente
evaluados en results/ding_percentiles.tsv.

Tambien reporta, para cada candidata, la distancia de frecuencia a su vecino
mas cercano y mas lejano dentro del fondo apareado, que es lo que permite
juzgar si 300 vecinos constituyen un apareamiento aceptable.

Correr: python3 scripts/15_auditoria.py
Salida: results/auditoria_filtro.tsv y results/auditoria_filtro.txt
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
N_VECINOS = 300

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def main():
    f_pct = RES / "ding_percentiles.tsv"
    f_mcps = PROC / "mcps_frecuencias.tsv"
    if not f_pct.exists() or not f_mcps.exists():
        sys.exit("Faltan ding_percentiles.tsv o mcps_frecuencias.tsv.")

    pct = pd.read_csv(f_pct, sep="\t")
    mc = pd.read_csv(f_mcps, sep="\t", low_memory=False)
    af_a, af_b = f"AF_IMX_{VERSION}", f"AF_EUR_{VERSION}"

    di("=" * 76)
    di("AUDITORIA DEL FILTRO DE FRECUENCIA Y DEL DENOMINADOR")
    di("=" * 76)

    # --- cuantos loci entraron realmente al analisis ----------------------
    di(f"\nLoci presentes en ding_percentiles.tsv: {len(pct)}")
    for est, g in pct.groupby("estrato"):
        di(f"  estrato {est}: {len(g)}")

    # --- las dos definiciones de frecuencia menor alelica media ----------
    p1, p2 = pct["AF_IMX"].values, pct["AF_EUR"].values
    defA = (np.minimum(p1, 1 - p1) + np.minimum(p2, 1 - p2)) / 2
    pbar = (p1 + p2) / 2
    defB = np.minimum(pbar, 1 - pbar)

    pct = pct.assign(maf_defA=np.round(defA, 5), maf_defB=np.round(defB, 5))
    pct["coincide_con_maf_media"] = np.isclose(pct["maf_media"], pct["maf_defA"],
                                               atol=1e-4)
    pct["excluido_defA"] = pct["maf_defA"] < MAF_MIN
    pct["excluido_defB"] = pct["maf_defB"] < MAF_MIN
    pct["mayoritario_cambia"] = (p1 > 0.5) != (p2 > 0.5)

    di(f"\nLa columna maf_media de la salida coincide con la definicion A en "
       f"{int(pct['coincide_con_maf_media'].sum())} de {len(pct)} loci.")

    di(f"\n{'gen':<9}{'rsid':<13}{'AF_IMX':>8}{'AF_EUR':>8}"
       f"{'maf_A':>8}{'maf_B':>8}{'exc_A':>7}{'exc_B':>7}{'invierte':>10}")
    for _, r in pct.sort_values("maf_defA").iterrows():
        di(f"{r['gen']:<9}{r['rsid']:<13}{r['AF_IMX']:>8.4f}{r['AF_EUR']:>8.4f}"
           f"{r['maf_defA']:>8.4f}{r['maf_defB']:>8.4f}"
           f"{'si' if r['excluido_defA'] else 'no':>7}"
           f"{'si' if r['excluido_defB'] else 'no':>7}"
           f"{'si' if r['mayoritario_cambia'] else 'no':>10}")

    nA = int(pct["excluido_defA"].sum())
    nB = int(pct["excluido_defB"].sum())
    di(f"\nExcluidos bajo definicion A (la implementada): {nA}")
    di(f"Excluidos bajo definicion B:                    {nB}")
    di(f"Loci con alelo mayoritario distinto entre ancestrias: "
       f"{int(pct['mayoritario_cambia'].sum())}")

    di("\nConclusion del punto 1: el denominador del analisis dirigido es "
       f"{len(pct)} loci.")
    if nA == 0:
        di("Ningun locus fue excluido por el filtro de frecuencia. La afirmacion")
        di("del manuscrito sobre cinco loci excluidos es incorrecta y debe")
        di("retirarse; esos cinco loci SI entraron al analisis de percentiles.")

    # --- loci con ALT casi ausente en IMX, que es otra cosa ---------------
    bajos = pct[pct["AF_IMX"] < 0.002]
    di(f"\nLoci con AF_IMX < 0.002 (frecuencia baja del alelo ALT, no "
       f"exclusion): {len(bajos)}")
    for _, r in bajos.sort_values("AF_IMX").iterrows():
        di(f"  {r['gen']:<9}{r['rsid']:<13}AF_IMX {r['AF_IMX']:.4f}  "
           f"AF_EUR {r['AF_EUR']:.4f}  percentil {r['pct_FST']:.1f}")

    # --- calidad del apareamiento por frecuencia --------------------------
    di("\n" + "=" * 76)
    di("CALIDAD DEL APAREAMIENTO POR FRECUENCIA")
    di("=" * 76)

    mc = mc[(mc["reparto_default"] == 0) & (mc["es_indel"] == 0)]
    mc = mc.dropna(subset=[af_a, af_b])
    q1, q2 = mc[af_a].values, mc[af_b].values
    mc = mc.assign(maf_media=(np.minimum(q1, 1 - q1) + np.minimum(q2, 1 - q2)) / 2)
    mc = mc[mc["maf_media"] >= MAF_MIN]
    mc["estrato"] = np.where(mc["source"] == "WGS", "wgs", "alta")
    fondo = mc[mc["clase"] == "control"]

    di(f"\n{'gen':<9}{'rsid':<13}{'estrato':<7}{'maf':>8}"
       f"{'d_min':>9}{'d_max':>9}{'%fondo':>8}")
    filas = []
    for _, r in pct.iterrows():
        f = fondo[fondo["estrato"] == r["estrato"]]
        if not len(f):
            continue
        d = (f["maf_media"] - r["maf_defA"]).abs().nsmallest(
            min(N_VECINOS, len(f)))
        filas.append({
            "gen": r["gen"], "rsid": r["rsid"], "estrato": r["estrato"],
            "maf": r["maf_defA"], "d_min": round(float(d.min()), 5),
            "d_max": round(float(d.max()), 5),
            "pct_fondo_usado": round(100 * len(d) / len(f), 1),
        })
        di(f"{r['gen']:<9}{r['rsid']:<13}{r['estrato']:<7}{r['maf_defA']:>8.4f}"
           f"{d.min():>9.5f}{d.max():>9.5f}{100*len(d)/len(f):>7.1f}%")

    pd.DataFrame(filas).to_csv(RES / "auditoria_apareamiento.tsv",
                               sep="\t", index=False)
    cols = ["gen", "rsid", "variant_id", "estrato", "AF_IMX", "AF_EUR",
            "maf_media", "maf_defA", "maf_defB", "excluido_defA",
            "excluido_defB", "mayoritario_cambia", "FST", "pct_FST",
            "IC95_inf", "IC95_sup"]
    pct[[c for c in cols if c in pct.columns]].to_csv(
        RES / "auditoria_filtro.tsv", sep="\t", index=False)

    di("\nd_max es la mayor distancia de frecuencia entre la candidata y un")
    di("vecino del fondo. Valores altos indican apareamiento pobre y justifican")
    di("imponer una distancia maxima ademas del numero de vecinos.")

    (RES / "auditoria_filtro.txt").write_text("\n".join(salida))
    print(f"\n\nTabla  -> {RES / 'auditoria_filtro.tsv'}")
    print(f"Apareo -> {RES / 'auditoria_apareamiento.tsv'}")
    print(f"Resumen-> {RES / 'auditoria_filtro.txt'}")


if __name__ == "__main__":
    main()
