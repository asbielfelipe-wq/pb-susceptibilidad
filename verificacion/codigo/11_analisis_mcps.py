#!/usr/bin/env python3
"""
Analisis principal: diferenciacion IMX vs EUR de los loci de susceptibilidad al
plomo definidos por Ding et al. 2016, situada dentro de la distribucion
empirica de genes control.

CONJUNTO DE VARIANTES
---------------------
Los 27 loci en 10 genes provienen de Ding et al. 2016 (PLoS One 11:e0161472,
PMID 27584680), seleccionados por esos autores con base en funcionalidad
establecida en estudios biologicos y epidemiologicos de asociacion entre plomo
y variantes geneticas. Tres no son evaluables aqui: las deleciones de GSTT1 y
GSTM1 y el repetido (GT)n del promotor de HMOX1.

Correccion de identificador: Ding reporta rs1073581 como VDR FokI, pero ese
rsID corresponde en dbSNP build 157 a una variante intergenica de
chr18:25,864,897. FokI es rs10735810, fusionado en rs2228570
(chr12:47,879,112, GRCh38), que es el analizado aqui.

DISENO
------
Comparacion IMX vs EUR dentro de la misma cohorte MCPS, por inferencia de
ancestria local. Esto elimina el confusor de mezcla de cualquier comparacion
AMR vs NFE e iguala plataforma y llamado de variantes entre los grupos.

Estratos obligatorios por tamano de muestra efectivo:
  alta = Array + WES (n_ef IMX ~91,800)
  wgs  = WGS         (n_ef IMX ~6,550)

Metricas: diferencia absoluta y FST de Hudson (Bhatia et al. 2013), con AN como
conteo efectivo de cromosomas segun los Metodos de Ziyatdinov et al. 2023.

Inferencia por percentiles contra fondo de genes control apareado por
frecuencia, con intervalo por bootstrap. No se usan pruebas de hipotesis: los
AC son fraccionarios por la asignacion probabilistica de ancestria, y con estos
AN cualquier diferencia alcanzaria significancia.

Correr: python3 scripts/11_analisis_mcps.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
ENTRADA = RAIZ / "data" / "processed" / "mcps_frecuencias.tsv"
COORD = [RAIZ / "data" / "raw" / "ding_coordenadas.tsv",
         RAIZ / "data" / "raw" / "ding_coordenadas_rescate.tsv"]
LISTA = RAIZ / "data" / "raw" / "candidatos_ding2016.tsv"
RES = RAIZ / "results"

VERSION = "PLANA"
A, B = "IMX", "EUR"
N_VECINOS = 300
N_BOOT = 2000
MAF_MIN = 0.01

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def fst_hudson(p1, n1, p2, n2):
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return np.where(den > 0, num / den, np.nan)


def main():
    if not ENTRADA.exists():
        sys.exit(f"Falta {ENTRADA}. Correr antes el script 10.")

    partes = [pd.read_csv(f, sep="\t", dtype=str) for f in COORD if f.exists()]
    if not partes:
        sys.exit("Faltan las coordenadas de Ding. Correr los scripts 13 y 13b.")
    coords = pd.concat(partes, ignore_index=True).drop_duplicates(subset="rsid")

    lista = pd.read_csv(LISTA, sep="\t", dtype=str).fillna("")
    snv = lista[lista["tipo"] == "SNV"]
    coords = coords.merge(snv[["rsid", "nomenclatura"]], on="rsid", how="left")

    df = pd.read_csv(ENTRADA, sep="\t", low_memory=False)
    af_a, an_a = f"AF_{A}_{VERSION}", f"AN_{A}_{VERSION}"
    af_b, an_b = f"AF_{B}_{VERSION}", f"AN_{B}_{VERSION}"

    n0 = len(df)
    df = df[(df["reparto_default"] == 0) & (df["es_indel"] == 0)]
    df = df.dropna(subset=[af_a, af_b, an_a, an_b])
    df = df[(df[an_a] > 0) & (df[an_b] > 0)].copy()
    df["maf_media"] = ((np.minimum(df[af_a], 1 - df[af_a]) +
                        np.minimum(df[af_b], 1 - df[af_b])) / 2)
    df["estrato"] = np.where(df["source"] == "WGS", "wgs", "alta")
    df["dif_abs"] = (df[af_a] - df[af_b]).abs()
    df["FST"] = fst_hudson(df[af_a].values, df[an_a].values,
                           df[af_b].values, df[an_b].values)
    df = df.dropna(subset=["FST"])

    di("=" * 78)
    di("DIFERENCIACION IMX vs EUR EN LOS LOCI DE DING ET AL. 2016")
    di("=" * 78)
    di(f"\nMCPS, frecuencias por ancestria local, version {VERSION}")
    di(f"{n0:,} variantes -> {len(df):,} tras filtros")

    fondo_todo = df[(df["clase"] == "control") & (df["maf_media"] >= MAF_MIN)]
    di(f"\n{'estrato':<8}{'control (fondo)':>17}{'n_ef IMX mediano':>19}")
    for est, g in fondo_todo.groupby("estrato"):
        di(f"{est:<8}{len(g):>17,}{g['n_efectivo_IMX'].median():>19,.0f}")

    # --- candidatas por variant_id ---------------------------------------
    cand = df.merge(coords[["rsid", "variant_id", "nomenclatura", "gen"]],
                    on="variant_id", how="inner", suffixes=("", "_lit"))
    cand["gen_lit"] = cand["gen_lit"].fillna(cand["gen"])

    di(f"\n{len(cand)} de {len(snv)} loci SNV localizados en MCPS")
    faltan = set(coords["rsid"]) - set(cand["rsid"])
    if faltan:
        di(f"No evaluables ({len(faltan)}): {', '.join(sorted(faltan))}")
        di("  (reparto por defecto, indel, o sin frecuencias)")

    bajas = cand[cand["maf_media"] < MAF_MIN]
    cand = cand[cand["maf_media"] >= MAF_MIN]
    if len(bajas):
        di(f"\n{len(bajas)} con MAF media < {MAF_MIN}, excluidas de la "
           f"comparacion de percentiles:")
        for _, r in bajas.sort_values("gen_lit").iterrows():
            di(f"  {r['gen_lit']:<8}{r['rsid']:<13}{r['nomenclatura']:<22}"
               f"AF_IMX {r[af_a]:.4f}  AF_EUR {r[af_b]:.4f}")
        di("  Estas variantes estan practicamente ausentes del componente")
        di("  indigena; su ausencia es el dato, no su diferenciacion.")

    # --- percentiles ------------------------------------------------------
    rng = np.random.default_rng(20261004)
    filas = []
    for _, r in cand.iterrows():
        fondo = fondo_todo[fondo_todo["estrato"] == r["estrato"]]
        if len(fondo) < 50:
            continue
        d = (fondo["maf_media"] - r["maf_media"]).abs()
        vec = fondo.loc[d.nsmallest(min(N_VECINOS, len(fondo))).index]
        fst_v = vec["FST"].values
        pct = (fst_v < r["FST"]).mean() * 100
        idx = rng.integers(0, len(vec), size=(N_BOOT, len(vec)))
        pcts = (fst_v[idx] < r["FST"]).mean(axis=1) * 100
        lo, hi = np.percentile(pcts, [2.5, 97.5])
        filas.append({
            "gen": r["gen_lit"], "rsid": r["rsid"],
            "nomenclatura": r["nomenclatura"],
            "variant_id": r["variant_id"],
            "estrato": r["estrato"], "source": r["source"],
            "maf_media": round(r["maf_media"], 4),
            f"AF_{A}": round(r[af_a], 4), f"AF_{B}": round(r[af_b], 4),
            "dif_abs": round(r["dif_abs"], 4),
            "FST": round(r["FST"], 4),
            "pct_FST": round(pct, 1),
            "IC95_inf": round(lo, 1), "IC95_sup": round(hi, 1),
            "n_fondo": len(vec),
            "FST_fondo_p95": round(float(np.percentile(fst_v, 95)), 4),
        })

    t = pd.DataFrame(filas).sort_values("pct_FST", ascending=False)
    RES.mkdir(parents=True, exist_ok=True)
    t.to_csv(RES / "ding_percentiles.tsv", sep="\t", index=False)
    cand.to_csv(RES / "ding_candidatas.tsv", sep="\t", index=False)

    di("\n" + "=" * 78)
    di("PERCENTIL DE FST CONTRA FONDO CONTROL APAREADO POR FRECUENCIA")
    di("=" * 78)
    for est in ["alta", "wgs"]:
        s = t[t["estrato"] == est]
        if not len(s):
            continue
        di(f"\nEstrato {est} ({len(s)} loci):")
        di(f"{'gen':<8}{'rsid':<13}{'nomencl.':<16}{'AF_IMX':>8}{'AF_EUR':>8}"
           f"{'dif':>8}{'FST':>8}{'pct':>7}{'IC95':>12}")
        for _, r in s.iterrows():
            ic = f"[{r['IC95_inf']:.0f}-{r['IC95_sup']:.0f}]"
            nom = (r["nomenclatura"] or "")[:15]
            di(f"{r['gen']:<8}{r['rsid']:<13}{nom:<16}{r['AF_IMX']:>8.3f}"
               f"{r['AF_EUR']:>8.3f}{r['dif_abs']:>8.3f}{r['FST']:>8.4f}"
               f"{r['pct_FST']:>7.1f}{ic:>12}")

    di("\n" + "=" * 78)
    di("RESUMEN")
    di("=" * 78)
    di(f"\nLoci evaluados: {len(t)} de {len(snv)} SNV del conjunto")
    di(f"Percentil >= 95:                    {(t['pct_FST'] >= 95).sum()}")
    di(f"Con IC95 completo por encima de 95: {(t['IC95_inf'] >= 95).sum()}")
    if len(t):
        di(f"\nRango de percentiles: {t['pct_FST'].min():.1f} a "
           f"{t['pct_FST'].max():.1f}   mediana {t['pct_FST'].median():.1f}")
    di("\nPercentiles repartidos sin concentracion en el extremo superior")
    di("indican que la diferenciacion de estos loci no excede la de variantes")
    di("de frecuencia comparable en genes sin relacion con metales.")

    (RES / "ding_resumen.txt").write_text("\n".join(salida))
    print(f"\n\nPercentiles -> {RES / 'ding_percentiles.tsv'}")
    print(f"Resumen     -> {RES / 'ding_resumen.txt'}")


if __name__ == "__main__":
    main()
