#!/usr/bin/env python3
"""
Situa la diferenciacion AMR vs NFE de las variantes candidatas dentro de la
distribucion empirica de todas las variantes de las mismas regiones genicas.

Motivacion: con AN de decenas de miles, practicamente cualquier variante alcanza
significancia al comparar poblaciones. La pregunta util no es si difiere, sino
si difiere mas de lo que difiere una variante cualquiera de frecuencia
comparable. Eso exige una distribucion de fondo.

El fondo sale de los archivos que ya estan en disco: las respuestas crudas de
las consultas por gen y por ventana contienen todas las variantes de esas
regiones, no solo las candidatas.

Metricas por variante:
  dif_abs   |AF_amr - AF_nfe|
  razon_AF  AF_amr / AF_nfe
  FST       estimador de Hudson, insesgado respecto al tamano de muestra

El FST de Hudson es la metrica principal: no se infla con frecuencias bajas
como la razon, y corrige por tamano de muestra.

Entradas:
  data/raw/gnomad_respuestas.json
  data/raw/gnomad_cache/*.json
  results/tabla1_frecuencias.tsv

Salidas:
  results/tabla3_fondo.tsv        todas las variantes de fondo con sus metricas
  results/tabla4_percentiles.tsv  cada candidata con su percentil

Correr: python3 scripts/05_fondo_y_percentiles.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
RAW = RAIZ / "data" / "raw"
RES = RAIZ / "results"

POB_A, POB_B = "amr", "nfe"

# Umbrales para que una variante entre al fondo. Con AN bajo las frecuencias
# son inestables y ensucian la distribucion.
AN_MINIMO_A = 10000
AN_MINIMO_B = 20000
AC_MINIMO = 10          # excluye singletons y casi-singletons

# Intervalos de frecuencia media para estratificar. La diferenciacion depende
# fuertemente de la frecuencia, asi que comparar sin estratificar es enganoso.
BINS = [0, 0.01, 0.05, 0.10, 0.20, 0.35, 0.50]


def fst_hudson(p1, n1, p2, n2):
    """Estimador de Hudson (Bhatia et al. 2013). Puede dar valores negativos
    cuando la diferenciacion real es nula; no se truncan, porque truncar sesga
    la distribucion de fondo hacia arriba."""
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return np.where(den > 0, num / den, np.nan)


def pobs(variante, bloque, pob):
    b = variante.get(bloque)
    if not b:
        return 0, 0
    for p in b.get("populations") or []:
        if p["id"] == pob:
            return p["ac"] or 0, p["an"] or 0
    return 0, 0


def conteos(variante, pob):
    ac_g, an_g = pobs(variante, "genome", pob)
    ac_e, an_e = pobs(variante, "exome", pob)
    return ac_g + ac_e, an_g + an_e


def recolectar():
    """Devuelve todas las variantes encontradas en los archivos en disco."""
    vistas, filas = set(), []

    def agregar(v, origen):
        vid = v.get("variant_id")
        if not vid or vid in vistas:
            return
        vistas.add(vid)
        ac_a, an_a = conteos(v, POB_A)
        ac_b, an_b = conteos(v, POB_B)
        filas.append({
            "variant_id": vid,
            "rsids": ";".join(v.get("rsids") or []),
            "consecuencia": v.get("consequence") or "",
            "origen": origen,
            "AC_a": ac_a, "AN_a": an_a,
            "AC_b": ac_b, "AN_b": an_b,
        })

    # Respuestas de la consulta por gen
    f = RAW / "gnomad_respuestas.json"
    if f.exists():
        d = json.loads(f.read_text())
        for simbolo, gen in (d.get("respuestas") or {}).items():
            if not gen:
                continue
            for v in gen.get("variants") or []:
                agregar(v, simbolo)
    else:
        print(f"aviso: falta {f.name}", file=sys.stderr)

    # Caches de region y de ventana
    cache = RAW / "gnomad_cache"
    if cache.exists():
        for f in sorted(cache.glob("*.json")):
            d = json.loads(f.read_text())
            simbolo = d.get("simbolo", f.stem)
            vs = (d.get("region") or {}).get("variants") or d.get("variants") or []
            for v in vs:
                agregar(v, simbolo)

    if not filas:
        sys.exit("No se encontraron variantes en disco. "
                 "Revisar data/raw/gnomad_respuestas.json y data/raw/gnomad_cache/")
    return pd.DataFrame(filas)


def main():
    df = recolectar()
    print(f"{len(df)} variantes unicas recuperadas de los archivos en disco")

    # --- filtros de calidad ----------------------------------------------
    antes = len(df)
    df = df[(df["AN_a"] >= AN_MINIMO_A) & (df["AN_b"] >= AN_MINIMO_B)]
    df = df[(df["AC_a"] + df["AC_b"]) >= AC_MINIMO]
    print(f"{len(df)} pasan los filtros de AN y AC ({antes - len(df)} excluidas)")

    df["AF_a"] = df["AC_a"] / df["AN_a"]
    df["AF_b"] = df["AC_b"] / df["AN_b"]
    df = df[(df["AF_a"] > 0) | (df["AF_b"] > 0)]

    df["dif_abs"] = (df["AF_a"] - df["AF_b"]).abs()
    df["razon_AF"] = np.where(df["AF_b"] > 0, df["AF_a"] / df["AF_b"], np.nan)
    df["FST"] = fst_hudson(df["AF_a"].values, df["AN_a"].values,
                           df["AF_b"].values, df["AN_b"].values)

    # MAF media de las dos poblaciones, para estratificar
    maf = pd.concat([df["AF_a"], 1 - df["AF_a"]], axis=1).min(axis=1)
    maf_b = pd.concat([df["AF_b"], 1 - df["AF_b"]], axis=1).min(axis=1)
    df["MAF_media"] = (maf + maf_b) / 2
    df["bin_MAF"] = pd.cut(df["MAF_media"], bins=BINS, include_lowest=True)

    # --- candidatas -------------------------------------------------------
    t1 = RES / "tabla1_frecuencias.tsv"
    if not t1.exists():
        sys.exit(f"Falta {t1}. Correr antes el script 04.")
    cand = pd.read_csv(t1, sep="\t")
    ids_cand = set(cand["variant_id"])

    df["es_candidata"] = df["variant_id"].isin(ids_cand)
    fondo = df[~df["es_candidata"]].copy()

    print(f"\nFondo: {len(fondo)} variantes")
    print(f"Candidatas localizadas en el fondo: "
          f"{df['es_candidata'].sum()} de {len(ids_cand)}")

    print("\nDistribucion del fondo por intervalo de MAF:")
    print(f"{'intervalo':<16}{'n':>8}{'FST med':>10}{'FST p95':>10}{'dif med':>10}")
    for b, g in fondo.groupby("bin_MAF", observed=True):
        if len(g) < 20:
            continue
        print(f"{str(b):<16}{len(g):>8}{g['FST'].median():>10.4f}"
              f"{g['FST'].quantile(0.95):>10.4f}{g['dif_abs'].median():>10.4f}")

    # --- percentil de cada candidata --------------------------------------
    filas = []
    for _, r in df[df["es_candidata"]].iterrows():
        b = r["bin_MAF"]
        pares = fondo[fondo["bin_MAF"] == b]
        if len(pares) < 20:
            pares = fondo  # intervalo escaso: usar el fondo completo
            nota = f"intervalo con n<20; percentil contra el fondo completo"
        else:
            nota = ""
        pct_fst = (pares["FST"] < r["FST"]).mean() * 100
        pct_dif = (pares["dif_abs"] < r["dif_abs"]).mean() * 100
        meta = cand[cand["variant_id"] == r["variant_id"]]
        filas.append({
            "gen": meta["gen"].iloc[0] if len(meta) else "",
            "rsid": r["rsids"],
            "variant_id": r["variant_id"],
            "nomenclatura": (meta["nomenclatura_literatura"].iloc[0]
                             if len(meta) and "nomenclatura_literatura" in meta
                             else ""),
            f"AF_{POB_A}": round(r["AF_a"], 4),
            f"AF_{POB_B}": round(r["AF_b"], 4),
            "dif_abs": round(r["dif_abs"], 4),
            "razon_AF": round(r["razon_AF"], 3) if pd.notna(r["razon_AF"]) else None,
            "FST": round(r["FST"], 4),
            "bin_MAF": str(b),
            "n_fondo": len(pares),
            "pct_FST": round(pct_fst, 1),
            "pct_dif_abs": round(pct_dif, 1),
            "nota": nota,
        })

    t4 = pd.DataFrame(filas).sort_values("pct_FST", ascending=False)

    RES.mkdir(parents=True, exist_ok=True)
    cols_f = ["variant_id", "rsids", "consecuencia", "origen", "AF_a", "AF_b",
              "AN_a", "AN_b", "dif_abs", "razon_AF", "FST", "MAF_media", "bin_MAF"]
    fondo[cols_f].to_csv(RES / "tabla3_fondo.tsv", sep="\t", index=False)
    t4.to_csv(RES / "tabla4_percentiles.tsv", sep="\t", index=False)

    print(f"\nTabla 3 -> {RES / 'tabla3_fondo.tsv'}")
    print(f"Tabla 4 -> {RES / 'tabla4_percentiles.tsv'}")

    print(f"\nCandidatas ordenadas por percentil de FST dentro de su "
          f"intervalo de frecuencia:")
    print(f"{'gen':<10}{'rsid':<12}{'AF_amr':>8}{'AF_nfe':>8}"
          f"{'dif':>8}{'FST':>8}{'pctFST':>8}{'n fondo':>9}")
    for _, r in t4.iterrows():
        print(f"{r['gen']:<10}{r['rsid']:<12}{r['AF_amr']:>8.3f}"
              f"{r['AF_nfe']:>8.3f}{r['dif_abs']:>8.3f}{r['FST']:>8.4f}"
              f"{r['pct_FST']:>8.1f}{r['n_fondo']:>9}")

    altos = t4[t4["pct_FST"] >= 95]
    print(f"\n{len(altos)} de {len(t4)} candidatas en el percentil 95 o superior "
          f"de diferenciacion para su intervalo de frecuencia.")


if __name__ == "__main__":
    main()
