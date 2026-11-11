#!/usr/bin/env python3
"""
Compara la tasa de BLOQUES atipicos entre genes candidatos y genes control.

POR QUE
-------
El script 12 compara proporciones de VARIANTES atipicas. Eso sobreestima el
exceso cuando los atipicos se concentran en pocos haplotipos: 19 variantes
contiguas en AGT son una senal, no diecinueve. La comparacion honesta es por
bloques, y normalizada por el tamano de la region evaluada, porque el numero de
bloques posibles depende de cuanta secuencia se examina.

METODO
------
Se agrupan las variantes atipicas de cada grupo en bloques (separacion menor a
VENTANA_LD en el mismo cromosoma) y se compara:
  - bloques por Mb de region evaluada
  - bloques por cada 1000 variantes evaluadas
en candidatos y en control, dentro de cada estrato.

Tambien se reporta el tamano de cada bloque, para distinguir una senal
haplotipica extendida de atipicos dispersos.

LIMITACION: el agrupamiento es por distancia fisica, no por r2. Sin genotipos
individuales no se puede calcular desequilibrio de ligamiento directamente.

Correr: python3 scripts/12b_bloques_comparados.py
Salidas:
  results/bloques_candidatos.tsv
  results/bloques_control.tsv
  results/bloques_resumen.txt
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
TODAS = RAIZ / "results" / "barrido_todas.tsv"
BED = RAIZ / "data" / "raw" / "regiones.bed"
RES = RAIZ / "results"

PCT_UMBRAL = 99.0
VENTANA_LD = 50000

salida = []


def di(t=""):
    print(t)
    salida.append(t)


def bloques(df):
    """Agrupa por proximidad fisica. Devuelve el df con columna 'bloque'."""
    if not len(df):
        return df.assign(bloque=pd.Series(dtype=int))
    d = df.sort_values(["chrom", "pos"]).copy()
    b, actual, pc, pp = [], 0, None, None
    for c, p in zip(d["chrom"].astype(str), d["pos"]):
        if pc != c or (p - pp) > VENTANA_LD:
            actual += 1
        b.append(actual)
        pc, pp = c, p
    d["bloque"] = b
    return d


def mb_por_clase():
    """Mb evaluadas por clase de gen, desde el BED."""
    if not BED.exists():
        return {}
    tot = {}
    with open(BED, encoding="utf-8") as fh:
        for linea in fh:
            if not linea.strip():
                continue
            c, i, f, nombre = linea.split("\t")[:4]
            _, _, clase = nombre.strip().rpartition("_")
            tot[clase] = tot.get(clase, 0) + (int(f) - int(i))
    return {k: v / 1e6 for k, v in tot.items()}


def main():
    if not TODAS.exists():
        sys.exit(f"Falta {TODAS}. Correr antes el script 12.")

    df = pd.read_csv(TODAS, sep="\t", low_memory=False)
    mb = mb_por_clase()

    di("=" * 76)
    di("TASA DE BLOQUES ATIPICOS: CANDIDATOS FRENTE A CONTROL")
    di("=" * 76)
    di(f"\nUmbral: percentil >= {PCT_UMBRAL:.0f}   "
       f"Agrupamiento: < {VENTANA_LD/1000:.0f} kb")
    if mb:
        di(f"\nRegion evaluada: " +
           ", ".join(f"{k} {v:.2f} Mb" for k, v in sorted(mb.items())))

    filas = []
    for est in sorted(df["estrato"].unique()):
        for grp in ["candidato", "control"]:
            s = df[(df["estrato"] == est) & (df["grupo"] == grp)]
            if not len(s):
                continue
            at = s[s["pct"] >= PCT_UMBRAL]
            bl = bloques(at)
            n_bloques = bl["bloque"].nunique() if len(bl) else 0
            mbs = mb.get(grp, np.nan)
            filas.append({
                "estrato": est, "grupo": grp,
                "n_variantes": len(s),
                "n_atipicas": len(at),
                "n_bloques": n_bloques,
                "Mb": round(mbs, 2) if mbs == mbs else "",
                "bloques_por_Mb": round(n_bloques / mbs, 3) if mbs == mbs and mbs else "",
                "bloques_por_1000var": round(1000 * n_bloques / len(s), 3),
                "tam_bloque_max": int(bl.groupby("bloque").size().max()) if n_bloques else 0,
                "tam_bloque_mediano": float(bl.groupby("bloque").size().median()) if n_bloques else 0,
            })
            if grp == "candidato" and n_bloques:
                bl.to_csv(RES / f"bloques_candidatos_{est}.tsv", sep="\t", index=False)
            elif grp == "control" and n_bloques:
                bl.to_csv(RES / f"bloques_control_{est}.tsv", sep="\t", index=False)

    t = pd.DataFrame(filas)
    t.to_csv(RES / "bloques_comparacion.tsv", sep="\t", index=False)

    di(f"\n{'estrato':<8}{'grupo':<11}{'n var':>8}{'atipicas':>10}"
       f"{'bloques':>9}{'Mb':>7}{'blq/Mb':>9}{'blq/1000v':>11}{'tam max':>9}")
    for _, r in t.iterrows():
        di(f"{r['estrato']:<8}{r['grupo']:<11}{r['n_variantes']:>8,}"
           f"{r['n_atipicas']:>10}{r['n_bloques']:>9}{str(r['Mb']):>7}"
           f"{str(r['bloques_por_Mb']):>9}{r['bloques_por_1000var']:>11.3f}"
           f"{r['tam_bloque_max']:>9}")

    di("\nLa columna atipicas sobreestima la senal cuando los atipicos son")
    di("contiguos: un haplotipo diferenciado produce decenas de variantes")
    di("atipicas que son una sola observacion. Las columnas de bloques")
    di("corrigen eso.")

    # --- interpretacion ---------------------------------------------------
    di("\n" + "=" * 76)
    di("RAZON CANDIDATO / CONTROL")
    di("=" * 76)
    di(f"\n{'estrato':<10}{'por variante':>15}{'por bloque/Mb':>16}"
       f"{'por bloque/1000v':>19}")
    for est in sorted(df["estrato"].unique()):
        c = t[(t["estrato"] == est) & (t["grupo"] == "candidato")]
        k = t[(t["estrato"] == est) & (t["grupo"] == "control")]
        if not len(c) or not len(k):
            continue
        c, k = c.iloc[0], k.iloc[0]
        r_var = ((c["n_atipicas"] / c["n_variantes"]) /
                 (k["n_atipicas"] / k["n_variantes"])) if k["n_atipicas"] else np.nan
        r_mb = (c["bloques_por_Mb"] / k["bloques_por_Mb"]) \
            if isinstance(c["bloques_por_Mb"], float) and k["bloques_por_Mb"] else np.nan
        r_1k = (c["bloques_por_1000var"] / k["bloques_por_1000var"]) \
            if k["bloques_por_1000var"] else np.nan
        def f(x):
            return f"{x:.2f}" if x == x else "-"
        di(f"{est:<10}{f(r_var):>15}{f(r_mb):>16}{f(r_1k):>19}")

    di("\nUna razon cercana a 1 por bloque indica que el aparente exceso de")
    di("variantes atipicas se explica por agrupamiento haplotipico y no por")
    di("una mayor densidad de senales independientes en genes candidatos.")

    # --- detalle de bloques candidatos -------------------------------------
    di("\n" + "=" * 76)
    di("BLOQUES EN GENES CANDIDATOS")
    di("=" * 76)
    for est in sorted(df["estrato"].unique()):
        s = df[(df["estrato"] == est) & (df["grupo"] == "candidato")]
        at = s[s["pct"] >= PCT_UMBRAL]
        if not len(at):
            continue
        bl = bloques(at)
        di(f"\nEstrato {est}:")
        di(f"{'gen':<9}{'variante pico':<22}{'n var':>7}{'extension kb':>14}"
           f"{'FST max':>9}{'AF_IMX':>8}{'AF_EUR':>8}")
        col_a = [c for c in bl.columns if c.startswith("AF_IMX")][0]
        col_b = [c for c in bl.columns if c.startswith("AF_EUR")][0]
        for b, g in bl.groupby("bloque"):
            pico = g.loc[g["FST"].idxmax()]
            ext = (g["pos"].max() - g["pos"].min()) / 1000
            di(f"{pico['gen']:<9}{pico['variant_id']:<22}{len(g):>7}"
               f"{ext:>14.1f}{pico['FST']:>9.4f}"
               f"{pico[col_a]:>8.3f}{pico[col_b]:>8.3f}")

    (RES / "bloques_resumen.txt").write_text("\n".join(salida))
    print(f"\n\nComparacion -> {RES / 'bloques_comparacion.tsv'}")
    print(f"Resumen     -> {RES / 'bloques_resumen.txt'}")


if __name__ == "__main__":
    main()
